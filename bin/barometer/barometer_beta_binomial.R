#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(glmmTMB))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: barometer_beta_binomial.R COUNTS.tsv RESULTS.tsv")
}

counts <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
required <- c("row_id", "group", "sample", "replicate", "successes", "trials")
if (!all(required %in% names(counts))) {
    stop("Input is missing required columns: ", paste(setdiff(required, names(counts)), collapse = ", "))
}

result_rows <- list()
append_result <- function(row_id, test, group1 = NA_character_, group2 = NA_character_,
                          statistic = NA_real_, estimate = NA_real_, mean1 = NA_real_,
                          mean2 = NA_real_, p_value = NA_real_, status = "ok",
                          dispersion = NA_real_) {
    result_rows[[length(result_rows) + 1]] <<- data.frame(
        row_id = as.character(row_id), test = test, group1 = group1, group2 = group2,
        statistic = statistic, estimate = estimate, mean1 = mean1, mean2 = mean2,
        p_value = p_value, status = status, dispersion = dispersion,
        stringsAsFactors = FALSE
    )
}

for (row_id in unique(counts$row_id)) {
    data <- counts[counts$row_id == row_id, , drop = FALSE]
    data$group <- as.character(data$group)

    data <- aggregate(
        cbind(successes, trials) ~ group + sample,
        data = data,
        FUN = sum
    )
    data <- data[data$trials > 0 & data$successes >= 0 & data$successes <= data$trials, , drop = FALSE]
    groups <- sort(unique(data$group))
    group_sizes <- table(factor(data$group, levels = groups))

    if (length(groups) < 2 || any(group_sizes < 2)) {
        append_result(row_id, "global", status = "insufficient_biological_replicates")
        next
    }

    if (all(data$successes == 0) || all(data$successes == data$trials)) {
        append_result(row_id, "global", status = "no_outcome_variation")
        next
    }

    data$group <- factor(data$group, levels = groups)
    fit_status <- "ok"
    fit_result <- tryCatch({
        fit_formula <- cbind(successes, trials - successes) ~ group
        full_fit <- glmmTMB(
            fit_formula,
            family = betabinomial(link = "logit"),
            data = data,
            control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
        )

        if (full_fit$fit$convergence != 0) {
            full_fit <- glmmTMB(
                fit_formula,
                family = betabinomial(link = "logit"),
                data = data,
                control = glmmTMBControl(
                    optimizer = optim,
                    optArgs = list(method = "BFGS"),
                    optCtrl = list(maxit = 5000, reltol = 1e-10)
                )
            )
        }

        if (full_fit$fit$convergence != 0) {
            stop(paste0(
                "model optimizer did not converge (code ", full_fit$fit$convergence,
                "): ", full_fit$fit$message
            ))
        }
        if (!isTRUE(full_fit$sdr$pdHess)) {
            fit_status <<- "converged_with_dispersion_hessian_warning"
        }

        design <- model.matrix(~ group, data = data.frame(group = factor(groups, levels = groups)))
        coefficients <- glmmTMB::fixef(full_fit)$cond
        covariance <- as.matrix(vcov(full_fit)$cond)
        if (any(!is.finite(covariance))) {
            stop("fixed-effect covariance contains non-finite values")
        }
        condition_coefficients <- coefficients[-1]
        condition_covariance <- covariance[-1, -1, drop = FALSE]
        statistic <- as.numeric(
            t(condition_coefficients) %*% solve(condition_covariance, condition_coefficients)
        )
        global_p <- pchisq(statistic, df = length(condition_coefficients), lower.tail = FALSE)
        means <- plogis(as.numeric(design %*% coefficients))

        list(
            statistic = statistic,
            global_p = global_p,
            design = design,
            covariance = covariance,
            coefficients = coefficients,
            means = means,
            dispersion = as.numeric(sigma(full_fit))
        )
    }, error = function(error) {
        fit_status <<- paste0("fit_failed: ", conditionMessage(error))
        NULL
    })

    if (is.null(fit_result)) {
        append_result(row_id, "global", status = fit_status)
        next
    }

    append_result(
        row_id, "global", statistic = fit_result$statistic,
        p_value = fit_result$global_p, status = fit_status,
        dispersion = fit_result$dispersion
    )

    for (i in seq_along(groups)) {
        append_result(
            row_id, "mean", group1 = groups[[i]], mean1 = fit_result$means[[i]],
            status = fit_status, dispersion = fit_result$dispersion
        )
    }

    if (length(groups) > 1) {
        for (i in seq_len(length(groups) - 1)) {
            for (j in seq.int(i + 1, length(groups))) {
                contrast <- fit_result$design[i, ] - fit_result$design[j, ]
                estimate <- sum(contrast * fit_result$coefficients)
                standard_error <- sqrt(as.numeric(contrast %*% fit_result$covariance %*% contrast))
                z_value <- estimate / standard_error
                p_value <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
                append_result(
                    row_id, "pairwise", group1 = groups[[i]], group2 = groups[[j]],
                    estimate = estimate, mean1 = fit_result$means[[i]],
                    mean2 = fit_result$means[[j]], p_value = p_value,
                    status = fit_status, dispersion = fit_result$dispersion
                )
            }
        }
    }
}

if (length(result_rows) == 0) {
    output <- data.frame(
        row_id = character(), test = character(), group1 = character(), group2 = character(),
        statistic = numeric(), estimate = numeric(), mean1 = numeric(), mean2 = numeric(),
        p_value = numeric(), status = character(), dispersion = numeric()
    )
} else {
    output <- do.call(rbind, result_rows)
}

write.table(output, args[[2]], sep = "\t", quote = FALSE, row.names = FALSE, na = "")