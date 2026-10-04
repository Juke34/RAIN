FROM condaforge/miniforge3:24.7.1-2

COPY env_barometer.yml .
RUN conda env update --name base --file env_barometer.yml \
    && Rscript -e 'library(glmmTMB)' \
    && conda clean --all --yes

COPY --from=rain-source barometer_beta_binomial.R /opt/rain/barometer_beta_binomial.R
ENV RAIN_BETA_BINOMIAL_SCRIPT=/opt/rain/barometer_beta_binomial.R

CMD ["bash"]