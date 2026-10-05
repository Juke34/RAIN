FROM ghcr.io/prefix-dev/pixi:0.81.0

WORKDIR /app

# Tool-specific dependencies
COPY barometer/pixi.toml barometer/pixi.lock ./
RUN pixi install --locked \
    && ln -sf /app/.pixi/envs/default/bin/* /usr/local/bin/ \
    && Rscript -e 'library(glmmTMB)'

COPY --from=rain-source barometer_beta_binomial.R /opt/rain/barometer_beta_binomial.R
ENV RAIN_BETA_BINOMIAL_SCRIPT=/opt/rain/barometer_beta_binomial.R

CMD ["bash"]