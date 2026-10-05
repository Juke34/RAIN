FROM ghcr.io/prefix-dev/pixi:0.81.0

WORKDIR /app

# Tool-specific dependencies
COPY drip/pixi.toml drip/pixi.lock ./
RUN pixi install --locked \
    && ln -sf /app/.pixi/envs/default/bin/* /usr/local/bin/

CMD ["bash"]
