FROM ghcr.io/prefix-dev/pixi:0.81.0

WORKDIR /app

# Tool-specific dependencies
COPY pluviometer/pixi.toml pluviometer/pixi.lock ./
RUN pixi install --locked \
    && ln -sf /app/.pixi/envs/default/bin/* /usr/local/bin/

CMD ["bash"]
