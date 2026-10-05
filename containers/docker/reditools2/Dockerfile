FROM ghcr.io/prefix-dev/pixi:0.81.0

WORKDIR /app

# Tool-specific dependencies
COPY reditools2/pixi.toml reditools2/pixi.lock ./
RUN pixi install --locked \
    && ln -sf /app/.pixi/envs/default/bin/* /usr/local/bin/

# Get and install tools
RUN apt-get update && apt-get install -y --no-install-recommends git \
    && rm -rf /var/lib/apt/lists/*
RUN git clone --depth=1 https://github.com/BioinfoUNIBA/REDItools2.git reditools2
RUN chmod +x reditools2/src/cineca/reditools.py
RUN ln -s /app/reditools2/src/cineca/reditools.py /usr/local/bin/reditools.py
RUN ln -s /app/reditools2/src/cineca/parallel_reditools.py /usr/local/bin/parallel_reditools.py

# Patch hostname retrieval in Reditools2
RUN sed -i.bak '849c\
addresses = netifaces.ifaddresses(interface)\
if netifaces.AF_INET in addresses:\
    ip_address = addresses[netifaces.AF_INET][0]["addr"]\
    hostname = socket.gethostbyaddr(ip_address)[0]\
else:\
    hostname = socket.gethostname()' reditools2/src/cineca/reditools.py

CMD ["bash"]
