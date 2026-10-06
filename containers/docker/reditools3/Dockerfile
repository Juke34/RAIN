FROM python:slim-bullseye

RUN apt update

# Install procps to provide ps, required by Nextflow
RUN apt install -y procps

RUN pip install REDItools3==3.7

# Patch reditools.logger: on HPC compute nodes the hostname is often not
# resolvable in DNS, so socket.gethostbyname() raises socket.gaierror and
# crashes `python -m reditools analyze` at startup. Wrap it in try/except,
# mirroring the fix in the REDItools3 main branch (not yet released on PyPI).
RUN python - <<'EOF'
import glob

for path in glob.glob("/usr/local/lib/python3*/site-packages/reditools/logger.py"):
    with open(path) as f:
        src = f.read()
    src = src.replace(
        "        ip_addr = socket.gethostbyname(hostname)\n",
        "        try:\n"
        "            ip_addr = socket.gethostbyname(hostname)\n"
        "        except socket.gaierror:\n"
        "            ip_addr = hostname\n",
    )
    with open(path, "w") as f:
        f.write(src)
    print(f"patched {path}")
EOF

CMD [ "bash" ]
