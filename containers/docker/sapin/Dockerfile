FROM python:slim-bullseye

RUN apt update

# Install procps to provide ps, required by Nextflow
RUN apt install -y procps

# Install Git to clone the repo
RUN apt install -y git
RUN git clone --depth=1 https://github.com/Juke34/SAPiN.git

# Use a pip constraint file to pin dependencies for reproducibility
COPY constraints.txt SAPiN/

RUN cd SAPiN && \
    pip install --constraint constraints.txt .

CMD [ "bash" ]
