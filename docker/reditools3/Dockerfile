FROM python:slim-bullseye

RUN apt update

# Install procps to provide ps, required by Nextflow
RUN apt install -y procps

RUN pip install REDItools3==3.6

CMD [ "bash" ]
