FROM condaforge/miniforge3:24.7.1-2

COPY env_drip.yml .
RUN conda env update --name base --file env_drip.yml \
    && conda clean --all --yes

CMD ["bash"]
