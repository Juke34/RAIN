FROM condaforge/miniforge3:24.7.1-2

# Tool-specific dependencies
COPY env_pluviometer.yml .
RUN conda env update --file env_pluviometer.yml

CMD ["bash"]
