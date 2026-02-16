FROM condaforge/miniforge3:24.7.1-2

# Tool-specific dependencies
COPY env_reditools2.yml .
RUN conda env create -n myenv --file env_reditools2.yml \
    && conda clean --all -f -y \
    && echo "conda activate myenv" >> ~/.bashrc # Auto-activating Conda setup for Python 2.7
    
# Get and install tools
RUN git clone --depth=1 https://github.com/BioinfoUNIBA/REDItools2.git reditools2
RUN chmod +x reditools2/src/cineca/reditools.py
RUN ln -n reditools2/src/cineca/reditools.py /usr/bin/reditools.py
RUN ln -n reditools2/src/cineca/parallel_reditools.py /usr/bin/parallel_reditools.py

ENTRYPOINT ["/bin/bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && conda activate myenv && exec \"$@\"", "--"]

CMD ["bash"]
