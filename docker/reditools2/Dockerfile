FROM condaforge/miniforge3:24.7.1-2
WORKDIR /home/work/

COPY env_reditools2.yml .
RUN conda init ${SHELL} \
    && . ~/.bashrc \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda update conda -y \
    && conda create -n myenv \
    && conda env update -n myenv --file env_reditools2.yml

RUN git clone --depth=1 https://github.com/BioinfoUNIBA/REDItools2.git reditools2
RUN chmod +x reditools2/src/cineca/reditools.py
RUN ln -n reditools2/src/cineca/reditools.py /usr/bin/reditools.py
RUN ln -n reditools2/src/cineca/parallel_reditools.py /usr/bin/parallel_reditools.py

RUN echo "conda activate myenv" >> ~/.bashrc

COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

CMD ["bash"]
