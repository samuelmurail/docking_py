FROM continuumio/miniconda3

COPY .conda.yml /tmp/environment.yml
COPY setup.py setup.py
COPY docking_py docking_py
RUN conda env create -f /tmp/environment.yml -n docking_py

RUN echo "conda activate docking_py" >> ~/.bashrc
RUN bash -c 'python setup.py install'

ENV PATH /opt/conda/envs/docking_py/bin:$PATH
ENV CONDA_DEFAULT_ENV docking_py