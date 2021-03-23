FROM continuumio/miniconda3

# Trick to ensure git clone is done if the deposit has been changed
ADD https://api.github.com/repos/samuelmurail/gromacs_py/git/refs/heads/master version.json

COPY .conda.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml -n docking_py
RUN echo "conda activate docking_py" >> ~/.bashrc
COPY setup.py setup.py
COPY README.rst README.rst
COPY HISTORY.rst HISTORY.rst
COPY docking_py docking_py
RUN conda run -n docking_py python setup.py install

ENV PATH /opt/conda/envs/docking_py/bin:$PATH
ENV CONDA_DEFAULT_ENV docking_py
