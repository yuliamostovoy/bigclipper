FROM continuumio/miniconda3

LABEL maintainer="Yulia Mostovoy"

# install conda packages
COPY docker/environment.yml /
ENV CONDA_PLUGINS_AUTO_ACCEPT_TOS=yes
RUN conda env create -f /environment.yml && conda clean -a -y
RUN echo "conda activate bigclipper" >> ~/.bashrc
ENV PATH=/opt/conda/envs/bigclipper/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gcloud (used to stream BAMs from GCS)
RUN apt-get update && apt-get install -y curl git-lfs parallel time && rm -rf /var/lib/apt/lists/*
RUN curl -sSL https://sdk.cloud.google.com | bash -s -- --disable-prompts --install-dir=/root

# copy BigClipper scripts from the build context (repository root)
COPY scripts/bigclipper_processbam.py scripts/bigclipper_getclusters.py scripts/bigclipper_verifystream.py scripts/bigclipper_cohort.py /opt/bigclipper/
