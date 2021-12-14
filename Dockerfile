FROM nfcore/base:1.14
LABEL authors="Sharada Baindoor" \
    description="Docker container for nextflow workshop"

WORKDIR ./
COPY environment.yml ./
RUN conda env create -f environment.yml && conda clean -a
ENV PATH /opt/conda/envs/test_env/bin:$PATH
