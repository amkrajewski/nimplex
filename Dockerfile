FROM mcr.microsoft.com/devcontainers/miniconda

RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

RUN conda install -y -c conda-forge nim && \
    conda install -y python=3.11 liblapack jupyter numpy pandas plotly && \
    nimble install -y arraymancer nimpy && \
    pip install pqam-rmsadtandoc2023