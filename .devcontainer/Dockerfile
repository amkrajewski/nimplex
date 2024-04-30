FROM continuumio/miniconda3

RUN apt-get update \
    && apt-get install -y gcc \
    && conda install -y -c conda-forge nim \
    && conda install -y python=3.11 liblapack jupyter numpy pandas plotly \
    && conda install -y -c conda-forge pycalphad \
    && nimble install -y arraymancer nimpy \
    && pip install pqam-rmsadtandoc2023 pathfinding \
    && nim c --d:release --out:examples/nimplex nimplex \
    && nim c --d:release --threads:on --app:lib --out:examples/nimplex.so nimplex \
    && nim c --d:release --threads:on --app:lib --out:examples/utils/plotting.so utils/plotting