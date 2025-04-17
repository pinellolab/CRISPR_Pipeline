FROM continuumio/miniconda3

# Set environment variables
ENV PATH="/usr/local/bin:/opt/conda/envs/nfenv/bin:$PATH"
# Set NUMBA cache directory to /tmp
ENV NUMBA_CACHE_DIR=/tmp 

# Set working directory
WORKDIR /app

# Copy configuration files
COPY nextflow.yaml /app/
COPY ./bin /usr/local/bin

# Update system and create conda environment
RUN apt update \
    && conda env create -f nextflow.yaml \
    && conda init bash \
    && conda create -n nfenv \
    && echo "conda activate nfenv" >> ~/.bashrc \
    && chmod +x /usr/local/bin/*  

# Set the default command
CMD ["/bin/bash", "--login"]