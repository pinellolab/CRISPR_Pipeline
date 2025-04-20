FROM condaforge/miniforge3:24.11.3-2

# Set environment variables
ENV PATH="/usr/local/bin:$PATH"
# Set NUMBA cache directory to /tmp
ENV NUMBA_CACHE_DIR=/tmp 

# Set working directory
WORKDIR /app

# Copy configuration files
COPY nextflow.yaml /app/
COPY ./bin /usr/local/bin

# Update system and create mamba environment
RUN apt update 
SHELL ["mamba", "run", "--no-capture-output", "-n", "base", "/bin/bash", "-c"]
RUN mamba env create -f nextflow.yaml
RUN mamba run -n nfenv pip install "pyseq-align>=1.0.0"
RUN mamba run -n nfenv pip install kb-python==0.28.2

# RUN echo "mamba activate nfenv" >> ~/.bash_profile

# Make binary file executable
RUN chmod +x /usr/local/bin/*

ENV PATH="/opt/conda/envs/nfenv/bin:$PATH"

# Set the default command
CMD ["/bin/bash", "--login"]

# And the entrypoint
ENTRYPOINT ["mamba", "run", "--no-capture-output", "-n", "nfenv"]
