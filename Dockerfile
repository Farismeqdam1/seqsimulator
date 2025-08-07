FROM ubuntu:22.04

LABEL maintainer="SeqSim Development Team"
LABEL description="SeqSimulator - Sequencing Data Simulation Tool"

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    wget \
    curl \
    git \
    build-essential \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install bioinformatics tools
WORKDIR /tmp

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar -xjf samtools-1.17.tar.bz2 && \
    cd samtools-1.17 && \
    ./configure && make && make install && \
    cd .. && rm -rf samtools-1.17*

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 && \
    tar -xjf bcftools-1.17.tar.bz2 && \
    cd bcftools-1.17 && \
    ./configure && make && make install && \
    cd .. && rm -rf bcftools-1.17*

# Install BWA
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make && \
    cp bwa /usr/local/bin/ && \
    cd .. && rm -rf bwa

# Install minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 && \
    tar -xjf minimap2-2.26_x64-linux.tar.bz2 && \
    cp minimap2-2.26_x64-linux/minimap2 /usr/local/bin/ && \
    rm -rf minimap2-2.26*

# Install ART
RUN wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz && \
    tar -xzf artbinmountrainier2016.06.05linux64.tgz && \
    cp art_bin_MountRainier/art_illumina /usr/local/bin/ && \
    rm -rf art*

# Install DWGSIM
RUN git clone --recursive https://github.com/nh13/DWGSIM.git && \
    cd DWGSIM && \
    make && \
    cp dwgsim /usr/local/bin/ && \
    cd .. && rm -rf DWGSIM

# Install Badread
RUN pip3 install badread

# Install PBSIM3
RUN wget https://github.com/yukiteruono/pbsim3/releases/download/v3.0.0/pbsim3-3.0.0-Linux.tar.gz && \
    tar -xzf pbsim3-3.0.0-Linux.tar.gz && \
    cp pbsim3-3.0.0-Linux/pbsim3 /usr/local/bin/ && \
    rm -rf pbsim3*

# Create working directory
WORKDIR /app

# Copy SeqSimulator files
COPY seqsimulator.py /app/
COPY requirements.txt /app/
COPY setup.py /app/
COPY README.md /app/

# Install Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt

# Make the script executable
RUN chmod +x seqsimulator.py

# Create data directory
RUN mkdir -p /data

# Set the entrypoint
ENTRYPOINT ["python3", "/app/seqsimulator.py"]
CMD ["--help"]
