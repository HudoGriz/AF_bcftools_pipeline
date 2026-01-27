FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# -----------------------------
# 1) System dependencies
# -----------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    ca-certificates \
    curl \
    wget \
    git \
    unzip \
    bzip2 \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libgsl-dev \
    perl \
    libperl-dev \
    openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------
# 2) Build htslib 1.19
# -----------------------------
RUN wget -q https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2 \
    && tar -xjf htslib-1.19.tar.bz2 \
    && cd htslib-1.19 \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && make install \
    && cd / \
    && rm -rf /htslib-1.19 /htslib-1.19.tar.bz2

# -----------------------------
# 3) Build bcftools 1.19 (GSL + perl-filters)
# -----------------------------
RUN wget -q https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 \
    && tar -xjf bcftools-1.19.tar.bz2 \
    && cd bcftools-1.19 \
    && ./configure --prefix=/usr/local --enable-libgsl --enable-perl-filters \
    && make -j"$(nproc)" \
    && make install \
    && cd / \
    && rm -rf /bcftools-1.19 /bcftools-1.19.tar.bz2

ENV BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools

RUN ldconfig

RUN bcftools --version && bcftools plugin -l

# -----------------------------
# 4) Install Nextflow
# -----------------------------
ENV NXF_HOME=/opt/nextflow
RUN mkdir -p "${NXF_HOME}" \
    && curl -s https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/nextflow \
    && chmod +x /usr/local/bin/nextflow \
    && nextflow -version

# -----------------------------
# 5) Copy pipeline into image
# -----------------------------
WORKDIR /pipeline
COPY main.nf /pipeline/main.nf
COPY nextflow.config /pipeline/nextflow.config

# sceVCF must exist in build context
COPY sceVCF /usr/local/bin/sceVCF
RUN chmod +x /usr/local/bin/sceVCF

# -----------------------------
# 6) Default runtime settings
# -----------------------------
WORKDIR /work
# keep the container from dying 
# CMD ["sleep","356d"] 