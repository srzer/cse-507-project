# docker build -t dreal-pyright Dockerfile
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        curl \
        wget \
        gnupg \
        ca-certificates \
        software-properties-common \
        build-essential \
        python3 \
        python3-pip \
        python3-setuptools && \
    apt-get clean

RUN curl -fsSL https://raw.githubusercontent.com/dreal/dreal4/master/setup/ubuntu/22.04/install.sh | bash -

# python lsp
RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - && \
    apt-get install -y nodejs

RUN npm install -g pyright

# FIXME: latest version released, probably won't change.
ENV DREAL_VERSION=4.21.06.2
ENV PATH="/opt/dreal/${DREAL_VERSION}/bin:${PATH}"

WORKDIR /workspace
