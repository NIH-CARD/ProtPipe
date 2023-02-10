ARG VARIANT="16-buster"
FROM mcr.microsoft.com/vscode/devcontainers/javascript-node:0-${VARIANT}

# [Optional] Uncomment if you want to install an additional version of node using nvm
# ARG EXTRA_NODE_VERSION=10
# RUN su node -c "source /usr/local/share/nvm/nvm.sh && nvm install ${EXTRA_NODE_VERSION}"

# [Optional] Uncomment if you want to install more global node modules
# RUN su node -c "npm install -g <your-package-list-here>"

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/* \
 && export GOVERSION=1.17.3 OS=linux ARCH=amd64   \
 && rm -rf /usr/local/go \
 && wget -O /tmp/go${GOVERSION}.${OS}-${ARCH}.tar.gz https://dl.google.com/go/go${GOVERSION}.${OS}-${ARCH}.tar.gz \
 && tar -C /usr/local -xzf /tmp/go${GOVERSION}.${OS}-${ARCH}.tar.gz \
 && echo 'export PATH=$PATH:/usr/local/go/bin' >> ~/.bashrc \
 && . ~/.bashrc \
 && git clone https://github.com/hpcng/singularity.git \
 && cd singularity \
 && git checkout v3.8.4 \
 && ./mconfig \
 && cd ./builddir \
 && make \
 && make install
 