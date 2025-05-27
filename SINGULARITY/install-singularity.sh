#!/usr/bin/env bash


# Ensure repositories are up-to-date
sudo DEBIAN_FRONTEND=noninteractive apt-get update

# Install debian packages for dependencies
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
    autoconf \
    automake \
    cryptsetup \
    fuse2fs \
    git \
    fuse \
    libfuse-dev \
    libseccomp-dev \
    libtool \
    pkg-config \
    runc \
    squashfs-tools \
    squashfs-tools-ng \
    uidmap \
    wget \
    zlib1g-dev


# Remove previous GO installation; reinstall updated GO
sudo rm -rf /usr/local/go

export VERSION=1.24.3 OS=linux ARCH=amd64  # change this as you need

wget -O /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz \
  https://dl.google.com/go/go${VERSION}.${OS}-${ARCH}.tar.gz
sudo tar -C /usr/local -xzf /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz

echo 'export PATH=$PATH:/usr/local/go/bin' >> ~/.bashrc
source ~/.bashrc


# Clone singularity repo
git clone --recurse-submodules https://github.com/sylabs/singularity.git
cd singularity

git submodule update --init

git checkout --recurse-submodules v4.3.1

./mconfig --without-libsubid
make -C builddir
sudo make -C builddir install