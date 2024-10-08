FROM rocker/r-ver:4.3.2

RUN apt-get -y update && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
        dirmngr \
        gnupg \
        wget \
        apt-transport-https \
        ca-certificates \
        software-properties-common \
        gcc \
        make \
        libbz2-dev \
        zlib1g-dev \
        libncurses5-dev  \
        libncursesw5-dev \
        liblzma-dev \
        tabix \
        imagemagick \
        libmagick++-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        g++ \
        curl \
        libxml2 \
        libxml2-dev \
        xauth

RUN R --no-echo -e "install.packages('data.table', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('R.utils', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('ggplot2', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('foreach', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('doMC', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('ggthemes', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('ggrepel', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('optparse', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('umap', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('corrplot', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('reshape2', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('magick', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('ggbeeswarm', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('ggdendro', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('ecodist', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('eulerr', repos='http://cran.us.r-project.org')" \
    && R --no-echo -e "install.packages('BiocManager')" \
    && R --no-echo -e "BiocManager::install('org.Hs.eg.db')" \
    && R --no-echo -e "BiocManager::install('clusterProfiler')" \
    && R --no-echo -e "BiocManager::install('STRINGdb')"

# Install python3.8 and pip3
RUN apt-get clean \
    && apt-get update \
    && apt-get install -y python3-pip

# Install mhcflurry
RUN pip3 install mhcflurry tensorflow \
    && unset DEBIAN_FRONTEND
