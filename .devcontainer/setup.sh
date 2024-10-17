DEBIAN_FRONTEND=noninteractive
apt-get -y update
apt-get -y install libzmq3-dev jupyter

R -e 'install.packages("data.table")'
R -e 'install.packages("foreach")'
R -e 'install.packages("ggthemes")'
R -e 'install.packages("ggplot2")'
R -e 'install.packages("data.table")'
R -e 'install.packages("optparse")'
R -e 'install.packages("corrplot")'
R -e 'install.packages("umap")'
R -e 'install.packages("magick")'
R -e 'install.packages("ggdendro")'
R -e 'install.packages("ecodist")'
R -e 'install.packages("ggbeeswarm")'
R -e 'install.packages("ggrepel")'
R -e 'install.packages("ggthemes")'
R -e 'install.packages("foreach")'
R -e 'install.packages("reshape2")'
R -e 'install.packages("org.Hs.eg.db")'
R -e 'install.packages("clusterProfiler")'
R -e 'install.packages("pheatmap")'
R -e 'install.packages("cowplot")'
R -e 'install.packages("IRkernel")'
R -e 'IRkernel::installspec(user = FALSE)'

git submodule update --init --recursive

unset DEBIAN_FRONTEND