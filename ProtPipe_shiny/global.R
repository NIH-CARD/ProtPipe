devtools::load_all()
library(shiny)
library(bslib)
library(clusterProfiler)
library(grid)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(org.Dm.eg.db)
library(org.Ce.eg.db)

# Define static organism mapping to use throughout the app
organism_map <- list(
  "Human"     = list(OrgDb = org.Hs.eg.db, kegg = "hsa"),
  "Mouse"     = list(OrgDb = org.Mm.eg.db, kegg = "mmu"),
  "Rat"       = list(OrgDb = org.Rn.eg.db, kegg = "rno"),
  "Fly"       = list(OrgDb = org.Dm.eg.db, kegg = "dme"),
  "Nematode"  = list(OrgDb = org.Ce.eg.db, kegg = "cel")
)

# increase file size limit
options(shiny.maxRequestSize = 30*1024^2)

package_list = c('ggplot2', 'data.table', 'corrplot', 'umap', 'magick', 'ggdendro', 'ecodist','ggbeeswarm',
                 'ggrepel', 'ggthemes', 'foreach','reshape2','org.Hs.eg.db','clusterProfiler','pheatmap')
if(all((lapply(package_list, require, character.only=TRUE)))) {
  cat("INFO: All packages successfully loaded\n")
} else {
  cat("ERROR: One or more packages not available. Are you running this within the container?\n")
}
