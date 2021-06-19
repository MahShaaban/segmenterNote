FROM bioconductor/bioconductor_docker
ADD ./chromhmmData ./chromhmmData
ADD ./segmenter ./segmenter
RUN Rscript -e "devtools::install('./chromhmmData')"
RUN Rscript -e "devtools::install('./segmenter')"
RUN Rscript -e "BiocManager::install('Gviz')"
RUN Rscript -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg18.knownGene')"
RUN Rscript -e "update.packages(ask = FALSE)"
RUN Rscript -e "install.packages('rmarkdown')"
RUN Rscript -e "install.packages('tinytex')"
