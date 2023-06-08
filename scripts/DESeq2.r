if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
library(DESeq2)
library(ggplot2)
library(tidyr)
library(dplyr, quietly = TRUE)
library(data.table, quietly = TRUE)

design=read.table("/groups/dog/stage/enora/isoformswitch/Data/design_lncrnaResist_allSamples.csv", header=TRUE, sep=',')
allCancer_samples=read.table("/groups/dog/nanopore/lncrna_resist_cgo/secondary/1_annexa/allCancer_samples.txt")
GTF="/groups/dog/data/hg38_GRCh38/annotation/Ensembl108/Homo_sapiens.GRCh38.108.chr.UCSCformat.gtf"
FASTA="/groups/dog/stage/enora/isoformswitch/Data/transcriptome_lncrna-resist.fa"
counts_matrice=read.table("/groups/dog/stage/enora/isoformswitch/Data/counts_transcript_filter.full.txt", header=TRUE)
exon_annot="/groups/dog/stage/enora/isoformswitch/Data/extended_annotations_filter.full.gtf"

design=as.data.frame.matrix(design)
design <- unite(design, condition, sample_type, condition, sep = "_")

names(design)[names(design) == "sample_id"] <- "sampleID"
samples=(c(design[,1]))
samples=sort(samples)
samples[1:2] = rev(samples[1:2])
counts_matrice=as.data.frame.matrix(counts_matrice)
counts_matrice = subset(counts_matrice,select = -c(GENEID))
colnames(counts_matrice)=c("",samples)

target = colnames(counts_matrice)[2:25]
design <- design[match(target, design$sampleID),]
rownames(design) <- NULL
rownames(counts_matrice) <- NULL


for (cancer in c("U251", "501Mel", "ADCA72", "PC3")) {
    nrow(counts_matrice[, grepl(cancer, names(counts_matrice))])
    df <- cbind(ensgene = counts_matrice[,1], counts_matrice[, grepl(cancer, names(counts_matrice))])
    df <- df[rowSums(df[,2:7])>0,]
    row.names(df) <- df[,1]
    df <- df[,-1]
    subdesign <- design[design$sampleID %like% cancer, ]
    dds <- DESeqDataSetFromMatrix(
        countData = round(df),
        colData = subdesign,
        design = sampleID ~ condition
    )
    saveRDS(dds, file = paste0(cancer,"_all_DESeq_lncrna_resist.Rds"))
    dds <- DESeq(dds)
    saveRDS(dds, file = paste0(cancer,"DESeq_lncrna_resist.Rds"))
}

