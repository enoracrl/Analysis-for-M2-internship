# ---------------------------------------------------------------------------- #
#                             IsoformSwitchAnalyzeR                            #
# ---------------------------------------------------------------------------- #




# Setup R Environment

### Load library
suppressMessages(library(IsoformSwitchAnalyzeR))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))
packageVersion("dplyr")
packageVersion("tidyr")
packageVersion("data.table")
packageVersion("IsoformSwitchAnalyzeR")

##########################################################################################################

#args            <- commandArgs(trailingOnly = TRUE)
#exon_annot      <- args[1]
#transcriptome   <- args[2]
#design          <- args[3]
#counts          <- args[4]
#comparisons     <- args[5]

# Experimental Data
design      <- as.data.frame.matrix(read.table("/groups/dog/stage/enora/isoformswitch/Data/design_lncrnaResist_allSamples.csv", header=TRUE, sep=','))

message(paste("Design"))
print(head(design))

types       <- unique(design$sample_type)
design      <- unite(data=design, condition, sample_type, condition, sep = "_")
design      <- subset(design, select = c(sample_id, condition))
names(design)[names(design) == "sample_id"] <- "sampleID"

# Header : samples names to match between the design file and the counts file
samples         <- (c(design[,1]))
samples         <- sort(samples)
samples[1:2]    <- rev(samples[1:2])

# Isoforms counts

counts              <- as.data.frame.matrix(read.table("/groups/dog/stage/enora/isoformswitch/Data/counts_transcript_filter.full.txt", header=TRUE))
counts              <- subset(counts, select = -c(GENEID))
colnames(counts)    <- c("isoform_id",samples)

message(paste("Count Matrix"))
print(head(counts))

# Comparisons table
comparisons <- data.frame(read.table("/groups/dog/stage/enora/isoformswitch/Data/comparisons.txt", header=TRUE))

message(paste("Comparisons"))
print(head(comparisons))

message(paste("Transcriptome path"))
transcriptome="/groups/dog/stage/enora/isoformswitch/Data/transcriptome_lncrna-resist.fa"
print(head(transcriptome))

message(paste("Annotation path"))
exon_annot="/groups/dog/stage/enora/isoformswitch/Data/extended_annotations_filter.full.gtf"
print(head(exon_annot))

#txdb <- makeTxDbFromGFF("/groups/dog/stage/enora/isoformswitch/Data/extended_annotations_filter.full.gtf", format="gtf")
#tx <- data.frame(trAanscripts(txdb))
#names(tx)[names(tx) == 'tx_name'] <- 'isoform_id'
#tx <- tx[,-5]

mainDir <- "/groups/dog/stage/enora/isoformswitch/"
subDir  <- "outputISA"

dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

# ##################### #
# STEP 2 - LOADING DATA #
# ##################### #

message('Importing data')
SwitchList <- importRdata(
    ### Core arguments
    isoformCountMatrix = counts,
    isoformRepExpression = NULL,
    detectUnwantedEffects = FALSE,
    designMatrix = design,
    isoformExonAnnoation = exon_annot,
    isoformNtFasta = transcriptome,
    comparisonsToMake = comparisons,
    showProgress = FALSE,
    quiet = TRUE
)
saveRDS(SwitchList, file = "SwitchList.Rds")


for (condition in types) {

    message(paste0(condition, 'subseting'))
    subsetSwitch <- subsetSwitchAnalyzeRlist(
        switchAnalyzeRlist=SwitchList,
        SwitchList$isoformFeatures$condition_1 %in% c(
            paste0(condition, "_sensitive"),
            paste0(condition, "_resistant")
        ) &
        SwitchList$isoformFeatures$condition_2 %in% c(
            paste0(condition, "_sensitive"),
            paste0(condition, "_resistant")
        )
    )

    message(paste0(condition, 'filtering'))
    subsetSwitch <- preFilter(
        switchAnalyzeRlist= subsetSwitch,
        geneExpressionCutoff = 1,
        isoformExpressionCutoff = 0,
        quiet=TRUE
    )
    saveRDS(subsetSwitch, file = paste0(condition,"_lncrna_resist.Rds"))

    
    subsetSwitchDE <- isoformSwitchTestDEXSeq(
        ### Core arguments
        switchAnalyzeRlist= subsetSwitch,
        alpha = 0.05,
        dIFcutoff = 0.1,
        quiet=TRUE
    )

    SwitchDE <- subsetSwitchDE$isoformFeatures[abs(
        subsetSwitchDE$isoformFeatures$dIF) > 0.1  & 
        subsetSwitchDE$isoformFeatures$isoform_switch_q_value < 0.05, ]
    #SwitchDE <- SwitchDE %>% mutate(
    #    isoform_annot = case_when(
    #        startsWith(isoform_id, "BambuTx") ~ "NOVEL",
    #        startsWith(isoform_id, "ENST") ~ "KNOWN"
    #        )) %>%
    #    mutate(gene_annot = case_when(
    #        startsWith(gene_id, "BambuGene") ~ "NOVEL",
    #        startsWith(gene_id, "ENSG") ~ "KNOWN"
    #    ))
    #SwitchDE <- merge(SwitchDE, tx, by = "isoform_id")
    saveRDS(SwitchDE, file = paste0(condition, "isoformFeatures.Rds"))
    saveRDS(subsetSwitchDE, file = paste0(condition,"_DEXSeq_lncrna_resist.Rds"))
}

# All

SwitchList <- preFilter(
    switchAnalyzeRlist= SwitchList,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    quiet=TRUE
)
saveRDS(SwitchList, file = "SwitchList_filtered_lncrna_resist.Rds")

SwitchList <- isoformSwitchTestDEXSeq(
    ### Core arguments
    switchAnalyzeRlist= SwitchList,
    alpha = 0.05, 
    dIFcutoff = 0.1,
    quiet=TRUE)
saveRDS(SwitchList, file = "SwitchList_DEXSeq_lncrna_resist.Rds")
