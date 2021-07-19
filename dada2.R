#reproduce run with seed
set.seed("1234")

#initialize libraries
library(dada2)
library(ShortRead)
library(Biostrings) 

#read filtered dataset
filtFs<- sort(list.files(getwd(), pattern = "_16S_F_filt.fastq.gz", full.names = TRUE))
filtRs<- sort(list.files(getwd(), pattern = "_16S_R_filt.fastq.gz", full.names = TRUE))

#build error profiles and save profiled objects. Will be useful to start here at a later date
errF_randomizeFalse <- learnErrors(filtFs, multithread = TRUE, randomize=FALSE) 
errR_randomizeFalse <- learnErrors(filtRs, multithread = TRUE, randomize=FALSE) 
errF_randomizeTrue <- learnErrors(filtFs, multithread = TRUE, randomize=TRUE) 
errR_randomizeTrue <- learnErrors(filtRs, multithread = TRUE, randomize=TRUE) 

saveRDS(errF_randomizeFalse, "errF_randomizeFalse.rds")
saveRDS(errR_randomizeFalse, "errR_randomizeFalse.rds")
saveRDS(errF_randomizeTrue, "errF_randomizeTrue.rds")
saveRDS(errR_randomizeTrue, "errR_randomizeTrue.rds")

png(file="errF_randomizeTrue.png", width=10, height=10, units="in", res=300)
plotErrors(errF_randomizeTrue, nominalQ = TRUE)
dev.off()

png(file="errF_randomizeFalse.png", width=10, height=10, units="in", res=300)
plotErrors(errF_randomizeFalse, nominalQ = TRUE)
dev.off()


png(file="errR_randomizeTrue.png", width=10, height=10, units="in", res=300)
plotErrors(errR_randomizeTrue, nominalQ = TRUE)
dev.off()

png(file="errR_randomizeFalse.png", width=10, height=10, units="in", res=300)
plotErrors(errR_randomizeFalse, nominalQ = TRUE)
dev.off()

#dereplicate using the error profiles from above, merge to generate a sequence table.
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
dadaFs_randomizeFalse <- dada(derepFs, err = errF_randomizeFalse, multithread = TRUE)
dadaFs_randomizeTrue <- dada(derepFs, err = errF_randomizeTrue, multithread = TRUE)
dadaRs_randomizeTrue <- dada(derepRs, err = errR_randomizeTrue, multithread = TRUE)
dadaRs_randomizeFalse <- dada(derepRs, err = errR_randomizeFalse, multithread = TRUE)

mergers_randomizeFalse <- mergePairs(dadaFs_randomizeFalse, derepFs, dadaRs_randomizeFalse, derepRs, verbose=TRUE)
mergers_randomizeTrue <- mergePairs(dadaFs_randomizeTrue, derepFs, dadaRs_randomizeTrue, derepRs, verbose=TRUE)

seqtab_randomizeTrue <- makeSequenceTable(mergers_randomizeTrue)
dim(seqtab_randomizeTrue)

seqtab_randomizeFalse <- makeSequenceTable(mergers_randomizeFalse)
dim(seqtab_randomizeFalse)

#remove chimeric sequences using multiple approches. Will be useful to base taxonomy on all these.
consensus_seqtab_randomizeFalse_nochim <- removeBimeraDenovo(seqtab_randomizeFalse, method="consensus", multithread=TRUE, verbose=TRUE)
sum(consensus_seqtab_randomizeFalse_nochim)/sum(seqtab_randomizeFalse)

consensus_seqtab_randomizeTrue_nochim <- removeBimeraDenovo(seqtab_randomizeTrue, method="consensus", multithread=TRUE, verbose=TRUE)
sum(consensus_seqtab_randomizeTrue_nochim)/sum(seqtab_randomizeTrue)

pooled_seqtab_randomizeFalse_nochim <- removeBimeraDenovo(seqtab_randomizeFalse, method="pooled", multithread=TRUE, verbose=TRUE)
sum(pooled_seqtab_randomizeFalse_nochim)/sum(seqtab_randomizeFalse)

pooled_seqtab_randomizeTrue_nochim <- removeBimeraDenovo(seqtab_randomizeTrue, method="pooled", multithread=TRUE, verbose=TRUE)
sum(pooled_seqtab_randomizeTrue_nochim)/sum(seqtab_randomizeTrue)

sample_seqtab_randomizeFalse_nochim <- removeBimeraDenovo(seqtab_randomizeFalse, method="per-sample", multithread=TRUE, verbose=TRUE)
sum(sample_seqtab_randomizeFalse_nochim)/sum(seqtab_randomizeFalse)

sample_seqtab_randomizeTrue_nochim <- removeBimeraDenovo(seqtab_randomizeTrue, method="per-sample", multithread=TRUE, verbose=TRUE)
sum(sample_seqtab_randomizeTrue_nochim)/sum(seqtab_randomizeTrue)

saveRDS(seqtab_randomizeTrue, "seqtab_randomizeTrue.rds")
saveRDS(seqtab_randomizeFalse, "seqtab_randomizeFalse.rds")
saveRDS(consensus_seqtab_randomizeFalse_nochim, "consensus_seqtab_randomizeFalse_nochim.rds")
saveRDS(consensus_seqtab_randomizeTrue_nochim, "consensus_seqtab_randomizeTrue_nochim.rds")
saveRDS(sample_seqtab_randomizeFalse_nochim, "sample_seqtab_randomizeFalse_nochim.rds")
saveRDS(sample_seqtab_randomizeTrue_nochim, "sample_seqtab_randomizeTrue_nochim.rds")
saveRDS(pooled_seqtab_randomizeFalse_nochim, "pooled_seqtab_randomizeFalse_nochim.rds")
saveRDS(pooled_seqtab_randomizeTrue_nochim, "pooled_seqtab_randomizeTrue_nochim.rds")

#Taxonomy assignment with assignTaxonomy in dada2
taxa1 <- assignTaxonomy(consensus_seqtab_randomizeFalse_nochim, "GTDB_bac120_arc122_ssu_r95_Genus.fa.gz", multithread=TRUE, tryRC=TRUE)
saveRDS(taxa1, "consensus_tax_randomizeFalse_nochim_GTDB.rds")

taxa2 <- assignTaxonomy(consensus_seqtab_randomizeFalse_nochim, "RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz", multithread=TRUE, tryRC=TRUE)
saveRDS(taxa2, "consensus_tax_randomizeFalse_nochim_RefSeq.rds")

taxa3 <- assignTaxonomy(consensus_seqtab_randomizeTrue_nochim, "GTDB_bac120_arc122_ssu_r95_Genus.fa.gz", multithread=TRUE, tryRC=TRUE)
saveRDS(taxa3, "consensus_tax_randomizeTrue_nochim_GTDB.rds")

taxa4 <- assignTaxonomy(consensus_seqtab_randomizeTrue_nochim, "RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz", multithread=TRUE, tryRC=TRUE)
saveRDS(taxa4, "consensus_tax_randomizeTrue_nochim_RefSeq.rds")

taxa5 <- assignTaxonomy(pooled_seqtab_randomizeFalse_nochim, "GTDB_bac120_arc122_ssu_r95_Genus.fa.gz", multithread=TRUE, tryRC=TRUE)
saveRDS(taxa5, "pooled_tax_randomizeFalse_nochim_GTDB.rds")

taxa6 <- assignTaxonomy(pooled_seqtab_randomizeFalse_nochim, "RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz", multithread=TRUE, tryRC=TRUE)
saveRDS(taxa6, "pooled_tax_randomizeFalse_nochim_RefSeq.rds")

taxa7 <- assignTaxonomy(pooled_seqtab_randomizeTrue_nochim, "GTDB_bac120_arc122_ssu_r95_Genus.fa.gz", multithread=TRUE, tryRC=TRUE)
saveRDS(taxa7, "pooled_tax_randomizeTrue_nochim_GTDB.rds")

taxa8 <- assignTaxonomy(pooled_seqtab_randomizeTrue_nochim, "RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz", multithread=TRUE, tryRC=TRUE)
saveRDS(taxa8, "pooled_tax_randomizeTrue_nochim_RefSeq.rds")

taxa9 <- assignTaxonomy(sample_seqtab_randomizeFalse_nochim, "GTDB_bac120_arc122_ssu_r95_Genus.fa.gz", multithread=TRUE, tryRC=TRUE)
saveRDS(taxa9, "sample_tax_randomizeFalse_nochim_GTDB.rds")

taxa10 <- assignTaxonomy(sample_seqtab_randomizeFalse_nochim, "RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz", multithread=TRUE, tryRC=TRUE)
saveRDS(taxa10, "sample_tax_randomizeFalse_nochim_RefSeq.rds")

taxa11 <- assignTaxonomy(sample_seqtab_randomizeTrue_nochim, "GTDB_bac120_arc122_ssu_r95_Genus.fa.gz", multithread=TRUE, tryRC=TRUE)
saveRDS(taxa11, "sample_tax_randomizeTrue_nochim_GTDB.rds")

taxa12 <- assignTaxonomy(sample_seqtab_randomizeTrue_nochim, "RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz", multithread=TRUE, tryRC=TRUE)
saveRDS(taxa12, "sample_tax_randomizeTrue_nochim_RefSeq.rds")

#Taxonomy assignment with DECIPHER
library(DECIPHER)
load("GTDB_r95-mod_August2020.RData")
tax_info1 <- IdTaxa(test=DNAStringSet(getSequences(consensus_seqtab_randomizeFalse_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info1, "consensus_tax_randomizeFalse_nochim_gtdb.rds")
tax_info2 <- IdTaxa(test=DNAStringSet(getSequences(consensus_seqtab_randomizeTrue_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info2, "consensus_tax_randomizeTrue_nochim_gtdb.rds")

load("RDP_v18_July2020.RData")
tax_info3 <- IdTaxa(test=DNAStringSet(getSequences(consensus_seqtab_randomizeTrue_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info3, "consensus_tax_randomizeTrue_nochim_RDP.rds")
tax_info4 <- IdTaxa(test=DNAStringSet(getSequences(consensus_seqtab_randomizeFalse_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info4, "consensus_tax_randomizeFalse_nochim_RDP.rds")

load("SILVA_SSU_r138_2019.RData")
tax_info5 <- IdTaxa(test=DNAStringSet(getSequences(consensus_seqtab_randomizeTrue_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info5, "consensus_tax_randomizeTrue_nochim_SILVA.rds")
tax_info6 <- IdTaxa(test=DNAStringSet(getSequences(consensus_seqtab_randomizeFalse_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info6, "consensus_tax_randomizeFalse_nochim_SILVA.rds")

load("GTDB_r95-mod_August2020.RData")
tax_info7 <- IdTaxa(test=DNAStringSet(getSequences(pooled_seqtab_randomizeFalse_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info7, "pooled_tax_randomizeFalse_nochim_gtdb.rds")
tax_info8 <- IdTaxa(test=DNAStringSet(getSequences(pooled_seqtab_randomizeTrue_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info8, "pooled_tax_randomizeTrue_nochim_gtdb.rds")

load("RDP_v18_July2020.RData")
tax_info9 <- IdTaxa(test=DNAStringSet(getSequences(pooled_seqtab_randomizeTrue_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info9, "pooled_tax_randomizeTrue_nochim_RDP.rds")
tax_info10 <- IdTaxa(test=DNAStringSet(getSequences(pooled_seqtab_randomizeFalse_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info10, "pooled_tax_randomizeFalse_nochim_RDP.rds")

load("SILVA_SSU_r138_2019.RData")
tax_info11 <- IdTaxa(test=DNAStringSet(getSequences(pooled_seqtab_randomizeTrue_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info11, "pooled_tax_randomizeTrue_nochim_SILVA.rds")
tax_info12 <- IdTaxa(test=DNAStringSet(getSequences(pooled_seqtab_randomizeFalse_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info12, "pooled_tax_randomizeFalse_nochim_SILVA.rds")

load("GTDB_r95-mod_August2020.RData")
tax_info13 <- IdTaxa(test=DNAStringSet(getSequences(sample_seqtab_randomizeFalse_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info13, "sample_tax_randomizeFalse_nochim_gtdb.rds")
tax_info14 <- IdTaxa(test=DNAStringSet(getSequences(sample_seqtab_randomizeTrue_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info14, "sample_tax_randomizeTrue_nochim_gtdb.rds")

load("RDP_v18_July2020.RData")
tax_info15 <- IdTaxa(test=DNAStringSet(getSequences(sample_seqtab_randomizeTrue_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info15, "sample_tax_randomizeTrue_nochim_RDP.rds")
tax_info16 <- IdTaxa(test=DNAStringSet(getSequences(sample_seqtab_randomizeFalse_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info16, "sample_tax_randomizeFalse_nochim_RDP.rds")

load("SILVA_SSU_r138_2019.RData")
tax_info17 <- IdTaxa(test=DNAStringSet(getSequences(sample_seqtab_randomizeTrue_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info17, "sample_tax_randomizeTrue_nochim_SILVA.rds")
tax_info18 <- IdTaxa(test=DNAStringSet(getSequences(sample_seqtab_randomizeFalse_nochim)), trainingSet=trainingSet, strand="both", processors=NULL)
saveRDS(tax_info18, "sample_tax_randomizeFalse_nochim_SILVA.rds")

##################sessionInfo########
sessionInfo()

