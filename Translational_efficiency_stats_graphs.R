#### Translational_efficiency_stats_graphs.R
## 
## Maria Mikedis mikedis@wi.mit.edu 
## Run as: Rscript Translational_efficiency_stats_graphs.R








library(ggplot2)
library(car) ## to assess log linear regression model later using vif() function
library(DESeq2)
library(RColorBrewer)
library(tximport)
library(dplyr)
library(MASS) ### for stepAIC() function

####################
## "Theme" for figures
####################

theme_mk <- theme(text=element_text(family="Helvetica",face="plain",size=6),
                  axis.text=element_text(size=6),
                  axis.title=element_text(size=6),
                  axis.line.y=element_blank(),
                  axis.line.x=element_line(size=0.25),
                  axis.ticks=element_line(size=0.25),
                  legend.title=element_blank(),
                  legend.background=element_rect(fill="transparent"),
                  legend.text=element_text(size=6), plot.title=element_text(size=6, face="plain"),
                  panel.background = element_blank(),
                  legend.key=element_blank())



### obtain counts from kallisto
### normalize counts via DESeq
### find transcripts that are enriched in Ribotag sample over wildtype sample
### use kallisto to convert normalized counts to TPMs
### proceed wtih analysis



###############
### remove everything except protein coding transcripts (i.e., if a gene has both protein coding and noncoding isoforms, remove noncoding isoforms) and re-normalize TPM to 1e6
###############

### list of transcripts to keep
genes.to.keep = read.table("/lab/solexa_page/maria/genome_files/kallisto/for_ASPeak_analysis/mm10_GRCm38_refGene_mRNA_names_no.ncRNAs.txt", stringsAsFactors = FALSE)


for (folder in c("../kallisto/AGCGATAG-ATAGAGGC.trimmed", 
                "../kallisto/ATTACTCG-ATAGAGGC.trimmed", 
                "../kallisto/ATTCAGAA-ATAGAGGC.trimmed", 
                "../kallisto/CGCTCATT-ATAGAGGC.trimmed", 
                "../kallisto/CGCTCATT-ATAGAGGC.trimmed", 
                "../kallisto/CGGCTATG-ATAGAGGC.trimmed", 
                "../kallisto/CTGAAGCT-ATAGAGGC.trimmed", 
                "../kallisto/GAATTCGT-ATAGAGGC.trimmed", 
                "../kallisto/GAGATTCC-ATAGAGGC.trimmed", 
                "../kallisto/TAATGCGC-ATAGAGGC.trimmed", 
                "../kallisto/TCCGCGAA-ATAGAGGC.trimmed", 
                "../kallisto/TCCGGAGA-ATAGAGGC.trimmed",
                "../kallisto/TCTCGCGC-ATAGAGGC.trimmed")
      ) {
   
  ### load transcript files
  file = read.delim(paste(folder, "/abundance.tsv", sep = "", collapse = NULL), stringsAsFactors = FALSE)

  ### assign rownames
  row.names(file) = file$target_id
  file = file[!names(file) %in% c("target_id")]


  ## select only those transcripts in our "to keep" list
  ## trim to TPM data

  file.trim = file[row.names(file) %in% genes.to.keep$V1 , ]
  file.trim2 = file.trim[names(file.trim) %in% c("target_id", "tpm")]


  ### renormalize TPM data
  file.trim.norm = file.trim2/colSums(file.trim2)*1000000

  ### remove old TPM column; select only first three columns
  file = file[,1:3]

  ### add renormalized TPM column
  file.new = merge(file, file.trim.norm, by="row.names")

  ### rename Row.names to target_id
  colnames(file.new)[1] = "target_id"

  ### write files
  write.table(file.new, file=paste(folder,"/abundance_coding.mRNA.only.tsv", sep = "", collapse = NULL), sep="\t", quote=F, col.names = T, row.names=F)

}



annofile <- read.delim("/lab/solexa_page/maria/genome_files/kallisto/for_ASPeak_analysis/mm10_GRCm38_refGene_mRNA_names_no.ncRNAs.txt", stringsAsFactors = FALSE)
tx2gene <- data.frame("TXNAME"=annofile$X.name, "GENEID"=annofile$name2)
dir <- "/lab/solexa_page/maria/Ribotag_undiff_gonia/180807_kallisto/kallisto_output_180904"

### for control input samples only
samples_all <- read.table("/lab/solexa_page/maria/Ribotag_undiff_gonia/180807_kallisto/kallisto_output_180904/MM291_samplenames_180807",stringsAsFactors = FALSE)
files <- file.path(dir,"kallisto",samples_all$V1,"abundance_coding.mRNA.only.tsv")
names(files) = c(samples_all$V2)
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene)




counts = txi.kallisto.tsv$counts 

synch.counts = read.table("/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/20180730_3S_RNAseq_kallisto/190622_kallisto/kallisto/mRNA_coding_only/MM257_3S.undiff.gonia_kallisto_mm10_GRCm38_coding.mRNA.only_counts_180821", header=TRUE, sep="\t", row.names=1, stringsAsFactors = FALSE)
colnames(synch.counts)[1] = "MM257.1b_3S"
colnames(synch.counts)[2] = "MM257.1a_3S"


### merge counts with synch.counts
counts_data = merge(counts, synch.counts, by="row.names", all=TRUE) 
row.names(counts_data) = counts_data$Row.names 
counts_data = counts_data[ , !(names(counts_data) %in% c("Row.names"))] ## to drop extra column with row.names



################
### Identification of transcripts that are enriched in Ribotag relative to control
### using Deseq for normalization
################

### using count output from tximport; rounding estimated counts to nearest integer
data <- round(counts_data, digits=0)


batch <- c("1", "1",
            "1", "1",
            "1", "1",
            "1", "1",
            "1", "1",
            "1", "1",
            "0", "0")


Ribotag.specific <- c("0", "1",
                      "0", "0",
                      "0", "1",
                      "0", "0",
                      "0", "1",
                      "0", "0",
                      "0", "0")

Ribotag.nonspecific <- c("0", "1",
                      "0", "1",
                      "0", "1",
                      "0", "1",
                      "0", "1",
                      "0", "1",
                      "0", "0")

sample <- c("1", "1",
            "2", "2",
            "3", "3",
            "4", "4",
            "5", "5",
            "6", "6",
            "7", "8")

germ.cell.specific <- c("0", "1",
                      "0", "0",
                      "0", "1",
                      "0", "0",
                      "0", "1",
                      "0", "0",
                      "1", "1")

### I can't include batch in model because I get error message "Model matrix not full rank"; this is a limitation of the model

meta <- data.frame(Ribotag.specific=Ribotag.specific, Ribotag.nonspecific=Ribotag.nonspecific, germ.cell.specific=germ.cell.specific)
row.names(meta) = colnames(data)
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~Ribotag.specific+Ribotag.nonspecific+germ.cell.specific)
	### Adding batch or sample to the design creates a model matrix that is not full rank because batch and sample are linear combinations of other variables

dds <- DESeq(dds)
colnames(dds) = colnames(data)




### obtain normalized counts
norm.counts = counts(dds, normalized=TRUE)
write.table(norm.counts, "MM291.trans.eff_MM257.3S.undiff.gonia_DESeq_counts.norm_190312", sep="\t", quote=FALSE, row.names=TRUE)


dds.Ribo = results(dds, contrast=c("Ribotag.specific", "1", "0"))
write.table(dds.Ribo, "MM291.trans.eff_MM257.3S.undiff.gonia_DESeq_Ribotag.specific_190312", sep="\t", quote=FALSE, row.names=TRUE)
	## translational efficiency is equivalent to the log2FoldChange for variable Ribotag.specific

dds.Ribo = as.data.frame(dds.Ribo)



### calculate TE via TPM only

abundance = as.data.frame(txi.kallisto.tsv$abundance)

synch.TPM = read.table("/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/20180730_3S_RNAseq_kallisto/190622_kallisto/kallisto/mRNA_coding_only/MM257_3S.undiff.gonia_kallisto_mm10_GRCm38_coding.mRNA.only_TPMs_180821", header=TRUE, sep="\t", row.names=1, stringsAsFactors = FALSE)
colnames(synch.counts)[1] = "MM257.1b_3S"
colnames(synch.counts)[2] = "MM257.1a_3S"


### merge counts with synch.counts
TPM_data = merge(abundance, synch.TPM, by="row.names", all=TRUE) 
row.names(TPM_data) = TPM_data$Row.names 
TPM_data = TPM_data[ , !(names(TPM_data) %in% c("Row.names"))] ## to drop extra column with row.names

write.table(TPM_data, "MM291.trans.eff_MM257.3S.undiff.gonia_TPMs_not.filtered_190718", sep="\t", quote=FALSE, row.names=TRUE)



#################
### Obtain a master data file for analysis
#################
### merge trans.eff with cds.length, CAI, 3S TPM, utr3.length


cds.length = read.table("CDS/mm10_GRCm38_cds_length.txt", sep="\t", stringsAsFactors = FALSE, header=T) ### file containing length of coding region (nt) for all coding transcripts
### modify header for mergining later
colnames(cds.length)[1]="transcript.id"



CAI = read.table("CDS/CDS_seqs_expressed_CAI", sep="\t")

### TPM file contains TPMs for coding mRNA transcripts only, renormalized after pseudoaligning to all transcripts
### file only contains the most robustly isoform per gene
synch.TPM = read.table("/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/20180730_3S_RNAseq_kallisto/190622_kallisto/kallisto/mRNA_coding_only/MM257_3S.undiff.gonia_kallisto_mm10_GRCm38_coding.mRNA.only_TPMs_180821_genes.transcripts_min1_most.expressed.isoform.per.gene", header=TRUE, sep="\t", stringsAsFactors = FALSE)


### file of 3' UTR lengths
utr3.length = read.table("/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/20180730_ASPeak/ASPeak_IgG_RNAseq_gapno0_CL.site.only_w.intergenic/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/utr3_metaplot/3UTRs_length.txt", header=T)
## modify headers
colnames(utr3.length) = c("transcript.id", "utr3_length")




CAI$transcript.id = row.names(CAI)


### add column for merging
dds.Ribo$gene.id = row.names(dds.Ribo)

norm.counts2 = as.data.frame(norm.counts)  ### need to do this to avoid error when creating a new column out of row.names
norm.counts2$gene.id = row.names(norm.counts2)




### merging Deseq output with TPM file, coding region length file, and CAI file
temp1 = merge(dds.Ribo, synch.TPM, by="gene.id")
temp2 = merge(temp1, CAI, by="transcript.id")
temp3 = merge(temp2, cds.length, by="transcript.id")
temp4 = merge(temp3, norm.counts2, by="gene.id")
temp5= merge(temp4, utr3.length, by="transcript.id")

trans.eff =unique(na.omit(temp5))

## Filter out genes that cannot be reliably detected via RNAseq and via Ribotag IP-seq (RiboTag IP samples must add up to 15 normalized counts; 3S RNAseq samples must add up to 10 counts)
  ## Ribotag filter
Ribotag_sums = (rowSums(trans.eff[, c("MM291.1.Ribotag.IP", "MM291.4.Ribotag.IP", "MM291.7.Ribotag.IP")]) >= 15)
filtering = trans.eff[Ribotag_sums,]

  ## 3S RNAseq counts filter
X3S_sums = (rowSums(filtering[, c("MM257.1b_3S", "MM257.1a_3S")]) >= (10))
filtering2 = filtering[X3S_sums,]


  ## 3S RNAseq TPM filter
X3S_TPM = (filtering2$mean.gene.TPM >=1)
filtering3 = filtering2[X3S_TPM,]

## NM_001099920 (Gm15093) and NM_001122734 (Gm15085) each have 2 CDS lengths: 1345 and 1350. Only 1350 is divisible by 3 so select this length for each isform.
## make sure all othe rtranscripts are divisible by 3

filtering4 = filtering3[order(filtering3$transcript.id, -filtering3$CDS_length),]
filtering5 = filtering4[match(unique(filtering4$transcript.id), filtering4$transcript.id),]




#### final dataset of genes with translational efficiency values
trans.eff.filtered=filtering5

### mark transcripts as Dazl targets or nontargets
Dazl_targets = read.table("utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_gene.ID.only")

Dazl_data = trans.eff.filtered[trans.eff.filtered$gene.id %in% Dazl_targets$V1 ,]
Dazl_background_data = trans.eff.filtered[!trans.eff.filtered$gene.id  %in% Dazl_targets$V1 ,]

Dazl_data$Dazl_target = "DAZL target"
Dazl_background_data$Dazl_target = "Nontarget"

total.data = rbind(Dazl_background_data, Dazl_data)

write.table(total.data, "MM291.trans.eff_MM257.3S.undiff.gonia_DESeq_Ribotag.specific_count.TPM.filtered_190618", sep="\t", quote=FALSE, row.names=F)




## extract only columns for model; write to file
mycol <- c("gene.id", "transcript.id", "log2FoldChange", "mean.gene.TPM", "CDS_length", "CAI", "utr3_length","Dazl_target")
temp4 <- total.data[, mycol]

## add DESeq-based TE as a non-log2 value
temp4$trans.eff.Ribotag.DESeq = 2^(temp4$log2FoldChange)

mycol2 <- c("gene.id", "transcript.id", "trans.eff.Ribotag.DESeq", "mean.gene.TPM", "CDS_length", "CAI", "utr3_length","Dazl_target")
total.data.filtered <- temp4[, mycol2]

write.table(total.data.filtered, "TE.from.DESeq_counts_cds.length_CAI_Dazl.target_190618.txt", sep="\t", quote=FALSE, row.names=F)










### Does my data show low translation for histones, which are known to be translated in a replication-dependent manner?
#  For histone analysis, the gene list from GO terms0006335 “DNA replication-dependent nucleosome assembly” was downloaded from AmiGO 2 v.2.5.12, last updated 2019-07-02 (DOI 10.5281//zenodo.3267438) and was filtered for histone genes. 
### boxplot comparing histones TE to other TEs calculated from DESeq

### mark transcripts as histones or non histones
histones_rep.dep = read.table("GO-0006335_DNA_replication-dependent_nucleosome_assembly_histones.only", sep="\t", header=F, stringsAsFactors = FALSE)




histone_data_rep.dep = total.data.filtered[total.data.filtered$gene.id %in% histones_rep.dep$V2 ,]
histone_background_data = total.data.filtered[!total.data.filtered$gene.id %in% histones_rep.dep$V2 ,]

histone_data_rep.dep$histone_target = "Replication-dependent histones"
histone_background_data$histone_target = "Other"

total.histone.data = rbind(histone_background_data, histone_data_rep.dep)
  ## order data
total.histone.data$histone_target = factor(total.histone.data$histone_target, c("Other", "Replication-dependent histones"))

### are different variables that affect translation different between Dazl targets and nontargets?
# Function to produce summary statistics (median and interquartile range)
data_summary <- function(x) {
   m <- median(x)
   ymin <- unname(quantile(x, 0.25))
   ymax <- unname(quantile(x, 0.75))
   return(c(y=m,ymin=ymin,ymax=ymax))
}

ggplot(total.histone.data, aes(histone_target, (log2(trans.eff.Ribotag.DESeq)),fill=histone_target, colour=histone_target)) +
  geom_violin() +
  stat_summary(fun.data=data_summary, colour="black", size=0.05) +
  ylab("Translational efficiency (log2)") +
  scale_fill_manual(values=c("gray", "#9e9ac8")) +
  scale_colour_manual(values=c("gray", "#9e9ac8")) +
  coord_cartesian(ylim=c(-10, 8)) +
  scale_y_continuous(breaks=c(-10, -8, -6, -4, -2, 0,2,4,6,8)) + 
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))

ggsave("FIGURE_trans.eff_DESeq_boxplot_histones.v.others_190313.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)






## test for statistical difference between histones and other transcripts using the 2 different TEs
my.stats= paste("Statistical difference between TE calculated via DESeq\n")
my.stats = append(my.stats, capture.output(ks.test(histone_data_rep.dep$trans.eff.Ribotag.DESeq, histone_background_data$trans.eff.Ribotag.DESeq)))
my.stats = append(my.stats, capture.output(wilcox.test(histone_data_rep.dep$trans.eff.Ribotag.DESeq, histone_background_data$trans.eff.Ribotag.DESeq, alternative="less")))


writeLines(my.stats, "TE_DESeq_histones.v.others_ks.test_wilcox.test_190313")






### boxplot comparing ribosome TEs to other TEs calculated from DESeq

### mark transcripts as ribosomes or non ribosomes
ribosomes = read.table("GO_term003735_structural_constituent_of_ribosome_20190404_gene.id.only.txt")

ribosome_data = total.data.filtered[total.data.filtered$gene.id %in% ribosomes$V1 ,]
ribosome_background_data = total.data.filtered[!total.data.filtered$gene.id %in% ribosomes$V1 ,]

ribosome_data$ribosome_target = "Ribosome"
ribosome_background_data$ribosome_target = "Other"

total.ribosome.data = rbind(ribosome_background_data, ribosome_data)

  ## order data
total.ribosome.data$ribosome_target = factor(total.ribosome.data$ribosome_target, c("Other", "Ribosome"))



ggplot(total.ribosome.data, aes(ribosome_target, 2^(trans.eff.Ribotag.DESeq),fill=ribosome_target)) +
  geom_boxplot(outlier.size=0, outlier.stroke=0, lwd=0.2) +
  ylab("Translational efficiency (DESeq)") +
  scale_fill_manual(values=c("gray", "#9e9ac8")) +
  coord_cartesian(ylim=c(0, 12)) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9, 10, 11, 12)) + 
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))

ggsave("FIGURE_trans.eff_DESeq_boxplot_ribosomes.v.others_190313.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)





## test for statistical difference between ribosomes and other transcripts using the 2 different TEs
my.stats= paste("Statistical difference between TE calculated via DESeq\n")
my.stats = append(my.stats, capture.output(ks.test(ribosome_data$trans.eff.Ribotag.DESeq, ribosome_background_data$trans.eff.Ribotag.DESeq)))


writeLines(my.stats, "TE_DESeq_ribosomes.v.others_ks.test_wilcox.test_190313")



### Does my data show that more robustly expressed transcripts are more highly translated than less robustly expressed transcripts? (sanity check)
### How do the other variables (CDS length, CAI) correlate with TE?
### Do the correlations looks different between DAZL targets and nontargets?

    ### assign Dazl targets 1, nontargets 0 for point biserial correlation (for a categorical vs. conitinous correlation comparison)

mm = total.data.filtered
i=1 
while (i<=nrow(mm)) {
  if (mm$Dazl_target[i]=="Nontarget") {
    mm$Dazl[i]=as.numeric(0)
  }
  if (mm$Dazl_target[i]=="DAZL target") {
    mm$Dazl[i]=as.numeric(1)
  }
  i=i+1
}

## check for correlation of each variable to TE with Pearson's 
out1 = capture.output(cor.test(log2(total.data.filtered$mean.gene.TPM) , log2(total.data.filtered$trans.eff.Ribotag.DESeq)))
out2 = capture.output(cor.test(log2(total.data.filtered$CDS_length) , log2(total.data.filtered$trans.eff.Ribotag.DESeq)))
out3 = capture.output(cor.test(log2(total.data.filtered$CAI) , log2(total.data.filtered$trans.eff.Ribotag.DESeq)))  
out4 = capture.output(cor.test(log2(total.data.filtered$utr3_length) , log2(total.data.filtered$trans.eff.Ribotag.DESeq)))
out5 = capture.output(cor.test(as.numeric(mm$Dazl) , log2(mm$trans.eff.Ribotag.DESeq))) ### point biserial coefficient
out = append(out1, out2)
out = append(out, out3)
out = append(out, out4)
out = append(out, out5)

writeLines(out, "TE_correlations_independent.variables_Dazl.targets.and.nontargets")


## check for correlation of each variable to TE
out1 = capture.output(cor.test(log2(total.data.filtered$mean.gene.TPM) , log2(total.data.filtered$trans.eff.Ribotag.DESeq),  method = c("spearman")))
out2 = capture.output(cor.test(log2(total.data.filtered$CDS_length) , log2(total.data.filtered$trans.eff.Ribotag.DESeq),  method = c("spearman")))
out3 = capture.output(cor.test(log2(total.data.filtered$CAI) , log2(total.data.filtered$trans.eff.Ribotag.DESeq),  method = c("spearman")))  
out4 = capture.output(cor.test(log2(total.data.filtered$utr3_length) , log2(total.data.filtered$trans.eff.Ribotag.DESeq),  method = c("spearman")))
out5 = capture.output(cor.test(as.numeric(mm$Dazl) , log2(mm$trans.eff.Ribotag.DESeq))) ### point biserial coefficient
out = append(out1, out2)
out = append(out, out3)
out = append(out, out4)
out = append(out, out5)

writeLines(out, "TE_correlations_independent.variables_Dazl.targets.and.nontargets_Spearmans")



## are TPM, CDS length, CAI, and 3' UTR length statistically different between DAZL targets and nontargets?
writeLines(capture.output(ks.test(Dazl_data$CDS_length, Dazl_background_data$CDS_length)), "DAZL.targets.v.all.nontargets_CDS_length_kstst")
writeLines(capture.output(ks.test(Dazl_data$CAI, Dazl_background_data$CAI)), "DAZL.targets.v.all.nontargets_CAI_kstst")
writeLines(capture.output(ks.test(Dazl_data$mean.gene.TPM, Dazl_background_data$mean.gene.TPM)), "DAZL.targets.v.all.nontargets_TPM_kstst")
writeLines(capture.output(ks.test(Dazl_data$utr3_length, Dazl_background_data$utr3_length)), "DAZL.targets.v.all.nontargets_utr3.length_kstst")

writeLines(capture.output(wilcox.test(Dazl_data$CDS_length, Dazl_background_data$CDS_length)), "DAZL.targets.v.all.nontargets_CDS_length_wilcoxtst")
writeLines(capture.output(wilcox.test(Dazl_data$CAI, Dazl_background_data$CAI)), "DAZL.targets.v.all.nontargets_CAI_wilcoxtst")
writeLines(capture.output(wilcox.test(Dazl_data$mean.gene.TPM, Dazl_background_data$mean.gene.TPM)), "DAZL.targets.v.all.nontargets_TPM_wilcoxtst")
writeLines(capture.output(wilcox.test(Dazl_data$utr3_length, Dazl_background_data$utr3_length)), "DAZL.targets.v.all.nontargets_utr3.length_wilcoxtst")

### TE vs TPM scatterplot, all targets
ggplot(total.data.filtered, aes(x=log2(mean.gene.TPM), y=log2(trans.eff.Ribotag.DESeq), color="mediumpurple")) +
  geom_point(size = 0.1, color = "mediumpurple", alpha =  0.15) + 
  xlab("Transcript abundance (log2 TPM)") + 
  ylab("Translational efficiency (log2 RiboTag enrichment)") +
  scale_y_continuous(breaks=c(-8, -6, -4, -2, 0, 2,4,6,8)) +
  scale_colour_manual(values=c("mediumpurple")) +
  coord_cartesian(ylim=c(-8, 8)) +
  theme_mk + 
  theme(#axis.title.x=element_blank(),
        #axis.ticks.x=element_blank(),
                 axis.text=element_text(size=6, colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))

ggsave("FIGURE_trans.eff.DESeq_all.transcripts_scatterplot_TPM_190313.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)






### TE vs CDS scatterplot
ggplot(total.data.filtered, aes(x=log2(CDS_length), y=log2(trans.eff.Ribotag.DESeq))) +
  geom_point(size = 0.1, color = "mediumpurple", alpha =  0.15) + 
  xlab("CDS length (log2 nt)") + 
  ylab("Translational efficiency (log2 RiboTag enrichment)") +
  scale_y_continuous(breaks=c(-8, -6, -4, -2, 0, 2,4,6,8)) +
  scale_colour_manual(values=c("mediumpurple")) +
  coord_cartesian(ylim=c(-8, 8)) +
  theme_mk + 
  theme(#axis.title.x=element_blank(),
        #axis.ticks.x=element_blank(),
         axis.text=element_text(size=6, colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))

ggsave("FIGURE_trans.eff.DESeq_all.transcripts_scatterplot_CDS_length_190313.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)





### TE vs CAI
ggplot(total.data.filtered, aes(x=log2(CAI), y=log2(trans.eff.Ribotag.DESeq))) +
  geom_point(size = 0.1, color = "mediumpurple", alpha =  0.15) + 
  xlab("Codon usage (log2 CAI)") + 
  ylab("Translational efficiency (log2 RiboTag enrichment)") +
  scale_y_continuous(breaks=c(-8, -6, -4, -2, 0, 2,4,6,8)) +
  scale_colour_manual(values=c("mediumpurple")) +
  coord_cartesian(ylim=c(-8, 8), xlim=c(-0.55, -0.15)) +
  theme_mk + 
  theme(#axis.title.x=element_blank(),
        #axis.ticks.x=element_blank(),
         axis.text=element_text(size=6, colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))

ggsave("FIGURE_trans.eff.DESeq_all.transcripts_scatterplot_CAI_190313.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)


### TE vs utr3 scatterplot
ggplot(total.data.filtered, aes(x=log2(utr3_length), y=log2(trans.eff.Ribotag.DESeq))) +
  geom_point(size = 0.1, color = "mediumpurple", alpha =  0.15) + 
  xlab("3' UTR length (log2 nt))") + 
  ylab("Translational efficiency (log2 RiboTag enrichment)") +
  scale_x_continuous(breaks =c(5, 10, 15)) +
  scale_y_continuous(breaks=c(-8, -6, -4, -2, 0, 2,4,6,8)) +
  scale_colour_manual(values=c("mediumpurple")) +
  coord_cartesian(ylim=c(-8, 8), xlim=c(2.5, 15)) +
  theme_mk + 
  theme(#axis.title.x=element_blank(),
        #axis.ticks.x=element_blank(),
         axis.text=element_text(size=6, colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))
ggsave("FIGURE_trans.eff.DESeq_all.transcripts_scatterplot_utr3_190313.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)








###################
### does Dazl binding enhance translational efficiency
###################
### Test for difference between DAZL targets and nontargets using DESeq-based TE



### as a first pass, are DAZL targets different from nontargets?
### are different variables that affect translation different between Dazl targets and nontargets?
# Function to produce summary statistics (median and interquartile range)
data_summary <- function(x) {
   m <- median(x)
   ymin <- unname(quantile(x, 0.25))
   ymax <- unname(quantile(x, 0.75))
   return(c(y=m,ymin=ymin,ymax=ymax))
}

### relabel target v. nontarget
nontargets = total.data.filtered[total.data.filtered$Dazl_target %in% "Nontarget",]
Dazl_targets = total.data.filtered[total.data.filtered$Dazl_target %in% "DAZL target",]

nontargets$Dazl_target = paste("Nontarget\n(",nrow(nontargets),")", sep="")
Dazl_targets$Dazl_target = paste("DAZL target\n(",nrow(Dazl_targets),")", sep="")



mm_data = rbind(nontargets[,c("trans.eff.Ribotag.DESeq", "Dazl_target")], 
            Dazl_targets[,c("trans.eff.Ribotag.DESeq", "Dazl_target")])

  ## order data
mm_data$Dazl_target = factor(mm_data$Dazl_target,c(paste("DAZL target\n(",nrow(Dazl_targets),")", sep=""), 
                                                paste("Nontarget\n(",nrow(nontargets),")", sep="")
                                              ))

### violin plot TE, dazl target v. nontarget
ggplot(mm_data, aes(Dazl_target, trans.eff.Ribotag.DESeq,fill=Dazl_target)) +
  geom_boxplot(outlier.size=0, outlier.stroke=0, lwd=0.2) +
  ylab("Translational efficiency") +
  scale_fill_manual(values=c("#80b1d3", "gray")) +
  coord_cartesian(ylim=c(0, 7)) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7)) + 
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))
ggsave("FIGURE_trans.eff.DESeq_target.v.nontarget_boxplot_total.data_191003.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)


  ## test TE difference between DAZL targets and nontargets
writeLines(capture.output(ks.test(nontargets$trans.eff.Ribotag.DESeq, Dazl_targets$trans.eff.Ribotag.DESeq)), "DAZL.targets.v.all.nontargets_undiff.gonia.TE_kstst")
writeLines(capture.output(wilcox.test(nontargets$trans.eff.Ribotag.DESeq, Dazl_targets$trans.eff.Ribotag.DESeq)), "DAZL.targets.v.all.nontargets_undiff.gonia.TE_wilcox.test")






### are different variables that affect translation different between Dazl targets and nontargets?
# Function to produce summary statistics (median and interquartile range)
data_summary <- function(x) {
   m <- median(x)
   ymin <- unname(quantile(x, 0.25))
   ymax <- unname(quantile(x, 0.75))
   return(c(y=m,ymin=ymin,ymax=ymax))
}

  ## order data
total.data.filtered$Dazl_target = factor(total.data.filtered$Dazl_target,c("Nontarget", "DAZL target"))

### violin plot TPM, dazl target v. nontarget
ggplot(total.data.filtered, aes(Dazl_target, log2(mean.gene.TPM),fill=Dazl_target, colour=Dazl_target)) +
  geom_violin() +
  stat_summary(fun.data=data_summary, colour="black", size=0.05) +
  ylab("Transcript abundance (log2 TPM)") +
  scale_fill_manual(values=c("gray", "#80b1d3")) +
   scale_colour_manual(values=c("gray", "#80b1d3")) +
  scale_y_continuous(breaks=c(0,5, 10, 15)) + 
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))
ggsave("FIGURE_TPM_target.v.nontarget_violinplot_190409.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)




### violin plot CAI, dazl target v. nontarget
ggplot(total.data.filtered, aes(Dazl_target, CAI,fill=Dazl_target, colour=Dazl_target)) +
  geom_violin() +
  stat_summary(fun.data=data_summary, colour="black", size=0.05) +
  ylab("Codon usage (CAI)") +
  scale_fill_manual(values=c("gray", "#80b1d3")) +
   scale_colour_manual(values=c("gray", "#80b1d3")) +
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))
ggsave("FIGURE_CAI_target.v.nontarget_violinplot_190409.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)

### violin plot CDS, dazl target v. nontarget
ggplot(total.data.filtered, aes(Dazl_target, log2(CDS_length),fill=Dazl_target, colour=Dazl_target)) +
  geom_violin() +
  stat_summary(fun.data=data_summary, colour="black", size=0.05) +
  ylab("CDS length (log2 nt)") +
  scale_fill_manual(values=c("gray", "#80b1d3")) +
   scale_colour_manual(values=c("gray", "#80b1d3")) + 
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))
ggsave("FIGURE_CDS.length_target.v.nontarget_violinplot_190409.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)


### violin plot utr3_length, dazl target v. nontarget
ggplot(total.data.filtered, aes(Dazl_target, log2(utr3_length),fill=Dazl_target, colour=Dazl_target)) +
  geom_violin() +
  stat_summary(fun.data=data_summary, colour="black", size=0.05) +
  ylab("3' UTR length (log2 nt)") +
  scale_fill_manual(values=c("gray", "#80b1d3")) +
   scale_colour_manual(values=c("gray", "#80b1d3")) +
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))
ggsave("FIGURE_utr3.length_target.v.nontarget_violinplot_190409.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)











### resubset data from the total.data.filtered (the smaller dataset that includes are the necessary datapoints)
Dazl_targets = read.table("utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_gene.ID.only")

Dazl_data = total.data.filtered[total.data.filtered$gene.id %in% Dazl_targets$V1 ,]
Dazl_background_data = total.data.filtered[!total.data.filtered$gene.id  %in% Dazl_targets$V1 ,]






#function to compare compare an dependent variable between two groups (foreground and background) after re-sampling to match an attribute (att, aka independent variable) of background to that of foreground
#attribute is the name of the column that needs to be controlled for
## match_type is the variable that is being matched, i.e., counts, CDS length, etc. Will be used to name output files
sampleMatchedPcts <- function(foreground_att,background_att,attribute_column_number, match_type) {

foreground_att$bin = cut(as.numeric(foreground_att[,attribute_column_number]), breaks=quantile(as.numeric(foreground_att[,attribute_column_number]),seq(0,1,0.05), na.rm=TRUE)) ### bin foreground
background_att$bin = cut(as.numeric(background_att[,attribute_column_number]), breaks=quantile(as.numeric(foreground_att[,attribute_column_number]),seq(0,1,0.05), na.rm=TRUE))  ### bin background based on foreground bins

matched_bckgrnd_genes = c()
matched_pctvec = c()
for (ibin in unique(foreground_att$bin)) { # for each foreground bin
  sampsize = nrow(subset(foreground_att, bin==ibin))
  sampsize_bckgnd = nrow(subset(background_att, bin==ibin))
  #if (nrow(subset(foreground_att, bin==ibin)) > nrow(subset(background_att, bin == ibin))){
  #  sampsize = floor((nrow(subset(foreground_att, bin==ibin))/nrow(foreground_att))*nrow(background_att))
  #}
  if (sampsize>=sampsize_bckgnd) {
    matched_bckgrnd_genes = c(matched_bckgrnd_genes,subset(background_att, bin == ibin )[,1])
  }
  if (sampsize<sampsize_bckgnd) {
    matched_bckgrnd_genes = c(matched_bckgrnd_genes,sample(x=subset(background_att, bin == ibin )[,1],size=sampsize,replace = FALSE))
  }
}

matched_bckgrnd_genes_pcts = background_att[background_att[,1] %in% matched_bckgrnd_genes ,]

### print stats comparing matched background set to actual set
writeLines(capture.output(ks.test(matched_bckgrnd_genes_pcts[,attribute_column_number], 
  foreground_att[,attribute_column_number])), 
  paste("DAZL.targets.v.", as.character(match_type), ".matched.bckgnd_",  as.character(match_type), "_kstest", sep=""))

writeLines(capture.output(wilcox.test(matched_bckgrnd_genes_pcts[,attribute_column_number], 
  foreground_att[,attribute_column_number])), 
  paste("DAZL.targets.v.", as.character(match_type), ".matched.bckgnd_",  as.character(match_type), "_wilcoxtest", sep=""))

return(matched_bckgrnd_genes_pcts)
}







#####################
### generate background dataset that matched the 3S RNA-seq TPM distribution of the DAZL dataset; graph via boxplot
#####################
TPM="TPM"
TPM_matched = sampleMatchedPcts(Dazl_data, Dazl_background_data,4, TPM)
write.table(TPM_matched, file="MM291_translational.efficiency.DESeq_3S.undiff.gonia_190313_TPM.matched.background.set", sep="\t", quote=F, col.names = NA)
  ##ks.test and wilcox.test results saved as separate file


  ## test TE difference between DAZL targets and TPM matched background
writeLines(capture.output(ks.test(Dazl_data$trans.eff.Ribotag.DESeq, TPM_matched$trans.eff.Ribotag.DESeq)), "DAZL.targets.v.TPM.matched_undiff.gonia.TE_kstst")
writeLines(capture.output(wilcox.test(Dazl_data$trans.eff.Ribotag.DESeq, TPM_matched$trans.eff.Ribotag.DESeq,  alternative = c("greater"))), "DAZL.targets.v.TPM.matched_undiff.gonia.TE_wilcoxtst")




data = rbind(Dazl_data[,c("trans.eff.Ribotag.DESeq", "Dazl_target")], TPM_matched[,c("trans.eff.Ribotag.DESeq", "Dazl_target")])

  ### translational efficiency graph
ggplot(data, aes(Dazl_target, log2(trans.eff.Ribotag.DESeq),fill=Dazl_target)) +
  geom_hline(yintercept=0, colour="gray", lwd=0.2) +
  geom_boxplot(outlier.size=0, outlier.stroke=0, lwd=0.2) +
  ylab("log2(translational efficiency)") +
  scale_fill_manual(values=c("Other expressed genes"="lightslategray","STRA8 targets"="lightskyblue1")) +
  coord_cartesian(ylim=c(-10, 10)) +
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))

ggsave("FIGURE_trans.eff.DESeq_boxplot_DAZL.targets.v.TPM.matched.bckgnd_190313.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)



#####################
### generate background dataset that maches the CDS length distribution of the DAZL dataset; graph via bloxplot
#####################
CDS="CDS_length"
CDS_matched = sampleMatchedPcts(Dazl_data, Dazl_background_data,5,CDS)
write.table(CDS_matched, file="MM291_translational.efficiency.DESeq_3S.undiff.gonia_130313_CDS.matched.background.set", sep="\t", quote=F, col.names = NA)
  ##ks.test and wilcox.test results saved as separate file

  ## test TE difference between DAZL targets and CDS matched background
writeLines(capture.output(ks.test(Dazl_data$trans.eff.Ribotag.DESeq, CDS_matched$trans.eff.Ribotag.DESeq)), "DAZL.targets.v.CDS.matched_undiff.gonia.TE_kstst")
writeLines(capture.output(wilcox.test(Dazl_data$trans.eff.Ribotag.DESeq, CDS_matched$trans.eff.Ribotag.DESeq,  alternative = c("greater"))), "DAZL.targets.v.CDS.matched_undiff.gonia.TE_wilcoxtst")




data = rbind(Dazl_data[,c("trans.eff.Ribotag.DESeq","Dazl_target")], CDS_matched[,c("trans.eff.Ribotag.DESeq", "Dazl_target")])

### translational efficiency graph
ggplot(data, aes(Dazl_target, log2(trans.eff.Ribotag.DESeq),fill=Dazl_target)) +
  geom_hline(yintercept=0, colour="gray", lwd=0.2) +
  geom_boxplot(outlier.size=0, outlier.stroke=0, lwd=0.2) +
  ylab("log2(translational efficiency)") +
  scale_fill_manual(values=c("Other expressed genes"="lightslategray","STRA8 targets"="lightskyblue1")) +
  coord_cartesian(ylim=c(-10, 10)) +
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))

ggsave("FIGURE_trans.eff.DESeq_boxplot_DAZL.targets.v.CDS.matched.bckgnd_190313.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)






#########
## Create a single graph with DAZL targets and TPM_matched and  CDS_matched backgrounds
#########


CDS_matched$Dazl_target = paste("CDS length (",nrow(CDS_matched),")", sep="")
TPM_matched$Dazl_target = paste("Transcript abundance (",nrow(TPM_matched),")", sep="")

Dazl_data$Dazl_target = paste("DAZL targets (",nrow(Dazl_data),")", sep="")

data = rbind(Dazl_data[,c("trans.eff.Ribotag.DESeq", "Dazl_target")], 
            CDS_matched[,c("trans.eff.Ribotag.DESeq", "Dazl_target")], 
            TPM_matched[,c("trans.eff.Ribotag.DESeq","Dazl_target")])

  ## order data
data$Dazl_target = factor(data$Dazl_target,c(paste("DAZL targets (",nrow(Dazl_data),")", sep=""), 
                                              paste("Transcript abundance (",nrow(TPM_matched),")", sep=""),
                                              paste("CDS length (",nrow(CDS_matched),")", sep="")
                                              ))


ggplot(data, aes(Dazl_target, trans.eff.Ribotag.DESeq,fill=Dazl_target)) +
  geom_boxplot(outlier.size=0, outlier.stroke=0, lwd=0.2) +
  ylab("Translational efficiency") +
  scale_fill_manual(values=c("#80b1d3", "gray","gray")) +
  coord_cartesian(ylim=c(0, 7)) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7)) + 
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))

ggsave("FIGURE_trans.eff.DESeq_boxplot_DAZL.targets.v.TPM.matched.v.CDS.matched_190313.pdf", scale = 1, width = 55, height=45, units = c("mm"), dpi = 300)



















#### add no. binding sites to dataset

no.sites = read.table("utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_peaks.per.gene", header=T)


Dazl.peaks.data.temp= total.data.filtered[total.data.filtered$gene.id %in% no.sites$gene.id  ,]
Dazl.peaks.data.temp2 = merge(Dazl.peaks.data.temp, no.sites, by="gene.id", all=F)

other.peaks.data= total.data.filtered[!total.data.filtered$gene.id %in% no.sites$gene.id  ,]
other.peaks.data$no.peaks=0

Dazl.peaks.data=rbind(Dazl.peaks.data.temp2, other.peaks.data)
## remove NAs
Dazl.peaks.data = na.omit(Dazl.peaks.data)


## adjust column name
colnames(Dazl.peaks.data)[9] = "no.Dazl.peaks"

write.table(Dazl.peaks.data, "TE.DESeq_TPM_CDS_length_utr3.length_CAI_Dazl.peaks.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names = NA)

#Dazl.peaks.data = read.table("TE.DESeq_TPM_CDS_length_utr3.length_CAI_Dazl.peaks.txt", sep="\t",header=T, row.names=1)











####################
### Does DAZL enhance translational efficiency?
###################
######
## create linear model with Dazl targets and nontargets
#####



# Check linear model using subsets of variables first
model1 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ Dazl.peaks.data$Dazl_target)
model2 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$mean.gene.TPM))
model3 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$CDS_length))
model4 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$CAI))
model5 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$utr3_length))
model6 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$mean.gene.TPM)+log2(Dazl.peaks.data$CDS_length))
model7 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$CDS_length)+log2(Dazl.peaks.data$CAI))
model8 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$mean.gene.TPM)+log2(Dazl.peaks.data$CAI))
model9 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$mean.gene.TPM)+log2(Dazl.peaks.data$CDS_length)+log2(Dazl.peaks.data$utr3_length))
model10 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$mean.gene.TPM)+log2(Dazl.peaks.data$CDS_length)+log2(Dazl.peaks.data$utr3_length)+log2(Dazl.peaks.data$CAI))
model11 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$mean.gene.TPM)+log2(Dazl.peaks.data$CDS_length)+log2(Dazl.peaks.data$CAI))

summaries = append(capture.output(summary(model1)),capture.output(summary(model2)))
summaries = append(summaries,capture.output(summary(model3)))
summaries = append(summaries,capture.output(summary(model4)))
summaries = append(summaries,capture.output(summary(model5)))
summaries = append(summaries,capture.output(summary(model6)))
summaries = append(summaries,capture.output(summary(model7)))
summaries = append(summaries,capture.output(summary(model8)))
summaries = append(summaries,capture.output(summary(model9)))
summaries = append(summaries,capture.output(summary(model10)))
summaries = append(summaries,capture.output(summary(model11)))


writeLines(summaries, "Linear.model_alternates_DAZL.binding_TPM_CDS.length_CAI_utr3.length")


## Full linear model with TPM, CDS, CAI, and 3' UTR
# verify that assumptions are met: http://r-statistics.co/Assumptions-of-Linear-Regression.html
model = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ Dazl.peaks.data$Dazl_target+log2(Dazl.peaks.data$mean.gene.TPM)+log2(Dazl.peaks.data$CDS_length)+log2(Dazl.peaks.data$utr3_length)+log2(Dazl.peaks.data$CAI))


out1 =  capture.output(summary(model))  ## Summary of model
out2 =  capture.output(confint(model))  ## Confidence interval for each variable's coefficient

## Assumption 1: regression model is linear in parameters
## Assumption 2: Mean of residuals is 0
out3 = paste("mean of residuals", capture.output(mean(model$residuals)), sep="\t") ## Mean of residuals

## Assumption 3: homoscedasticity of residuals (plot of fitted values vs. regression)
## Assumption 10: normality of residuals (qqnorm() plot)
pdf("Linear.model_Dazl.peaks_TPM_CDS.length_utr3.length_CAI_graphs.pdf")
plot(model)
dev.off()

## Assumption 4: no autocorrelation of residuals
pdf("Linear.model_DAZL.peaks_TPM_CDS.length_utr3.length_CAI_autocorrleation.graph.pdf")
acf(model$residuals) 
dev.off()

## Assumption 5: X variables and residuals are uncorrelated 
#out4 = capture.output(cor.test(log2(Dazl.peaks.data$mean.gene.TPM) , model$residuals))
#out5 = capture.output(cor.test(log2(Dazl.peaks.data$CDS_length) , model$residuals))
#out6 = capture.output(cor.test(log2(Dazl.peaks.data$CAI) , model$residuals))
#out7 = capture.output(cor.test(Dazl.peaks.data$Dazl_target , model$residuals))


## Assumption 6: no. observations must be greater than no X variables
## Assumption 7: variability in X variable values is positive
out8 = paste("TPM variance", capture.output(var(log2(Dazl.peaks.data$mean.gene.TPM))), sep="\t")
out9 = paste("CDS length variance", capture.output(log2(var(Dazl.peaks.data$CDS_length))), sep="\t")
out10 = paste("CAI variance", capture.output(var(log2(Dazl.peaks.data$CAI))), sep="\t")
out11 = paste("3' UTR variance", capture.output(var(log2(Dazl.peaks.data$utr3_length))), sep="\t")
out12 = paste("No. DAZL peaks variance", capture.output(var(as.numeric(Dazl.peaks.data$Dazl_target))), sep="\t")


## Assumption 8: regression model is correctly specified

## Assumption 9: no perfect multicollinearity
## VIF should be less then 4; less than 2 is an even more stringent cutoff
out13 = paste("Variance Inflation Factor (VIF)", capture.output(vif(model)), sep="\t") 



## Write model and assumption tests to file
output = append(out1, out2)
output = append(output, out3)
#output = append(output, out4)
#output = append(output, out5)
#output = append(output, out6)
#output = append(output, out7)
output = append(output, out8)
output = append(output, out9)
output = append(output, out10)
output = append(output, out11)
output = append(output, out12)
output = append(output, out13)

writeLines(output, "Linear.model_DAZL.peaks_TPM_CDS.length_utr3.length_CAI")







## Full linear model with TPM, CDS, CAI
# verify that assumptions are met: http://r-statistics.co/Assumptions-of-Linear-Regression.html
model2 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ Dazl.peaks.data$Dazl_target+log2(Dazl.peaks.data$mean.gene.TPM)+log2(Dazl.peaks.data$CDS_length)+log2(Dazl.peaks.data$CAI))


out1 =  capture.output(summary(model2))  ## Summary of model2
out2 =  capture.output(confint(model2))  ## Confidence interval for each variable's coefficient

## Assumption 1: regression model2 is linear in parameters
## Assumption 2: Mean of residuals is 0
out3 = paste("mean of residuals", capture.output(mean(model2$residuals)), sep="\t") ## Mean of residuals

## Assumption 3: homoscedasticity of residuals (plot of fitted values vs. regression)
## Assumption 10: normality of residuals (qqnorm() plot)
pdf("Linear.model_Dazl.peaks_TPM_CDS.length_CAI_graphs.pdf")
plot(model2)
dev.off()

## Assumption 4: no autocorrelation of residuals
pdf("Linear.model_DAZL.peaks_TPM_CDS.length_CAI_autocorrleation.graph.pdf")
acf(model2$residuals) 
dev.off()

## Assumption 5: X variables and residuals are uncorrelated 
#out4 = capture.output(cor.test(log2(Dazl.peaks.data$mean.gene.TPM) , model2$residuals))
#out5 = capture.output(cor.test(log2(Dazl.peaks.data$CDS_length) , model2$residuals))
#out6 = capture.output(cor.test(log2(Dazl.peaks.data$CAI) , model2$residuals))
#out7 = capture.output(cor.test(Dazl.peaks.data$Dazl_target , model2$residuals))


## Assumption 6: no. observations must be greater than no X variables
## Assumption 7: variability in X variable values is positive
out8 = paste("TPM variance", capture.output(var(log2(Dazl.peaks.data$mean.gene.TPM))), sep="\t")
out9 = paste("CDS length variance", capture.output(log2(var(Dazl.peaks.data$CDS_length))), sep="\t")
out10 = paste("CAI variance", capture.output(var(log2(Dazl.peaks.data$CAI))), sep="\t")
out11 = paste("3' UTR variance", capture.output(var(log2(Dazl.peaks.data$utr3_length))), sep="\t")



## Assumption 8: regression model2 is correctly specified

## Assumption 9: no perfect multicollinearity
## VIF should be less then 4; less than 2 is an even more stringent cutoff
out13 = paste("Variance Inflation Factor (VIF)", capture.output(vif(model2)), sep="\t") 



## Write model2 and assumption tests to file
output = append(out1, out2)
output = append(output, out3)
#output = append(output, out4)
#output = append(output, out5)
#output = append(output, out6)
#output = append(output, out7)
output = append(output, out8)
output = append(output, out9)
output = append(output, out10)
output = append(output, out11)
output = append(output, out12)
output = append(output, out13)

writeLines(output, "Linear.model_DAZL.peaks_TPM_CDS.length_CAI")




## Full linear model with TPM, CDS
# verify that assumptions are met: http://r-statistics.co/Assumptions-of-Linear-Regression.html
model3 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ Dazl.peaks.data$Dazl_target+log2(Dazl.peaks.data$mean.gene.TPM)+log2(Dazl.peaks.data$CDS_length))


out1 =  capture.output(summary(model3))  ## Summary of model3
out2 =  capture.output(confint(model3))  ## Confidence interval for each variable's coefficient

## Assumption 1: regression model3 is linear in parameters
## Assumption 2: Mean of residuals is 0
out3 = paste("mean of residuals", capture.output(mean(model3$residuals)), sep="\t") ## Mean of residuals

## Assumption 3: homoscedasticity of residuals (plot of fitted values vs. regression)
## Assumption 10: normality of residuals (qqnorm() plot)
pdf("Linear.model_Dazl.peaks_TPM_CDS.length_graphs.pdf")
plot(model3)
dev.off()

## Assumption 4: no autocorrelation of residuals
pdf("Linear.model_DAZL.peaks_TPM_CDS.length_autocorrleation.graph.pdf")
acf(model3$residuals) 
dev.off()

## Assumption 5: X variables and residuals are uncorrelated 
#out4 = capture.output(cor.test(log2(Dazl.peaks.data$mean.gene.TPM) , model3$residuals))
#out5 = capture.output(cor.test(log2(Dazl.peaks.data$CDS_length) , model3$residuals))
#out6 = capture.output(cor.test(log2(Dazl.peaks.data$CAI) , model3$residuals))
#out7 = capture.output(cor.test(Dazl.peaks.data$Dazl_target , model3$residuals))


## Assumption 6: no. observations must be greater than no X variables
## Assumption 7: variability in X variable values is positive
out8 = paste("TPM variance", capture.output(var(log2(Dazl.peaks.data$mean.gene.TPM))), sep="\t")
out9 = paste("CDS length variance", capture.output(log2(var(Dazl.peaks.data$CDS_length))), sep="\t")
out10 = paste("CAI variance", capture.output(var(log2(Dazl.peaks.data$CAI))), sep="\t")


## Assumption 8: regression model3 is correctly specified

## Assumption 9: no perfect multicollinearity
## VIF should be less then 4; less than 2 is an even more stringent cutoff
out13 = paste("Variance Inflation Factor (VIF)", capture.output(vif(model3)), sep="\t") 



## Write model3 and assumption tests to file
output = append(out1, out2)
output = append(output, out3)
#output = append(output, out4)
#output = append(output, out5)
#output = append(output, out6)
#output = append(output, out7)
output = append(output, out8)
output = append(output, out9)
output = append(output, out10)
output = append(output, out11)
output = append(output, out12)
output = append(output, out13)

writeLines(output, "Linear.model_DAZL.peaks_TPM_CDS.length")





##########
### Does adding DAZL binding to the model statistically improve the model's fit to the data?
### Does adding other variables to the model statistically improve the model's fit to the data?
#########


### set up cumulative linear models

alt.model1 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$mean.gene.TPM))
alt.model2 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$mean.gene.TPM)+log2(Dazl.peaks.data$CDS_length))
alt.model3 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$mean.gene.TPM)+log2(Dazl.peaks.data$CDS_length)+log2(Dazl.peaks.data$utr3_length))
alt.model4 = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$mean.gene.TPM)+log2(Dazl.peaks.data$CDS_length)+log2(Dazl.peaks.data$utr3_length)+log2(Dazl.peaks.data$CAI))
model.full = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$mean.gene.TPM)+log2(Dazl.peaks.data$CDS_length)+log2(Dazl.peaks.data$utr3_length)+log2(Dazl.peaks.data$CAI)+Dazl.peaks.data$Dazl_target)


# calculate the log likelihood for each fit...
# this is the log of the exact probability of seeing your data
# given your model (this quantity is central to AIC); i.e., 
# it's a measure of how well the model fits the data; bigger
# numbers (i.e., less negative numbers) mean a better fit; 
# L1 will be bigger than L0
L.alt1 <- as.numeric(logLik(alt.model1))
L.alt2 <- as.numeric(logLik(alt.model2))
L.alt3 <- as.numeric(logLik(alt.model3))
L.alt4 <- as.numeric(logLik(alt.model4))
L.full <- as.numeric(logLik(model.full))

# now we are looking to see if the log likelihoods are more
# different than we would expect by chance for two models that
# are the same except one has an extra variable;
# it turns out that -2 x the difference in log likelihoods
# is a value that follows a chi-squared distribution with one
# degree of freedom; we use one degree of freedom because 
# the second model has one extra variable; we can use this 
# to calculate a p-value...
pval.full <- pchisq(-2 * (L.alt4 -L.full), df=1, lower.tail = FALSE)
pval.alt4 <- pchisq(-2 * (L.alt3 - L.alt4), df=1, lower.tail = FALSE)
pval.alt3 <- pchisq(-2 * (L.alt2 - L.alt3), df=1, lower.tail = FALSE)
pval.alt2 <- pchisq(-2 * (L.alt1 - L.alt2), df=1, lower.tail = FALSE)


### output this analysis to file
log.likelihood.output = "Does adding DAZL binding to the linear model statistically improve the fit to the data"
log.likelihood.output = append(log.likelihood.output,capture.output(alt.model4))
log.likelihood.output = append(log.likelihood.output,capture.output(model.full))
log.likelihood.output = append(log.likelihood.output, "Likelihood ratio test pval:")
log.likelihood.output = append(log.likelihood.output, capture.output(pval.full))
log.likelihood.output = append(log.likelihood.output, "\n\n")

log.likelihood.output = append(log.likelihood.output,"Does adding codon usage (CAI) to the linear model statistically improve the fit to the data")
log.likelihood.output = append(log.likelihood.output,capture.output(alt.model4))
log.likelihood.output = append(log.likelihood.output,capture.output(alt.model3))
log.likelihood.output = append(log.likelihood.output, "Likelihood ratio test pval:")
log.likelihood.output = append(log.likelihood.output, capture.output(pval.alt4))
log.likelihood.output = append(log.likelihood.output, "\n\n")

log.likelihood.output = append(log.likelihood.output,"Does adding 3' UTR length to the linear model statistically improve the fit to the data")
log.likelihood.output = append(log.likelihood.output,capture.output(alt.model3))
log.likelihood.output = append(log.likelihood.output,capture.output(alt.model2))
log.likelihood.output = append(log.likelihood.output, "Likelihood ratio test pval:")
log.likelihood.output = append(log.likelihood.output, capture.output(pval.alt3))
log.likelihood.output = append(log.likelihood.output, "\n\n")

log.likelihood.output = append(log.likelihood.output,"Does adding CDS length to the linear model statistically improve the fit to the data")
log.likelihood.output = append(log.likelihood.output,capture.output(alt.model2))
log.likelihood.output = append(log.likelihood.output,capture.output(alt.model1))
log.likelihood.output = append(log.likelihood.output, "Likelihood ratio test pval:")
log.likelihood.output = append(log.likelihood.output, capture.output(pval.alt2))
log.likelihood.output = append(log.likelihood.output, "\n\n")

writeLines(log.likelihood.output, "Linear.model_DAZL.peaks_TPM_CDS.length_utr3.length_CAI_likelihood.ratio.test.txt")






##############
### use model to compare adjusted TE in DAZL targets vs nontargets
##############

model.full = lm(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) ~ log2(Dazl.peaks.data$mean.gene.TPM)+log2(Dazl.peaks.data$CDS_length)+log2(Dazl.peaks.data$utr3_length)+log2(Dazl.peaks.data$CAI)+Dazl.peaks.data$Dazl_target)


intercept = model.full$coefficients[1]
coeff.TPM = model.full$coefficients[2]
coeff.CDS_length = model.full$coefficients[3]
coeff.utr3_length = model.full$coefficients[4]
coeff.CAI = model.full$coefficients[5]



### calculate adjusted trans eff that controls for TPM, CDS length, 3' UTR length, and CAI
Dazl.peaks.data$adj.trans.eff= 2^(log2(Dazl.peaks.data$trans.eff.Ribotag.DESeq) - 
              (coeff.TPM*log2(Dazl.peaks.data$mean.gene.TPM)+coeff.CDS_length*log2(Dazl.peaks.data$CDS_length)+coeff.utr3_length*log2(Dazl.peaks.data$utr3_length)+coeff.CAI*log2(Dazl.peaks.data$CAI)))

## write data

write.table(Dazl.peaks.data, "MM291.trans.eff_MM257.3S.undiff.gonia_DESeq_Ribotag.specific_count.TPM.filtered_adjTE", sep="\t", quote=FALSE, row.names=F)

### test difference in adj. TE between DAZL targets and nontargets
out.Dazl = capture.output(summary(Dazl.peaks.data[Dazl.peaks.data$Dazl_target %in% "DAZL target",]))
out.nontargets = capture.output(summary(Dazl.peaks.data[Dazl.peaks.data$Dazl_target %in% "Nontarget",]))
out.ks = capture.output(ks.test(Dazl.peaks.data[Dazl.peaks.data$Dazl_target %in% "DAZL target",]$adj.trans.eff, Dazl.peaks.data[Dazl.peaks.data$Dazl_target %in% "Nontarget",]$adj.trans.eff))
out.wil = capture.output(wilcox.test(Dazl.peaks.data[Dazl.peaks.data$Dazl_target %in% "DAZL target",]$adj.trans.eff, Dazl.peaks.data[Dazl.peaks.data$Dazl_target %in% "Nontarget",]$adj.trans.eff, alternative = c("greater")))
output = paste("Summary of Dazl targets")
output = append(output,out.Dazl)
output = append(output, paste("Summary of nontargets"))
output = append(output, out.nontargets)
output = append(output, out.ks)
output = append(output, out.wil)

writeLines(capture.output(output), "DAZL.targets.v.all.nontargets_undiff.gonia.adjTE_kstest_wilcoxtest.txt")


### violin plot adj TE, dazl target v. nontarget
ggplot(Dazl.peaks.data, aes(Dazl_target, adj.trans.eff,fill=Dazl_target)) +
  geom_boxplot(outlier.size=0, outlier.stroke=0, lwd=0.2) +
  ylab("Adjusted translational efficiency") +
  #scale_x_discrete(limits=c("STRA8 targets","Other expressed genes"),labels=c("Expressed\nSTRA8\ntargets","Other\nexpressed\ngenes")) +
  scale_fill_manual(values=c("#80b1d3", "gray")) +
                    #labels=c("Expressed STRA8\ntargets","Other expressed\ngenes"),
                    #breaks=c("STRA8 targets","Other expressed genes")) +
  coord_cartesian(ylim=c(0, 3)) +
  #scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7)) + 
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))
ggsave("FIGURE_adj.trans.eff_target.v.nontarget_boxplot_total.data_191003.pdf", scale = 1, width = 30, height=45, units = c("mm"), dpi = 300)






######
## Graph of cumulative linear R2
######

linear.R2 <- read.table("Cumulative_linear_R2.txt",stringsAsFactors = FALSE, header=T)
  ## order data
linear.R2$Variable = factor(linear.R2$Variable,c( "DAZL_binding", "CAI"   , "utr3_length" ,  "CDS_length", "TPM"))


ggplot(linear.R2, aes(x=Variable, y=Cum_linear_R2, fill=Variable) )+
  geom_bar(stat="identity", width=0.8) +
  coord_flip() +  ### horizontal graph
  #labs(y= expression(paste("-log"[10], "(", italic("P"), "-value)",sep="")), x="") +
  theme_mk +
  scale_fill_manual(values=c("gray", "gray", "gray", "gray", "gray")) + 
  theme(#axis.title.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0))) 

ggsave("FIGURE_Cumulative_linear_R2_180901.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)


writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

