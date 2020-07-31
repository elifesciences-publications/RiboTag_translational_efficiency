#### extract_CDS_calculate_CAI.R
## 
## Maria Mikedis mikedis@wi.mit.edu	
## Run as: Rscript extract_CDS_calculate_CAI.R
## Extract the CDS sequences from mm10 using UCSC table "refGene"
## Obtain the CAI (codon adaptation index for a list of transcripts)

## Output 5 files
## File 1: mm10_GRCm38_cds.fa of sequences
## File 2: mm10_GRCm38_cds_length.txt containing lengths of CDS
	## format: Refseq_id 	length_CDS
## File 3: CDS_seqs_expressed.fa containing sequences of transcripts expressed in undiff gonia (for most expressed isoform per gene; see synch.TPM variable)
## File 4: CDS_seqs_expressed_CAI containing CAI for transcripts expressed in undiff gonia (based on sequence for most expressed isoform per gene)
## File 5: CDS_seqs_expressed_no.CDS.no.CAI is a list of transcript isoform for which there was no sequence in the table downloaded here and therefore are not included in the CAI table
## Dependencies: GenomicFeatures, BSgenome.Mmusculus.UCSC.mm10


args <- commandArgs()
print(args)

#input <- args[6]  ### filename is passed as argument #6
#output <- args[7]

library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
library(coRdon)

## upload TPM data; data already filtered for min TPM of 1
## file contains coding transcripts only (pseudoaligned to everything, then pulled out coding transcripts on and re-normalized TPM)
synch.TPM = read.table("/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/20180730_3S_RNAseq_kallisto/190622_kallisto/kallisto/mRNA_coding_only/MM257_3S.undiff.gonia_kallisto_mm10_GRCm38_coding.mRNA.only_TPMs_180821_genes.transcripts_min1_most.expressed.isoform.per.gene", header=TRUE, sep="\t", stringsAsFactors = FALSE)



#########
## Extract CDS sequences for desired transcripts
#########

## Create gene model of mm10 genome from UCSC table ncbiRefSeqCurated
txdb <- makeTxDbFromUCSC(genome="mm10", tablename="refGene")

## Load mm10 genome as DNAstring object
genome <- BSgenome.Mmusculus.UCSC.mm10

## Extract the transcript coordinates from this gene model
transcripts(txdb)

## Create a map between gene and transcript identifiers by outing both
txbygene <- transcriptsBy(txdb, "gene")

## Extract the CDS coordinates; use.names option keeps the RefSeq ID with the coordinate
CDS <- cdsBy(txdb, use.names=TRUE)

## Extract the CDS sequences from the BSgenome
CDS_seqs <- extractTranscriptSeqs(genome, CDS)


filename="mm10_GRCm38_cds"
## Export sequences as fasta file
writeXStringSet(CDS_seqs, file=paste(filename, ".fa",sep = ""))

## Create dataframe containing transcript name and width of 3' UTR
df=data.frame(names(CDS_seqs), width(CDS_seqs)) 
colnames(df)=c("id", "CDS_length")

## For some reason, there are duplicate lines; remove duplicate lines.

df = unique(df)


## write this dataframe to a file
write.table(df, file=paste(filename, "_length.txt",sep = ""), quote=F, sep="\t", row.names=F, col.names=T)






### Upload file with translational efficiency
total.data.filtered = read.table("../MM291.trans.eff_MM257.3S.undiff.gonia_DESeq_Ribotag.specific_count.TPM.filtered_190618", header=TRUE, sep="\t", stringsAsFactors = FALSE)

### Take the top 5% of genes in translational efficiency
total.data.filtered.top5perc.TE= total.data.filtered[total.data.filtered$log2FoldChange > quantile(total.data.filtered$log2FoldChange,prob=1-5/100),]






## Extract sequences for transcripts in subset list
i=1
CDS_seqs_subset = DNAStringSet()
while (i<=length(total.data.filtered.top5perc.TE$transcript.id)) {
	CDS_seqs_subset = append(CDS_seqs_subset, CDS_seqs[as.character(total.data.filtered.top5perc.TE$transcript.id[i])])
	i=i+1
}




## Calculate codon usage
codon.usage <- codonTable(CDS_seqs_subset)
cc <- t(codonCounts(codon.usage)) ## also transpose so codons are in row names

## get sums
cc.sums = data.frame(row.names = row.names(cc),total.no=rowSums(cc))

## adjust ATG to not count for start codon
cc.sums["ATG" ,] = cc.sums["ATG" ,] - length(codon.usage)

## Load codon table
table = read.table("codon.table", row.names = 1, sep="\t", stringsAsFactors = FALSE)
colnames(table) = "AA"

## merge codon number with amino acid info 
data = cbind(cc.sums, table[, "AA"][match(rownames(cc.sums), rownames(table))])
colnames(data) = c("no", "amino.acid")
data$amino.acid = as.character(data$amino.acid)


########
## Calculate relative adaptives for each codon
########

## calculate relative adaptiveness per codon
	## Calculate fraction of codon relative to synonymous codons
	## For each line, relative adaptiveness equals codon fraction over
	## the max codon fraction for the subset with the same amino acid
i=1
while (i<= nrow(data)) {
	sum.per.AA = sum(subset(data, amino.acid %in% data$amino.acid[i])$no)
	data$fraction[i] = data$no[i]/sum.per.AA
	i=i+1
}
i=1
while (i<= nrow(data)) {
	data$rel.adapt[i] = data$fraction[i]/(max(subset(data, amino.acid %in% data$amino.acid[i])$fraction))
	i=i+1
}



## Create relative adaptiveness w table

w = as.matrix(data$rel.adapt)
row.names(w) = row.names(data)
w = as.table(w)



## get sequences of expressed genes; get sequence for those genes
expressed = synch.TPM$transcript.id

i=1
CDS_seqs_subset = DNAStringSet()
no_sequence = c() ## for transcript isoform that have no sequence in this downloaded table and will therefore be skipped
while (i<=length(expressed)) {
	if (expressed[i] %in% names(CDS_seqs)) { ### check if transcript is in the dataset of CDS sequences; if not, can skip this transcript
		CDS_seqs_subset = append(CDS_seqs_subset, CDS_seqs[expressed[i]])
	}
	else {
		no_sequence = append(no_sequence, expressed[i])
	}
	i=i+1
}


## Export sequences as fasta file, reimport using read.fasta tool from seqinr
writeXStringSet(CDS_seqs_subset, file=paste("CDS_seqs_expressed.fa",sep = ""))
input2 <- read.fasta(file = "CDS_seqs_expressed.fa")


########
## Calculate CAI
########
cai.res <- sapply(input2, cai, w = w)

cai.res.matrix = as.matrix(cai.res)

colnames(cai.res.matrix) = c("CAI")


## export to file
write.table(cai.res.matrix, file="CDS_seqs_expressed_CAI", quote=F, sep="\t", row.names=T, col.names=T)
write.table(no_sequence, file="CDS_seqs_expressed_no.CDS.no.CAI", quote=F, sep="\t")


