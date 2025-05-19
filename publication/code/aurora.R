rm(list = ls())
graphics.off()
set.seed(010124)
#
setwd("~/KarakachLab Dropbox/Research/BreastCancer/Methylation")
dataDir <- paste0(getwd(),"/AURORA/DNAm/")
resultsDir <- paste0(getwd(),"/AURORA/results/")
load.libs <- c("ArrayExpress","arrayQualityMetrics","Biobase","biomaRt","BSgenome.Hsapiens.UCSC.hg19","caret", "ChAMP","chipseq","circlize","clusterProfiler",
               "ComplexHeatmap", "cowplot","data.table","DESeq2","DMRcate","DOSE","doSNOW","dplyr","edgeR","ELMER","factoextra","FactoMineR","foreign","forestplot","forestploter",
               "GEOquery","genefilter","geneplotter", "GenomicRanges","gplots","ggplot2","ggpubr","ggrepel","GO.db","Glimma","grid","gridExtra","GSEABase","glmnet",
               "gProfileR","Gviz","oligo", "here","Homo.sapiens","hugene10sttranscriptcluster.db","IlluminaHumanMethylation450kanno.ilmn12.hg19","IlluminaHumanMethylation450kmanifest",
               "IMA","IRanges","kableExtra","kniter","limma","matrixStats","MEDIPS","missMethyl",
               "minfi","minfiData","oligoClasses","org.Hs.eg.db",  "org.Mm.eg.db","org.Rn.eg.db","openxlsx","parallel","pd.hugene.1.0.st.v1","pheatmap","PMA","RCy3","RColorBrewer",
               "sesame", "sesameData","ReactomePA","readxl","rmeta","RTCGA","RTCGA.clinical","RTCGA.mutations","RTCGA.rnaseq","rtracklayer","rWikiPathways","scales", "stringi",
               "stringr","SummarizedExperiment","survival","survminer", "TCGAbiolinks", "tidyr","tidyverse","topGO")

#
if (!require("pacman")) install.packages("pacman"); library(pacman)
#
pacman::p_load(load.libs, update = FALSE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
  print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}
cat("\014")

#
# create directories that may not exist
if (!file.exists(dataDir)) dir.create(dataDir)
if (!file.exists(resultsDir)) dir.create(resultsDir)

#
metadataFile<-"metaData_aurora.csv" # File containing metadata
exp_des<-as.data.frame(read.csv(metadataFile,stringsAsFactors = F,header = T))
# Select only the columns of interest from the metadata

idx <- which(colnames(exp_des) %in%
               c("BCR.Portion.barcode", "Sample.Type","Tissue.primary.type.treatment","Race","Anatomic.Site.Simplified",
                 "Days.primary.diagnosis","Days.Metastasisastais.diagnosis","PAM50.Call","TNBC_Subtype","TMB"))
#
metaData <- data.frame((exp_des)[,idx])
rownames(metaData) <- exp_des$BCR.Portion.barcode
# replace "." with "-" in rownames
rownames(metaData) <- gsub("\\.","-",rownames(metaData))

# Load the data
load("combined_aurora.RData")
#
hub = AnnotationHub()
query(hub, c("EnsDb", "Homo sapiens", "97"))
edb = hub[["AH73881"]]
#keytypes(edb)
#columns(edb)
keys = keys(edb, "GENENAME")
columns =  c("GENEID", "ENTREZID", "GENEBIOTYPE","DESCRIPTION")

tbl <-ensembldb::select(edb, keys, columns, keytype = "GENENAME") %>%
  as_tibble()
#
#filter = ~ gene_name %in% keys & gene_biotype == "protein_coding"

#tbl = ensembldb::select(edb, filter, columns) %>%
# as_tibble()
coding_Annotation<-as.data.frame(tbl)

#
# get the 450k annotation data
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
#
# subset the data to keep rows where GeneIDs and coding_genes match
#

# get the least variable methylation probes
varProbes <- rownames(X.meth)[order(rowVars(X.meth) , decreasing=FALSE)]
varProbes <- varProbes[1:round(length(varProbes)*0.2)]
#
varGenes <-rownames(X.rna)[order(rowVars(X.rna) , decreasing=FALSE)]
varGenes <- varGenes[1:round(length(varGenes)*0.2)]
#
HVMeth <- X.meth[varProbes,]
HVRNA <- X.rna[varGenes,]



#
# combine the two datasets
#
combinedData <- rbind(HVMeth, HVRNA)

# now combine the metadata corresponding to the two datasets
samp_annot<- merge(meta.meth, meta.rna, by="row.names")
rownames(samp_annot)<-samp_annot$Row.names
samp_annot$Row.names<-NULL
# merge the annottations from columns Sample.Type and Sample_Group
samp_annot$Sample.Type<-paste(samp_annot$Sample.Type,samp_annot$Sample_Group,sep="_")
#
# remove values before "_" in the Sample.Type column
samp_annot$Sample.Type<-gsub(".*_","",samp_annot$Sample.Type)

mData<-metaData[match(rownames(samp_annot),rownames(metaData)),]
# add the TMB column to the samp_annot
samp_annot<-cbind(samp_annot,mData$TMB)
colnames(samp_annot)[ncol(samp_annot)]<-"TMB"

# save # save X.rna, X.meth and samp_annot
dim(X.rna)
dim(X.meth)
dim(samp_annot)

# drop last two columns of X.rna and X.meth
X.rna<-X.rna[,-c(ncol(X.rna)-1,ncol(X.rna))]
X.meth<-X.meth[,-c(ncol(X.meth)-1,ncol(X.meth))]

# save X.rna, X.meth and samp_annot in one file in rds
saveRDS(list(X.rna, X.meth, samp_annot), file = paste0("aurora_data.rds"))

# remove Normal samples from the dataset
ii<-which(samp_annot$Sample.Type=="Normal tissue")
rmv<-rownames(samp_annot)[ii]

sample_annot<-samp_annot[-which(rownames(samp_annot)%in%rmv),]
jj<-which(colnames(combinedData)%in%rmv)
combinedData<-combinedData[,-jj]
#
samp_annot<-sample_annot[colnames(combinedData),]

nsmpl<-min(table(samp_annot$Sample.Type))
#
#balance primary vs metastatic samples
subset_primary<-sample(which(samp_annot$Sample.Type =="Primary tumor"),nsmpl)
subset_metastatic<-sample(which(samp_annot$Sample.Type=="Metastatic tumor"),nsmpl)
#
metaData<-samp_annot[c(subset_primary,subset_metastatic),]
balanced_data<-combinedData[,c(subset_primary,subset_metastatic)]
#
pca_res<-prcomp(t(balanced_data))
pca_res<-as.data.frame(pca_res$x)
pca_res$Sample.Type<-metaData$Sample.Type

res_pca<-ggplot(pca_res, aes(x=PC2, y=PC1, color=Sample.Type)) +  geom_point(size = 3)+ theme_minimal()
#
graphics.off()
res_pca
library(KPPACart)
#plot the pca elbow plot if i=1 then pick the corresponding number of pca components to project on
#(nfeat<-ceil(dim(balanced_data)[1]*0.01)) # use 10% of the features
nfeat=450
pca <- prcomp(t(balanced_data[sample(rownames(balanced_data))[1:nfeat],]), scale=TRUE)
eigval <- pca$sdev^2
elbow_plot <- data.frame(PC = seq_along(eigval), Eigenvalue = eigval)

ggplot(elbow_plot, aes(x = PC, y = Eigenvalue)) +
  geom_point(size = 3) +
  geom_line() +
  ggtitle("Elbow Plot of PCA Eigenvalues") +
  xlab("Principal Component") +
  ylab("Eigenvalue") +
  scale_x_continuous(breaks = seq(0, length(eigval), by = 10)) +  # Increase x-axis resolution
  theme_minimal()

#make an empty list to store results of kppa
aurora_kppa_results<-list()
# run a loop to get the results of kppa for each of the 4 subtypes randomly selected from the dataset 100 time
nboot<-1
for (i in 1:nboot){
  print(paste("iteration number:",i))
  #
  xx<-as.data.frame(balanced_data)
  # now run KPPACart
  aurora_kppa_results[[i]]<-KPPACart(xx, n_features=nfeat, n_iterations=1000, k_dim=10, exp_clusters=2, n_cores=16)
}
i=1
#genes<-rowData(balanced_data)
scores<-as.data.frame(aurora_kppa_results[[i]]$T)
all_scores<-cbind(scores,metaData)
colnames(all_scores)[1:2]<-c("Proj1","Proj2")
# use ggplot to plot the scores labelled by the PAM50 subtype
#
graphics.off()
kppa_res<-ggplot(all_scores, aes(x=Proj1, y=Proj2, color=Sample.Type)) + geom_point(size = 3) + #geom_density_2d(bins=25)+
  theme_bw()
kppa_res#|res_pca
#
xmat<- aurora_kppa_results[[i]]$BestData#
#
#remove NAs from xmat
#
xmat<-xmat[complete.cases(xmat),]
#
reportDir<-"~/KarakachLab Dropbox/Research/BreastCancer/MultiOmics/Aurora/Results/spreadsheets"
#
cpgs<-grep("cg",rownames(xmat))
#
marks<-ann450k[rownames(xmat)[cpgs],]
mrnas<-rownames(xmat)[-cpgs]
annot<-tbl[which(tbl$GENENAME %in%mrnas),]

write.csv(annot,file=paste(reportDir,"KPPA_genes.csv",sep="/"))
write.csv(marks,file=paste(reportDir,"KPPA_probes.csv",sep="/"))

#show rownames that do not start with cg



sample_annot<-metaData
# replace black or african american with black
sample_annot$Race<-gsub("Black or African American","Black",sample_annot$Race)
# replace not reported with unknown
sample_annot$Race<-gsub("Other","Unknown",sample_annot$Race)
#
# convert age from days to years
#
#colnames(sample_annot)[3]<-“age”
#sample_annot$age<-sample_annot$age/365
base_mean = rowMeans(xmat)
mat_scaled = t(apply(xmat, 1, scale))
l = rownames(xmat)[cpgs]
desc<-ann450k[l,c("chr","pos","strand","UCSC_RefGene_Name","UCSC_RefGene_Accession","Relation_to_Island","Islands_Name","Probe_maf")]
#remove rownames that are NA
desc<-desc[-which(is.na(rownames(desc))),]
#collapse USC_RefGene_Name column entries to remove anything after semicolon
desc$UCSC_RefGene_Name<-gsub(";.*","",desc$UCSC_RefGene_Name)
#create a new column combining UCSC_RefGene_Name and rownames of desc seperated by semi colon
desc$gene<-paste0(desc$UCSC_RefGene_Name,";",rownames(desc))

#genes_for_show<-paste(genes_to_show)
genes_to_show<-c(rownames(desc),mrnas[1:5])#top_genes[1:ceil(dim(xmat)[1]*0.1)]#show 10% of the most important genes
geneLabels<-c(desc$gene,mrnas[1:5])
mark_at = which(rownames(xmat) %in% genes_to_show)

ha1 = rowAnnotation(genelabels = anno_mark(at = mark_at, labels = geneLabels))

library(ComplexHeatmap)
library(circlize)
#---------------------------------------------
#  Now make the heatmap annotation colors
#----------------------------------------
TNBCtype<-as.factor(sample_annot$TNBC_Subtype)
PAM50<-as.factor(sample_annot$PAM50.Call)
#stage<-as.factor(sample_annot$stage)
Race<-as.factor(sample_annot$Race)
#
samp_annot<-sample_annot[,c("Sample.Type","PAM50.Call","TNBC_Subtype","Race","TMB")]
#------------------------
#
column_tree = hclust(dist(t(mat_scaled)))
column_order = column_tree$order
graphics.off()
#met_col = c(“Metastasis, NOS” = “#FF0099”, “No Metastasis” = “#66FF00”)
ha = HeatmapAnnotation(TMB = anno_points(sample_annot[[16]],
                                         gp = gpar(col = ifelse(sample_annot[[16]] > 50, "red", "green")),
                                         height = unit(3, "cm")), #dkfz_cluster = phenotype[[1]],
                       SampleType = sample_annot[[6]],
                       PAM50 = sample_annot[[14]],
                       TNBCSubtype = sample_annot[[15]],
                       Race = sample_annot[[10]],
                       col = list(PAM50 = structure(names = c("Basal","Claudin", "Her2","LumA", "LumB", "Normal"), brewer.pal(6, "Dark2")),
                                  Race = structure(names = c("White", "Asian", "Black", "Unknown"),
                                                   c("white", "#FFFF33", "#2F4F4F", "grey")),
                                  SampleType=c("Primary tumor"="#5AB4AC","Metastatic tumor"="#FF0099"),
                                  TNBCSubtype = structure(names = c("BL1", "BL2", "IM", "notTNBC", "LAR","M","MSL","UNS"), brewer.pal(8, "Set3"))),
                       annotation_name_side = "left",
                       na_col = "grey", border = TRUE,
                       show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
                       show_annotation_name = FALSE,
                       annotation_legend_param = list(PAM50 = list(title = "PAM50 Classification"),
                                                      Race= list(title="Race"),
                                                      SampleType = list(title = "Sample Type"),
                                                      TNBCSubtype = list(title = "TNBC Subtype")))
#
col_fun = colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C"))
ht_list<-Heatmap(mat_scaled, col = col_fun, name = "Scaled \nRNAseq counts",
                 #column_order = column_order,
                 clustering_distance_columns = "spearman",#"spearman",
                 show_row_dend = FALSE,
                 show_column_dend = TRUE,
                 show_column_names = FALSE,
                 show_row_names = F,
                 bottom_annotation = ha,
                 right_annotation = ha1,
                 column_title = paste0("BRCA samples (n = ", ncol(mat_scaled), ")"),
                 #row_split = factor(cgi2, levels = c("Island", "Shore", "Shelf", "OpenSea")),
                 row_title_gp = gpar(col = "#FFFFFF00"))

# draw the heatmap
draw(ht_list, row_title = paste0("Top GeneEpression Features (n = ", nrow(mat_scaled), ")"),
     annotation_legend_side = "left", heatmap_legend_side = "left")

annotation_titles = c(TNBCSubtype = "TNBC Subtype",
                      Race = "Race",
                      PAM50 = "PAM50",
                      SampleType = "Sample Type")

for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    grid.text(annotation_titles[an], unit(-2, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}
#
decorate_annotation("TMB", {
  grid.text("TMB", unit(-10, "mm"), rot = 90, just = "left")
  grid.rect(gp = gpar(fill = NA, col = "black"))
  grid.lines(unit(c(0, 1), "npc"), unit(c(50, 50), "native"), gp = gpar(lty = 20))
})
#
# Now load TCGA data
#
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Locations)

dat<-readRDS("~/KarakachLab Dropbox/Research/BreastCancer/Genomics/TCGA_BRCA_methylation_data.rds")
mat<-assay(dat)
tcga_annot<-dat@colData
tcga_annot$sample_type<-gsub("Primary Tumor","Tumor",tcga_annot$sample_type)
#replace Solid Tissue Normal with Normal
tcga_annot$sample_type<-gsub("Solid Tissue Normal","Normal",tcga_annot$sample_type)
tcga_annot$sample_type<-gsub("Metastatic","Metastatic",tcga_annot$sample_type)
#
tcga_annot$CNV_Clusters<-tcga_annot$`paper_CNV Clusters`
tcga_annot$DNA_Methylation_Clusters<-tcga_annot$`paper_DNA.Methylation Clusters`
tcga_annot$mRNA_Clusters<-tcga_annot$`paper_mRNA Clusters`
#
tcga_annot<-as_tibble(tcga_annot) %>% dplyr::select(barcode, sample_type,race,age_at_index,paper_BRCA_Subtype_PAM50, metastasis_at_diagnosis,CNV_Clusters, DNA_Methylation_Clusters,paper_BRCA_Pathology,mRNA_Clusters)
#
tcga_annot$race<-ifelse(tcga_annot$race=="black or african american","black",tcga_annot$race)
#rename american indian or alaska native to native
tcga_annot$race<-ifelse(tcga_annot$race=="american indian or alaska native","native",tcga_annot$race)
#rename some column names in samp_annot
colnames(tcga_annot)<-c("barcode","sample_type","race","age","subtype","metastasis","CNV_Clusters","DNA_Methylation_Clusters","BRCA_Pathology","mRNA_Clusters")

data(SNPs.137CommonSingle)
data(Islands.UCSC)

l = rownames(xmat)[cpgs]#Locations$chr %in% paste0("chr", 1) & is.na(SNPs.137CommonSingle$Probe_rs)
desc<-ann450k[l,c("chr","pos","strand","UCSC_RefGene_Name","UCSC_RefGene_Accession","Relation_to_Island","Islands_Name","Probe_maf")]
#collapse USC_RefGene_Name column entries to remove anything after semicolon
desc$UCSC_RefGene_Name<-gsub(";.*","",desc$UCSC_RefGene_Name)
#create a new column combining UCSC_RefGene_Name and rownames of desc seperated by semi colon
desc$gene<-paste0(desc$UCSC_RefGene_Name,";",rownames(desc))
#mat = mat[l, ]
#Get subsets for locations of probes and the annotation to CpG Islands accordingly.
cgi = Islands.UCSC$Relation_to_Island[l]
loc = Locations[l, ]
#

mat<-mat[,tcga_annot$barcode]
# Separate the matrix into a matrix for tumor samples and a matrix for normal samples.
# Also modify column names for the tumor samples to be consistent with the phenotype data which we will read later.
#
mat1 = as.matrix(mat[, grep("Tumor", tcga_annot$sample_type)])   # tumor samples
mat2 = as.matrix(mat[, grep("Normal", tcga_annot$sample_type)])  # normal samples
#
#filter barcodes in samp_annot to match the column names in mat1
phenotype<-tcga_annot[which(tcga_annot$barcode %in% colnames(mat1)),]
#
#Extract the top 8000 probes with most variable methylation in the tumor samples,
# and also subset other information correspondingly.
#
ind = order(rowVars(mat1, na.rm = TRUE), decreasing = FALSE)[1:8000]
m1 = mat1[ind, ]
m2 = mat2[ind, ]
cgi2 = cgi[ind]
cgi2 = ifelse(grepl("Shore", cgi2), "Shore", cgi2)
cgi2 = ifelse(grepl("Shelf", cgi2), "Shelf", cgi2)
loc = loc[ind, ]
#
#For each probe, find the distance to the closest TSS. pc_tx_tss.bed contains
# positions of TSS from protein coding genes.

gr = GRanges(loc[, 1], ranges = IRanges(loc[, 2], loc[, 2]+1))
tss = read.table("data/pc_tx_tss.bed", stringsAsFactors = FALSE)
tss = GRanges(tss[[1]], ranges = IRanges(tss[, 2], tss[, 3]))

tss_dist = distanceToNearest(gr, tss)
tss_dist = tss_dist@elementMetadata$distance
#
m1[is.na(m1)] = 0.5
m2[is.na(m2)] = 0.5
# get z-scores for the methylation data
m1_scaled = t(scale(t(m1)))
m2_scaled = t(scale(t(m2)))
#
ii<-which(gr@seqnames=="chr1")
jj<-which(gr[ii,]@ranges@start >= 2136845 & gr[ii,]@ranges@start <= 2169373)
ii<-which(rownames(m1) %in% rownames(xmat))
mark_at<-ii[jj]
#mark_at<-which(gr@ranges@start >= 2136845 & gr@ranges@start <= 2169373) # promoter region for Rap1GAP is: 86980530 and 86980552
rap1cpgs<-rownames(m1)[mark_at]
ha1 = rowAnnotation(CpGs = anno_mark(at = mark_at, labels = rap1cpgs))


