library('edgeR')

##read in files
samples <- read.delim("samples.txt", stringsAsFactors = FALSE)
samples=samples[samples$sex=='F',]
counts = readDGE(samples$file)$counts

#filter by coverage
noint = rownames(counts) %in% c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
keep = rowSums(counts)>=10 & !noint
counts = counts[keep,]

#normalize
d = DGEList(counts=counts, group=samples$group)
d = calcNormFactors(d)

#estimate tagwise dispersion
d=estimateCommonDisp(d)
d=estimateTagwiseDisp(d)

#test for differential expression
de=exactTest(d,pair=c('A997b','GenomeC'))

tt = topTags(de, n = nrow(d))

#####################optional######################################

nc = cpm(d, normalized.lib.sizes=TRUE)
rn = rownames(tt$table)

#subset output to only target genes or only Or genes of interest
DmojID=read.csv("DmojID.csv", header=TRUE,stringsAsFactors = FALSE) #list of target genes
tt$table$FBgn <- rownames(tt$table)
rownames(tt$table) <- NULL
final=merge(DmojID,tt$table,by.x="Final_verdict",by.y="FBgn")
names(final)[6]='allgenesFDR'
targetgenesFDR=p.adjust(final$PValue, method="fdr")
final=cbind(final,targetgenesFDR)
DmojOr=read.csv("Or_list.csv", header=TRUE,stringsAsFactors = FALSE) #list of Or genes
final=merge(DmojOr,final,by.x="Or_list",by.y="Final_verdict")
OrgenesFDR=p.adjust(final$PValue, method="fdr")
final=cbind(final,OrgenesFDR)

######################generate volcano plot#############################

colnames(final)[8]="padj" #column of adjusted p-value (FDR)
with(final, plot(logFC, -log10(padj), pch=20))
with(subset(final, padj<.05), points(logFC, -log10(padj), pch=20, col="green"))
with(subset(final, padj<.01), points(logFC, -log10(padj), pch=20, col="orange"))
with(subset(final, padj<.001), points(logFC, -log10(padj), pch=20, col="red"))

library(calibrate)
with(subset(final, padj<.05), textxy(logFC, -log10(padj), labs=Dmoj_trans, cex=.8))
