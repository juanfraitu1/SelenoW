


# Package names
packages <- c("ggrepel", "pheatmap", "DOSE", "clusterProfiler", "ensembldb", "annotables", "apeglm", "RColorBrewer", "vsn", "DESeq2", "ggplot2", "reshape2", "gage", "gageData", "AnnotationDbi", "tidyr", "dplyr", "org.Mm.eg.db", "gplots", "GO.db", "GOstats", "AnnotationHub", "pathview", "plyr", "tidyverse", "purrr", "GSEABase", "igraph", "vissE", "patchwork", "msigdb", "tibble", "msigdbr","fgsea", "GEOquery","EnhancedVolcano","ggvenn","readxl","msigdbr","Biobase","DEGreport","viridis")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))









tabler<-function(res){
  res=as.data.frame(res)
  res=rownames_to_column(res,var="gene")
  sig <- res[which(res[,7]<0.05 & abs(res[,3])>2),]
  #Remove unannotated genes ending in RIK
  sig_red=sig[grep("*Rik", sig$gene, invert=TRUE),]
  #Remove unnanotated genes starting with Gm
  sig_red=sig_red[grep("Gm*", sig$gene, invert=TRUE),]
  return(sig)
}


norm_sig <- function(nc, vec){
  require(tibble)
  require(dplyr)
  nor=dplyr::filter(nc,nc$gene %in% vec$gene)
  nor=data.frame(nor) 
  nor=column_to_rownames(nor,var = "gene") 
  
  order=c("WT0_1","WT0_2","WT0_3","KO0_1","KO0_2","KO0_3","WT4_1","WT4_2","WT4_3","KO4_1","KO4_2","KO4_3","WT8_1","WT8_2","WT8_3","KO8_1","KO8_2","KO8_3","WT24_1","WT24_2","WT24_3","KO24_1","KO24_2")
  
  
  nor=nor[,order]
  return(nor)
}



fc2list=function(nor, res){
  ids<-bitr(rownames(nor), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Mm.eg.db)
  dedup_ids = ids[!duplicated(ids[c("ENTREZID")]),]
  df2 = res[res$gene %in% dedup_ids$SYMBOL,]
  df2$Y = dedup_ids$ENTREZID
  # Name vector with ENTREZ ids
  lis <- df2$log2FoldChange
  names(lis) <- df2$Y
  # omit any NA values 
  lis<-na.omit(lis)
  # sort the list in decreasing order (required for clusterProfiler)
  lis = sort(lis, decreasing = TRUE)
  return(lis)
}




gse <- function(kegg_gene_list){
  gseGO(gene = kegg_gene_list, 
        keyType = "ENTREZID",
        OrgDb = org.Mm.eg.db, 
        ont = "BP",     pAdjustMethod = "BH",eps = 0)
}




GSEAa=function(kegg_gene_list){
  GSEA(kegg_gene_list,
       TERM2GENE = term2gene,
       TERM2NAME = term2name,
       pvalueCutoff = 1.00,
       minGSSize = 15,
       maxGSSize = 500)
}


go <- function(sig_genes){
  enrichGO(gene = names(sig_genes), 
           universe = as.character(grcm38$entrez),
           keyType = "ENTREZID",
           OrgDb = org.Mm.eg.db, 
           ont = "BP", 
           pAdjustMethod = "BH", 
           qvalueCutoff = 0.05, 
           readable = TRUE)}



kgg <- function(kegg_gene_list){
  gseKEGG(geneList= kegg_gene_list, organism= "mmu", minGSSize    = 3,maxGSSize    = 1800, pvalueCutoff = 0.05,pAdjustMethod = "BH",by = "fgsea")
}

reactom=function(kegg_gene_list){
  enrichPathway(gene=names(kegg_gene_list),pvalueCutoff=0.05, readable=T,organism = "mouse")
}










term2gene <- msigdbr(species = "Mus musculus", category = "H") %>%
  select(gs_name, entrez_gene)
term2name <- msigdbr(species = "Mus musculus", category = "H") %>%
  select(gs_name, gs_description) %>%
  distinct()





countData <- read_excel("Gene_count_matrix_DEseq2_V4_NOT normalized.xlsx")
colnames(countData)=c("gene_id","WT0_1","WT0_2","WT0_3","WT4_1","WT4_2","WT4_3","WT8_1","WT8_2","WT8_3","WT24_1","WT24_2","WT24_3","KO0_1","KO0_2","KO0_3","KO4_1","KO4_2","KO4_3","KO8_1","KO8_2","KO8_3","KO24_1","KO24_2")
countData=countData[,c(1:4,14:16)]
countData=aggregate(. ~gene_id , data = countData, sum)
countData=column_to_rownames(countData,'gene_id')
genotype=c(rep("WT",3),rep("KO",3))
#time=c(rep("0",3),rep("0",3))
colData=as_tibble(cbind(colnames(countData),genotype))

colData$genotype=relevel(factor(colData$genotype),ref = "WT")


colnames(colData)=c("samplename","genotype")
dds=DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~genotype)
dds=DESeq(dds)
dds =dds[rowSums(counts(dds))>5,]
sizeFactors(dds)
#pca=plotPCA(rld, intgroup=c("group"))
resultsNames(dds)


normalized_counts_c=counts(dds, normalized=TRUE)


normalized_counts_c <- normalized_counts_c %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()


res_base=results(dds)

ko_base_tb=tabler(res_base)

nko_base=norm_sig(normalized_counts_c,ko_base_tb)

ko_base_genelist=fc2list(nor = nko_base,res = ko_base_tb )




gse_base <- gse(ko_base_genelist)

##dotplot(gse_base)

go_base=go(ko_base_genelist)

dotplot(go_base)

kgg_base=kgg(ko_base_genelist)

##dotplot(kgg_base)

reactom_base=reactom(ko_base_genelist)

dotplot(reactom_base)







  
  
  
countData <- read_excel("Gene_count_matrix_DEseq2_V4_NOT normalized.xlsx")
colnames(countData)=c("gene_id","WT0_1","WT0_2","WT0_3","WT4_1","WT4_2","WT4_3","WT8_1","WT8_2","WT8_3","WT24_1","WT24_2","WT24_3","KO0_1","KO0_2","KO0_3","KO4_1","KO4_2","KO4_3","KO8_1","KO8_2","KO8_3","KO24_1","KO24_2")
countData=aggregate(. ~gene_id , data = countData, sum)
countData=column_to_rownames(countData,'gene_id')
genotype=c(rep("WT",12),rep("KO",11))
time=c(rep("0",3),rep("4",3),rep("8",3),rep("24",3),rep("0",3),rep("4",3),rep("8",3),rep("24",2))
colData=as_tibble(cbind(colnames(countData),genotype,time))

colData$genotype=relevel(factor(colData$genotype),ref = "WT")

colData$time=relevel(factor(colData$time), ref = "0")

colnames(colData)=c("samplename","genotype","time")
dds=DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~genotype+time+time:genotype)
dds=DESeq(dds)
dds =dds[rowSums(counts(dds))>5,]
sizeFactors(dds)
#pca=plotPCA(rld, intgroup=c("group"))
resultsNames(dds)



normalized_counts_c=counts(dds, normalized=TRUE)


normalized_counts_c <- normalized_counts_c %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()






# KO_4h_vs_0h : 4h vs 0h in KO.
# It compares the expression levels between KO and WT at 4 hours relative to time 0
# The ko_4h_0h  takes into account the change in gene expression from time 0 to 4 hours for both genotypes
# and then estimates the differential expression between KO and WT specifically at 4 hours, relative to time 0.
ko_4h_0h=results(dds, contrast=list(c("time_4_vs_0","genotypeKO.time4")), test="Wald")
write.csv(as.data.frame(ko_4h_0h),"KO_time_4h_vs_0h.csv")

ko_4h_0h_tb=tabler(ko_4h_0h)

nko_4h_0h=norm_sig(normalized_counts_c,ko_4h_0h_tb)

ko_4h_0_genelist=fc2list(nor = nko_4h_0h,res = ko_4h_0h_tb )


GSEAa_4h_0=GSEAa(ko_4h_0_genelist)

GSEAa_4h_0_<- GSEAa_4h_0
GSEAa_4h_0_@result$Description <- GSEAa_4h_0_@result$ID
dotplot(GSEAa_4h_0_)


gse4h_0 <- gse(ko_4h_0_genelist)

go4h_0=go(ko_4h_0_genelist)


kgg_4h_0=kgg(ko_4h_0_genelist)

reactom_4h_0=reactom(ko_4h_0_genelist)


# 
# 
# mydf <- data.frame(Entrez=names(ko_4h_0_genelist),FC=as.numeric(ko_4h_0_genelist))
# mydf <- mydf[abs(mydf$FC) > 2,]
# mydf$group <- "upregulated"
# mydf$group[mydf$FC < 0] <- "downregulated"
# mydf$othergroup <- "C"
# mydf$othergroup[abs(mydf$FC) > 2] <- "A"
# 
# formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichGO", OrgDb = org.Mm.eg.db)
#   
# 
# 
# head(formula_res)
# 
# dotplot(formula_res, x="group") + facet_grid(~othergroup)
# 
# 




# KO_8h_vs_0h : 8h vs 0h in KO.
# It compares the expression levels between KO and WT at 8 hours relative to time 0
ko_8h_0h=results(dds, contrast=list(c("time_8_vs_0","genotypeKO.time8")), test="Wald")
write.csv(ko_8h_0h,"KO_time_8h_vs_0h.csv")


ko_8h_0h_tb=tabler(ko_8h_0h)

nko8h_0h=norm_sig(normalized_counts_c,ko_8h_0h_tb)







ko_8h_0_genelist=fc2list(nor = nko8h_0h,res = ko_8h_0h_tb )




gse8h_0 <- gse(ko_8h_0_genelist)




GSEAa_8h_0=GSEAa(ko_8h_0_genelist)

GSEAa_8h_0_<- GSEAa_8h_0
GSEAa_8h_0_@result$Description <- GSEAa_8h_0_@result$ID
dotplot(GSEAa_8h_0_)



go8h_0=go(ko_8h_0_genelist)



kgg_8h_0=kgg(ko_8h_0_genelist)

reactom_8h_0=reactom(ko_8h_0_genelist)


# 
# 
# mydf2 <- data.frame(Entrez=names(ko_8h_0_genelist),FC=as.numeric(ko_8h_0_genelist))
# mydf2 <- mydf[abs(mydf$FC) > 2,]
# mydf2$group <- "upregulated"
# mydf2$group[mydf$FC < 0] <- "downregulated"
# mydf2$othergroup <- "A"
# mydf2$othergroup[abs(mydf$FC) > 2] <- "B"
# 
# mydf3=rbind(mydf,mydf2)
# 
# formula_res <- compareCluster(Entrez~group+othergroup, data=mydf2, fun="gseGO", OrgDb = org.Mm.eg.db)
#   
# head(formula_res)
# 
# 
# 
# dotplot(formula_res, x="group") + facet_grid(~othergroup)
# 



# KO_24h_vs_0h : 24h vs 0h in KO.
# It compares the expression levels between KO and WT at 24 hours relative to time 0
ko_24h_0h=results(dds, contrast=list(c("time_24_vs_0","genotypeKO.time24")), test="Wald")
write.csv(ko_24h_0h,"KO_time_24h_vs_0h.csv")

ko_24h_0h_tb=tabler(ko_24h_0h)

nko_24h_0h=norm_sig(normalized_counts_c,ko_24h_0h_tb)


ko_24h_0_genelist=fc2list(nor = nko_24h_0h,res = ko_24h_0h_tb )






gse24h_0 <- gse(ko_24h_0_genelist)


GSEAa_24h_0=GSEAa(ko_24h_0_genelist)
GSEAa_24h_0_<- GSEAa_24h_0
GSEAa_24h_0_@result$Description <- GSEAa_24h_0_@result$ID





go24h_0=go(ko_24h_0_genelist)

kgg_24_0=kgg(ko_24h_0_genelist)

reactom_24h_0=reactom(ko_24h_0_genelist)

# 
# 
# mydf4 <- data.frame(Entrez=names(ko_24h_0_genelist),FC=as.numeric(ko_24h_0_genelist))
# mydf4 <- mydf[abs(mydf$FC) > 2,]
# mydf4$group <- "upregulated"
# mydf4$group[mydf$FC < 0] <- "downregulated"
# mydf4$othergroup <- "D"
# mydf4$othergroup[abs(mydf$FC) > 2] <- "E"
# 
# mydf5=rbind(mydf,mydf2,mydf4)
# 
# formula_res <- compareCluster(Entrez~othergroup:group, data=mydf5, fun="enrichGO", organism="mmu")
#   
# head(formula_res)
# 
# 
# 
# dotplot(formula_res, x="group") + facet_grid(~othergroup)
# 
# 








# 
# 
# gl1=list("4_0"=names(ko_4h_0_genelist),"8_0"=names(ko_8h_0_genelist),"24_0"=names(ko_24h_0_genelist))
# 
# ck1 <- compareCluster(geneCluster = gl1, fun = enrichGO, OrgDb = org.Mm.eg.db)
# 
# 
# dotplot(ck1, showCategory=25)
# 
# #ck2k <- compareCluster(geneCluster = gl1, fun = gseKEGG, organism="mmu")
# 
# ck2 <- compareCluster(geneCluster = gl1, fun = enrichPathway, organism="mouse")
# 
# 
# dotplot(ck2,showCategory=25)





go_mer=merge_result(list("4"=go4h_0,"8"=go8h_0,"24"=go24h_0))

dotplot(mer, showCategory=25, font.size=5)+facet_grid(~Cluster)

gsea_mer=merge_result(list("4"=gse4h_0,"8"=gse8h_0,"24"=gse24h_0))

dotplot(gsea_mer, showCategory=25, font.size=5)+facet_grid(~Cluster)


GSEA_mer=merge_result(list("4"=GSEAa_4h_0_,"8"=GSEAa_8h_0_,"24"=GSEAa_24h_0_))

dotplot(GSEA_mer, showCategory=25, x="NES", font.size=5)+facet_grid(~Cluster)



kgg_mer=merge_result(list("4"=kgg_4h_0,"8"=kgg_8h_0,"24"=kgg_24_0))

dotplot(mer, showCategory=25, x="NES", font.size=5)+facet_grid(~Cluster)


reac_mer=merge_result(list("4"=reactom_4h_0,"8"=reactom_8h_0,"24"=reactom_24h_0))


dotplot(reac_mer, showCategory=25, font.size=5)+facet_grid(~Cluster)
















# Genotype specific effect at 4h : Main effect + genotypeKO.time4
# It compares the expression levels between KO and WT directly at 4 hours, considering the genotype-specific effect.
# It identifies the genes that are differentially expressed between the KO and WT genotypes specifically at 4 hours,
# without considering the change from time 0. Instead, it focuses on the difference between KO and WT directly at 4 hours,
# incorporating both the main effect of genotype and the interaction effect at 4 hours.
ko_4h=results(dds, contrast=list(c("genotype_KO_vs_WT","genotypeKO.time4")), test="Wald")
write.csv(ko_4h,"KO_4h.csv")


ko_4h_tb=tabler(ko_4h)

nko_4h=norm_sig(normalized_counts_c,ko_4h_tb)


ko_4h_genelist=fc2list(nor = nko_4h,res = ko_4h_tb )


go4h=go(ko_4h_genelist)

gse4h <- gse(ko_4h_genelist)


GSEAa_4h=GSEAa(ko_4h_genelist)

GSEAa_4h_<- GSEAa_4h
GSEAa_4h_@result$Description <- GSEAa_4h_@result$ID
##dotplot(GSEAa_4h_)



kgg_4h=kgg(ko_4h_genelist)

reactom_4h=reactom(ko_4h_genelist)






# Genotype specific effect at 8h : Main effect + genotypeKO.time8
ko_8h=results(dds, contrast=list(c("genotype_KO_vs_WT","genotypeKO.time8")), test="Wald")
write.csv(ko_8h,"KO_8h.csv")


ko_8h_tb=tabler(ko_8h)

nko_8h=norm_sig(normalized_counts_c,ko_8h_tb)


ko_8h_genelist=fc2list(nor = nko_8h,res = ko_8h_tb )

go8h=go(ko_8h_genelist)

gse8h <- gse(ko_8h_genelist)


GSEAa_8h=GSEAa(ko_8h_genelist)

GSEAa_8h_<- GSEAa_8h
GSEAa_8h_@result$Description <- GSEAa_8h_@result$ID
##dotplot(GSEAa_8h_)


kgg_8h=kgg(ko_8h_genelist)

reactom_8h=reactom(ko_8h_genelist)



# Genotype specific effect at 24h : Main effect + genotypeKO.time24
ko_24h=results(dds, contrast=list(c("genotype_KO_vs_WT","genotypeKO.time24")), test="Wald")
write.csv(ko_24h,"KO_24h.csv")

ko_24h_tb=tabler(ko_24h)

nko_24h=norm_sig(normalized_counts_c,ko_24h_tb)

ko_24h_genelist=fc2list(nor = nko_24h,res = ko_24h_tb )


gse24h <- gse(ko_24h_genelist)


GSEAa_24h=GSEAa(ko_24h_genelist)

GSEAa_24h_<- GSEAa_24h
GSEAa_24h_@result$Description <- GSEAa_24h_@result$ID
##dotplot(GSEAa_24h_)


go24h=go(ko_24h_genelist)


kgg_24h=kgg(ko_24h_genelist)

reactom_24h=reactom(ko_24h_genelist)










# 
# 
# gl1=list("4_0"=names(ko_4h_0_genelist),"8_0"=names(ko_8h_0_genelist),"24_0"=names(ko_24h_0_genelist))
# 
# ck1 <- compareCluster(geneCluster = gl1, fun = enrichGO, OrgDb = org.Mm.eg.db)
# 
# 
# dotplot(ck1, showCategory=25)
# 
# #ck2k <- compareCluster(geneCluster = gl1, fun = gseKEGG, organism="mmu")
# 
# ck2 <- compareCluster(geneCluster = gl1, fun = enrichPathway, organism="mouse")
# 
# 
# dotplot(ck2,showCategory=25)







go_mer=merge_result(list("4"=go4h,"8"=go8h,"24"=go24h))

dotplot(mer, showCategory=25, font.size=5)+facet_grid(~Cluster)

gsea_mer=merge_result(list("4"=gse4h,"8"=gse8h,"24"=gse24h))

kgg_mer=merge_result(list("4"=kgg_4h,"8"=kgg_8h,"24"=kgg_24h))

dotplot(mer, showCategory=25, x="NES", font.size=5)+facet_grid(~Cluster)


reac_mer=merge_result(list("4"=reactom_4h,"8"=reactom_8h,"24"=reactom_24h))


dotplot(reac_mer, showCategory=50, font.size=5)+facet_grid(~Cluster)












dds=DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~genotype+time+time:genotype)
dds=DESeq(dds, test = "LRT",reduced = ~time)
dds =dds[rowSums(counts(dds))>5,]
sizeFactors(dds)


normalized_counts_c=counts(dds, normalized=TRUE)


normalized_counts_c <- normalized_counts_c %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



res_table=as.data.frame(results(dds,contrast=list(c("genotype_KO_vs_WT","time_24_vs_0","time_4_vs_0","time_8_vs_0"))))




LRT_tb=tabler(res_table)

nLRT=norm_sig(normalized_counts_c,LRT_tb)


LRT_genelist=fc2list(nor = nLRT,res = LRT_tb )


goLRT=go(LRT_genelist)

dotplot(goLRT, showCategory=25, font.size=5)

gseLRT <- gse(LRT_genelist)


dotplot(gseLRT, showCategory=25, font.size=5, x="NES")


GSEAa_LRT=GSEAa(LRT_genelist)

GSEAa_LRT_<- GSEAa_LRT
GSEAa_LRT_@result$Description <- GSEAa_LRT_@result$ID
dotplot(GSEAa_LRT_, showCategory=25)




kggLRT_=kgg(LRT_genelist)

dotplot(kggLRT_, showCategory=25, x="NES")

reactomLRT_=reactom(LRT_genelist)

dotplot(reactomLRT_, showCategory=25)






###pdf("PCA.pdf")
###print(pca)
# dev.off()
# 
# 
# pdf("dispest.pdf")
# print(plotDispEsts(dds))
# dev.off()

















pheatmap(nLRT, scale = "row", clustering_distance_rows = "correlation", fontsize = 3, color = viridis(255), cellwidth = 12.0, main = "All times", cluster_cols = TRUE)


# a=sort(cutree(phea$tree_row, k=2))
# 
# 
# write.csv(names(a[a==1]),"WT_only.csv")
# 
# 
# write.csv(names(a[a==2]),"KO_only.csv")












# 
# idx = grcm38$symbol %in% sig$gene
# 
# ids <- grcm38[idx, ]
# 
# non_duplicates <- which(duplicated(ids$symbol) == FALSE)
# 
# ids <- ids[non_duplicates, ] 
# 
# res_ids <- inner_join(sig, ids, by=c("gene"="symbol"))     
# 
# all_genes <- as.character(res_ids$entrez)
# 
# 
# 
# 
# sig_genes <- as.character(res_ids$entrez)
# 
# 
# res_entrez <- filter(res_ids, entrez != "NA")
# res_entrez <- res_entrez[which(duplicated(res_entrez$entrez) == F), ]
# all_foldchanges <- res_entrez$log2FoldChange
# names(all_foldchanges) <- res_entrez$entrez
# all_foldchanges <- sort(all_foldchanges, decreasing = TRUE)
# 
# symbol_foldchanges=sig$log2FoldChange
# 
# names(symbol_foldchanges)=sig$gene
# 
# 
# all_foldchanges=sort(all_foldchanges,decreasing=TRUE)
# 
# 
# 
# go <- enrichGO(gene = sig_genes, 
#                 universe = as.character(grcm38$entrez),
#                 keyType = "ENTREZID",
#                 OrgDb = org.Mm.eg.db, 
#                 ont = "BP", 
#                 pAdjustMethod = "BH", 
#                 qvalueCutoff = 0.05, 
#                 readable = TRUE)
# 
# 
# 
# 
# 
# 
# go_dp=dotplot(go, showCategory=60)
# 
# 
# up=sig[sig[,3]>0,]
# 
# 
# down=sig[sig[,3]<0,]
# 
# up_entrez <- as.character(up$entrez)
# 
# down_entrez =  as.character(down$entrez)
# 
# upg=res_entrez[res_entrez$entrez %in% up_entrez,]
# 
# up_genes=upg$log2FoldChange
# names(up_genes) <- upg$entrez
# up_genes <- sort(up_genes, decreasing = TRUE)
# 
# downg=res_entrez[res_entrez$entrez %in% down_entrez,]
# 
# down_genes=downg$log2FoldChange
# 
# names(down_genes) <- downg$entrez
# 
# down_genes <- sort(down_genes, decreasing = TRUE)
# 
# 
# 
# 
# go_up <- enrichGO(gene = up_entrez, 
#                 universe = all_genes,
#                 keyType = "ENTREZID",
#                 OrgDb = org.Mm.eg.db, 
#                 ont = "BP", 
#                 pAdjustMethod = "BH", 
#                 qvalueCutoff = 0.05, 
#                 readable = TRUE,)
# 
# 
# 
# 
# 
# go_down <- enrichGO(gene = down_entrez, 
#                 universe = all_genes,
#                 keyType = "ENTREZID",
#                 OrgDb = org.Mm.eg.db, 
#                 ont = "BP", 
#                 pAdjustMethod = "BH", 
#                 qvalueCutoff = 0.05, 
#                 readable = TRUE)
# 
# 
# 
# 
# 
# 
# 
# 
# netplot=cnetplot(go, foldchange=all_foldchanges,fixed=F)
# 
# 
# ids<-bitr(rownames(norm_sig), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Mm.eg.db)
# 
# 
# dedup_ids = ids[!duplicated(ids[c("ENTREZID")]),]
# 
# 
# df2 = res_table_tb[res_table_tb$gene %in% dedup_ids$SYMBOL,]
# 
# df2$Y = dedup_ids$ENTREZID
# 
# 
# # Name vector with ENTREZ ids
# kegg_gene_list <- df2$log2FoldChange
# 
# names(kegg_gene_list) <- df2$Y
# 
# # omit any NA values 
# kegg_gene_list<-na.omit(kegg_gene_list)
# 
# # sort the list in decreasing order (required for clusterProfiler)
# kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
# 
# 
# 
# gse <- gseGO(gene = kegg_gene_list, 
#                 keyType = "ENTREZID",
#                 OrgDb = org.Mm.eg.db, 
#                 ont = c("BP","MF"),     pAdjustMethod = "BH",eps = 0)
# 
# gsea=dotplot(gse, showCategory=50, split=".sign") + facet_grid(.~.sign)
# 
# gsea1=dotplot(gse, showCategory=25, split=".sign", x="NES", decreasing = TRUE,font.size=6)
# 
# 
# 
# #lapply(1:length(gse$Description), function(i) {gseaplot(gse, by = "all", title = gse$Description[i], geneSetID = i)})
# 
# 
# 
# #pdf("gseas.pdf",width = 10)
# #print(lapply(1:length(gse$Description), function(i) {gseaplot(gse, by = "all", title = gse$Description[i], geneSetID = i)}))
# #dev.off()
# 
# 
# kgg <- gseKEGG(geneList= kegg_gene_list, organism= "mmu", minGSSize    = 3,maxGSSize    = 1800, pvalueCutoff = 0.05,pAdjustMethod = "none",by = "fgsea")
# 
# 
# 
# dotplot(kgg , showCategory=25, split=".sign", x="NES", decreasing = TRUE,font.size=6)
# 
# 
# 
# 
# 
# #a=lapply(1:10, function(x) {
#   #pathview(gene.data = foldchanges, pathway.id = kgg$ID[x], species = "mmu")})
# 
# 
# 
# pdf("GO_dotplot.pdf", width=10, height=33)
# print(go_dp)
# dev.off()
# 
# 
# 
# pdf("gsea_dotplot.pdf", width=10, height=33)
# print(gsea)
# dev.off()
# 
# 
# pdf("kegg_dotplot.pdf", width=10, height=35)
# print(kgg)
# dev.off()
# 
# 
# pdf("network.pdf", width=27, height=20)
# print(netplot)
# dev.off()
# 








# 
# normalized_counts=rlog(counts(dds))
# 
# 
# 
# 
# 
# normalized_counts <- normalized_counts %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>%
#   as_tibble()
# 
# normalized_counts$entrez=mapIds(org.Mm.eg.db,key,keys = normalized_counts$gene,column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
# 
# normalized_counts$gene=normalized_counts$entrez
# 
# normalized_counts=normalized_counts[,1:24]
# 
# normalized_counts=aggregate(. ~gene , data = normalized_counts, sum)
# 
# normalized_counts <- normalized_counts %>%
#   data.frame() %>%
#  column_to_rownames(var = "gene")
# 
# 
# require(Biobase)
# normalized_counts1<-new("ExpressionSet", exprs=as.matrix(normalized_counts))
# 
# 
# normalized_counts1$time= colData$group
# 
# 
# 
# 
# pathwaysDF <- msigdbr("mouse", category=c("H"))
# pathways <- split(as.character(pathwaysDF$entrez_gene), pathwaysDF$gs_name)
# 
# 
# set.seed(1)
# gesecaRes <- geseca(pathways, exprs(normalized_counts1), minSize = 10, maxSize = 1000, eps = 0,nPermSimple = 10000)
# 
# geseres=plotGesecaTable(gesecaRes |> head(5), pathways,E=exprs(normalized_counts1),)
# 
# IFNG=plotCoregulationProfile(pathway=pathways[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
#                         E=exprs(normalized_counts1), conditions=normalized_counts1$time)
# 
# 
# INF=plotCoregulationProfile(pathway=pathways[["HALLMARK_INFLAMMATORY_RESPONSE"]],
#                         E=exprs(normalized_counts1), conditions=normalized_counts1$time)
# 
# 
# TNF_NFKB=plotCoregulationProfile(pathway=pathways[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
#                         E=exprs(normalized_counts1), conditions=normalized_counts1$time)
# 
# 
# TNFA=plotCoregulationProfile(pathway=pathways[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
#                         E=exprs(normalized_counts1), conditions=normalized_counts1$time)
# 
# E2F=plotCoregulationProfile(pathway=pathways[["HALLMARK_E2F_TARGETS"]],
#                         E=exprs(normalized_counts1), conditions=normalized_counts1$time)
# 
# 
# 
# 
# 
# 
# 
# fgseaRes <- fgsea(pathways = pathways, 
#                   stats    = all_foldchanges,)
# 
# 
# fgseaResTidy <- fgseaRes %>%
#   as_tibble() %>%
#   arrange(desc(NES))
# 
# # Show in a nice table:
# fgseaResTidy %>% 
#   dplyr::select(-leadingEdge, -ES) %>% 
#   arrange(padj) %>% 
#   DT::datatable()
# 
# 
# 
# fgsearesgraph=ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
#   geom_col(color=fgseaResTidy$NES>0) +
#   coord_flip() +
#   labs(x="Pathway", y="Normalized Enrichment Score",
#        title="Hallmark pathways NES from GSEA") + 
#   theme_minimal()
# 
# 
# 
# pdf("IFNG.pdf",width = 10)
# print(IFNG)
# dev.off()
# 
# 
# pdf("INF.pdf",width = 10)
# print(INF)
# dev.off()
# 
# 
# pdf("TNF_NFKB.pdf",width = 10)
# print(TNF_NFKB)
# dev.off()
# 
# 
# pdf("TNFA.pdf",width = 10)
# print(TNFA)
# dev.off()
# 
# 
# pdf("E2F.pdf",width = 10)
# print(E2F)
# dev.off()
# 
# 
# pdf("geseca.pdf",width = 10)
# print(geseres)
# dev.off()
# 
# 
# 
# pdf("fgsearesgraph.pdf",width = 10)
# print(fgsearesgraph)
# dev.off()
# 
# 
# 







