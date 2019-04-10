BiocManager::install("affy")
library(affy)
setwd("~/Downloads/Moliner_CELfiles/")
files = list.files(".",
full.names = TRUE)
affy.data = ReadAffy(filenames = files)
eset.rma <- rma(affy.data)
colnames(exprSet.nologs)
colnames(eset.rma)
eset.rma@annotation
eset.rma@phenoData
eset.rma@phenoData@data
pData(eset.rma)
groups <- factor(rep("None", 9), rep("ES", 3), rep("MEF",1), "revES")
groups
groups <- factor(c(rep("None", 9), rep("ES", 3), rep("MEF",1), "revES"))
groups
dMat <- model.matrix(~0 + groups)
colnames(dMat) <- levels(groups)
makeContrasts
cMat <- makeContrasts(levels = colnames(dMat), growth=MEF-ES)
cMat
fit.ls <- lmFit(eset.rma, dMat, method = "ls")
fit.ls
fit.ls <- contrasts.fit(fit.ls, cMat)
eb.ls <- eBayes(fit.ls, proportion = 0.01)
eb.ls
tG <- topTable(eb.ls, coef = "growth", number = 1000, resort = "logFC", genelist = rownames(eb.ls))
tG
readr::write_tsv(tG, "~/Code/SignacBenchMark/Moliner_1000DEG.tsv")
eset.rma
BiocManager::install("mgu74av2.db")
library(mgu74av2.db)
columns(mgu74av2.db)
mgu74av2.db[["ENSEMBL"]]
mgu74av2.db::mgu74av2ENSEMBL
as.data.frame(mgu74av2.db::mgu74av2ENSEMBL)
mgu74av2ENSEMBL <- as.data.frame(mgu74av2.db::mgu74av2ENSEMBL)
mgu74av2ENSEMBL
match(rownames(tG), mgu74av2ENSEMBL$probe_id)
table(is.na(match(rownames(tG), mgu74av2ENSEMBL$probe_id)))
tG <- topTable(eb.ls, coef = "growth", number = 1100, resort = "logFC", genelist = rownames(eb.ls))
table(is.na(match(rownames(tG), mgu74av2ENSEMBL$probe_id)))
tG <- topTable(eb.ls, coef = "growth", number = 1050, resort = "logFC", genelist = rownames(eb.ls))
table(is.na(match(rownames(tG), mgu74av2ENSEMBL$probe_id)))
tG <- topTable(eb.ls, coef = "growth", number = 1050, resort = "logFC", genelist = rownames(eb.ls))
50
tG <- topTable(eb.ls, coef = "growth", number = 1040, resort = "logFC", genelist = rownames(eb.ls))
50
tG <- topTable(eb.ls, coef = "growth", number = 1050, resort = "logFC", genelist = rownames(eb.ls))
table(is.na(match(rownames(tG), mgu74av2ENSEMBL$probe_id)))
tG <- topTable(eb.ls, coef = "growth", number = 1040, resort = "logFC", genelist = rownames(eb.ls))
table(is.na(match(rownames(tG), mgu74av2ENSEMBL$probe_id)))
mgu74av2ENSEMBL <- as.data.frame(mgu74av2.db::mgu74av2GENENAME)
table(is.na(match(rownames(tG), mgu74av2ENSEMBL$probe_id)))
tG <- topTable(eb.ls, coef = "growth", number = 1000, resort = "logFC", genelist = rownames(eb.ls))
table(is.na(match(rownames(tG), mgu74av2ENSEMBL$probe_id)))
mgu74av2ENSEMBL <- as.data.frame(mgu74av2.db::mgu74av2ENTREZID)
table(is.na(match(rownames(tG), mgu74av2ENSEMBL$probe_id)))
mgu74av2ENSEMBL <- as.data.frame(mgu74av2.db::mgu74av2SYMBOL)
table(is.na(match(rownames(tG), mgu74av2ENSEMBL$probe_id)))
mgu74av2ENSEMBL <- as.data.frame(mgu74av2.db::mgu74av2ENSEMBL)
table(is.na(match(rownames(tG), mgu74av2ENSEMBL$probe_id)))
tG <- topTable(eb.ls, coef = "growth", number = 1040, resort = "logFC", genelist = rownames(eb.ls))
table(is.na(match(rownames(tG), mgu74av2ENSEMBL$probe_id)))
!is.na(match(rownames(tG), mgu74av2ENSEMBL$probe_id))
mgu74av2ENSEMBL$ensembl_id[(match(rownames(tG), mgu74av2ENSEMBL$probe_id))]
rownames(tG) <- make.unique(mgu74av2ENSEMBL$ensembl_id[(match(rownames(tG), mgu74av2ENSEMBL$probe_id))])
rownames(tG) <- as.character(make.unique(mgu74av2ENSEMBL$ensembl_id[(match(rownames(tG), mgu74av2ENSEMBL$probe_id))]))
r.names <- make.unique(mgu74av2ENSEMBL$ensembl_id[(match(rownames(tG), mgu74av2ENSEMBL$probe_id))])
r.names[is.na(r.names)] <- "NA"
rownames(tG) <- r.names
tg
tG
readr::write_tsv(tG, "~/Code/SignacBenchMark/Moliner_1000DEG.tsv")
savehistory("~/Downloads/Moliner_CELfiles/process_Moliner.R")
