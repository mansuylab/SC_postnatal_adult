data(exprs)
## divide exprs into two parts corresponding to condition 1 
## (exprs.1) and condition 2 (exprs.2) respectively
expGenes<-rownames(exprs)
exprs<-exprs[1:100,]
exprs.1<-exprs[1:100,1:16]
exprs.2<-exprs[1:100,17:63]

data(tf2target)
DCp.res<-DCp(exprs.1,exprs.2,
             link.method = 'qth',cutoff=0.25)
DCe.res<-DCe(exprs.1,exprs.2,
             link.method = 'qth',cutoff=0.25,nbins=10,p=0.1)
DCsum.res<-DCsum(DCp.res,DCe.res,DCpcutoff=0.25,DCecutoff=0.4)

## rank all the potential TFs
data(tf2target)
DRrank.TED.res<-DRrank(DCsum.res$DCGs, DCsum.res$DCLs, 
                       tf2target, expGenes,
                       rank.method=c('TED','TDD')[1],
                       Nperm=0)

DRrank.TED.res[1:3,]

DRrank.TDD.res<-DRrank(DCsum.res$DCGs, DCsum.res$DCLs, 
                       tf2target, expGenes,
                       rank.method=c('TED','TDD')[2],
                       Nperm=0)

DRrank.TDD.res[1:3,]


DCGL::DRplot(DCsum.res$DCGs, DCsum.res$DCLs, 
             tf2target, expGenes)




load("/mnt/IM/anno/regulon_weighted.RData")


convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = g , mart = human,
                   attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  # humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  # print(head(humanx))
  return(genesV2)
}

genes <- convertHumanGeneList(g)


t1 <- lapply(regulon, data.frame)

tf.df <- NULL
for(i in 1:length(t1)){
  tab <- data.frame(HGNC.symbol = rownames(t1[[i]]), TF = names(t1)[i],
                    tfmode = t1[[i]]$tfmode, likelihood = t1[[i]]$likelihood, 
                    stringsAsFactors = F)
  tf.df <- rbind.data.frame(tf.df, tab, stringsAsFactors = F)
}

tf.df.mm <- merge(genes, tf.df)


test <- DRsort(DCGs = pnd8.pnd15$dca$DCGs, DCLs = pnd8.pnd15$dca$DCLs, tf2target = tf.df.mm[,c(3,2)], expGenes = pnd8.pnd15$dca$DCGs$Gene)

test1 <- DRrank(DCGs = pnd8.pnd15$dca$DCGs, DCLs = pnd8.pnd15$dca$DCLs, tf2target = tf.df.mm[,c(3,2)], expGenes = pnd8.pnd15$dca$DCGs$Gene)
