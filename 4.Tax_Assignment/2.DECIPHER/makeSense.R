phytab <- function(x){phylum<-sapply(x,
                 function(x) {
                   w <- which(x$rank=="phylum")
                   if (length(w) != 1) {
                     "unknown"
                   } else {
                     x$taxon[w]
                   }
                 })
mytab<-table(phylum)
data.frame(mytab)
}


taxas<- sort(list.files(getwd(), pattern = "rds", full.names = FALSE))
allrds<-lapply(taxas,function(x){y<-readRDS(x);df<-phytab(y);df$iter<-x;return(df)})
alldf<-data.table::rbindlist(allrds)


#Silva, Ribosomal Database Project  and GreenGenes
alldf[str_detect(alldf$iter,'SILVA'),'db']<-"SILVA"
alldf[str_detect(alldf$iter,'gtdb'),'db']<-"gtdb"
alldf[str_detect(alldf$iter,'RDP'),'db']<-"RDP"

#pooling allows information to be shared across samples, which makes it easier to resolve rare variants that were present as singletons or doubletone in one sample but were present many times across samples.
alldf[str_detect(alldf$iter,'consensus'),'chimeraremoval']<-"consensus"
alldf[str_detect(alldf$iter,'sample'),'chimeraremoval']<-"sample"
alldf[str_detect(alldf$iter,'pooled'),'chimeraremoval']<-"pooled"


# If FALSE, samples are read in the provided order until enough reads are obtained. If TRUE, samples are picked at random from those provided.
alldf[str_detect(alldf$iter,'randomizeTrue'),'randomize']<-"true"
alldf[str_detect(alldf$iter,'randomizeFalse'),'randomize']<-"false"

ggplot(alldf,aes(phylum,Freq))+geom_col()+facet_wrap(~db)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))