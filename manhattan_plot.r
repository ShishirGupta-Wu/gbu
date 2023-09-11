
library(tidyverse)
library(dplyr)
vignette('qqman')
ara11<-read.csv("~/lastversion_gwas_analysis/ara11_50.csv")
setwd("/home/s417377/trans_genes/AT2G18735/maf_filtered/")
Y <- list()
csvfiles<-list.files(pattern = "*.csv")
Y <- lapply(csvfiles, function(x) {
  file <- read.csv(x)
  return(file)
})
for (i in 1:length((csvfiles))){
  Y[[i]]<-mutate(Y[[i]], SNP=paste(CHR, sep="- ", POS))
}
for (i in 1:length((csvfiles))){
  colnames(Y[[i]])[4]<-"chr"
  colnames(Y[[i]])[5]<-"pos"
  colnames(Y[[i]])[6]<-"pvalue"
}
#"AT2G29710_AT2G24600"
#make annotation factor
ann<-rep(1, length(Y[[1]]$pvalue))
ann[with(Y[[1]], chr==2  &  pos>=12648573 &  pos<= 12750343)]<-2
ann[with(Y[[1]], chr==2  &  pos==9137974)]<-3
ann[with(Y[[1]], chr==2  &  pos==9303977)]<-4
ann[with(Y[[1]], chr==2  &  pos==17833766)]<-5
ann[with(Y[[1]], chr==4  &  pos==840051)]<-6
ann<-factor(ann, levels=1:6, labels=c("","-", "9137974","9303977","17833766","840051"))
#draw plot with annotation
manhattan.plot(Y[[1]]$chr, Y[[1]]$pos, Y[[1]]$pvalue, annotate=ann, ylim=c(0,50))
#"AT2G24600_AT2G29710"
#make annotation factor
ann<-rep(1, length(Y[[3]]$pvalue))
ann[with(Y[[3]], chr==2  &  pos>=10402034 &  pos<= 10504593)]<-2
ann[with(Y[[3]], chr==2  &  pos==9137974)]<-3
ann[with(Y[[3]], chr==2  &  pos==9303977)]<-4
ann[with(Y[[3]], chr==2  &  pos==17833769)]<-5
ann[with(Y[[3]], chr==4  &  pos==840051)]<-6
ann<-factor(ann, levels=1:6, labels=c("","-", "9137974","9303977","17833769","840051"))
#draw plot with annotation
manhattan.plot(Y[[3]]$chr, Y[[3]]$pos, Y[[3]]$pvalue, annotate=ann, ylim=c(0,50))
#"AT2G24600_AT4G00970.csv"
#make annotation factor
ann<-rep(1, length(Y[[2]]$pvalue))
ann[with(Y[[2]], chr==4  &  pos>=368327 &  pos<= 471885)]<-2
ann[with(Y[[2]], chr==2  &  pos==8122148)]<-3
ann[with(Y[[2]], chr==2  &  pos==9345189)]<-4
ann[with(Y[[2]], chr==2  &  pos==17833784)]<-5
ann[with(Y[[2]], chr==4  &  pos==839000)]<-6
ann<-factor(ann, levels=1:6, labels=c("","-", "8122148","9345189","17833784","839000"))
#draw plot with annotation
manhattan.plot(Y[[2]]$chr, Y[[2]]$pos, Y[[2]]$pvalue, annotate=ann, ylim=c(0,50))
#"AT4G00970_AT2G24600.csv"
#make annotation factor
ann<-rep(1, length(Y[[6]]$pvalue))
ann[with(Y[[6]], chr==2  &  pos>=10402034 &  pos<= 10504593)]<-2
ann[with(Y[[6]], chr==2  &  pos==8122148)]<-3
ann[with(Y[[6]], chr==2  &  pos==9345189)]<-4
ann[with(Y[[6]], chr==2  &  pos==17833784)]<-5
ann[with(Y[[6]], chr==4  &  pos==839000)]<-6
ann<-factor(ann, levels=1:6, labels=c("","-", "8122148","9345189","17833784","839000"))
#draw plot with annotation
manhattan.plot(Y[[6]]$chr, Y[[6]]$pos, Y[[6]]$pvalue, annotate=ann, ylim=c(0,50))
#"AT3G12020_AT3G19040.csv" 
ann<-rep(1, length(Y[[4]]$pvalue))
ann[with(Y[[4]], chr==3  &  pos>=6517021 &  pos<= 6517021)]<-2
ann[with(Y[[4]], chr==2  &  pos==8122723)]<-3
ann<-factor(ann, levels=1:3, labels=c("","-", "8122723"))
#draw plot with annotation
manhattan.plot(Y[[4]]$chr, Y[[4]]$pos, Y[[4]]$pvalue, annotate=ann, ylim=c(0,50))
#"AT3G19040_AT3G12020.csv" 
#make annotation factor
ann<-rep(1, length(Y[[5]]$pvalue))
ann[with(Y[[5]], chr==3  &  pos>=3775980 &  pos<= 3885146)]<-2
ann[with(Y[[5]], chr==2  &  pos==8122723)]<-3
ann<-factor(ann, levels=1:3, labels=c("","-", "8122723"))
#draw plot with annotation
manhattan.plot(Y[[5]]$chr, Y[[5]]$pos, Y[[5]]$pvalue, annotate=ann, ylim=c(0,50))