#originally created on: September 19th, 2014
#by: Nick Sard

#ABOUT: This script was simulate some data for a potential rockfish kinship analysis

#loading in libraries
library(related)
library(dplyr)
library(ggplot2)

#loading my own functions
source("C:/Users/Nick.Sard/Documents/Research/R/Source scripts/share.R")

#setting working directory and reading in data
setwd("C:/Users/Nick.Sard/Documents/Research/R/Data analysis/2015/Rockfish/Simulations/")
list.files()


##################################
### generating breeding matrix ###
##################################

#making the breeding matrix and RS distribution (poisson)
breeding.matrix <- data.frame(moms = c(1,1,2),
                              dads = c(1,2,3))

#putting  reproductive success values in breeding matrix
breeding.matrix$rs <- 2
breeding.matrix

#checking the RS
sum(breeding.matrix$rs)

###################################
### generate allele frequencies ###
###################################

#making some gts get allele frequencies for
#use information from excel spreadsheet
tmp <- sim.pop.gts(n = 1000,loci = 16,pops = 1,k = 20,k.sd = 5)
head(tmp)

#calcuatealf freqs and preping to simulate data
alfs <- alf.freq(tmp,one.pop = F, opt = 2)
head(alfs)

##################################################################
### Generate genotypes for offspring and unrelated individuals ###
##################################################################

df <- sim.gt.data2(alfs = alfs,breeding.matrix = breeding.matrix,error = 0,Nunrelated = 540-sum(breeding.matrix$rs),
                   write2file = F)
head(df)
tail(df)

#getting just the gts from offspring and unrelated individuals
off <- df[grepl("Offspring", df$IDs),-1]
unrel <- df[grepl("Ind.off", df$IDs),-1]

#combing back together
tmp <- rbind(off,unrel)
head(tmp)
dim(tmp)

#adding an underscore at the end here
tmp$IDs <- gsub("Offspring ","Offspring_",tmp$IDs)
tmp$IDs <- gsub("Ind.off ","Ind.off_",tmp$IDs)

#writing to file
write.table(tmp,"sample.540.1.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = F)

###################################
### Now calculating relatedness ###
###################################

head(tmp)
outfile <- coancestry(tmp[1:20,],trioml = 1,wang = 0,lynchli = 0,lynchrd = 0,ritland = 0,quellergt = 0,dyadml = 0,working.directory = getwd(),output.file = FALSE)
str(outfile)
tmp1 <- outfile$relatedness
head(tmp1)

#writing to file
write.table(tmp1,"relatedness.estimates.540.1.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = F)

#########################################
### now processing the related output ###
#########################################
tmp <- tmp1[,c(2,3,5)]

tmp <- tbl_df(tmp)
colnames(tmp) <- c("Ind1","Ind2","TrioEst")
tmp

#now getting mom and dad IDs from names of Ind1 and Ind2
head(tmp)
tail(tmp)
#first offspring parent
tmp$par1 <- tmp$Ind1
tmp$par1 <- gsub(pattern = "Offspring_",replacement = "",x = tmp$par1)
tmp$par1 <- gsub(pattern = "Ind.off_",replacement = "",x = tmp$par1)
tmp$par1[grepl(pattern = "Ind.off_",x = tmp$Ind1)==T] <- "UK.UK.UK"
x <- data.frame(matrix(unlist(strsplit(tmp$par1, ".",fixed = T)), ncol=3, byrow=TRUE), stringsAsFactors=F)
colnames(x) <- c("mom1","dad1","IDkid1")
tmp$par1 <- NULL
tmp <- cbind(tmp,x)
head(tmp)

#second offspring's parents
#first offsprings parents
tmp$par2 <- tmp$Ind2
tmp$par2 <- gsub(pattern = "Offspring_",replacement = "",x = tmp$par2)
tmp$par2 <- gsub(pattern = "Ind.off_",replacement = "",x = tmp$par2)
tmp$par2[grepl(pattern = "Ind.off_",x = tmp$Ind2)==T] <- "UK.UK.UK"
x <- data.frame(matrix(unlist(strsplit(tmp$par2, ".",fixed = T)), ncol=3, byrow=TRUE), stringsAsFactors=F)
colnames(x) <- c("mom2","dad2","IDkid2")
tmp$par2 <- NULL
tmp <- cbind(tmp,x)
head(tmp)

#removing x
x <- NULL

#looking at if they had the same moms and dads
tmp$moms <- ifelse(tmp$mom1==tmp$mom2, "same","diff")
tmp$dads <- ifelse(tmp$dad1==tmp$dad2, "same","diff")

#making the kids with UK for the parents 
tmp$moms[tmp$mom1 == "UK" & tmp$mom2 == "UK"] <- "UK"
tmp$dads[tmp$dad1 == "UK" & tmp$dad2 == "UK"] <- "UK"

#now identifing maternal half-sibs and paternal half-sibs & full sibs
tmp$rel[tmp$moms == "same" & tmp$dads == "same"] <- 0.5
tmp$rel[tmp$moms != "same" & tmp$dads == "same"] <- 0.25
tmp$rel[tmp$moms == "same" & tmp$dads != "same"] <- 0.25
tmp$rel[tmp$moms != "same" & tmp$dads != "same"] <- 0

##now identifing maternal half-sibs and paternal half-sibs & full sibs (for graphics)
tmp$sibs[tmp$moms == "same" & tmp$dads == "same"] <- "Full"
tmp$sibs[tmp$moms != "same" & tmp$dads == "same"] <- "Paternal"
tmp$sibs[tmp$moms == "same" & tmp$dads != "same"] <- "Maternal"
tmp$sibs[tmp$moms != "same" & tmp$dads != "same"] <- "Not"
#tmp

tmp$sibs <- as.factor(tmp$sibs)
#tmp$sibs <- relevel(tmp$sibs,ref="Paternal")
tmp$sibs <- relevel(tmp$sibs,ref="Maternal")
tmp$sibs <- relevel(tmp$sibs,ref="Not")

#range of relatedness values used for cut-off
a <- c(0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50)

#removing tmp1 to save on RAM
tmp1 <- NULL

#trying to clean up tmp so that I dont max out my RAM
tmp <- select(tmp, sibs, TrioEst)
head(tmp)

#a little graph
ggplot(tmp, aes(x=sibs, y=TrioEst))+
  geom_boxplot()

#generating summary tatistics table: mean sdss counts, and counts at various relatedness cut-offs as well as the percentage of total captured
tmp1 <- tmp %>%
  group_by(sibs) %>%
  summarize(
    means = round(mean(TrioEst),2),
    sds = round(sd(TrioEst),2),
    counts = n())
tmp1

#writing summary to file
write.table(tmp, "relatedness.summary.stats.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

#getting some counts
i<- NULL
j <- NULL
x<- NULL

for(j in 1:length(a)){
  
  #first filtering data on correlation values in variable a
  tmp2 <-filter(tmp, TrioEst > a[j])
  tmp2
  
  #now getting counts for each type of estimater and relatedness 
  x<- NULL
  for(i in 1:nrow(tmp1)){  
    #    print(paste(j,i))
    x1 <- nrow(filter(tmp2, sibs == tmp1$sibs[i] & sibs ==  tmp1$sibs[i]))
    x1
    if((x1>0)==F){x1<-0}
    x <- c(x,x1)
  }
  
  #adding them to the original tmp1
  tmp1 <- cbind(tmp1,x)
}

#fixing names
colnames(tmp1) <- c(colnames(tmp1)[1:4],paste("cor",a,sep="_"))
tmp1

#getting percent
tmp1$cor.15.p <- round(tmp1$cor_0.15/tmp1$counts,4)*100
tmp1$cor.2.p <- round(tmp1$cor_0.2/tmp1$counts,4)*100
tmp1$cor.25.p <- round(tmp1$cor_0.25/tmp1$counts,4)*100
tmp1$cor.3.p <- round(tmp1$cor_0.3/tmp1$counts,4)*100
tmp1$cor.35.p <- round(tmp1$cor_0.35/tmp1$counts,4)*100
tmp1$cor.4.p <- round(tmp1$cor_0.4/tmp1$counts,4)*100
tmp1$cor.45.p <- round(tmp1$cor_0.45/tmp1$counts,4)*100
tmp1$cor.5.p <- round(tmp1$cor_0.5/tmp1$counts,4)*100
tmp1

#writing output to file
write.table(tmp1, "relatedness.summary.detailed.txt", sep="\t",append = F, quote=F,sep = "\t",row.names=F,col.names = T)

#fin!