#originally created on: September 2nd, 2014
#by: Nick Sard

#ABOUT: This script was written to get genotypes from the specific reporting groups and assess their genetic variation
#adapted from CROOS analysis September 19th, 2014

#install.packages(c("ggplot2","dplyr","reshape2"))

#loading in libraries
library(ggplot2)
library(reshape2)
library(dplyr)

#loading my own functions
source("C:/Users/Nick.Sard/Documents/Research/R/Source scripts/share.R")

#setting working directory and reading in data
setwd("C:/Users/Nick.Sard/Documents/Research/R/Misc/Rockfish kinship/graphs/")
list.files("Input/")

#reading in gts
tmp <- read.table("Input/relatedness.cutoff.analysis.txt",header=T,sep="\t",stringsAsFactors=F)
head(tmp)
tail(tmp)

tmp1 <- melt(tmp, id.vars = c("Run", "r.type","relatedness", "mean.r","sd.r"))
unique(tmp1$variable)
head(tmp1)
unique(tmp1$Run)
unique(tmp1$relatedness)

#looking at just counts
counts <- tmp1[grepl("count",tmp1[,6])==T,]
head(counts)

#getting a variable to graph
counts$variable1 <- gsub(pattern = "count.","",counts$variable)
counts$variable1 <- gsub(pattern = "total.count",0,counts$variable1)
counts$variable1 <- as.numeric(counts$variable1)/100

#doing some releveling
counts$relatedness <- as.factor(counts$relatedness)
levels(counts$relatedness)
counts$relatedness <- factor(counts$relatedness, levels = c("Full", "Maternal", "Paternal","Not"))

table(counts$relatedness,counts$Run)

head(counts)

#now the plot
ggplot(counts, aes(x=variable1,y=value))+
  facet_grid(relatedness~Run, scales = "free_y")+
  geom_smooth(color="Black",size=1)+
  labs(x="Relatedness cutoff", y= "Total number of individuals")+
  theme_bw()+
  theme(axis.title = element_text(size = 32),
        legend.title = element_text(size = 22),
        legend.text =  element_text(size = 24),
        strip.text = element_text(size = 24),
        axis.text = element_text(size = 24, colour = "Black"))

#looking at just percents
percents <- tmp1[grepl("percent",tmp1[,6])==T,]
p2 <- tmp1[grepl("total.count",tmp1[,6])==T,]
p2$variable <- gsub("total.count","percent.00",p2$variable)
p2$value <- 100
head(percents)
head(p2)
tail(percents)
tail(p2)
percents <- rbind(percents,p2)

#getting a variable to graph
percents$variable1 <- gsub(pattern = "percent.","",percents$variable)
percents$variable1 <- as.numeric(percents$variable1)/100
table(percents$variable1)

#doing some releveling
percents$relatedness <- as.factor(percents$relatedness)
levels(percents$relatedness)
percents$relatedness <- factor(percents$relatedness, levels = c("Full", "Maternal", "Paternal","Not"))


head(percents)
table(percents$Run)

#getting some columns so that I can use facet grid
percents$run[grepl(pattern = "_05_",x = percents$Run)==T]<- "5"
percents$run[grepl(pattern = "_10_",x = percents$Run)==T]<- "10"
percents$run[grepl(pattern = "_25_",x = percents$Run)==T]<- "25"
percents$run[grepl(pattern = "_50_",x = percents$Run)==T]<- "50"
percents$run[grepl(pattern = "_75_",x = percents$Run)==T]<- "75"

head(percents)
table(percents$run,percents$Run)
percents$run <- as.numeric(percents$run)

percents$loci[grepl(pattern = "_20L",x = percents$Run)==T]<- "20 loci"
percents$loci[grepl(pattern = "_10L",x = percents$Run)==T]<- "10 loci"
table(percents$loci,percents$run)
head(percents)


#now the plot
ggplot(percents, aes(x=variable1,y=value))+
  facet_grid(relatedness~loci, scales = "free_y",)+
  geom_smooth(color="Black",size=1)+
  labs(x="Relatedness cutoff", y= "Percent total")+
  theme_bw()+
  theme(axis.title = element_text(size = 32),
        legend.title = element_text(size = 22),
        legend.text =  element_text(size = 24),
        strip.text = element_text(size = 24),
        axis.text = element_text(size = 24, colour = "Black"))

#now making maternal and parental halfsibs into one category and removing the not related indviduals
table(percents$relatedness, percents$loci)
percent1 <- filter(percents, as.character(relatedness) != "Not")
table(percent1$relatedness, percent1$loci)

percent1$relatedness <- as.character(percent1$relatedness)
percent1$relatedness <- ifelse(percent1$relatedness == "Full", "Full-siblings retained", "Half-siblings retained")
table(percent1$relatedness, percent1$loci)

#now the plot
ggplot(percent1, aes(x=variable1,y=value))+
  facet_grid(relatedness~loci, scales = "free_y",)+
  geom_smooth(color="Black",size=1)+
  labs(x="Relatedness cutoff", y= "Percent total")+
  theme_bw()+
  theme(axis.title = element_text(size = 32),
        legend.title = element_text(size = 22),
        legend.text =  element_text(size = 24),
        strip.text = element_text(size = 24),
        axis.text = element_text(size = 24, colour = "Black"))

#now percent not related falsely considered related
head(tmp1)

tmp2 <- data.frame(run = rep(x = unique(tmp$Run),each = length(unique(tmp$r.type))),
           run.type = rep(x = unique(tmp$r.type),times = length(unique(tmp$Run))),
           stringsAsFactors=F)
tmp2 <- select(tmp1, Run, r.type, relatedness, variable) %.%
  filter(relatedness == "Not") %.%
  select(-relatedness)
head(tmp2)
tmp2 <- tmp2[grepl("count",tmp2$variable)==T,]
head(tmp2)

#making a for loop to calculate the percent total related individuals 
i<-1
i <- NULL
for(i in 1:nrow(tmp2)){
  x <- filter(tmp1, Run == tmp2$Run[i] & r.type == tmp2$r.type[i] & variable == tmp2$variable[i])
  x <- sum(x[x$relatedness %in% c("Maternal","Paternal","Full"),"value"])/sum(x$value)
  tmp2$percent.total[i] <- x
}

head(tmp2)
str(tmp2)

#fixing the variable column
tmp2$variable1 <- gsub(pattern = "count.","",as.character(tmp2$variable))
tmp2$variable1 <- gsub(pattern = "total.count","0",as.character(tmp2$variable1))
tmp2$variable1 <- as.numeric(tmp2$variable1)/100

head(tmp2)
tail(tmp2)

#getting some columns so that I can use facet grid
tmp2$run[grepl(pattern = "_05_",x = tmp2$Run)==T]<- "5"
tmp2$run[grepl(pattern = "_10_",x = tmp2$Run)==T]<- "10"
tmp2$run[grepl(pattern = "_25_",x = tmp2$Run)==T]<- "25"
tmp2$run[grepl(pattern = "_50_",x = tmp2$Run)==T]<- "50"
tmp2$run[grepl(pattern = "_75_",x = tmp2$Run)==T]<- "75"
head(tmp2)
tmp2$run <- as.numeric(tmp2$run)
table(tmp2$run,tmp2$Run)
tmp2$loci[grepl(pattern = "_10",x = tmp2$Run)==T]<- "10 loci"
tmp2$loci[grepl(pattern = "_20",x = tmp2$Run)==T]<- "20 loci"

table(tmp2$Run,tmp2$loci)

#fixing the names
head(tmp2)                
str(tmp2)

#now for a plot
ggplot(tmp2, aes(x=variable1, y=percent.total*100))+
  facet_wrap(~loci,nrow=1)+
  geom_smooth(color="black",size=1)+
  labs(x="Relatedness cutoff", y= "Perent total of related individuals")+
  theme_bw()+
  theme(axis.title = element_text(size = 32),
        legend.title = element_text(size = 22),
        legend.text =  element_text(size = 24),
        strip.text = element_text(size = 24),
        axis.text = element_text(size = 24, colour = "Black"))

#now adding to the percent1
head(percent1)
head(tmp2)

#changing a bunch of stuff with tmp2 so that it can be rbinded with percent1
tmp2$relatedness <- "True Related Individuals"
tmp2 <- move.me(data = tmp2,tomove = "relatedness",where = "before",ba = "variable")
head(tmp2)
head(percent1)
percent1$mean.r <- NULL
percent1$sd.r <- NULL
names(tmp2)[5] <- "value"
tmp2$value <- round(tmp2$value*100,2)
head(tmp2)
head(percent1)

#double checking the names are right and then rbinding
table(colnames(tmp2)==colnames(percent1))
tmp3 <- rbind(tmp2,percent1)
head(tmp3)

#making somem reordering in relatedness
tmp3$relatedness <- as.factor(tmp3$relatedness)
levels(tmp3$relatedness)
tmp3$relatedness <- relevel(tmp3$relatedness, ref="Half-siblings retained")
tmp3$relatedness <- relevel(tmp3$relatedness, ref="Full-siblings retained")
tmp3$relatedness <- relevel(tmp3$relatedness, ref="True Related Individuals")

table(tmp3$run)
table(tmp3$loci,tmp3$run)
table(tmp3$loci,tmp3$Run)
tmp3$run <- as.factor(tmp3$run)
levels(tmp3$run)


#now the graph
ggplot(tmp3, aes(x=variable1,y=value, fill=run))+
  facet_grid(relatedness~loci)+
  geom_smooth(color="Black",size=1,method="loess")+
  scale_y_continuous(limits = c(0, 100))+
  labs(x="Relatedness cutoff", y= "Percent total", fill="Percent of\nrelated Individuals")+
  theme_bw()+
  theme(axis.title = element_text(size = 32),
        legend.title = element_text(size = 22),
        legend.text =  element_text(size = 24),
        strip.text = element_text(size = 22),
        axis.text = element_text(size = 24, colour = "Black"))

table(is.na(tmp3$run))
table(is.na(tmp3$relatedness))
table(is.na(tmp3$variable1))
table(is.na(tmp3$value))
table(is.na(tmp3$loci))

