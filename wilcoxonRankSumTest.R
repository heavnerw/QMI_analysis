#Whitney E. Heavner
#Wilcoxon Rank Sum of PiSCES (Co-IP pairs) assigned to different tSNE clusters
#To identify which PiSCES are the most different between clusters
#Updated 20191213
#Updated 20200923

#install necessary packages

install.packages("dplyr")
install.packages("magrittr")

library("dplyr")
library("magrittr")

#Let's load average log2Fold Change values for all interactions in all sets labeled by their tSNE-assigned group

setwd("~/Dropbox/")

#Input formatted as .csv
#Column 1 = Condition
#Column 2 = Group (tSNE Group)
#Column 3-Column 3+N_PiSCES = PiSCES (IP_Probe)
#Rows=Median Fluorescent Intensity (MFI)

groupData<-read.csv("full_tSNE_input_wGroups.csv")  

#Now, let's look at the structure of the input data

head(groupData)
attach(groupData)
names(groupData)
class(Interaction)
class(Group)
class(Shank3_SAPAP)
levels(Group)

#Let's also look at boxplots of the relationship between select interactions and group assignment

boxplot(Homer1_panShank ~ Group)
boxplot(GluR2_GluR2 ~ Group)
boxplot(PSD95_PSD95 ~ Group)
boxplot(PSD95_GluR2 ~ Group)
boxplot(Ube3a_Shank3 ~ Group)
boxplot(Homer1_mGluR5 ~ Group)
boxplot(Fyn_Homer1 ~ Group)
boxplot(Fyn_Fyn ~ Group)
boxplot(Ube3a_Fyn ~ Group)
boxplot(Ube3a_Ube3a ~ Group)

#Now let's calculate the wilcoxon p-value for a single interaction between two groups
#First we need to load the data with two groups you want to compare
#Input formatted same as above; contains only the groups you want to compare.

groups<-read.csv("full_tSNE_input_groups_3_4.csv") 

#Check the structure of the input, attach, and check column names

head(groups)
attach(groups)
names(groups)

#To test the probability that a single interaction is different between the two groups, run wilcox.test for that interaction

wilcox.test(Homer1_panShank ~ Group, mu=0, alt="two.sided", conf.int=T, conf.level=0.95, paired=F, exact=F, correct=T)

#you can also designate the column number rather than the interaction name if you're a masochist (not recommended!)
#define the test for the PiSCES in the column as "model" or someting

model <- wilcox.test(groups[ ,98] ~ Group, mu=0, alt="two.sided", conf.int=T, conf.level=0.95, data = groups, paired=F, exact=F, correct=T)
model

#or we can run the comparison for all interactions if you're lazy (recommended!)
#define the test for all PiSCES in columns 3:3+N PiSCES as "modelList"

modelList <- list()
for(i in 3:380){
  fmla <- formula(paste(names(groups)[i], "~ Group"))
  modelList[[i]]<-wilcox.test(fmla, data = groups, paired=F, exact=F, correct=T)
}

#print all the results for each PiSCES to the screen
#if you're a super masochist (not recommended)

modelList

#or just return a vector of the pvals (recommended!)

pvals <- unlist(lapply(modelList, '[[',3))

#save all the the pvals to a .csv file

write.csv(pvals, "pvals_forGroups_3_4.csv")
