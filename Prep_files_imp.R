#####this script makes the files needed for marker imputation####################

##loading libraries
library(tidyverse)
library(tibble)
library(dplyr)
library(NAM)
require(ggpubr)
require(tidyverse)
library(genetics)
library(NAM)


####cluster1
Geno<- read.csv('/class/datamine/corporate/bayer/students/team_3/combined_progenygenotypeCL1.csv') ##this file was obtaining using python script "coca_csv_files.py"
Geno<-Geno[,-c(1,2)]
##leaving the one column with the file name as C1.##
filename<- Geno$filename
filename<-gsub(".combinedgeno.csv","", x=filename, perl = TRUE)
#filename<-gsub("\\.","x", x=filename, perl = TRUE)

##removing unused columns
Geno<-Geno[,c(1,2913,2:2912)]
Geno<-Geno[,-2]
##joiniing filename and data
Geno<- cbind(filename,Geno)

###extracting the markers information- chromosome and centimorgans data
Mark_info<- Geno[Geno$LINE=="",] ##to check that is the same info in all files
Mark_info<- Mark_info[1:2,-c(1:3)]
Mark_info<- t(Mark_info)
Mark_info<- Mark_info[-2912,]
Mark_info<- as.data.frame(Mark_info)
colnames(Mark_info)<-c("Chromosome","Position")
Mark_info$Chromosome<- as.numeric(Mark_info$Chromosome)
Mark_info$Position<- as.numeric(Mark_info$Position)
markers<- rownames(Mark_info)
Mark_info<- cbind(markers,Mark_info)
Mark_info$markers<- as.factor(as.character(Mark_info$markers))

write.csv(Mark_info, file="Markers_infoC1.csv") ##create the file for being used later


#now cleaning the genotype file for no require lines for imputation (remove parents rows and rows with marker information)

Geno<- Geno[-(grep("^PID",Geno$LINE)),] ##remove parents rows
Geno<- Geno[-which(Geno$LINE == ""), ] ##removes rows with marker information 
Geno<- distinct(Geno)
write.csv(Geno,file='GenoCL1forImp.csv') ## In this file I already filter parents information and no needed lines. it's the file use for imputation

MInfo<-read.csv("/Markers_infoC1.csv") ##this file was created on the script Descriptive_statistics line 233
MInfo<- MInfo[,-1]

#mrks<- mrks[,-c(1,2)]
##conver from -1,0,1 to 0,1,2 by adding 1 
#mrks<- mrks + 1
#mrks <- mrks[, colSums(is.na(mrks)) == 0 ]
##conver from -1,0,1 to 0,1,2 by adding 1then run alencar LD 
#ld<- NAM::LD(mrks)

Geno$filename<- as.factor(as.character(Geno$filename))


####loop for estimating the average recombination rates among markers which is also needed for marker imputation
b<- unique(MInfo$Chromosome)
final_r<-data.frame(matrix(ncol = 4, nrow = 0))
for (i in 1:length(b)){
  dat<- MInfo[MInfo$Chromosome== i, ]
  pair<-as.data.frame(permutations(n= length(dat$markers), r=2, dat$markers))
  V1<- pair$V1
  V2<- pair$V2
  r_vec<- c()
  for (k in 1:length(V1)){
  m<- abs(dat[dat$markers== V1[[k]],3]-dat[dat$markers== V2[k],3])
   if (m == 0){
     m<- 0.1
     r= (exp(m/25)-1)/(2*(exp(m/25)+1))
     r_vec<- c(r_vec,r)
   }else{
  r= (exp(m/25)-1)/(2*(exp(m/25)+1))
  r_vec<- c(r_vec,r)
  }}
  chr<-rep(i, length(pair$V1))
  pair_dat<-cbind(pair,chr,r_vec)
  final_r<-rbind(final_r,pair_dat) }

write.csv(final_r, file="c_estimationsCL1.csv") ##file with the recombinations among marker 

####cluster2##################################################################
Geno<- read.csv('/combined_progenygenotypeCL2.csv') ##this file was obtaining using python script "coca_csv_files.py"
Geno<-Geno[,-c(1,2)]
##leaving the one column with the file name as C1.##
filename<- Geno$filename
filename<-gsub(".combinedgeno.csv","", x=filename, perl = TRUE)
#filename<-gsub("\\.","x", x=filename, perl = TRUE)

##removing unused columns
Geno<-Geno[,c(1,2913,2:2912)]
Geno<-Geno[,-2]
##joiniing filename and data
Geno<- cbind(filename,Geno)

###extracting the markers information- chromosome and centimorgans data
Mark_info<- Geno[Geno$LINE=="",] ##to check that is the same info in all files
Mark_info<- Mark_info[1:2,-c(1:3)]
Mark_info<- t(Mark_info)
Mark_info<- Mark_info[-2912,]
Mark_info<- as.data.frame(Mark_info)
colnames(Mark_info)<-c("Chromosome","Position")
Mark_info$Chromosome<- as.integer(as.character(Mark_info$Chromosome))
Mark_info$Position<- as.numeric(as.character(Mark_info$Position))
markers<- rownames(Mark_info)
Mark_info<- cbind(markers,Mark_info)
Mark_info$markers<- as.factor(as.character(Mark_info$markers))

write.csv(Mark_info, file="Markers_infoC2.csv") ##create the file for being used later


#now cleaning the genotype file for no require lines for imputation (remove parents rows and rows with marker information)

Geno<- Geno[-(grep("^PID",Geno$LINE)),] ##remove parents rows
Geno<- Geno[-which(Geno$LINE == ""), ] ##removes rows with marker information 
Geno<- distinct(Geno)
write.csv(Geno,file='GenoCL2forImp.csv') ## In this file I already filter parents information and no needed lines. it's the file use for imputation

MInfo<-read.csv("/Markers_infoC1.csv") ##this file 
#mrks<- mrks[,-c(1,2)]
##conver from -1,0,1 to 0,1,2 by adding 1 
#mrks<- mrks + 1
#mrks <- mrks[, colSums(is.na(mrks)) == 0 ]
##conver from -1,0,1 to 0,1,2 by adding 1then run alencar LD 
#ld<- NAM::LD(mrks)

Geno$filename<- as.factor(as.character(Geno$filename))


####loop for estimating the average recombination rates among markers which is also needed for marker imputation
b<- unique(MInfo$Chromosome)
final_r<-data.frame(matrix(ncol = 4, nrow = 0))
for (i in 1:length(b)){
  dat<- MInfo[MInfo$Chromosome== i, ]
  pair<-as.data.frame(permutations(n= length(dat$markers), r=2, dat$markers))
  V1<- pair$V1
  V2<- pair$V2
  r_vec<- c()
  for (k in 1:length(V1)){
    m<- abs(dat[dat$markers== V1[[k]],3]-dat[dat$markers== V2[k],3])
    if (m == 0){
      m<- 0.1
      r= (exp(m/25)-1)/(2*(exp(m/25)+1))
      r_vec<- c(r_vec,r)
    }else{
      r= (exp(m/25)-1)/(2*(exp(m/25)+1))
      r_vec<- c(r_vec,r)
    }}
  chr<-rep(i, length(pair$V1))
  pair_dat<-cbind(pair,chr,r_vec)
  final_r<-rbind(final_r,pair_dat) }

write.csv(final_r, file="c_estimationsCL2.csv") ##file with the recombinations among marker 


