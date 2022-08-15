# FILENAME: server1.R ##this files was created to impute genotype files from a large dataset and we did 20 families at the time to reduce the computation time
###### installing required packages and loading them
install.packages('readxl',repos = "http://ftp.ussg.iu.edu/CRAN/")
library(readxl)
install.packages('plyr',repos = "http://ftp.ussg.iu.edu/CRAN/")
library(plyr)
install.packages('dplyr',repos = "http://ftp.ussg.iu.edu/CRAN/")
library(dplyr)
install.packages('tidyverse',repos = "http://ftp.ussg.iu.edu/CRAN/")
library(tidyverse)
install.packages('parallel',repos = "http://ftp.ussg.iu.edu/CRAN/")
library(parallel)
install.packages('MASS',repos = "http://ftp.ussg.iu.edu/CRAN/")
library(MASS)


###Data Import######################
geno_path <- "/GenoCL1forImp.csv" ##genotype file to be imputation
minfo_path <- "/Markers_infoC1.csv" ##marker positions in cM
final_path <- "c_estimationsCL1.csv" ##recombination between markers

##genotype file
Geno <- read.csv(geno_path)
Geno<- Geno[,-1]
car<- Geno[,c(1,2)]
num<- seq(3:2913)
Geno[,num] <- lapply(Geno[,num],as.numeric)
Geno<- Geno[,-c(1,2)]
Geno<- cbind(car,Geno)

###filename to identify by families
Geno$filename<- as.factor(as.character(Geno$filename)) ##


##file with the marker information, position and chromosome
MInfo <- read.csv(minfo_path)
MInfo<-MInfo[,-1]
##file that has the recombination information
final_r<- read.csv(final_path)
final_r<-final_r[,-1]



###creating the list of elements to run lapply, I subset them by families by picking them using the information below, 
families0<- unique(Geno$filename)
families<- families0[1:40]


families ## I did run the first 40 families.  It's important thing to check that the total of families you are running match with the number of cores request at the en of the script in lmapply and in the submison file 
#families<-droplevels(families)

##creating the list of markers
chr<- unique(MInfo$Chromosome)

funct<- function(fam) { ##this loop is going through the families
  
  ###obtaining the subset of markers by family
  mrks1<- Geno[which(Geno$filename == fam), ]
  Line<-mrks1[,2]
  filename<-mrks1[,1]
  impute<-matrix(ncol =0, nrow = nrow(mrks1))
  for (m in chr){
    markers<- MInfo[MInfo$Chromosome== m,1]  ##extracting the markers from each chromosome 
    mrks<- mrks1[,names(mrks1) %in% (markers)] ##subsetting from each family the markers of the chromosome
    imp<-mrks
    for (j in 1:ncol(mrks)){
      for (k in 1:nrow(mrks)){
        Pi<- mrks[k,j]
        if(is.na(Pi)){
          dat1<- mrks[k,]
          a<- which(!is.na(dat1))
          index<-j-a
          right<-min(a[index<0])
          left<-max(a[index>0])
          if(right!= Inf && right>0 && left>0 &&left != -Inf && length(a)>=2){
            right<- dat1[right]
            right<- colnames(right)
            left<-dat1[left]
            left<- colnames(left)
            c_r<-final_r[final_r$V1 == left & final_r$V2 == right,4]
            r1<- final_r[final_r$V1 == left & final_r$V2 == colnames(mrks)[j], 4 ]
            r1<-as.numeric(r1)
            r2<-final_r[final_r$V1 == colnames(mrks)[j] & final_r$V2 == right,4 ]
            r2<-as.numeric(r2)
            Ml<- mrks[k,left]
            MR<-mrks[k,right]
            if (Ml+MR == 2 && !is.na(r1) && !is.na(r2)){
              mprob<- (1-c_r)^2/4
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*(1-r1)^2*(1-r2)^2
                G1<- 1/2*r1*r2*(1-r1)*(1-r2)
                G0<-1/4*r1^2*r2^2
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                print("recombination of 1")}
              
            }else if (Ml+MR == -2 && !is.na(r1) && !is.na(r2))  {
              mprob<- (1-c_r)^2/4  
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r1^2*r2^2
                G1<- 1/2*r1*r2*(1-r1)*(1-r2)
                G0<-1/4*(1-r1)^2*(1-r2)^2
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if (Ml== 0 && MR == 0 && !is.na(r1) && !is.na(r2))  {
              mprob<-((1-c_r)^2 + c_r^2)/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r1*r2*(1-r1)*(1-r2)
                G1<- (1/2*(1-(2*r1)+(2*r1^2))/(1-(2*r2)+(2*r2^2)))
                G0<-1/4*r1*r2*(1-r1)*(1-r2)
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if (Ml==-1 && MR== 1 && !is.na(r1) && !is.na(r2)){
              mprob<-(c_r)^2/4
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4 * (r1)^2* (1-r2)^2
                G1<- 1/2*r1*r2*(1-r1)*(1-r2)
                G0<-1/4*(1-r1)^2*r2  
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if(Ml==1 && MR== -1 && !is.na(r1) && !is.na(r2)) {
              mprob<-(c_r)^2/4
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*(1-r1)^2*(r2)^2 
                G1<- 1/2*r1*r2*(1-r1)*(1-r2)
                G0<-1/4 * (r1)^2* (1-r2)^2 
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if (Ml== 1 && MR==0 && !is.na(r1) && !is.na(r2)){
              mprob<- (c_r*(1-c_r))/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r2*(1-r1)^2*(1-r2)
                G1<- 1/2*r1*(1-r1)*(1-(2*r2)+(2*(r2)^2))
                G0<-1/4 *(r1)^2*r2* (1-r2)^2 
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if (Ml== 0 && MR==1 && !is.na(r1) && !is.na(r2)){
              mprob<- (c_r*(1-c_r))/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r1*(1-r1)*(1-r2)^2
                G1<- 1/2*r2*(1-r2)*(1-(2*r1)+(2*(r1)))^2
                G0<-1/4 *r1* (1-r1)*(r2)^2 
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
              
            }else if (Ml== 0 && MR==-1 && !is.na(r1) && !is.na(r2)){
              mprob<- (c_r*(1-c_r))/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r1*(1-r1)*(r2)^2
                G1<- 1/2*r2*(1-r2)*(1-(2*r1)-(2*(r1)^2))
                G0<-1/4 *r1* (1-r1)*(1-r2)^2 
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
              
            }else if (Ml== -1 && MR==0 && !is.na(r1) && !is.na(r2)){
              mprob<- (c_r*(1-c_r))/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*(r1)^2*(r2)*(1-r2)
                G1<- 1/2*r1*(1-r1)*(1-(2*r2)+(2*(r2)^2))
                G0<-1/4 *r2* (1-r1)^2*(1-r2)
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else{
              print("error, no recombination available")
            }
            
          }else if (right== Inf && left>0 &&left != -Inf && length(a)>=2) {
            right<-max(a[index>0])
            a<-a[-which(a== right)]
            index<-j-a
            left<- max(a[index>0])
            right<- dat1[right]
            right<- colnames(right)
            left<-dat1[left]
            left<- colnames(left)
            c_r<-final_r[final_r$V1 == left & final_r$V2 == right,4]
            r1<- final_r[final_r$V1 == left & final_r$V2 == colnames(mrks)[j], 4 ]
            r1<-as.numeric(r1)
            r2<-final_r[final_r$V1 == colnames(mrks)[j] & final_r$V2 == right,4 ]
            r2<-as.numeric(r2)
            Ml<- mrks[k,left]
            MR<-mrks[k,right]
            if (Ml+MR == 2 && !is.na(r1) && !is.na(r2)){
              mprob<- (1-c_r)^2/4
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*(1-r1)^2*(1-r2)^2
                G1<- 1/2*r1*r2*(1-r1)*(1-r2)
                G0<-1/4*r1^2*r2^2
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"}
              
            }else if (Ml+MR == -2 && !is.na(r1) && !is.na(r2))  {
              mprob<- (1-c_r)^2/4  
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r1^2*r2^2
                G1<- 1/2*r1*r2*(1-r1)*(1-r2)
                G0<-1/4*(1-r1)^2*(1-r2)^2
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if (Ml== 0 && MR == 0 && !is.na(r1) && !is.na(r2))  {
              mprob<-((1-c_r)^2 + c_r^2)/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r1*r2*(1-r1)*(1-r2)
                G1<- (1/2*(1-(2*r1)+(2*r1^2))/(1-(2*r2)+(2*r2^2)))
                G0<-1/4*r1*r2*(1-r1)*(1-r2)
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if (Ml==-1 && MR== 1 && !is.na(r1) && !is.na(r2)){
              mprob<-(c_r)^2/4
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4 * (r1)^2* (1-r2)^2
                G1<- 1/2*r1*r2*(1-r1)*(1-r2)
                G0<-1/4*(1-r1)^2*r2  
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if(Ml==1 && MR== -1 && !is.na(r1) && !is.na(r2)) {
              mprob<-(c_r)^2/4
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*(1-r1)^2*(r2)^2 
                G1<- 1/2*r1*r2*(1-r1)*(1-r2)
                G0<-1/4 * (r1)^2* (1-r2)^2 
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if (Ml== 1 && MR==0 && !is.na(r1) && !is.na(r2)){
              mprob<- (c_r*(1-c_r))/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r2*(1-r1)^2*(1-r2)
                G1<- 1/2*r1*(1-r1)*(1-(2*r2)+(2*(r2)^2))
                G0<-1/4 *(r1)^2*r2* (1-r2)^2 
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if (Ml== 0 && MR==1 && !is.na(r1) && !is.na(r2)){
              mprob<- (c_r*(1-c_r))/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r1*(1-r1)*(1-r2)^2
                G1<- 1/2*r2*(1-r2)*(1-(2*r1)+(2*(r1)))^2
                G0<-1/4 *r1* (1-r1)*(r2)^2 
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
              
            }else if (Ml== 0 && MR==-1 && !is.na(r1) && !is.na(r2)){
              mprob<- (c_r*(1-c_r))/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r1*(1-r1)*(r2)^2
                G1<- 1/2*r2*(1-r2)*(1-(2*r1)-(2*(r1)^2))
                G0<-1/4 *r1* (1-r1)*(1-r2)^2 
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
              
            }else if (Ml== -1 && MR==0 && !is.na(r1) && !is.na(r2)){
              mprob<- (c_r*(1-c_r))/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*(r1)^2*(r2)*(1-r2)
                G1<- 1/2*r1*(1-r1)*(1-(2*r2)+(2*(r2)^2))
                G0<-1/4 *r2* (1-r1)^2*(1-r2)
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else{
              print("error, no recombination available")}
            
          }else if (right!= Inf && right>0 && left == -Inf && length(a)>=2) {
            a<- which(!is.na(dat1))
            left<-min(a[index<0])
            a<-a[-which(a== left)]
            index<-j-a
            right<-min(a[index<0])
            right<- dat1[right]
            right<- colnames(right)
            left<-dat1[left]
            left<- colnames(left)
            c_r<-final_r[final_r$V1 == left & final_r$V2 == right,4]
            r1<- final_r[final_r$V1 == left & final_r$V2 == colnames(mrks)[j], 4 ]
            r1<-as.numeric(r1)
            r2<-final_r[final_r$V1 == colnames(mrks)[j] & final_r$V2 == right,4 ]
            r2<-as.numeric(r2)
            Ml<- mrks[k,left]
            MR<-mrks[k,right]
            if (Ml+MR == 2 && !is.na(r1) && !is.na(r2)){
              mprob<- (1-c_r)^2/4
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*(1-r1)^2*(1-r2)^2
                G1<- 1/2*r1*r2*(1-r1)*(1-r2)
                G0<-1/4*r1^2*r2^2
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"}
              
            }else if (Ml+MR == -2 && !is.na(r1) && !is.na(r2))  {
              mprob<- (1-c_r)^2/4  
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r1^2*r2^2
                G1<- 1/2*r1*r2*(1-r1)*(1-r2)
                G0<-1/4*(1-r1)^2*(1-r2)^2
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if (Ml== 0 && MR == 0 && !is.na(r1) && !is.na(r2))  {
              mprob<-((1-c_r)^2 + c_r^2)/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r1*r2*(1-r1)*(1-r2)
                G1<- (1/2*(1-(2*r1)+(2*r1^2))/(1-(2*r2)+(2*r2^2)))
                G0<-1/4*r1*r2*(1-r1)*(1-r2)
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if (Ml==-1 && MR== 1 && !is.na(r1) && !is.na(r2)){
              mprob<-(c_r)^2/4
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4 * (r1)^2* (1-r2)^2
                G1<- 1/2*r1*r2*(1-r1)*(1-r2)
                G0<-1/4*(1-r1)^2*r2  
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if(Ml==1 && MR== -1 && !is.na(r1) && !is.na(r2)) {
              mprob<-(c_r)^2/4
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*(1-r1)^2*(r2)^2 
                G1<- 1/2*r1*r2*(1-r1)*(1-r2)
                G0<-1/4 * (r1)^2* (1-r2)^2 
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if (Ml== 1 && MR==0 && !is.na(r1) && !is.na(r2)){
              mprob<- (c_r*(1-c_r))/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r2*(1-r1)^2*(1-r2)
                G1<- 1/2*r1*(1-r1)*(1-(2*r2)+(2*(r2)^2))
                G0<-1/4 *(r1)^2*r2* (1-r2)^2 
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
            }else if (Ml== 0 && MR==1 && !is.na(r1) && !is.na(r2)){
              mprob<- (c_r*(1-c_r))/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r1*(1-r1)*(1-r2)^2
                G1<- 1/2*r2*(1-r2)*(1-(2*r1)+(2*(r1)))^2
                G0<-1/4 *r1* (1-r1)*(r2)^2 
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
              
            }else if (Ml== 0 && MR==-1 && !is.na(r1) && !is.na(r2)){
              mprob<- (c_r*(1-c_r))/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*r1*(1-r1)*(r2)^2
                G1<- 1/2*r2*(1-r2)*(1-(2*r1)-(2*(r1)^2))
                G0<-1/4 *r1* (1-r1)*(1-r2)^2 
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                "recombination of 1"
              }
              
              
            }else if (Ml== -1 && MR==0 && !is.na(r1) && !is.na(r2)){
              mprob<- (c_r*(1-c_r))/2
              if(mprob>0 && r1>0 && r2>0){
                G2<- 1/4*(r1)^2*(r2)*(1-r2)
                G1<- 1/2*r1*(1-r1)*(1-(2*r2)+(2*(r2)^2))
                G0<-1/4 *r2* (1-r1)^2*(1-r2)
                PG2<-G2/mprob
                PG1<-G1/mprob
                PG0<-G0/mprob
                a<- c(PG2,PG1, PG0)
                b<-c(22,11,00)
                prob<-which.max(a)
                imp[k,j]<- b[prob]
                print("imputed")
              } else{
                print("recombination of 1")
              }
              
            }else{
              print("no any of the cases")
            }  
          }else{
            print("error, no recombination available")}
          
          
        }else{
          print("no NA")}
        
      }
    }
    impute<-cbind(impute, imp)
  }
  impute<-cbind(Line,filename,impute)    
  
}


results<- mclapply(families, funct, mc.cores =40)

save(results, file="results.RData") ##I have to re run some families that why here is results 19 so i just repeat them in the initial scripts. 