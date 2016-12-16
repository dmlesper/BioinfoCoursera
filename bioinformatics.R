library(Biostrings)
genome <-DNAString("TCAACTGTTATCGTCCGGATCGTCCGGGAAGGGGCATCGTCCGGATCGTCCGGCACGATCGCCACGATCGCGAAGGGGCCACGATCGCCCCTTTACAGAAGGGGCTCAACTGTTCACGATCGCTCAACTGTTTCAACTGTTTCAACTGTTCCCTTTACACACGATCGCATCGTCCGGCACGATCGCCACGATCGCGAAGGGGCATCGTCCGGCACGATCGCGAAGGGGCCACGATCGCGAAGGGGCATCGTCCGGTCAACTGTTGAAGGGGCGAAGGGGCATCGTCCGGATCGTCCGGATCGTCCGGGAAGGGGCCCCTTTACATCAACTGTTATCGTCCGGGAAGGGGCATCGTCCGGCACGATCGCCACGATCGCCCCTTTACACCCTTTACATCAACTGTTATCGTCCGGGAAGGGGCCACGATCGCCCCTTTACACCCTTTACAGAAGGGGCATCGTCCGGATCGTCCGGATCGTCCGGGAAGGGGCATCGTCCGGCACGATCGCCCCTTTACAATCGTCCGGGAAGGGGCATCGTCCGGCACGATCGCTCAACTGTTCCCTTTACAATCGTCCGGTCAACTGTTCACGATCGCGAAGGGGCCACGATCGCGAAGGGGCATCGTCCGGCACGATCGCTCAACTGTTCACGATCGCATCGTCCGGCACGATCGCGAAGGGGCCCCTTTACACCCTTTACAATCGTCCGGGAAGGGGCATCGTCCGGCACGATCGCTCAACTGTTCCCTTTACAGAAGGGGCCACGATCGCGAAGGGGCATCGTCCGGCACGATCGCGAAGGGGCCCCTTTACACACGATCGCGAAGGGGCCCCTTTACACCCTTTACATCAACTGTTTCAACTGTTGAAGGGGCCCCTTTACATCAACTGTTCACGATCGCCCCTTTACAGAAGGGGCTCAACTGTT")

#get most frequent 12-mer          
genome11<- oligonucleotideFrequency(genome,11)
genome11[genome11 == min(genome11)]
b<- genome11[genome11 == max(genome11)] 



#obtain reverse completement
install.packages('ape')
library(stringr)
getwd()
new <- readLines(file.choose())
genome <- DNAString(new) 
x<- reverseComplement(genome) #finds reverse DNA completement to genome
as.character(x)
summary(x)
save(x,file="answer.txt")
mat<-matchPattern("GAATAGAGA",new) #find the pattern matches going forward
tam <- matchPattern("AGAGATAAG",new) #find the pattern matches going backward
tam

#pattern matching problem
#set variables
genomeTxt <- readLines('temp.txt')
pattern <- 'CTTGATCAT'
pattern2 <- reverse(pattern)
forward <-matchPattern(pattern,genomeTxt)
reverse<- matchPattern(pattern2,genomeTxt)
answer <-forward@ranges@start
answer2 <- reverse@ranges@start
answer
answer2
save(answer,file="care.txt")
ans<- answer -1
answer2 -1 
t<- c(answer, answer2)
as.vector(rbind(answer,answer2))
k<- sort(t)
correct <- k-1
correct
correct
k
t
forward
ans
l <- 12
x[0:l]


##Clump Finding Problem
library(Biostrings)

#load in the genome and make it a DNA string
Genome <- readLines('clumpfindingproblem.txt')
Genome <- DNAString(Genome) 

#find clumps within the genome
clumpfinding <- function(Genome,k,L,t){
        clumpvector <- c(0) 
        End <- length(Genome)-L
        #read a substring of L and find k-mers within it
        #ensure k-mers appear t times
        for (i in 1:length(Genome)){
                if (i <= End){
                        sub <- substring(Genome,i,L)
                        clumps <- oligonucleotideFrequency(sub,k)
                        clumps<- clumps[clumps==t]
                        if (length(clumps)>0){
                                #don't use duplicates
                                if (new!=(names(clumps))){
                                        new <-names(clumps)
                                        clumpvector <- c(clumpvector,new)
                                        }
                                }
                        L<- L+1
                }
        }
        print(clumpvector)
}

clumpfinding(Genome,12,572,17)

