library(Biostrings)
###MOD 1
##get most frequent 12-mer 
genome <-DNAString("TCAACTGTTATCGTCCGGATCGTCCGGGAAGGGGCATCGTCCGGATCGTCCGGCACGATCGCCACGATCGCGAAGGGGCCACGATCGCCCCTTTACAGAAGGGGCTCAACTGTTCACGATCGCTCAACTGTTTCAACTGTTTCAACTGTTCCCTTTACACACGATCGCATCGTCCGGCACGATCGCCACGATCGCGAAGGGGCATCGTCCGGCACGATCGCGAAGGGGCCACGATCGCGAAGGGGCATCGTCCGGTCAACTGTTGAAGGGGCGAAGGGGCATCGTCCGGATCGTCCGGATCGTCCGGGAAGGGGCCCCTTTACATCAACTGTTATCGTCCGGGAAGGGGCATCGTCCGGCACGATCGCCACGATCGCCCCTTTACACCCTTTACATCAACTGTTATCGTCCGGGAAGGGGCCACGATCGCCCCTTTACACCCTTTACAGAAGGGGCATCGTCCGGATCGTCCGGATCGTCCGGGAAGGGGCATCGTCCGGCACGATCGCCCCTTTACAATCGTCCGGGAAGGGGCATCGTCCGGCACGATCGCTCAACTGTTCCCTTTACAATCGTCCGGTCAACTGTTCACGATCGCGAAGGGGCCACGATCGCGAAGGGGCATCGTCCGGCACGATCGCTCAACTGTTCACGATCGCATCGTCCGGCACGATCGCGAAGGGGCCCCTTTACACCCTTTACAATCGTCCGGGAAGGGGCATCGTCCGGCACGATCGCTCAACTGTTCCCTTTACAGAAGGGGCCACGATCGCGAAGGGGCATCGTCCGGCACGATCGCGAAGGGGCCCCTTTACACACGATCGCGAAGGGGCCCCTTTACACCCTTTACATCAACTGTTTCAACTGTTGAAGGGGCCCCTTTACATCAACTGTTCACGATCGCCCCTTTACAGAAGGGGCTCAACTGTT")
genome11<- oligonucleotideFrequency(genome,11)
genome11[genome11 == min(genome11)]
b<- genome11[genome11 == max(genome11)] 



##obtain reverse completement
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

##pattern matching problem
#set variables
genomeTxt <- readLines('temp.txt')
pattern <- 'CTTGATCAT'
pattern2 <- reverse(pattern)
forward <-matchPattern(pattern,genomeTxt) #find forward pattern match
reverse<- matchPattern(pattern2,genomeTxt) #find reverse pattern match
answer <-forward@ranges@start #how to subset when using matchPattern
answer2 <- reverse@ranges@start
answer
answer2
save(answer,file="care.txt")
#vectorize the answers for easier insert into Coursera
t<- c(answer, answer2)
as.vector(rbind(answer,answer2))
#answer and answer 2 needed to be sorted for clear answer
k<- sort(t)
correct <- k-1 #coursera's indexing requires a -1 to all values



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


###Mod 2
##finding ORI
#didn't create a fucntion as this, may have been easier
library(Biostrings)
genome <- readLines('genome.txt')
fin <- 0
vector <- c(0)
genelength <- nchar(genome)#how many characters exist
for (i in 1:genelength) {
        #is the ith spot = to c?
        if(substr(genome,i,i)=='C'){
               #c gives -1 value
                fin<- fin-1;
                vector <- c(vector,fin)
                
        } else if (substr(genome,i,i) == 'G') {
                #g gives +1 value
                fin <- fin +1;
                vector <- c(vector,fin)
                
        } else {
                #A and T don't provide values
                fin <- fin+0
                vector <- c(vector,fin)
        }
        
}
print(vector)
#again to provide a clear answer to Coursera
which(vector==min(vector))-1 

##Hamming Distance
library(Biostrings)
genes <- readLines("HammingDistanceProblem.txt")
Hammingdistance <- function(gene){
        gene1 <- gene[1] #read first string
        gene2 <- gene[2] #read second string
        gene1<- DNAString(gene1) #makes for easy length read
        gene2 <- DNAString(gene2) #allows for comparison
        HammingDistance <- 0
        for (i in 1:length(gene1)){
                if (substring(gene1,i,i) != substring(gene2,i,i)){ #check for mismatch
                        HammingDistance <- HammingDistance+1 #mismatch gives distance +1
                }
        }
        HammingDistance
}
Hammingdistance(genes)

##Approx Pattern Matching Problem
#very sloppy and non-function based
string <- readLines('approxpattern.txt')
string <- DNAString(string)
pattern <- 'TGTTAGGCCGTG'
k<- matchPattern(pattern,string,max.mismatch=5,fixed=FALSE) 
t<- as.vector(k@ranges@start,mode='any') #the @ranges@start is how you subset for start spot
t <- t-1 #Coursera uses different numbering 
cat(t) #remove inherited indexing to provide answer to Coursera 

##Approx Pattern Count
#Count the approximate pattern matches in a string
library(Biostrings)
string <- DNAString(readLines('aproxcount.txt'))
pattern <- 'CGATTTT'
PatMatch <- matchPattern(pattern,string,max.mismatch= 2,fixed=FALSE)
MatchVec <- as.vector(PatMatch@ranges@start,mode='any')
length(MatchVec)
length(PatMatch)

##Frequent Words with Mismatches
string <- DNAString('ACGTTGCATGTCGCATGATGCATGAGAGCT')
string <- DNAString(readLines('freqwordmismatch.txt'))
FWM <- function(string,k,d){
        kmer <- oligonucleotideFrequency(string,k)
        vector <- c()
        old <- 0
        for(i in 1:length(kmer)){
                kname <- names(kmer) #obtain all possible k-mers in existence 
                pattern <-kname[i]  #subset based on the ith name in the length of kmers
                PatMatch <- matchPattern(pattern,string,max.mismatch =d,fixed=FALSE) #find the number of pattern matches
                
                #if the same number of matches exist add the name to the vector
                if (length(PatMatch)==old){
                         vector <- c(vector,pattern)
                }
                #if more matches exist, replace vector with name of pattern        
                else if (length(PatMatch)>old){
                        vector <- c(pattern)
                        old <- length(PatMatch)
                } 
                #at the end print vector
                if (i==length(kmer)){
                        print(vector)
                }
        }
}
FWM(string,4,1)

##Frequent Word Mismatch and Rev Comp problem
string <- DNAString(readLines('fwmrcp.txt'))
FWMRCP <- function(string,k,d){
        kmer <- oligonucleotideFrequency(string,k)
        vector <- c()
        old <- 0
        for(i in 1:length(kmer)){
                kname <- names(kmer) #obtain all possible k-mers in existence 
                pattern <-kname[i]  #subset based on the ith name in the length of kmers
                #find the number of pattern matches
                PatMatch <- matchPattern(pattern,string,max.mismatch =d,fixed=FALSE)
                #find the reverse completement pattern and see how many matches exist in string
                RevPat <- as.character(reverseComplement(DNAString(pattern))) 
                RevPatMatch<- matchPattern(RevPat,string,max.mismatch =d,fixed=FALSE)
                
                together <- length(PatMatch)+ length(RevPatMatch)
                #if the same number of matches exist add the name to the vector
                if (together==old){
                        #see if RevPat and pattern are the same (take out any duplicates)
                        if (RevPat==pattern){ 
                                vector <- c(vector,pattern)
                        } else if (RevPat != pattern){
                                vector <- c(vector,pattern,RevPat)
                        }
                        #if more matches exist, replace vector with name of pattern        
                } else if (together>old){
                        if (RevPat==pattern){ 
                                vector <- c(pattern)
                        } else if (RevPat != pattern){
                                vector <- c(pattern,RevPat)
                        }
                        old <- together
                } 
                #at the end print vector
                if (i==length(kmer)){
                        print(vector)
                }
        }
}
FWMRCP(string,7,2)
