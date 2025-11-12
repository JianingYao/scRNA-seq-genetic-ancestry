param <- commandArgs(trailingOnly=T)

file  = eval(paste(text=param[1])) 
#file = "gbr.scSNPs"
error_rate  = as.numeric(eval(paste(text=param[2])))
# error_rate = 1-0.9222302

library(digest)

# Function to convert a string to a numeric hash suitable for set.seed
string_to_number <- function(input_string) {
    hash_value <- digest(input_string, algo = "md5", serialize = FALSE)
    int_value <- as.numeric(paste0("0x", substr(hash_value, 1, 8)))
    int_value <- as.integer(int_value %% .Machine$integer.max)
    return(int_value)
}

data = read.table(paste0(file,".tmp2.tped"))
snps = read.table(paste0(file,".tmp1.bim"))[,5:6]

for(i in 1:nrow(data)){
    IID = 5:ncol(data)
    torm = which(data[i,-(1:4)]==0)
    if(length(torm)>0){IID = IID[-torm]}
    n = length(IID)/2 #n is the number of genotypes
    set.seed(string_to_number(file)+i)
    error = IID[which(rbinom(n,1,error_rate)==1)*2]
    for (j in error){
        if (data[i,j]==snps[i,1]){
            data[i,j]=snps[i,2]
        } else {
            data[i,j]=snps[i,1]
        }
    }
}

write.table(data,file=paste0(file,".tmp2.tped"),quote=F,sep="\t",col.names=F,row.names=F)
