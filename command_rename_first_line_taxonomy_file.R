#!/bin/R

args = commandArgs(trailingOnly=TRUE)

print("changing tax table header...")

file <- args[1]
percentage <- args[2]


if ( !file.exists( file  ) ){
 print( paste( file, "not found by R script"))
 quit()}

#print(file)


dir_list <- unlist(strsplit( file, "/" ))
current_dir <- paste( dir_list[ 1: length(dir_list)-1], collapse="/" )

print("writing new table to to")
print( current_dir )
setwd( current_dir )


#ignore original col names
tab_ <- suppressWarnings(read.table( file , skip=1, comment.char="", sep="\t")) #, row.names=FALSE)
# tab_ <- read.table( file , skip=1, comment.char="", sep="\t") #, row.names=FALSE)

#replace fist line of taxonomy.tsv with "#OTUID       taxonomy        confidence"
colnames(tab_) <- c( "#OTUID", "taxonomy", "confidence")

#print(head(tab_))
 
write.table( tab_, file= paste(file, percentage ,"_renamed.tsv", sep=""), quote =FALSE, sep = "\t", col.names = TRUE, row.names=FALSE)

print("table header changed succesfully")


