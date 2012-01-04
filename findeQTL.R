# Function to find the eQTLs from an scanone object from R/QTL
# Returns an index of the markers of that meet a lod threshold for all
# phenotypes that have at least one QTL that meets threshold, the lodscore and
# the phenotype columns (lodcolumn)

# Written by Jeffrey Hsu 2010
require(qtl)

findeQTL <- function(scanone, genotype_file, threshold = 2.7, cisinterval =
20000000, type=c("all","chromosome","genome")){
    lodcolumn <- c()
    qtl <- c()
    lodscore <- c()
    marker_chromosome <- as.vector(genotype_file[,2])
    marker_chromosome <- ifelse(marker_chromosome =="X", 20, marker_chromosome)
    marker_chromosome <- ifelse(marker_chromosome == "Y", 21, marker_chromosome)
    print(paste("Marker length:",length(marker_chromosome)))
    marker_chromosome <- as.numeric(marker_chromosome)
    marker_chromosome <- as.factor(marker_chromosome)
    offset <- as.vector(cumsum(tapply(marker_chromosome, marker_chromosome, length)))
    offset <- append(offset,0,after=0)
    # :TODO Change to the number of chromosomes
    offset <- offset[1:20]
    tindex <- 1
    # Loads function to get the Max at each chromosome
    source('chrMax.R')
    for(lod in scanone[,3:length(scanone)]){
	    if(any(lod>=threshold)){
        if(type=="all"){
		      hits <- which(lod>threshold)
	      }else if (type=="chromosome"){
		      temp <- tapply(lod, marker_chromosome, chrMax, 
                         threshold = threshold)
		      hits <- temp + offset
		      hits <- as.vector(hits)
		      hits <- hits[!is.na(hits)]
	      }else if (type=="genome"){
		      hits <- which(lod==max(lod))
	      }
	    tindex_new <-
	    tindex*(sequence(length(hits))/sequence(length(hits)))
	    lodcolumn <- append(lodcolumn, tindex_new)
	    qtl <- append(qtl,hits)
	    lodscore <- append(lodscore, lod[hits])
	    }
	    tindex <- tindex + 1   
    }
    return(list(qtl_index = qtl, lodcolumns = lodcolumn, lodscores = lodscore))
}
