source('findeQTL.R')

################################################################
# Takes a scaone object and generates an eQTL table
################################################################

generateTable <- function(scanone, gf_file){
  gf <- read.table(gf_file, sep=",", header=T)
  qtl.temp <- findeQTL(scanone, gf, type="chromosome")
  lodcolumn <- qtl.temp$lodcolumns
  lodscore <- qtl.temp$lodscores
  qtlillmnresults <- matrix(colnames(scanone)[lodcolumn+2], ncol=1, byrow=TRUE)
  colnames(qtlillmnresults) <- c("Probe_ID")
  qtlillmnresults <- as.table(qtlillmnresults)
  qtlillmnresults<-cbind(qtlillmnresults, lodcolumn)
  qtl <- qtl.temp$qtl_index
  qtl_index <- qtl
  qtlillmnresults<-cbind(qtlillmnresults, qtl_index)
  ###############################
  # Handle the Marker Information
  ###############################
  qtl_chromosomes <- as.vector(gf[,2])[qtl]
  qtl_marker_pos <- gf[,3][qtl]
  qtl_names <- rownames(ecc.em)[qtl]
  
  qtlillmnresults<-cbind(qtlillmnresults, qtl_names)
  qtlillmnresults<-cbind(qtlillmnresults, lodscore)
   
  qtlillmnresults<-cbind(qtlillmnresults, qtl_chromosomes)
  qtlillmnresults<-cbind(qtlillmnresults, qtl_marker_pos)
  


  return(qtlillmnresults)
}

#########################################################
# Adds annotation information to reusult.  Right now 
# input is a specific file.  Need to handle other cases. 
#########################################################

annotate_eQTL <- function(result_table, annotation_table, qtl_interval = 20000000){
  
  probe_annot <- read.table(annotation_table, header=T, sep=",", fill=T, stringsAsFactors=F)
  index <- match(as.vector(result_table[,"Probe_ID"]), probe_annot[,"Probe_Id"])
  
  probe_chr <- probe_annot[index, "Chromosome"]
  probe_loc <- probe_annot[index, "Probe_Coordinates"]
  probe_loc <- sub("\\-.*", "", probe_loc)
  probe_loc <- as.numeric(probe_loc)
  probe_loc <- ifelse(probe_loc=="", NA, probe_loc)
  probe_chr <- ifelse(probe_chr=="", NA, probe_chr)
  probe_chr <- ifelse(probe_chr=="X", 20, probe_chr)
  probe_chr <- ifelse(probe_chr=="Y", 21, probe_chr)
  probe_chr <- ifelse(probe_chr=="Un", "", probe_chr)
  probe_chr <- ifelse(probe_chr=="", NA, probe_chr)
  ####################################################
  # Remove |NT_XXXXXXX identifies in annotation file
  ####################################################
  probe_chr <- sub("\\|.*", "", probe_chr)
  probe_chr <- as.numeric(probe_chr)
  
  gene_symbol <- as.character(probe_annot[index, "Symbol"])
  
  result_table <- cbind(result_table, gene_symbol)
  result_table <- cbind(result_table, probe_chr)
  result_table <- cbind(result_table, probe_loc)
  
  #######################################################
  # Calculates whether eQTL is cis or trans
  #
  #######################################################
  qtl_chr <- result_table[,"qtl_chromosomes"]
  qtl_chr <- ifelse(qtl_chr=='X', 20, qtl_chr)
  qtl_chr <- ifelse(qtl_chr=='Y', 21, qtl_chr)
  qtl_chromosomes <- as.numeric(qtl_chr)
  qtl_marker_pos <- as.numeric(result_table[,"qtl_marker_pos"])
  qtlpos <- (qtl_chromosomes-1) * 200000000 + qtl_marker_pos
  probepos <- (probe_chr - 1) * 200000000 + probe_loc
  cisstatus <- ifelse(abs(qtlpos-probepos) <= qtl_interval, "cis", "trans")
  
  result_table <- cbind(result_table, cisstatus)
  return(result_table)
}

  
#########################################################
#
# Annotate eQTL table with whether the probe contains a SNP 
# or not.  
#
#########################################################
  
SNPs_in_probe <- function(result_table, SNPS_in_probes){
  
}
  
  
##########################################################
# Plot qtls.  
#
##########################################################
  
plot_QTLs <- function(result_table, chrom="ALL"){
  # Needs to be run after annotate_qtl
  require(ggplot2)
  
  if (chrom=="ALL"){
    qtl_chr <- result_table[,"qtl_chromosomes"]
    qtl_chr <- ifelse(qtl_chr=='X', 20, qtl_chr)
    qtl_chr <- ifelse(qtl_chr=='Y', 21, qtl_chr)
    chroms <- c()
    qtl_chromosomes <- as.numeric(qtl_chr)
    qtl_marker_pos <- as.numeric(result_table[,"qtl_marker_pos"])
    qtl_endpoints <- cumsum(tapply(qtl_marker_pos, factor(qtl_chromosomes), max))
    qbreaks <- qtl_endpoints
    qtl_endpoints <-append(0, qtl_endpoints)
    qtl_endpoints <- qtl_endpoints[-length(qtl_endpoints)]
    plot_lines <- qtl_endpoints
    qbreaks <- (qbreaks-qtl_endpoints)/2
    qbreaks <- qbreaks + qtl_endpoints
    print(qbreaks)
    temp <- c()
    for (i in seq(1, length(qbreaks))){
      chroms <- append(chroms, paste('chr', i))
    }
    for (i in qtl_chromosomes){
    temp <- append(temp, qtl_endpoints[[i]])
    }
    qtlpos <-  temp + qtl_marker_pos
    probe_chr <- as.numeric(result_table[,"probe_chr"])
    probe_loc <- as.numeric(result_table[,"probe_loc"])
    probe_endpoints <- cumsum(tapply(probe_loc, factor(probe_chr), max))
    probe_endpoints <- probe_endpoints - probe_endpoints[[1]]
    temp <- c()
    for (i in probe_chr){
    if (is.na(i)){
      temp <- append(temp,"NA")
    }else{
      temp <- append(temp, probe_endpoints[[i]])
    }
    temp <- as.numeric(temp)
    }
    probepos <-  temp + probe_loc
    temp <- data.frame(Probe_pos=probepos, QTL_pos = qtlpos, LOD = result_table[,"lodscore"])
    angle = 45

  } else {
    chr.table <- result_table[result_table[,"qtl_chromosomes"] == chrom,]
    qtlpos <- as.numeric(chr.table[,"qtl_marker_pos"])
    #qbreaks <- (max(qtlpos)- min(qtlpos))/2
    qbreaks <- as.numeric(levels(as.factor(qtlpos)))
    plot_lines <- c()
    chroms = (levels(as.factor(qtlpos)))
    probe_chr <- as.numeric(chr.table[,"probe_chr"])
    probe_loc <- as.numeric(chr.table[,"probe_loc"])
    probe_endpoints <- cumsum(tapply(probe_loc, factor(probe_chr), max))
    probe_endpoints <- probe_endpoints - probe_endpoints[[1]]
    temp <- c()
    for (i in probe_chr){
    if (is.na(i)){
      temp <- append(temp,"NA")
    }else{
      temp <- append(temp, probe_endpoints[[i]])
    }
    temp <- as.numeric(temp)
    }
      probepos <-  temp + probe_loc
      temp <- data.frame(Probe_pos=probepos, QTL_pos = qtlpos, LOD = chr.table[,"lodscore"])
      angle = 90
  }
  g <- ggplot(temp, aes(QTL_pos, Probe_pos)) + geom_point()
  g <- g + geom_vline(xintercept=plot_lines) + 
    scale_x_continuous('QTL Positions', breaks=qbreaks, labels=chroms) + 
    scale_y_continuous('Probe Locations') + 
    opts(axis.text.y = theme_blank(),axis.text.x=theme_text(angle=angle, hjust=0.6), 
         title="ECC eQTL Hotspots")
  g
}

####################################################################################
# Compares 2 tables
#
####################################################################################
compareResults <- function(result.table, result.table2){
  
}

####################################################################################
# Convienence Function for grabbing eQTLs that fall within an interval
#
####################################################################################
eQTLs_intervals <- function(result.table, chr, start, end, type="None"){
  result.table <- result.table[result.table[,"qtl_chromosomes"] == chr,]
  result.table <- result.table[as.numeric(result.table[,"qtl_marker_pos"]) >= start,]
  result.table <- result.table[as.numeric(result.table[,"qtl_marker_pos"]) <= end,]
  
  if(type=="None"){
  } else {
    result.table <- result.table[result.table[,"cisstatus"]==type,]
  }
  return(result.table)
  
}
  
plot_number_eQTLs <- function(result.table){
  results <- c()
  for(i in seq(2,10,by=0.1)){
    results <- append(results, sum(as.numeric(result.table[,"lodscore"])>i))
  }
  t <- data.frame(threshold=seq(2,10,by=0.1), eQTLs=results)
  g <- ggplot(t, aes(threshold, eQTLs)) + geom_point()
  g
}

####################################################################################
# Convienence Function for running permutation tests for individual phenotypes.  Uses haley-knott
# instead of em.  Returns a new result table with the X% percentile genome-wide threshold
# These calculate the gene-specific threshold in order to deal with the correlation structure.
####################################################################################

run_perms <- function(cross, result.table, n.perm=1000){
  index <- match(as.vector(result.table[,1]), colnames(cross$pheno))
  nperall <- scanone(cross, method="hk", pheno.col=index, n.perm=1000, n.clust=8)
  n.quantile <- apply(nperall, MARGIN=1, FUN=quantile, probs=c(.95))
}
  