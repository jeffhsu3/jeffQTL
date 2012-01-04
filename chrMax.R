chrMax <- function(x, threshold){
    #Function to determine whether the max LOD marker on a chromosome meets
    # the LOD threshold and returns the marker row number with the highest LOD
    # score of that marker relative to its chromosome position.  Need to offset this value if you want to grab
    # rows form the aboslute row position.  
    if(any(x>=threshold)){
	intres <- match(max(x),x)
	return(intres)
	# return(which(x == max(x),arr.ind=T))
    }else{
	return(NA)
    }
}
