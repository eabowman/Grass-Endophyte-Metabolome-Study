varpart2.MEM <- 
	function(Y, X1, X2, is.MEM=NULL, method="hierarchical")
#
# Hierarchical and proportional variation partitioning of Y 
# with respect to X1 and X2.
#
# Usage:
# Y : Vector, matrix or data frame of response data.
# X1, X2 : Vectors, matrices or data frames of explanatory data.
# is.MEM=NULL : No matrix of MEM eigenfunctions, or only one. That matrix 
#       will receive no special treatment since no hierarchy of spatial 
#       submodels has to be set. Classical variation partitioning is used.
# is.MEM = c(1:2) : Matrices X1 and X2 contain MEM eigenfunctions. 
#		Place the broad-scale MEM model before the fine-scale model.
# method="hierarchical" : Hierarchical partitioning.
#       ="proportional" : Proportional apportioning of the shared fractions.
#       ="classical"    : classical variation partitioning.
# Methods can be abbreviated to any recognizable subset, e.g. "h", "p", "c".
#
# Note -- Classical variation partitioning is used if method="classical",
# if is.MEM=NULL, or if there is a single file of MEM (e.g. is.MEM=2).
#
# Author:: Pierre Legendre
# License: GPL-2
{
require(vegan)
method <- match.arg(method, c("hierarchical", "proportional","classical"))
x <- varpart(Y, X1, X2)   # The heavy work is done by varpart()
#
	if(method=="classical") return(x)
	if(any(is.MEM < 1) || any(is.MEM > 2)) stop("Error in is.MEM statement")
    if(!is.null(is.MEM)) {
    	if(length(is.MEM) == 1) {
    	   cat("Partitioning with no special treatment for a single MEM file\n")
    	   return(x)
    		} else {
    		are.MEM <- TRUE
    		cat("Files X1 and X2 contain MEM eigenfunctions\n")
    		}
    	} else { return(x) }
#
if(length(is.MEM) == 2) {    # Two matrices of MEM
### method="hierarchical"
	if(method=="hierarchical") {
		cat("Hierarchical partitioning of the shared fraction\n")
		x$part$indfract[1,3] <- sum(x$part$indfract[1:2,3])
		x$part$indfract[2,3] <- NA

		} else {
### method="proportional"
		cat("Proportional apportioning of the shared fraction\n")
		a <- x$part$indfract[1,3]
		b <- x$part$indfract[2,3]
		c <- x$part$indfract[3,3]
		b.to.a <- b * a/(a+c)
		b.to.c <- b * c/(a+c)
		x$part$indfract[1,3] <- a + b.to.a
		x$part$indfract[2,3] <- NA
		x$part$indfract[3,3] <- c + b.to.c
		}
	# Check results
	if(abs(sum(x$part$indfract[c(1,3),3])-x$part$fract[3,3]) > 1e-10) 
		cat("Error in fraction calculation\n")
	return(x)
	}
}