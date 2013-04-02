###
## Diagram triplex visualization
##
## Author: Kamil Rajdl, Jiri Hon
## Date: 2012/10/28
## Package: triplex
##


###
## This function computes necessary information to draw a triplex using one of
## function for triplex visualization.
## 
## triplex TriplexViews object
##
## RETURN: list containing the following components:
## DNAseq1  sequence of nucletides of the main (bent) thread of triplex
## DNAseq2  sequence of nucletides of the third thread of triplex
## Wbonds   matrix with numbers of nucleotides conected with W-C bonds, the 
##            matrix has two rows with corresponding numbers of nucleotides 
## Hbonds   matrix with numbers of nucleotides conected with Hoogsteen bonds, 
##            the matrix has two rows with corresponding numbers of nucleotides 
## bends    vector with informations about bending (sequence number of 
##            nucleotide beyond which starts bending and sequence number of 
##            nucleotide before which ends bending)
##
triplex.input <- function(view)
{
	triplex <- triplex.align(view)
	triplex.type <- type(view)
	
	m = nchar(triplex);
	triplex.vec <- toupper(substring(triplex, 1:m, 1:m));
	DNAseq.vec <- toupper(triplex.vec[-which((triplex.vec == "=")|(triplex.vec == "-"))]);
	DNAseq1 <- paste(DNAseq.vec, collapse = '')
	n <- length(DNAseq.vec);
	bends <- which(triplex.vec == "=");   
	gaps <- which(triplex.vec == "-"); 
	
	# Input conditions ====================
	
	# Input string conditions
	if (!all(triplex.vec %in% c("A","T","C","G","=","-"))) {
		stop("Input string contains illegal characters.");
	}
	if (length(bends) != 2) {
		stop("Wrong number of '=' characters (exactly two are needed).");
	}
	if ((bends[1] == 1) || (bends[2] == m)) {
		stop("Input string must not start or end with '=' character.");
	}
	if (length(triplex.vec[1:(bends[1]-1)]) != length(triplex.vec[(bends[2]+1):m])) {
		l.part <- length(triplex.vec[1:(bends[1]-1)]);
		r.part <- length(triplex.vec[(bends[2]+1):m]);
		stop(paste(c("Left and right part of input string are not the same length. (L=",l.part,",R=",r.part,")"), collapse = ""));
	}
	for (i in gaps) {
		if (i %in% c(1,m)) {
			stop("Input string must not start or end with '-' character.");
		}
		if ((i > bends[1]) && (i < bends[2])) {
			stop("Character '-' must not be between '=' characters.");
		}  
	}
	
	# triplex.type condition
	if (!(triplex.type %in% c(0,1,2,3,4,5,6,7))) {
		stop("wrong triplex.type");
	}
	
	# ======================================
	
	# creation of matrix of "pre-bonds"
	bonds <- numeric(0);
	ind1 = 1;
	ind2 = n;
	for (i in 1:(bends[1]-1)) {
		if ((triplex.vec[i] != "-") && (triplex.vec[m+1-i]) != "-") {
			bonds <- cbind(bonds, c(ind1, ind2));
		}
		if (triplex.vec[i] != "-") {
			ind1 <- ind1+1;
		} 
		if (triplex.vec[m+1-i] != "-") {
			ind2 <- ind2-1;
		}   
	} 
	# creation of matrices of W-C and Hoogsteen bonds and of the third thread
	Wbonds <- numeric(0);
	Hbonds <- numeric(0);
	DNAseq2 <- "";
	bonds.count <- dim(bonds)[2];
	
	if (triplex.type %in% c(0,3)) {
		for (i in 1:bonds.count) {
			Wbonds <- cbind(Wbonds, c(bonds[1,i], bonds.count-i+1));
			Hbonds <- cbind(Hbonds, c(bonds[2,i], bonds.count-i+1));
			DNAseq2 <- paste(complementary(DNAseq.vec[bonds[1,i]]), DNAseq2, sep="");   
		}
	}
	
	if (triplex.type %in% c(1,2)) {
		for (i in 1:bonds.count) {
			Wbonds <- cbind(Wbonds, c(bonds[2,i], i));
			Hbonds <- cbind(Hbonds, c(bonds[1,i], i));
			DNAseq2 <- paste(DNAseq2, complementary(DNAseq.vec[bonds[2,i]]), sep="");   
		}
	}
	
	if (triplex.type %in% c(4,7)) {
	Hbonds <- bonds;
		for (i in 1:bonds.count) {
			Wbonds <- cbind(Wbonds, c(bonds[2,i], i));
			DNAseq2 <- paste(DNAseq2, complementary(DNAseq.vec[bonds[2,i]]), sep="");   
		}
	}
	
	if (triplex.type %in% c(5,6)) {
	Hbonds <- bonds;
		for (i in 1:bonds.count) {
			Wbonds <- cbind(Wbonds, c(bonds[1,i], bonds.count-i+1));
			DNAseq2 <- paste(complementary(DNAseq.vec[bonds[1,i]]), DNAseq2, sep="");     
		}
	}
	
	# bends positions in DNAseq1
	gn <- length(gaps[gaps<bends[1]]);
	bends <- c(bends[1]-1-gn, bends[2]-1-gn);
	npl <- bends[1] - max(bonds[1,]);
	npr <- min(bonds[2,]) - bends[2];
	
	bends <- c(bends[1]-npl, bends[2]+npr);
	
	result <- list(DNAseq1, DNAseq2, Wbonds, Hbonds, bends, alignment.DNAStringSet(triplex, triplex.type));
	return(result);
}


###
## Get complementary nuclotide
##
complementary <- function(N)
{
	comp <- "A";

	if(N == "A") {
		comp <- "T";
	}
	if(N == "T") {
		comp <- "A";
	}
	if(N == "G") {
		comp <- "C";
	}
	if(N == "C") {
		comp <- "G";
	}
	return(comp);
}
