###
## Diagram triplex visualization
##
## Author: Kamil Rajdl
## Date: 2012/10/28
## Package: triplex
##


###
## This function draws a 2-D representation of a triplex 
##
## triplex      TriplexViews object
## circles      TRUE - nucleotides are drawn as characters in circles
##              FALSE - nucleotides are drawn just as characters
## wcbonds.lty  type of Watson-Crick bonds lines
## hbonds.lty   type of Hoogsteen bonds lines
## mbonds.lty   type of main (skelet) bonds lines
## wcbonds.lwd  width of Watson-Crick bonds lines
## hbonds.lwd   width of Hoogsteen bonds lines
## mbonds.lwd   width of main (skelet) bonds lines
## labels.cex   multiplier of size of labels of nucleotides
## circles.cex  multiplier of size of nucleotides
## margin       left and right margin of the picture
##
triplex.diagram <- function(
	triplex,
	circles = TRUE,
	mbonds.lty = 1,
	mbonds.lwd = 2.5,
	wcbonds.lty = 1,
	wcbonds.lwd = 1,
	hbonds.lty = 2,
	hbonds.lwd = 1,
	labels.cex = 1,
	circles.cex = 1,
	margin = 0.1,
	bonds.length = 0.07)
{
	# processing of inputs
	results <- triplex.input(triplex);
	triplex.type <- type(triplex);

	DNAseq1 <- results[[1]];
	DNAseq2 <- results[[2]];
	WCbonds <- results[[3]];
	Hbonds <- results[[4]];
	bends <- results[[5]];
	alignment <- results[[6]];

	# input condition
	if ((margin<0) || (margin>=0.5)) {
		stop("xmargin have to be from interval [0,0.5).");
	}

	# initial parameters and variables 
	bend.start <- bends[1];
	bend.end <- bends[2];       
	xsize <- 1-2*margin;
	n.bonds <- dim(WCbonds)[2];
	n.bend <- bend.end - bend.start - 1;
	n.seq1 <- nchar(DNAseq1);
	n.seq2 <- nchar(DNAseq2);
	dist2 <- bonds.length;                                # length of Watson-Crick and Hoogsteen bonds
	dist1 <- xsize/(n.bonds+(n.bend+2)^2/((n.bend+1)*pi)); # distance between nucleotides in sequence
	if (dist1 > 0.07) {
		dist1 = 0.07;
		xsize = dist1*(n.bonds+(n.bend+2)^2/((n.bend+1)*pi));
		margin = (1-xsize)/2;
		
	}

	# division of Hoogsteen bonds into DNAseq1-DNAseq1 bonds and DNAseq1-DNAseq2 bonds
	if (triplex.type %in% c("0","1","2","3")) {
		Hbonds1 <- Hbonds;
		Hbonds2 <- matrix(0,2,0);
	}
	if (triplex.type %in% c("4","5","6","7")) {
		Hbonds1 <- matrix(0,2,0);
		Hbonds2 <- Hbonds;
	}

	# processing of sequences into vector of characters
	DNAseq1.vec <- substring(DNAseq1, 1:n.seq1, 1:n.seq1);
	DNAseq2.vec <- substring(DNAseq2, 1:n.seq2, 1:n.seq2);

	# x-coordinates of paired nucleotides (the rest is set to NA for now)
	seq1.xcoords <- rep(NA, n.seq1);
	j <- 1;
	for (i in 1:bend.start) {
		if ((i %in% WCbonds[1,]) || (i %in% Hbonds1[1,]) ||  (i %in% Hbonds2)) {
			seq1.xcoords[i] <- margin+j*dist1;
			j <- j+1;   
		}
	}
	j <- 1;
	for (i in n.seq1:bend.end) {
		if ((i %in% WCbonds[1,]) || (i %in% Hbonds1[1,]) ||  (i %in% Hbonds2)) {
			seq1.xcoords[i] <- margin+j*dist1;
			j <- j+1;   
		}
	}

	# coordinates of direction labels (3', 5')
	dirlab.xcoords <- c(margin, margin, margin); 
	dirlab.ycoords <- c(0.5+dist2, 0.5, 0.5-dist2); 

	# the remaining coordinates of paired nucleotides and creation of direction labels
	if (triplex.type %in% c(0,3)) {
		seq1.ycoords <- c(rep(0.5+dist2, bend.start), rep(NA, n.bend), 
		                  rep(0.5-dist2, n.seq1-bend.end+1));    
		seq2.xcoords <- seq(from=margin+dist1*n.bonds, to=margin+dist1, by=-dist1);
		seq2.ycoords <- rep(0.5, n.bonds);

		dirlab <- textGrob(c("5'","3'","3'"), dirlab.xcoords, dirlab.ycoords);                   
	}
	if (triplex.type %in% c(1,2)) {
		seq1.ycoords <- c(rep(0.5+dist2, bend.start), rep(NA, n.bend), 
		                  rep(0.5-dist2, n.seq1-bend.end+1));  
		seq2.xcoords <- seq(from=margin+dist1, to=margin+dist1*n.bonds, by=dist1);
		seq2.ycoords <- rep(0.5, n.bonds);

		dirlab <- textGrob(c("5'","5'","3'"), dirlab.xcoords, dirlab.ycoords); 
	}
	if (triplex.type %in% c(4,7)) {
		seq1.ycoords <- c(rep(0.5-dist2, bend.start), rep(NA, n.bend), 
		                  rep(0.5, n.seq1-bend.end+1));     
		seq2.xcoords <- seq(from=margin+dist1, to=margin+dist1*n.bonds, by=dist1);
		seq2.ycoords <- rep(0.5+dist2, n.bonds);

		dirlab <- textGrob(c("5'","3'","5'"), dirlab.xcoords, dirlab.ycoords);  
	}
	if (triplex.type %in% c(5,6)) {
		seq1.ycoords <- c(rep(0.5, bend.start), rep(NA, n.bend), 
		                  rep(0.5+dist2, n.seq1-bend.end+1));    
		seq2.xcoords <- seq(from=margin+dist1*n.bonds, to=margin+dist1, by=-dist1);
		seq2.ycoords <- rep(0.5-dist2, n.bonds);

		dirlab <- textGrob(c("3'","5'","3'"), dirlab.xcoords, dirlab.ycoords);  
	}

	# coordinates of nucleotiodes in the bend
	if ((bend.end-bend.start)>1) {
		x1 <- c(seq1.xcoords[bend.start], seq1.ycoords[bend.start]);
		x2 <- c(seq1.xcoords[bend.end], seq1.ycoords[bend.end]);
		S <- c(x1[1]+dist1*(n.bend+1)/(2*pi), (x1[2]+x2[2])/2);
		y <- abs(x1[2]-x2[2])/2;
		r <- sqrt((x1[1]-S[1])^2+(x1[2]-S[2])^2);
		alpha1 <- asin(y/r);
		alpha2 <- (2*pi-2*alpha1)/(n.bend+1);
		for (i in (bend.start+1):(bend.end-1)) {
			if (triplex.type %in% c(0,1,2,3)) {
				seq1.xcoords[i] <- S[1] + r*cos(pi-alpha1-(i-bend.start)*alpha2);
				seq1.ycoords[i] <- S[2] + r*sin(pi-alpha1-(i-bend.start)*alpha2);
			}
			if (triplex.type %in% c(4,5,6,7)) {
				seq1.xcoords[i] <- S[1] + r*cos(pi+alpha1+(i-bend.start)*alpha2);
				seq1.ycoords[i] <- S[2] + r*sin(pi+alpha1+(i-bend.start)*alpha2);
			}
		}
	}

	# coordinates of nucleotides opposite the gaps
	for (i in 1:n.seq1) {
		if (is.na(seq1.xcoords[i])) {
			j<-1
			while (is.na(seq1.xcoords[i+j])) {
				j<-j+1;
			}
			x1 <- c(seq1.xcoords[i-1], seq1.ycoords[i-1]);
			x2 <- c(seq1.xcoords[i+j], seq1.ycoords[i+j]);
			y <- abs(x1[1]-x2[1])/2;
			z <- dist1*(j+1)/(2*pi);
			down = ((triplex.type %in% c(0,1,2,3)) && (i >= bend.end)) || 
			       ((triplex.type %in% c(4,7)) && (i <= bend.start));
			if (down) {
				S <- c((x1[1]+x2[1])/2, x1[2]-z);
			} else {
				S <- c((x1[1]+x2[1])/2, x1[2]+z); 
			}
			r <- sqrt((x1[1]-S[1])^2+(x1[2]-S[2])^2);
			alpha1 <- asin(y/r);
			alpha2 <- (2*pi-2*alpha1)/(j+1);
			if ((x1[1] < x2[1]) && (!down)) {
				for (k in i:(i+j-1)) {
					seq1.xcoords[k] <- S[1] + r*cos(-pi/2-alpha1-(k-i+1)*alpha2);
					seq1.ycoords[k] <- S[2] + r*sin(-pi/2-alpha1-(k-i+1)*alpha2);
				}
			}
			if ((x1[1] < x2[1]) && (down)) {
				for (k in i:(i+j-1)) {
					seq1.xcoords[k] <- S[1] + r*cos(pi/2+alpha1+(k-i+1)*alpha2);
					seq1.ycoords[k] <- S[2] + r*sin(pi/2+alpha1+(k-i+1)*alpha2);
				}
			}
			if ((x1[1] > x2[1]) && (!down)) {
				for (k in i:(i+j-1)) {
					seq1.xcoords[k] <- S[1] + r*cos(-pi/2+alpha1+(k-i+1)*alpha2);
					seq1.ycoords[k] <- S[2] + r*sin(-pi/2+alpha1+(k-i+1)*alpha2);
				}
			}
			if ((x1[1] > x2[1]) && (down)) {
				for (k in i:(i+j-1)) {
					seq1.xcoords[k] <- S[1] + r*cos(pi/2-alpha1-(k-i+1)*alpha2);
					seq1.ycoords[k] <- S[2] + r*sin(pi/2-alpha1-(k-i+1)*alpha2);
				}
			}             
		}
	}

	# creation of bonds in sequences
	bonds1 <- linesGrob(seq1.xcoords, seq1.ycoords, default.units = "npc",
	                    gp = gpar(lwd = mbonds.lwd, lty = mbonds.lty)); 
	bonds2 <- linesGrob(seq2.xcoords, seq2.ycoords, default.units = "npc",
	                    gp = gpar(lwd = mbonds.lwd, lty = mbonds.lty)); 

	# creation of Watson-Crick and Hoogsteen bonds
	WCbonds1.xcoords <- rep(NA, n.bonds);
	WCbonds1.ycoords <- rep(NA, n.bonds);
	WCbonds2.xcoords <- rep(NA, n.bonds);
	WCbonds2.ycoords <- rep(NA, n.bonds);
	for (i in 1:n.bonds) {
		WCbonds1.xcoords[i] <- seq1.xcoords[WCbonds[1,i]];
		WCbonds1.ycoords[i] <- seq1.ycoords[WCbonds[1,i]];
		WCbonds2.xcoords[i] <- seq2.xcoords[WCbonds[2,i]];
		WCbonds2.ycoords[i] <- seq2.ycoords[WCbonds[2,i]];
	} 
	bonds3 <- segmentsGrob(WCbonds1.xcoords, WCbonds1.ycoords, WCbonds2.xcoords,
	                       WCbonds2.ycoords, default.units = "npc",
	                       gp = gpar(lwd = wcbonds.lwd, lty = wcbonds.lty));

	Hbonds1.xcoords <- rep(NA, n.bonds);
	Hbonds1.ycoords <- rep(NA, n.bonds);
	Hbonds2.xcoords <- rep(NA, n.bonds);
	Hbonds2.ycoords <- rep(NA, n.bonds);
	if (dim(Hbonds1)[2] > 0) {
		for (i in 1:dim(Hbonds1)[2]) {
			Hbonds1.xcoords[i] <- seq1.xcoords[Hbonds[1,i]];
			Hbonds1.ycoords[i] <- seq1.ycoords[Hbonds[1,i]];
			Hbonds2.xcoords[i] <- seq2.xcoords[Hbonds[2,i]];
			Hbonds2.ycoords[i] <- seq2.ycoords[Hbonds[2,i]];
		} 
	}
	if (dim(Hbonds2)[2] > 0) {
		for (i in 1:dim(Hbonds2)[2]) {
			Hbonds1.xcoords[i] <- seq1.xcoords[Hbonds[1,i]];
			Hbonds1.ycoords[i] <- seq1.ycoords[Hbonds[1,i]];
			Hbonds2.xcoords[i] <- seq1.xcoords[Hbonds[2,i]];
			Hbonds2.ycoords[i] <- seq1.ycoords[Hbonds[2,i]];
		} 
	}
	bonds4 <- segmentsGrob(Hbonds1.xcoords, Hbonds1.ycoords, Hbonds2.xcoords,
	                       Hbonds2.ycoords, default.units = "npc",
	                       gp = gpar(lwd = hbonds.lwd, lty = hbonds.lty));        

	# creation of labels of nucleotides
	labels1 <- textGrob(DNAseq1.vec, seq1.xcoords, seq1.ycoords, 
	                    default.units = "npc", gp = gpar(cex = labels.cex));
	labels2 <- textGrob(DNAseq2.vec, seq2.xcoords, seq2.ycoords, 
	                    default.units = "npc", gp = gpar(cex = labels.cex));

	# creation of nucleotides  
	if (circles == TRUE) {
		color <- "black";
	} else {
		color <- "white";
	}            
	r <- (convertX(unit(get.gpar()$fontsize*labels.cex, "points"), "npc",
			valueOnly = TRUE)/2)*1.2*circles.cex;
	nucleotides1 <- circleGrob(seq1.xcoords, seq1.ycoords, rep(r, n.seq1), gp = 
	                           gpar(fill = "white", col = color));
	nucleotides2 <- circleGrob(seq2.xcoords, seq2.ycoords, rep(r, n.seq2), gp = 
	                           gpar(fill = "white", col = color));
	
	# creation of additional lines
	addlines <- segmentsGrob(c(margin+dist1/2, margin+dist1/2, margin+dist1/2), 
	                         c(0.5+dist2, 0.5, 0.5-dist2), 
	                         c(margin+dist1, margin+dist1, margin+dist1),
	                         c(0.5+dist2, 0.5, 0.5-dist2),
	                         default.units = "npc", gp = gpar(lwd = mbonds.lwd, 
	                         lty = mbonds.lty));
	# drawing
	grid.draw(dirlab);
	grid.draw(bonds1);
	grid.draw(bonds2);
	grid.draw(bonds3);
	grid.draw(bonds4);
	grid.draw(addlines);
	grid.draw(nucleotides1);
	grid.draw(nucleotides2);
	grid.draw(labels1);        
	grid.draw(labels2);
	
	return(alignment)
}
