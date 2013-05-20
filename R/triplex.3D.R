###
## 3D triplex visualization
##
## Author: Kamil Rajdl
## Date: 2012/10/28
## Package: triplex
##


###
## This function draws a 3-D model of a triplex 
##
## triplex   TriplexViews object
## opt       TRUE - structure of triplex will be optimalized
##           FALSE - structure will be drawn without optimalization
## A.col     color of Adine base
## T.col     color of Thymine base
## G.col     color of Guanine base
## C.col     color of Cytosine base
## bgr.col   color of background
## bbone.col color of backbone
## bbone.n   number of sides of backbone bonds
##
triplex.3D <- function(
	triplex,
	opt = TRUE,
	A.col = "red",
	T.col = "brown",
	G.col = "green",
	C.col = "blue",
	bbone.col = "violet",
	bgr.col = "white",
	bbone.n = 20)
{
	if ("rgl" %in% installed.packages()[,"Package"])
		library("rgl")
	else
		stop("Please install rgl package from CRAN to use this function.")
	
	# processing of inputs
	results <- triplex.input(triplex);
	triplex.type <- type(triplex);

	DNAseq1 <- results[[1]];
	DNAseq2 <- results[[2]];
	WCbonds <- results[[3]];
	Hbonds <- results[[4]];
	bends <- results[[5]];
	alignment <- results[[6]];

	bend.start <- bends[1];
	bend.end <- bends[2];   
	n.bonds <- dim(WCbonds)[2];
	n.bend <- bend.end - bend.start - 1;
	n.seq1 <- nchar(DNAseq1);
	n.seq2 <- nchar(DNAseq2);        

	# initial parameters and variables
	alpha.B <- 2.593247;
	dist.B <- 7.053836; 
	twist.B <- 36*pi/180;
	rise.B <- 3.4;
	r2 <- 0.5;               # radius of backbone

	# triplex must have length at least 2
	if (n.seq2 == 1) {
		stop("Triplex is too short, at least two triplets are necessary.");
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

	# coordinates of nucleotides
	seq1.xcoords <- rep(NA, n.seq1);
	seq1.ycoords <- rep(NA, n.seq1);
	seq1.zcoords <- rep(NA, n.seq1);
	seq2.xcoords <- rep(NA, n.seq2);
	seq2.ycoords <- rep(NA, n.seq2);
	seq2.zcoords <- rep(NA, n.seq2);

	# indexes of paired nucleotides
	ind1 <- rep(NA,n.seq2);
	ind2 <- rep(NA,n.seq2);
	ind3 <- rep(NA,n.seq2);

	# names of paired nucleotides
	nuc1 <- rep(NA,n.seq2);
	nuc2 <- rep(NA,n.seq2);
	nuc3 <- rep(NA,n.seq2);

	# angles and radii in levels 
	angles1 <- rep(NA,n.seq2);
	angles2 <- rep(NA,n.seq2);
	angles3 <- rep(NA,n.seq2);
	radii <- rep(NA,n.seq2); 

	i <- 1;
	j <- n.seq1;
	if (triplex.type %in% c("0","3","5","6")) {
		n1 <- n.seq2;
		n2 <- 1;
	}
	if (triplex.type %in% c("1","2","4","7")) {
		n1 <- 1;
		n2 <- n.seq2;
	}
	p <- 0;
	for (k in n1:n2) {
		while (!(i %in% WCbonds[1,]) && !(i %in% Hbonds1[1,]) &&  !(i %in% Hbonds2)) {
			i <- i+1;
		}
		while (!(j %in% WCbonds[1,]) && !(j %in% Hbonds1[1,]) &&  !(j %in% Hbonds2)) {
			j <- j-1;
		}
		ind1[p+1] <- i;
		ind2[p+1] <- j;
		ind3[p+1] <- k;
		nuc1[p+1] <- DNAseq1.vec[i];
		nuc2[p+1] <- DNAseq1.vec[j];
		nuc3[p+1] <- DNAseq2.vec[k];

		trpl <- toTriplet(DNAseq1.vec[i], DNAseq1.vec[j], DNAseq2.vec[k], triplex.type);
		tr <- getTR(trpl, triplex.type);
		t <- tr[1];
		r <- tr[2];  
		alpha <- pi-t;              
		angles <- getAngles(alpha, triplex.type);
		a1 <- angles[1];
		a2 <- angles[2];
		a3 <- angles[3];
		angles1[p+1] <- (a1 + p*twist.B);
		angles2[p+1] <- (a2 + p*twist.B);
		angles3[p+1] <- (a3 + p*twist.B);
		radii[p+1] <- r;
		i <- i+1;
		j <- j-1;
		p <- p+1;
	}

	# optimalization
	if ((opt) && (n.bonds>2)) {
		ret <- nlm(f=toMin, p=rep(0,n.seq2-1), alpha.opt = alpha.B, dist.opt = dist.B,
					s1a = angles1, s2a = angles2, s3a = angles3, r = radii, d = rise.B);
		angles <- ret[[2]];
	} else {
		angles <- rep(0,n.seq2-1);
	}

	r <- radii[1];
	seq1.xcoords[ind1[1]] <- r*cos(angles1[1]);
	seq1.ycoords[ind1[1]] <- r*sin(angles1[1]);
	seq1.zcoords[ind1[1]] <- 0;
	seq1.xcoords[ind2[1]] <- r*cos(angles2[1]);
	seq1.ycoords[ind2[1]] <- r*sin(angles2[1]);
	seq1.zcoords[ind2[1]] <- 0;
	seq2.xcoords[ind3[1]] <- r*cos(angles3[1]);
	seq2.ycoords[ind3[1]] <- r*sin(angles3[1]);
	seq2.zcoords[ind3[1]] <- 0;

	if (n.seq2 > 1) {
	for (p in 2:n.seq2) {
		i <- ind1[p];
		j <- ind2[p];
		k <- ind3[p];
		r <- radii[p];
		seq1.xcoords[i] <- r*cos(angles1[p] + angles[p-1]);
		seq1.ycoords[i] <- r*sin(angles1[p] + angles[p-1]);
		seq1.zcoords[i] <- (p-1)*rise.B;
		seq1.xcoords[j] <- r*cos(angles2[p] + angles[p-1]);
		seq1.ycoords[j] <- r*sin(angles2[p] + angles[p-1]);
		seq1.zcoords[j] <- (p-1)*rise.B;
		seq2.xcoords[k] <- r*cos(angles3[p] + angles[p-1]);
		seq2.ycoords[k] <- r*sin(angles3[p] + angles[p-1]);
		seq2.zcoords[k] <- (p-1)*rise.B;
	}                
	} 

	# coordinates of nucleotiodes in the bend
	if ((bend.end-bend.start)>1) {
		x1 <- c(seq1.xcoords[bend.start], seq1.ycoords[bend.start], seq1.zcoords[bend.start]);
		j <- 1;
		while (is.na(seq1.xcoords[bend.start-j])) {
			j <- j+1;
		}
		x2 <- c(seq1.xcoords[bend.start-j], seq1.ycoords[bend.start-j], seq1.zcoords[bend.start-j]);
		y1 <- c(seq1.xcoords[bend.end], seq1.ycoords[bend.end], seq1.zcoords[bend.end]);
		j <- 1;
		while (is.na(seq1.xcoords[bend.end+j])) {
			j <- j+1;  
		}
		y2 <- c(seq1.xcoords[bend.end+j], seq1.ycoords[bend.end+j], seq1.zcoords[bend.end+j]);
		m <- bend.end-bend.start-1; 
		coords <- countBend(x1,x2,y1,y2,dist.B,m);
		xcoords <- coords[[1]]; 
		ycoords <- coords[[2]]; 
		zcoords <- coords[[3]]; 
		
		for (i in 1:m) {
			seq1.xcoords[bend.start+i] <- xcoords[i];
			seq1.ycoords[bend.start+i] <- ycoords[i];
			seq1.zcoords[bend.start+i] <- zcoords[i];
		}
	}

	# coordinates of nucleotides opposite the gaps
	for (i in 1:n.seq1) {
		if (is.na(seq1.xcoords[i])) {
			j<-1
			while (is.na(seq1.xcoords[i+j])) {
				j<-j+1;
			}
			x1 <- c(seq1.xcoords[i-1], seq1.ycoords[i-1], seq1.zcoords[i-1]);
			y1 <- c(seq1.xcoords[i+j], seq1.ycoords[i+j], seq1.zcoords[i+j]);
			x2 <- c(0,0,x1[3]);
			y2 <- c(0,0,y1[3]);
			coords <- countBend(x1,x2,y1,y2,dist.B,j);
			xcoords <- coords[[1]]; 
			ycoords <- coords[[2]]; 
			zcoords <- coords[[3]]; 
			for (k in 1:j) {
				seq1.xcoords[i+k-1] <- xcoords[k];
				seq1.ycoords[i+k-1] <- ycoords[k];
				seq1.zcoords[i+k-1] <- zcoords[k];
			}                
		}
	}

	# coordinates of Watson-Crick and Hoogsteen bonds
	WCbonds1.xcoords <- rep(NA, n.bonds);
	WCbonds1.ycoords <- rep(NA, n.bonds);
	WCbonds1.zcoords <- rep(NA, n.bonds);
	WCbonds2.xcoords <- rep(NA, n.bonds);
	WCbonds2.ycoords <- rep(NA, n.bonds);
	WCbonds2.zcoords <- rep(NA, n.bonds);
	for (i in 1:n.bonds) {
		WCbonds1.xcoords[i] <- seq1.xcoords[WCbonds[1,i]];
		WCbonds1.ycoords[i] <- seq1.ycoords[WCbonds[1,i]];
		WCbonds1.zcoords[i] <- seq1.zcoords[WCbonds[1,i]];
		WCbonds2.xcoords[i] <- seq2.xcoords[WCbonds[2,i]];
		WCbonds2.ycoords[i] <- seq2.ycoords[WCbonds[2,i]];
		WCbonds2.zcoords[i] <- seq2.zcoords[WCbonds[2,i]];
	}


	Hbonds1.xcoords <- rep(NA, n.bonds);
	Hbonds1.ycoords <- rep(NA, n.bonds);
	Hbonds1.zcoords <- rep(NA, n.bonds);
	Hbonds2.xcoords <- rep(NA, n.bonds);
	Hbonds2.ycoords <- rep(NA, n.bonds);
	Hbonds2.zcoords <- rep(NA, n.bonds);
	if (dim(Hbonds1)[2] > 0) {
		for (i in 1:dim(Hbonds1)[2]) {
			Hbonds1.xcoords[i] <- seq1.xcoords[Hbonds[1,i]];
			Hbonds1.ycoords[i] <- seq1.ycoords[Hbonds[1,i]];
			Hbonds1.zcoords[i] <- seq1.zcoords[Hbonds[1,i]];
			Hbonds2.xcoords[i] <- seq2.xcoords[Hbonds[2,i]];
			Hbonds2.ycoords[i] <- seq2.ycoords[Hbonds[2,i]];
			Hbonds2.zcoords[i] <- seq2.zcoords[Hbonds[2,i]];
		} 
	}
	if (dim(Hbonds2)[2] > 0) {
		for (i in 1:dim(Hbonds2)[2]) {
			Hbonds1.xcoords[i] <- seq1.xcoords[Hbonds[1,i]];
			Hbonds1.ycoords[i] <- seq1.ycoords[Hbonds[1,i]];
			Hbonds1.zcoords[i] <- seq1.zcoords[Hbonds[1,i]];
			Hbonds2.xcoords[i] <- seq1.xcoords[Hbonds[2,i]];
			Hbonds2.ycoords[i] <- seq1.ycoords[Hbonds[2,i]];
			Hbonds2.zcoords[i] <- seq1.zcoords[Hbonds[2,i]];
		} 
	}

	# drawing
	rgl.bg(color = bgr.col);
	rgl.spheres(seq1.xcoords, seq1.ycoords, seq1.zcoords, rep(r2, n.seq1), col = bbone.col);
	rgl.spheres(seq2.xcoords, seq2.ycoords, seq2.zcoords, rep(r2, n.seq2), col = bbone.col);

	for (i in 1:(n.seq1-1)) {
		x <- c(seq1.xcoords[i], seq1.ycoords[i], seq1.zcoords[i]);
		y <- c(seq1.xcoords[i+1], seq1.ycoords[i+1], seq1.zcoords[i+1]);
		join(x,y,r2,bbone.col, bbone.n);
	}

	for (i in 1:(n.seq2-1)) {
		x <- c(seq2.xcoords[i], seq2.ycoords[i], seq2.zcoords[i]);
		y <- c(seq2.xcoords[i+1], seq2.ycoords[i+1], seq2.zcoords[i+1]);
		join(x,y,r2,bbone.col, bbone.n);
	}

	for (i in 1:n.bonds) {
	x1 <- c(WCbonds1.xcoords[i], WCbonds1.ycoords[i], WCbonds1.zcoords[i]);
	x2 <- c(WCbonds2.xcoords[i], WCbonds2.ycoords[i], WCbonds2.zcoords[i]);
	y1 <- c(Hbonds1.xcoords[i], Hbonds1.ycoords[i], Hbonds1.zcoords[i]);
	y2 <- c(Hbonds2.xcoords[i], Hbonds2.ycoords[i], Hbonds2.zcoords[i]);
	if (sum(abs(x1-y1)) < 10^(-3)) {
		x <- x2;
		y <- x1;
		z <- y2;
	}
	if (sum(abs(x2-y1)) < 10^(-3)) {
		x <- x1;
		y <- x2;
		z <- y2;
	} 
	if (sum(abs(x1-y2)) < 10^(-3)) {
		x <- x2;
		y <- x1;
		z <- y1;
	} 
	if (sum(abs(x2-y2)) < 10^(-3)) {
		x <- x1;
		y <- x2;
		z <- y1;
	}

	trpl <- toTriplet(nuc1[i], nuc2[i], nuc3[i], triplex.type);
	trpl <- substring(trpl, 1:3, 1:3);
	drawBonds(x, y, z, trpl[3], trpl[2], trpl[1], bbone.col, A.col, T.col, G.col, C.col);
	}
	return(alignment)
}


###
## Auxialiary functions
##

getColor <- function(nucleotide, A.col, T.col, G.col, C.col)
{
	color <- "black";
	if (nucleotide == "A") {
		color <- A.col;
	}
	if (nucleotide == "T") {
		color <- T.col;
	}
	if (nucleotide == "G") {
		color <- G.col;
	} 
	if (nucleotide == "C") {
		color <- C.col;
	}
	return(color);
}

join <- function(x,y,r,color,bbone.n)
{
	u <- x-y;
	a <- u[1];
	b <- u[2];
	cc <- u[3];
	if (abs(a) < 10^(-8)) {
		if (abs(b) < 10^(-8)) {
			v1 <- c(1,0,0);
			v2 <- c(0,1,0);
		} else {
			v1 <- c(1,0,0);
			v2 <- c(0,cc,-b);
		}
	} else {
		v1 <- c(b,-a,0);
		if ((abs(b) < 10^(-8)) || (abs(cc) < 10^(-8))) {
			v2 <- c(cc, 0, -a);
		} else {
			v2 <- c(b-(b^2+a^2)/b, -a, a*(b^2+a^2)/(b*cc));
		}   
	}
	v1 <- v1*r/sqrt(sum(v1*v1));
	v2 <- v2*r/sqrt(sum(v2*v2));
	n <- bbone.n;
	alpha <- 2*pi/n;
	for (i in 1:n) {
		p1 <- x + v1*cos((i-1)*alpha) + v2*sin((i-1)*alpha); 
		p2 <- x + v1*cos(i*alpha) + v2*sin(i*alpha);
		p3 <- y + v1*cos((i-1)*alpha) + v2*sin((i-1)*alpha); 
		p4 <- y + v1*cos(i*alpha) + v2*sin(i*alpha);
		rgl.quads(c(p1[1],p2[1],p4[1],p3[1]), c(p1[2],p2[2],p4[2],p3[2]), c(p1[3],p2[3],p4[3],p3[3]), col = color);
	}
}

getTR <- function(S, triplex.type)
{
	r = 10;
	t = 80;
	found = FALSE;
	if (triplex.type %in% c("4","5","6","7")) {
		if (S == "AAT") {
			r = 7.409632;
			t = 76.5364;
			found = TRUE;
		}
		if (S == "AGC") {
			r = 11.9128;
			t = 125.6857;
			found = TRUE;
		}
		if (S == "CAT") {
			r = 5.563016;
			t = 72.08654;
			found = TRUE;
		}
		if (S == "GGC") {
			r = 7.961633;
			t = 94.1352;
			found = TRUE;
		}
		if (S == "TAT") {
			r = 6.379035;
			t = 72.20278;
			found = TRUE;
		}
		if (S == "TCG") {
			r = 7.441224;
			t = 93.85905;
			found = TRUE;
		}
	}
	if (triplex.type %in% c("0","1","2","3")) {
		if (S == "CGC") {
			r = 8.897408;
			t = 108.8662;
			found = TRUE;
		}
		if (S == "GGC") {
			r = 7.161497;
			t = 74.82808;
			found = TRUE;
		}
		if (S == "GTA") {
			r = 12.28118;
			t = 125.5032;
			found = TRUE;
		}
		if (S == "TAT") {
			r = 7.896785;
			t = 103.7168;
			found = TRUE;
		}
		if (S == "TCG") {
			r = 6.469698;
			t = 70.47479;
			found = TRUE;
		}
		if (S == "TGC") {
			r = 6.630192;
			t = 77.89392;
			found = TRUE;
		}
	}
	if (!found) {
		warning("Triplex contains unknown triplet (its parameters were set to default).");
	}
	t = (t/180)*pi;
	return(c(t,r));
}

toTriplet <- function(nuc1, nuc2, nuc3, triplex.type)
{
	triplet <- paste(nuc1,nuc2,nuc3,sep="");
	if (triplex.type %in% c("0","3")) {
		triplet <- paste(nuc2,nuc3,nuc1,sep="");
	}
	if (triplex.type %in% c("1","2")) {
		triplet <- paste(nuc1,nuc3,nuc2,sep="");
	}
	if (triplex.type %in% c("4","7")) {
		triplet <- paste(nuc1,nuc2,nuc3,sep="");
	}
	if (triplex.type %in% c("5","6")) {
		triplet <- paste(nuc2,nuc1,nuc3,sep="");
	}
	return(triplet); 
}                                                                           

getAngles <- function(alpha, triplex.type)
{
	a1 <- 0;
	a2 <- 2*alpha;
	a3 <- alpha;
	
	if (triplex.type %in% c("0","3")) {
		a1 <- 0;
		a2 <- 2*alpha;
		a3 <- alpha;
	}
	if (triplex.type %in% c("1","2")) {
		a1 <- 2*alpha;
		a2 <- 0;
		a3 <- alpha;
	}
	if (triplex.type %in% c("4","7")) {
		a1 <- 2*alpha;
		a2 <- alpha;
		a3 <- 0;
	}
	if (triplex.type %in% c("5","6")) {
		a1 <- alpha;
		a2 <- 0;
		a3 <- 2*alpha;
	} 
	return(c(a1,a2,a3));
}

countBend <- function(x1,x2,y1,y2,bend.dist,m)
{
	xcoords <- rep(NA,m);
	ycoords <- rep(NA,m);
	zcoords <- rep(NA,m);   
	d1 <- sqrt(sum((x1-y1)^2));
	d2 <- bend.dist;
	z1 <- (x1+y1)/2;
	z2 <- (x2+y2)/2;                                         
	u <- y1-x1;
	v <- z1-z2;
	w <- c(u[2]*v[3]-u[3]*v[2],u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1]);
	direc <- c(w[2]*u[3]-w[3]*u[2],w[3]*u[1]-w[1]*u[3],w[1]*u[2]-w[2]*u[1]); 
		
	if ((m+1)*d2 <= d1) {                                # just line
		base1 <- u/(m+1);
		for (i in 1:m) {
			point <- x1 + base1*i;
			xcoords[i] <- point[1];
			ycoords[i] <- point[2];
			zcoords[i] <- point[3];
		}
	} else {     # circle

		if (d2 < d1*sin(pi/(2*(m+1)))) {                   # smaller than half circle
			r <- optimize(function(r) (2*r*sin((asin(d1/(2*r)))/(m+1))-d2)^2, 
											lower=d1/2, upper=(m+1)*d2)[[1]];
			S <- z1 - (direc/sqrt(sum(direc^2)))*sqrt(r^2-(d1^2)/4);
			base1 <- x1 - S;
			vv <- y1 - S;
			base2 <- base1 - (sum(base1*base1)/sum(base1*vv))*vv;
			if (sum(base1*vv) > 0) {
				base2 <- -base2*r/sqrt(sum(base2^2));
				} else {
				base2 <- base2*r/sqrt(sum(base2^2));
			}
			betta <- (2*asin(d1/(2*r)))/(m+1);
				
		} else if (d2 > d1*sin(pi/(2*(m+1)))) {             # bigger than half circle
			r <- optimize(function(r) (2*r*sin((pi-asin(d1/(2*r)))/(m+1))-d2)^2, 
											lower=d1/2, upper=(m+1)*d2)[[1]];
			S <- z1 + (direc/sqrt(sum(direc^2)))*sqrt(r^2-(d1^2)/4);
			base1 <- x1 - S;
			vv <- y1 - S;
			base2 <- base1 - (sum(base1*base1)/sum(base1*vv))*vv;
			if (sum(base1*vv) > 0) {
				base2 <- base2*r/sqrt(sum(base2^2));
				} else {
				base2 <- -base2*r/sqrt(sum(base2^2));
			}
			betta <- (2*pi - 2*asin(d1/(2*r)))/(m+1);
				
		} else {                                             # half circle
			r <- d1/2;
			S <- z1;
			base1 <- x1 - S;
			base2 <- (direc/sqrt(sum(direc^2)))*d1/2;
			betta <- pi/(m+1);
		}
			
		
		for (i in 1:m) {
			point <- S + base1*cos(i*betta)+base2*sin(i*betta);
			xcoords[i] <- point[1];
			ycoords[i] <- point[2];
			zcoords[i] <- point[3];    
		}
	}
	coords <- list(xcoords,ycoords,zcoords);
	return(coords);
}

drawBonds <- function(x, y, z, nuc1, nuc2, nuc3, bbone.col, A.col, T.col, G.col, C.col)
{
	S <- (x+y+z)/3;
	u <- S-x;
	v <- S-y;
	w <- S-z;    
	u <- (u/sqrt(sum(u^2)))*2.5;
	v <- (v/sqrt(sum(v^2)))*2.5;  
	w <- (w/sqrt(sum(w^2)))*2.5;
	joinRect(x, x+u, 0.2, bbone.col);
	joinRect(y, y+v, 0.2, bbone.col);
	joinRect(z, z+w, 0.2, bbone.col);
	drawBase(x, x+u, nuc1, A.col, T.col, G.col, C.col);
	drawBase(y, y+v, nuc2, A.col, T.col, G.col, C.col); 
	drawBase(z, z+w, nuc3, A.col, T.col, G.col, C.col);  
}

joinRect <- function(x,y,d,color)
{
	u <- y-x;
	v <- c(u[2],-u[1],0);
	v <- (v/sqrt(sum(v^2)))*(d/2);
	p1 <- x+v;
	p2 <- x-v;
	p3 <- y+v;
	p4 <- y-v;   
	rgl.quads(c(p1[1],p2[1],p4[1],p3[1]), c(p1[2],p2[2],p4[2],p3[2]), c(p1[3],p2[3],p4[3],p3[3]), col = color);
}

drawBase <- function(x,y,nuc, A.col, T.col, G.col, C.col)
{
	color = getColor(nuc, A.col, T.col, G.col, C.col);

	if (nuc %in% c("C","T")) {
		r <- 1.2;
		u <- y-x;
		S <- y + (u/sqrt(sum(u^2)))*r;
		base1 <- y-S;
		base2 <- c(base1[2],-base1[1],0); 
		xcoords <- c(NA,NA,S[1]);
		ycoords <- c(NA,NA,S[2]);
		zcoords <- c(NA,NA,S[3]);
		for (i in 1:6) {
			point1 <- S + base1*cos((i-1)*(pi/3))+base2*sin((i-1)*(pi/3));
			xcoords[1] <- point1[1];
			ycoords[1] <- point1[2];
			zcoords[1] <- point1[3];
			point2 <- S + base1*cos(i*(pi/3))+base2*sin(i*(pi/3));
			xcoords[2] <- point2[1];
			ycoords[2] <- point2[2];
			zcoords[2] <- point2[3];
			rgl.triangles(xcoords, ycoords, zcoords, col = color);    
		}
	}
	if (nuc %in% c("A","G")) {
		r <- 0.85*1.2;
		u <- y-x;
		S <- y + (u/sqrt(sum(u^2)))*r;
		base1 <- y-S;
		base2 <- c(base1[2],-base1[1],0); 
		if ((base1[1]*base2[2]-base1[2]*base2[1])<0) {
			base2 <- -base2;
		}
		xcoords <- c(NA,NA,S[1]);
		ycoords <- c(NA,NA,S[2]);
		zcoords <- c(NA,NA,S[3]);
		for (i in 1:5) {
			point1 <- S + base1*cos((i-1)*(2*pi/5)+4*pi/5)+base2*sin((i-1)*(2*pi/5)+4*pi/5);
			xcoords[1] <- point1[1];
			ycoords[1] <- point1[2];
			zcoords[1] <- point1[3];
			point2 <- S + base1*cos(i*(2*pi/5)+4*pi/5)+base2*sin(i*(2*pi/5)+4*pi/5);
			xcoords[2] <- point2[1];
			ycoords[2] <- point2[2];
			zcoords[2] <- point2[3];
			rgl.triangles(xcoords, ycoords, zcoords, col = color);    
		}
		r <- sqrt(sum((point1-point2)^2));
		z <- (point1+point2)/2;
		v <- z-S;
		v <- (v/sqrt(sum(v^2)))*(sqrt(3)*r/2);
		S <- z+v;
		base1 <- point1-S;
		base2 <- c(base1[2],-base1[1],0); 
		for (i in 1:6) {
			point1 <- S + base1*cos((i-1)*(pi/3))+base2*sin((i-1)*(pi/3));
			xcoords[1] <- point1[1];
			ycoords[1] <- point1[2];
			zcoords[1] <- point1[3];
			point2 <- S + base1*cos(i*(pi/3))+base2*sin(i*(pi/3));
			xcoords[2] <- point2[1];
			ycoords[2] <- point2[2];
			zcoords[2] <- point2[3];
			rgl.triangles(xcoords, ycoords, zcoords, col = color);    
		}
	}
}

toMin <- function(alpha, alpha.opt, dist.opt, s1a, s2a, s3a, r, d)
{
	n <- length(r);
	S <- 0;
	x1 <- c(r[1]*cos(s1a[1]), r[1]*sin(s1a[1]),-d);
	x2 <- c(r[2]*cos(s1a[2]+alpha[1]), r[2]*sin(s1a[2]+alpha[1]),0);
	x3 <- c(r[3]*cos(s1a[3]+alpha[2]), r[3]*sin(s1a[3]+alpha[2]),d);
	u <- x1-x2;
	v <- x3-x2;
	alpha.old <- acos(sum(u*v)/(sqrt(sum(u^2))*sqrt(sum(v^2))));
	S <- S+((alpha.old-alpha.opt)/(alpha.opt))^2 + ((sqrt(sum(u^2))-dist.opt)/(dist.opt))^2; 
	
	x1 <- c(r[1]*cos(s2a[1]), r[1]*sin(s2a[1]),-d);
	x2 <- c(r[2]*cos(s2a[2]+alpha[1]), r[2]*sin(s2a[2]+alpha[1]),0);
	x3 <- c(r[3]*cos(s2a[3]+alpha[2]), r[3]*sin(s2a[3]+alpha[2]),d);
	u <- x1-x2;
	v <- x3-x2;
	alpha.old <- acos(sum(u*v)/(sqrt(sum(u^2))*sqrt(sum(v^2))));
	S <- S+((alpha.old-alpha.opt)/(alpha.opt))^2 + ((sqrt(sum(u^2))-dist.opt)/(dist.opt))^2; 
		
	x1 <- c(r[1]*cos(s3a[1]), r[1]*sin(s3a[1]),-d);
	x2 <- c(r[2]*cos(s3a[2]+alpha[1]), r[2]*sin(s3a[2]+alpha[1]),0);
	x3 <- c(r[3]*cos(s3a[3]+alpha[2]), r[3]*sin(s3a[3]+alpha[2]),d);
	u <- x1-x2;
	v <- x3-x2;
	alpha.old <- acos(sum(u*v)/(sqrt(sum(u^2))*sqrt(sum(v^2))));
	S <- S+((alpha.old-alpha.opt)/(alpha.opt))^2 + ((sqrt(sum(u^2))-dist.opt)/(dist.opt))^2; 

	if (n>3) {
	for (i in 3:(n-1)) {
		x1 <- c(r[i-1]*cos(s1a[i-1]+alpha[i-2]), r[i-1]*sin(s1a[i-1]+alpha[i-2]),-d);
		x2 <- c(r[i]*cos(s1a[i]+alpha[i-1]), r[i]*sin(s1a[i]+alpha[i-1]),0);
		x3 <- c(r[i+1]*cos(s1a[i+1]+alpha[i]), r[i+1]*sin(s1a[i+1]+alpha[i]),d);
		u <- x1-x2;
		v <- x3-x2;
		alpha.old <- acos(sum(u*v)/(sqrt(sum(u^2))*sqrt(sum(v^2))));
		S <- S+((alpha.old-alpha.opt)/(alpha.opt))^2 + ((sqrt(sum(u^2))-dist.opt)/(dist.opt))^2; 
		
		x1 <- c(r[i-1]*cos(s2a[i-1]+alpha[i-2]), r[i-1]*sin(s2a[i-1]+alpha[i-2]),-d);
		x2 <- c(r[i]*cos(s2a[i]+alpha[i-1]), r[i]*sin(s2a[i]+alpha[i-1]),0);
		x3 <- c(r[i+1]*cos(s2a[i+1]+alpha[i]), r[i+1]*sin(s2a[i+1]+alpha[i]),d);
		u <- x1-x2;
		v <- x3-x2;
		alpha.old <- acos(sum(u*v)/(sqrt(sum(u^2))*sqrt(sum(v^2))));
		S <- S+((alpha.old-alpha.opt)/(alpha.opt))^2 + ((sqrt(sum(u^2))-dist.opt)/(dist.opt))^2;  
		
		x1 <- c(r[i-1]*cos(s3a[i-1]+alpha[i-2]), r[i-1]*sin(s3a[i-1]+alpha[i-2]),-d);
		x2 <- c(r[i]*cos(s3a[i]+alpha[i-1]), r[i]*sin(s3a[i]+alpha[i-1]),0);
		x3 <- c(r[i+1]*cos(s3a[i+1]+alpha[i]), r[i+1]*sin(s3a[i+1]+alpha[i]),d);
		u <- x1-x2;
		v <- x3-x2;
		alpha.old <- acos(sum(u*v)/(sqrt(sum(u^2))*sqrt(sum(v^2))));
		S <- S+((alpha.old-alpha.opt)/(alpha.opt))^2 + ((sqrt(sum(u^2))-dist.opt)/(dist.opt))^2; 
	}
	}
	return(S);
}
