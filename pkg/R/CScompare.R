# Project: CSFA
# 
# Author: lucp8394
###############################################################################

#' Compare CS Results.
#' 
#' After applying different CSanalysis on the same data, you can compare 2 different results of connectivity scores and gene scores
#' Unless the result came from a Zhang and Gant analysis, you choose from which component (factor, PC, bicluster) the scores should be derived.
#' 
#' @export
#' @param CSresult1 First result.
#' @param CSresult2 Second result.
#' @param component1.plot If you are using a non-Zhang&Gant result, specify the bicluster, factor or principal component which should be used to derive connectivity scores from for the \emph{first} result.
#' @param component2.plot If you are using a non-Zhang&Gant result, specify the bicluster, factor or principal component which should be used to derive connectivity scores from for the \emph{second} result.
#' @param which Choose one or both plots which should be created. 1: CS Comparison Plot, 2: GS Comparison Plot
#' @param color.columns Vector of colors for the reference and query columns (compounds). If \code{NULL}, blue will be used for reference and black for query. Use this option to highlight reference columns and query columns of interest.
#' @param gene.thresP Vector of length 2 containing the positive gene thresholds for \code{CSresult1} and \code{CSresult2}. Genes above the threshold will be colored. (e.g. \code{c(1,2)})
#' @param gene.thresN Vector of length 2 containing the negative gene thresholds for \code{CSresult1} and \code{CSresult2}. Genes below the threshold will be colored. (e.g. \code{c(-1,-2)})
#' @param thresP.col Vector of length 2 containing the colors for the high gene scores for \code{CSresult1} and \code{CSresult2} (e.g. \code{c("blue","light blue")}).
#' @param thresN.col Vector of length 2 containing the colors for the low gene scores for \code{CSresult1} and \code{CSresult2} (e.g. \code{c("red","pink")}).
#' @param legend.names Option to draw a legend (about the highlights in \code{color.columns}) in the CS plot. If \code{NULL}, only references are in the legend.
#' @param legend.cols Colors to be used for the \code{legend.names}.
#' @param legend.pos The location of the legend: \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"} and \code{"center"}.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param basefilename Basename of the graphs if saved in pdf files
#' @param threshold.pvalues If both CSresult1 and CSresult contain pvalues (and adjusted pvalues), this threshold will be used to compare the number of overlapping significant results. 
#' @return A list object with the vector with correlation between Connectivity Scores (and Gene Scores if available). If both CSresult contain p-values the other list slots are filled with some comparison between the number of significant p-values.
#' @examples
#' \dontrun{
#' data("dataSIM",package="CSFA")
#' Mat1 <- dataSIM[,c(1:6)]
#' Mat2 <- dataSIM[,-c(1:6)]
#' 
#' MFA_analysis <- CSanalysis(Mat1,Mat2,"CSmfa")
#' ZHANG_analysis <- CSanalysis(Mat1,Mat2,"CSzhang")
#' 
#' CScompare(MFA_analysis,ZHANG_analysis,1)
#' }
CScompare <- function(CSresult1,CSresult2,component1.plot,component2.plot,threshold.pvalues=0.05,which=c(1,2),color.columns=NULL,gene.thresP=NULL,gene.thresN=NULL,thresP.col=c("blue","light blue"),thresN.col=c("red","pink"),legend.names=NULL,legend.cols=NULL,legend.pos="topright",plot.type="device",basefilename="CScompare"){
	
	if(class(CSresult1)!="CSresult"){stop("CSresult1 is not of the correct class type")}
	if(class(CSresult2)!="CSresult"){stop("CSresult2 is not of the correct class type")}
	
	refdim1 <- dim(CSresult1@CS$CS.ref)[1]
	refdim2 <- dim(CSresult2@CS$CS.ref)[1]
	querdim1 <- dim(CSresult1@CS$CS.query)[1]
	querdim2 <- dim(CSresult2@CS$CS.query)[1]
	
	if(refdim1!=refdim2){stop("Using 2 different reference matrices")}
	if(querdim1!=querdim2){stop("Using 2 different query matrices")}
	
	
	
	CSGS1 <- get.CS.GS(CSresult1,component1.plot)
	CSGS2 <- get.CS.GS(CSresult2,component2.plot)
	
	loadings1 <- CSGS1$CS
	loadings2 <- CSGS2$CS
	scores1 <- CSGS1$GS
	scores2 <- CSGS2$GS
	name1 <- CSGS1$name
	name2 <- CSGS2$name
	axename1 <- CSGS1$axename
	axename2 <- CSGS2$axename
	
	if(is.null(color.columns)){color.columns <- c(rep("blue",refdim1),rep("black",querdim1))} else {color.columns <- color.columns}
	if(is.null(legend.names)){legend.names <- c("References")}
	if(is.null(legend.cols)){legend.cols <- c("blue")}
	
	
	if(1%in%which){cor.CS <- compare.CS.plot(loadings1=loadings1,loadings2=loadings2,name1=name1,name2=name2,axename1=axename1,axename2=axename2,nref=refdim1,color.columns=color.columns,legend.names=legend.names,legend.cols=legend.cols,legend.pos=legend.pos,plot.type=plot.type,basefilename=basefilename)}else{cor.CS <- NULL}
	
	if(!(is.null(scores1)|is.null(scores2)|!(2 %in% which))){
		
		cor.GS <- compare.GS.plot(scores1=scores1,scores2=scores2,name1=name1,name2=name2,axename1=axename1,axename2=axename2,nref=refdim1,gene.thresP=gene.thresP,gene.thresN=gene.thresN,thresP.col=thresP.col,thresN.col=thresN.col,plot.type=plot.type,basefilename=basefilename)
	
	}else{cor.GS <- NULL}
		
	
	if(!is.null(cor.GS) & !is.null(cor.CS)){
		cor.temp <- c(cor.CS,cor.GS)
		names(cor.temp) <- c("CS Correlation","GS Correlation")
	}
	else if(!is.null(cor.CS)){
		cor.temp <- cor.CS
		names(cor.temp) <- c("CS Correlation")
	}
	else{
		cor.temp <- cor.GS
		names(cor.temp) <- c("GS Correlation")
	}
#	return(cor.temp)
		
	
	if(!is.null(CSresult1@permutation.object) & !is.null(CSresult2@permutation.object)){
		
		list.pval.dataframe <- list(CSresult1@permutation.object$pval.dataframe,CSresult2@permutation.object$pval.dataframe)
		out_pval_compare <- pvalue2_compare(list.pval.dataframe,threshold=threshold.pvalues)
				
		return(list(correlation=cor.temp,compare.pvalues=out_pval_compare$compare.pvalues,compare.pvalues.adjusted=out_pval_compare$compare.pvalues.adjusted,pval.data1=out_pval_compare$list.pval.dataframe[[1]],pval.data2=out_pval_compare$list.pval.dataframe[[2]]))
		
	}
	else{
		return(list(correlation=cor.temp,compare.pvalues=NULL,compare.pvalues.adjusted=NULL,pval.data1=NULL,pval.data2=NULL))
	}
}



get.CS.GS <- function(CSresult,component.plot){
	type <- CSresult@type
	
	if(type=="CSfabia"){
		loadings <- CSresult@object@L[,component.plot]
		scores <- t(CSresult@object@Z)[,component.plot]
		
		return(list(CS=loadings,GS=scores,name="FABIA",axename=paste0("Fabia BC ",component.plot)))
		
	}
	else if(type=="CSmfa"){
		loadings <- CSresult@object$quanti.var$coord[,component.plot]		
		scores <- CSresult@object$ind$coord[,component.plot]	
		
		return(list(CS=loadings,GS=scores,name="MFA",axename=paste0("MFA Factor ",component.plot)))
	}
	else if(type =="CSpca"){
		loadings <- CSresult@object$var$coord[,component.plot]
		scores <- CSresult@object$ind$coord[,component.plot]
		
		return(list(CS=loadings,GS=scores,name="PCA",axename=paste0("PCA PC ",component.plot)))
	}
	else if(type =="CSsmfa"){
		loadings <- CSresult@object$loadings[,component.plot]
		scores <- CSresult@object$scores[,component.plot]
		
		return(list(CS=loadings,GS=scores,name="sMFA",axename=paste0("sMFA Factor ",component.plot)))
	}
	else if(type == "CSzhang"){
		loadings <- rbind(CSresult@CS$CS.ref,CSresult@CS$CS.query)[,1]
		names(loadings) <- c(rownames(CSresult@CS$CS.ref),rownames(CSresult@CS$CS.query))
		
		return(list(CS=loadings,GS=NULL,name="ZHANG",axename="Zhang CS"))
	}
	else{
		stop("Result type not recognised")
	}
}
	


compare.CS.plot <- function(loadings1,loadings2,name1,name2,axename1,axename2,nref,color.columns,legend.names,legend.cols,legend.pos,plot.type,basefilename){
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}	
	
	##
	if(!is.null(color.columns)){groupCol <- color.columns} else { groupCol <- "black"}
	##	
	
	
	minX <- min(-1,loadings1,na.rm=TRUE)
	maxX <- max(1,loadings1,na.rm=TRUE)
	minY <- min(-1,loadings2,na.rm=TRUE)
	maxY <- max(1,loadings2,na.rm=TRUE)
	
	plot.in(plot.type,paste0(basefilename,"_CS_",name1,"_VS_",name2,".pdf"))
	par(mfrow=c(1,1))
	plot(loadings1,loadings2,xlim=c(minX,maxX),ylim=c(minY,maxY),col=groupCol,bg="grey",main=paste0(name1," VS ",name2," CScores - ",nref," Ref Compound"),xlab=paste0(axename1," Connectivity Scores"),ylab=paste0(axename2," Connectivity Scores"),pch=21)
	text(loadings1,loadings2, names(loadings1),	pos=1,	cex=0.5,	col=groupCol)
	if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
	plot.out(plot.type)
	cor.CS <- cor(loadings1,loadings2,use="complete.obs")
	return(cor.CS)
}



compare.GS.plot <- function(scores1,scores2,name1,name2,axename1,axename2,nref,gene.thresP,gene.thresN,thresP.col,thresN.col,plot.type,basefilename){
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}	
	
	## Gene Scores Coloring function
	.give.gene.color <- function(data,P1,N1,P2,N2,col.P1,col.N1,col.P2,col.N2){
		data1 <- data[1]
		data2 <- data[2]
		
		if(data2>=P2){
			if(data1<=N1){return(c(col.P2,col.N1))}
			else if(data1>=P1){return(c(col.P2,col.P1))}
			else if(data1>N1 & data1<P1){return(c(col.P2,col.P2))}
			else {return(c("grey","grey"))}
		}
		else if(data2<=N2){
			if(data1<=N1){return(c(col.N2,col.N1))}
			else if(data1>=P1){return(c(col.N2,col.P1))}
			else if(data1>N1 & data1<P1){return(c(col.N2,col.N2))}
			else {return(c("grey","grey"))}
		}
		else if(data1<=N1){return(c(col.N1,col.N1))}
		else if(data1>=P1){return(c(col.P1,col.P1))}
		else {return(c("grey","grey"))}
	}
		
	minX <- min(-1,scores1,na.rm=TRUE)
	maxX <- max(1,scores1,na.rm=TRUE)
	minY <- min(-1,scores2,na.rm=TRUE)
	maxY <- max(1,scores2,na.rm=TRUE)
	
	# gene colors
	if(!is.null(gene.thresP) | !is.null(gene.thresN)){
		if(is.null(gene.thresP)){gene.thresP <- c(99999,99999)}
		if(is.null(gene.thresN)){gene.thresN <- c(-99999,-99999)}
		
		scores <- rbind(scores1,scores2)
		list.scores <- as.list(as.data.frame(scores))
		list.colors <- lapply(X=list.scores,FUN=.give.gene.color,P1=gene.thresP[1],P2=gene.thresP[2],N1=gene.thresN[1],N2=gene.thresN[2],col.P1=thresP.col[1],col.P2=thresP.col[2],col.N1=thresN.col[1],col.N2=thresN.col[2])
		list.colors <- t(as.data.frame(list.colors))
		colnames(list.colors) <- c("bg","col")
	}
	else{
		list.colors <- matrix("grey",ncol=2,nrow=length(scores2)) 
	}
	
	plot.in(plot.type,paste0(basefilename,"_GS_",name1,"_VS_",name2,".pdf"))
	par(mfrow=c(1,1))
	plot(scores1,scores2,xlim=c(minX,maxX),ylim=c(minY,maxY),col=list.colors[,2],bg=list.colors[,1],main=paste0(name1," VS ",name2," GScores - ",nref," Ref Compound"),xlab=paste0(axename1," Gene Scores"),ylab=paste0(axename2," Gene Scores"),pch=21)
	text(scores1,scores2, names(scores1),	pos=1,	cex=0.5,	col=list.colors[,2])
	if(!is.null(gene.thresP)){
		abline(v=gene.thresP[1],lty=3)
		abline(h=gene.thresP[2],lty=3)
	}
	if(!is.null(gene.thresN)){
		abline(v=gene.thresN[1],lty=3)
		abline(h=gene.thresN[2],lty=3)
	}
	plot.out(plot.type)
	cor.GS <- cor(scores1,scores2,use="complete.obs")
	return(cor.GS)

}