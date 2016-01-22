# Project: Connectivity
# 
# Author: Gebruiker
###############################################################################



#######################################################################################################################################
#######################################################################################################################################
###### ZHANG ANALYSIS
## which:	- 1: Zhang Scores

#### NOTES !!!
# MFA -> QUANTI.VAR$COR
# PCA -> VAR$COR (NOT COORD) 
# (coord are the standardized coefficients. For mfa the cor and coord are the same due to the standardizing of the data beforehand, namely the weighting)




analyse_zhang <- function(dataref,dataquery,nref=NULL,nquery=NULL,ord.query=TRUE,permute=FALSE,B=100000,ntop.pvalues=20,ntop.scores=20,
					basefilename="analyseZhang",
					#column.interest=NULL,
					colour.query=NULL,legend.names=NULL,legend.cols=unique(colour.query),legend.x=NULL,legend.y=NULL,
					result.available=NULL,plot.type="pdf",print.top=TRUE,
					which=c(1)){
	
	
	## Plot-in and -out functions			
	plot.in <- function(plot.type,name){
			if(plot.type=="pdf"){pdf(name)}
			if(plot.type=="device"){dev.new()}
			if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}

	
	## Doing Zhang Analysis
	if(is.null(result.available)){
		zhang_result <- zhangscore(dataref=dataref,dataquery=dataquery,nref=nref,nquery=nquery,ord.query=ord.query,permute=permute,B=B,ntop.pvalues=ntop.pvalues,ntop.scores=ntop.scores,pvalue.method="BH")
		add.pvalue.color <- FALSE
		
	}
	else{
		zhang_result <- result.available@object
		
		if(!is.null(result.available@permutation.object)){
			add.pvalue.color <- TRUE	
		}else{
			add.pvalue.color <- FALSE
		}
		 
	}
	
	## Printing Top Connections Scores + Pvalues
	if(print.top){print(zhang_result$Top)}
	if(permute==TRUE){print(zhang_result$Toppvalues)}
	
	##
	if(!is.null(colour.query)){groupCol <- colour.query} else { groupCol <- "black"}
	if(add.pvalue.color){
		pvaltype <- ifelse(result.available@permutation.object$extra.parameters$method.adjust=="none","pvalues","pvalues.adjusted")
		signCol <- ifelse(result.available@permutation.object$CS.pval.dataframe[,pvaltype]<=0.05,"purple","grey")
		signCol <- signCol[-c(1:dim(dataref)[2])]
	}
	else{
		signCol <- "grey"
	}
	
	
	ncol_q <- ncol(dataquery)
	##
	
	if(1 %in% which){
		## PLOT: Zhang Score Plot
		plot.in(plot.type,paste0(basefilename,"_zhangscore.pdf"))
		par(mfrow=c(1,1))
		plot(zhang_result$All[,1],xlim=c(1,ncol_q),ylim=c(-1,1),pch=21,bg=signCol,col=groupCol,main=paste0("Zhang Score"),xlab="Compound Index",ylab="Connection Score")
		text(c(1:length(zhang_result$All[,1])),zhang_result$All[,1],rownames(zhang_result$All),	pos=1,	cex=0.5,	col=groupCol)
		abline(0,0,lty=3)
		legend.bg <- c()
		if(add.pvalue.color){
			legend.names <- c(legend.names,paste0(result.available@permutation.object$extra.parameters$method.adjust," adj. p-value <= 0.05"))
			legend.bg <- c(rep("white",length(legend.names)-1),"purple")
			legend.cols <- c(legend.cols,"white")
		}
		else{
			legend.bg <- "white"
		}
#		if(is.null(legend.x)){legend.x <- (ncol_q-0.45*ncol_q)}
#		if(is.null(legend.y)){legend.y <- 1}
		if(length(legend.names)>0){legend("topright",legend.names,pch=21,col=legend.cols,pt.bg=legend.bg,bty="n")}
		plot.out(plot.type)
	}	

	return(zhang_result)
}

#######################################################################################################################################
#######################################################################################################################################
###### MFA Analysis  (For multiple ref compounds)
## which:	-1: Percentage Variance Explained by factors
##			-2: Loadings for reference compounds
##			-3: Factor 'factor.plot' VS Factor '?' : Loadings & Genes
##			-4: Loadings & Genes for Factor 'factor.plot'
##			-5: Compound Profiles (Select if necessary)
##			-6: CS Rank Scores for 'factor.plot'

analyse_MFA <- function(data,group,type=rep("s",length(group)),ind.sup=NULL,ncp=5,name.group=NULL,num.group.sup=NULL,graph=FALSE,weight.col.mfa=NULL,row.w=NULL,axes=c(1,2),tab.comp=NULL,
				basefilename="analyseMFA",
				factor.plot=1,column.interest=NULL,row.interest=NULL,gene.thresP=NULL,gene.thresN=NULL,
				colour.columns=NULL,legend.names=NULL,legend.cols=unique(colour.columns),thresP.col="blue",thresN.col="red",
				result.available=NULL,plot.type="pdf",
				CSrank.refplot=FALSE,gene.highlight=NULL,profile.type="gene",
				which=c(1,2,3,4,5)){
	
			
				
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	## Doing MFA Analysis
	if(is.null(result.available)){
		resMFA <- MFA(data,	group=group, type=type, ind.sup=ind.sup,ncp=ncp,name.group=name.group,num.group.sup=num.group.sup,graph=graph,weight.col.mfa=weight.col.mfa,row.w=row.w,axes=axes,tab.comp=tab.comp)
		add.pvalue.color <- FALSE
		
	}
	else{
		resMFA <- result.available@object
		if(!is.null(result.available@permutation.object)){
			add.pvalue.color <- TRUE
		}else{
			add.pvalue.color <- FALSE
		}
	}
	
	## Printing the Echoufier Rv Correlation
	cat("Echoufier Rv Correlation:\n")
	print(resMFA$group$RV)
	
	if(1 %in% which){
		## PLOT: Amount of variance explained by PC's
		perc.var <- resMFA$eig[,2]				
		plot.in(plot.type,paste0(basefilename,"_MFApercvar.pdf"))
		par(mfrow=c(1,1))
		plot(perc.var,main="Percentage of Variance Explained",xlab="Number of Components",ylab="Perc. Var. Explained")
		plot.out(plot.type)
	}
	
	##
	loadings <- resMFA$quanti.var$coord	# For compounds	
	scores <- resMFA$ind$coord			# For genes
	if(!is.null(colour.columns)){groupCol <- colour.columns} else { groupCol <- "black"}
	ref.index <- c(1:group[1])
	##
	
	
	if(2 %in% which){
		## PLOT: Loadings for reference compounds -> select which PC to use if factor.plot is NOT given
		plot.in(plot.type,paste0(basefilename,"_MFA_RefLoadings.pdf"))
		plot(0,0,type="n",main=paste0("Loadings for Ref ",paste(ref.index,collapse=",")),xlab="Factor Index",ylab="Loadings",ylim=c(min(loadings[ref.index,]),max(loadings[ref.index,])),xlim=c(1,dim(loadings)[2]))
		for(i.ref in ref.index){
			points(c(1:dim(loadings)[2]),loadings[i.ref,],col=i.ref)
		}
		abline(0,0,lty=3)
		legend(dim(loadings)[2]*(1-0.4),max(loadings[ref.index,]),colnames(data)[ref.index],col=ref.index,bty="n",pch=21)
		plot.out(plot.type)
	}
	
	
	## Selecting the factor.plot
	if(is.null(factor.plot) | (3 %in% which) | (4 %in% which) | (5 %in% which) | (is.null(column.interest)&(6 %in% which)) | ((6 %in% which) & ( (!is.null(gene.thresP))   | (!is.null(gene.thresN)) ) ) | (7 %in% which)){
		if(is.null(factor.plot)){
			if(plot.type=="pdf" | !(2 %in% which)){
				dev.new()
				plot(0,0,type="n",main=paste0("Loadings for Ref ",paste(ref.index,collapse=",")),xlab="Factor Index",ylab="Loadings",ylim=c(min(loadings[ref.index,]),max(loadings[ref.index,])),xlim=c(1,dim(loadings)[2]))
				for(i.ref in ref.index){
					points(c(1:dim(loadings)[2]),loadings[i.ref,],col=i.ref)
				}
				abline(0,0,lty=3)
				legend(dim(loadings)[2]*(1-0.4),max(loadings[ref.index,]),colnames(data)[ref.index],col=ref.index,bty="n",pch=21)
			}
			cat("Please select with left mousebutton which Factor should be investigated. (The Factor with highest loading for reference)\nIf multiple reference were used, click in the middle of the group of points.\n")
			y.mean <- apply(loadings[ref.index,,drop=FALSE],MARGIN=2,FUN=mean)
			factor.plot <- identify(x=c(1:dim(loadings)[2]),y=y.mean,plot=TRUE,n=1,tolerance=0.5)
		}
		
		##
		factor2.plot <- ifelse(factor.plot==1,2,1)
		##
	}
	

	if(3 %in% which){
		## PLOT: PC 'factor.plot'  vs PC 'factor2.plot' - Loadings
		plot.in(plot.type,paste0(basefilename,"_MFA_PC",factor.plot,"PC",factor2.plot,"_loadings.pdf"))
		par(mfrow=c(1,1))
		plot(loadings[,factor.plot],loadings[,factor2.plot],  type="p",xlim=c(min(loadings[,factor.plot]),max(loadings[,factor.plot])),ylim=c(min(loadings[,factor2.plot]),max(loadings[,factor2.plot])),
				xlab=paste0("Factor ",factor.plot," - Compound Loadings"), 
				ylab=paste0("Factor ",factor2.plot," - Compound Loadings"),
				pch=21,
				bg="grey",
				col=groupCol,
				cex=1,main=paste0("PCA weighted (MFA) - Factor ",factor.plot," vs ",factor2.plot)
		)
		text(loadings[,factor.plot],loadings[,factor2.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#		if(length(legend.names)>0){legend(0.8,0.8,legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
		if(length(legend.names)>0){legend("topright",legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
		
		plot.out(plot.type)
	
	
		## PLOT: PC `factor.plot' vs PC 'factor2.plot' - Factors
		plot.in(plot.type,paste0(basefilename,"_MFA_PC",factor.plot,"PC",factor2.plot,"_factors.pdf"))
		plot(scores[,factor.plot],scores[,factor2.plot],  type="p",
				xlab=paste0("Factor ",factor.plot,"- Gene Factor Scores"), 
				ylab=paste0("Factor ",factor2.plot," - Gene Factor Scores"),
				pch=21,
				bg="grey",
				col="grey",
				cex=1,main=paste0("PCA weighted (MFA) - Factor ",factor.plot," vs ",factor2.plot)
		)
		text(scores[,factor.plot],scores[,factor2.plot], rownames(data),	pos=1,	cex=0.5,	col="grey")
		plot.out(plot.type)
	}
	
	
	if(4 %in% which){
		## PLOT: Genes Factors of 'factor.plot'
		bg.temp <-  col.temp <- rep("grey",length(scores[,factor.plot]))
		if(!is.null(gene.thresP)){
			temp.boolean <- (scores[,factor.plot]>=gene.thresP)
			bg.temp[temp.boolean] <- thresP.col
			col.temp[temp.boolean] <- thresP.col
		}
		if(!is.null(gene.thresN)){
			temp.boolean <- (scores[,factor.plot]<=gene.thresN)
			bg.temp[temp.boolean] <- thresN.col
			col.temp[temp.boolean] <- thresN.col
		}
		
		# highlighting genes
		if(!is.null(gene.highlight)){
			if(class(gene.highlight)=="numeric" | class(gene.highlight)=="integer"){
				col.temp[gene.highlight] <- "green"
			}
			if(class(gene.highlight)=="list"){
				if(length(gene.highlight)>5){stop("Too many different gene.highlight",call.=FALSE)}
				col.pool <- c("green","deeppink","darkorchid4","gold3","tan1")
				for(i.list in 1:length(gene.highlight)){
					col.temp[gene.highlight[[i.list]]] <- col.pool[i.list]
				}
			}
		}
		
		plot.in(plot.type,paste0(basefilename,"_MFAFactors.pdf"))
		plot(c(1:length(scores[,factor.plot])),scores[,factor.plot],  type="p",
				xlab="Gene Index", 
				ylab="Gene Factor Scores",
				pch=21,
				bg=bg.temp,
				col=col.temp,
				cex=1,main=paste0("Factor ",factor.plot," - Gene Factor Scores")
		)
		text(c(1:length(scores[,factor.plot])),scores[,factor.plot], rownames(data),	pos=1,	cex=0.5,	col=col.temp)
		if(!is.null(gene.thresP)){abline(gene.thresP,0,lty=3)}
		if(!is.null(gene.thresN)){abline(gene.thresN,0,lty=3)}
		plot.out(plot.type)
	}
	
	## SELECTING THE GENES OF INTEREST (replot in device if necessary)
	genes_interest <- row.interest #make the object in case of a cmpd profile
	
	if(6 %in% which & profile.type=="gene"){
		if(is.null(row.interest)){
			if(plot.type=="pdf" | !(4 %in% which)){
				bg.temp <-  col.temp <- rep("grey",length(scores[,factor.plot]))
				if(!is.null(gene.thresP)){
					temp.boolean <- (scores[,factor.plot]>=gene.thresP)
					bg.temp[temp.boolean] <- thresP.col
					col.temp[temp.boolean] <- thresP.col
				}
				if(!is.null(gene.thresN)){
					temp.boolean <- (scores[,factor.plot]<=gene.thresN)
					bg.temp[temp.boolean] <- thresN.col
					col.temp[temp.boolean] <- thresN.col
				}
				
				# highlighting genes
				if(!is.null(gene.highlight)){
					if(class(gene.highlight)=="numeric" | class(gene.highlight)=="integer"){
						col.temp[gene.highlight] <- "green"
					}
					if(class(gene.highlight)=="list"){
						if(length(gene.highlight)>5){stop("Too many different gene.highlight",call.=FALSE)}
						col.pool <- c("green","deeppink","darkorchid4","gold3","tan1")
						for(i.list in 1:length(gene.highlight)){
							col.temp[gene.highlight[[i.list]]] <- col.pool[i.list]
						}
					}
				}
				dev.new()
				plot(c(1:length(scores[,factor.plot])),scores[,factor.plot],  type="p",
						xlab="Gene Index", 
						ylab="Gene Factor Scores",
						pch=21,
						bg=bg.temp,
						col=col.temp,
						cex=1,main=paste0("Factor ",factor.plot," - Gene Factor Scores")
				)
				text(c(1:length(scores[,factor.plot])),scores[,factor.plot], rownames(data),	pos=1,	cex=0.5,	col=col.temp)
				if(!is.null(gene.thresP)){abline(gene.thresP,0,lty=3)}
				if(!is.null(gene.thresN)){abline(gene.thresN,0,lty=3)}
			}
			cat("Select as many genes as desired with left mouse button. Right-click to end selection procedure.\n\n")
			genes_interest <- identify(c(1:length(scores[,factor.plot])),scores[,factor.plot],n=999,plot=TRUE,labels=rownames(data))
			
		}
		else{
			genes_interest <- row.interest
		}
	}
	
	if(5 %in% which){
		## PLOT: Compound Loadings of 'factor.plot'
		
		signCol <- "grey"
		
		if(add.pvalue.color){
			if(result.available@permutation.object$extra.parameters$mfa.factor==factor.plot){
				pvaltype <- ifelse(result.available@permutation.object$extra.parameters$method.adjust=="none","pvalues","pvalues.adjusted")
				signCol <- ifelse(result.available@permutation.object$CS.pval.dataframe[,pvaltype]<=0.05,"purple","grey")
				
			}
		}
		
		plot.in(plot.type,paste0(basefilename,"_MFALoadings.pdf"))
		par(mfrow=c(1,1))
		plot(c(1:length(loadings[,factor.plot])),loadings[,factor.plot],  type="p",
				xlab="Compound Index", 
				ylab="Compound Loadings",
				pch=21,
				bg=signCol,
				col=groupCol,
				cex=1,main=paste0("Factor ",factor.plot," - Compound Loadings")
		)
		text(c(1:length(loadings[,factor.plot])),loadings[,factor.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#		if(length(legend.names)>0){legend((length(loadings[,factor.plot])*(1-0.45)),max(loadings),legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
		
		legend.names.loadings <- legend.names
		legend.cols.loadings <- legend.cols
		
		legend.bg <- "white"
		if(add.pvalue.color){
			if(result.available@permutation.object$extra.parameters$mfa.factor==factor.plot){
				legend.names.loadings <- c(legend.names.loadings,paste0(result.available@permutation.object$extra.parameters$method.adjust," adj. p-value <= 0.05"))
				legend.bg <- c(rep("white",length(legend.names.loadings)-1),"purple")
				legend.cols.loadings <- c(legend.cols.loadings,"white")	
			}
		}
	
		if(length(legend.names.loadings)>0){legend("topright",legend.names.loadings,pch=21,col=legend.cols.loadings,pt.bg=legend.bg,bty="n")}
		
		plot.out(plot.type)
	}
	
	## SELECTING THE COMPOUNDS OF INTEREST (replot in device if necessary)
	if(6 %in% which){
		if(is.null(column.interest)){
			if(plot.type=="pdf" | !(5 %in% which)){
				## PLOT: Compound Loadings of 'factor.plot'
				dev.new()
				par(mfrow=c(1,1))
				plot(c(1:length(loadings[,factor.plot])),loadings[,factor.plot],  type="p",
						xlab="Compound Index", 
						ylab="Compound Loadings",
						pch=21,
						bg="grey",
						col=groupCol,
						cex=1,main=paste0("Factor ",factor.plot," - Compound Loadings")
				)
				text(c(1:length(loadings[,factor.plot])),loadings[,factor.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#				if(length(legend.names)>0){legend((length(loadings[,factor.plot])*(1-0.45)),max(loadings),legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
				if(length(legend.names)>0){legend("topright",legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
				
			}
			cat("Select as many compounds as desired with left mouse button. Right-click to end selection procedure.\n\n")
			cmpds_interest <- identify(c(1:length(loadings[,factor.plot])),loadings[,factor.plot],n=999,plot=TRUE,labels=colnames(data))
		}
		else{
			cmpds_interest <- column.interest
		}
	}
	
	if( 6 %in% which){
		## PLOT: Profiles Plot
		CSprofiles(data=data,ref_index=ref.index,gene.select=genes_interest,cmpd.select=cmpds_interest,profile.type=profile.type,
				cmpd.loadings=loadings,gene.scores=scores,component.plot=factor.plot,gene.thresP=gene.thresP,gene.thresN=gene.thresN,
				basefilename=basefilename,plot.type=plot.type,thresP.col=thresP.col,thresN.col=thresN.col,main.base=paste0("MFA Factor ",factor.plot))
			

	}
	
	
	if(7 %in% which){
		
		
		signCol <- "grey"
		if(add.pvalue.color){
			if(result.available@permutation.object$extra.parameters$mfa.factor==factor.plot){
				pvaltype <- ifelse(result.available@permutation.object$extra.parameters$method.adjust=="none","pvalues","pvalues.adjusted")
				signCol <- ifelse(result.available@permutation.object$CSRank.pval.dataframe[,pvaltype]<=0.05,"purple","grey")
			}
		}
		
		legend.bg <- "white"
		legend.names.csrank <- legend.names
		legend.cols.csrank <- legend.cols
	
		if(add.pvalue.color){
			if(result.available@permutation.object$extra.parameters$mfa.factor==factor.plot){
				legend.names.csrank <- c(legend.names.csrank,paste0(result.available@permutation.object$extra.parameters$method.adjust," adj. p-value <= 0.05"))
				legend.bg <- c(rep("white",length(legend.names.csrank)-1),"purple")
				legend.cols.csrank <- c(legend.cols.csrank,"white")	
			}
		}
		
		
		out_CS_rank <- list(CSrank2(loadings,ref.index,color.columns=groupCol,ref.plot=CSrank.refplot,loadings_names=colnames(data),component.plot=factor.plot,type.component="Factor",plot=TRUE,plot.type=plot.type,basefilename=basefilename,signCol=signCol,legend.bg=legend.bg,legend.names=legend.names.csrank,legend.cols=legend.cols.csrank))
				
		names(out_CS_rank) <- paste0("Factor",factor.plot)
		
	}
	else{
		out_CS_rank <- list(CSrank2(loadings,ref.index,color.columns=groupCol,ref.plot=CSrank.refplot,loadings_names=colnames(data),component.plot=factor.plot,type.component="Factor",plot=FALSE))
		names(out_CS_rank) <- paste0("Factor",factor.plot)
		
	}
	
	## Returning the MFA result for further use
	
	out <- list(factor.select=factor.plot,result=resMFA,CSRank=out_CS_rank)	
	return(out)
	

}
#######################################################################################################################################
#######################################################################################################################################
###### PCA ANALYSIS: 1 REFERENCE COMPOUND
## which:	-1: Percentage Variance Explained by PC's
##			-2: Loadings for reference compounds
##			-3: PC 'factor.plot' VS PC '?' : Loadings & Genes
##			-4: Loadings & Genes for PC 'factor.plot'
##			-5: Compound Profiles (Select if necessary)
##			-6: CS Rank Scores for 'factor.plot'



analyse_PCA <- function(data, scale.unit = TRUE, ncp = 5, ind.sup = NULL,
		quanti.sup = NULL, quali.sup = NULL, row.w = NULL,
		col.w = NULL, graph = FALSE, axes = c(1,2),
		basefilename="analysePCA",
		ref.index=1,factor.plot=NULL,column.interest=NULL,gene.thresP=NULL,gene.thresN=NULL,
		colour.columns=NULL,legend.names=NULL,legend.cols=unique(colour.columns),thresP.col="blue",thresN.col="red",
		CSrank.refplot=FALSE,gene.highlight=NULL,profile.type="gene",row.interest=NULL,
		result.available=NULL,plot.type="pdf",which=c(1,2,3,4,5,6)){
	
	## Checking reference index is correct
	if(length(ref.index)>1){stop("There can only be 1 reference compound.",call.=FALSE)}
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	
	## Doing PCA Analysis
	if(is.null(result.available)){
		resPCA <- PCA(X=data,scale.unit=scale.unit,ncp=ncp,ind.sup=ind.sup,quanti.sup=quanti.sup,quali.sup=quali.sup,row.w=row.w,col.w=col.w,graph=graph,axes=axes)
	}
	else{
		resPCA <- result.available@object
	}
	
	if(1 %in% which){
		## PLOT: Amount of variance explained by PC's
		perc.var <- resPCA$eig[,2]
		plot.in(plot.type,paste0(basefilename,"_PCApercvar.pdf"))
		plot(perc.var,ylab="Perc. Var. Explained",main="Percentage of Variance Explained for Components",xlab="Principal Components Index")
		plot.out(plot.type)
	}
	
	##
	loadings <- resPCA$var$cor
	scores <- resPCA$ind$coord
	if(!is.null(colour.columns)){groupCol <- colour.columns} else { groupCol <- "black"}
	
	##
	
	if(2 %in% which){
		## PLOT: Loadings for reference compound -> select which PC to use if factor.plot is NOT given
		plot.in(plot.type,paste0(basefilename,"_PCA_RefLoadings.pdf"))
		plot(loadings[ref.index,],ylab=paste0("Loading for Ref ",ref.index),xlab="Principal Components Index",main=paste0("Loadings for ",ref.index))
		plot.out(plot.type)
	}
	
	## Selecting the factor.plot
	if(is.null(factor.plot)|(3 %in% which) | (4 %in% which) | (5 %in% which) | (is.null(column.interest)&(6 %in% which)) | ((6 %in% which) & ( (!is.null(gene.thresP))   | (!is.null(gene.thresN)) ) ) | (7%in%which) ){
		if(is.null(factor.plot)){
			if(plot.type=="pdf" | !(2 %in% which)){
				dev.new()
				plot(loadings[ref.index,],ylab=paste0("Loading for Ref ",ref.index),xlab="Principal Components Index",main=paste0("Loadings for ",ref.index))
			}
			cat("Please select with left mousebutton which PC should be investigated. (The PC with highest loading for reference)\n\n")
			factor.plot <- identify(loadings[ref.index,],n=1,plot=FALSE)
		}
	
		##
		factor2.plot <- ifelse(factor.plot==1,2,1)
		##
	}
	
	if(3 %in% which){
		## PLOT: PC 'factor.plot'  vs PC 'factor2.plot' - Loadings
		plot.in(plot.type,paste0(basefilename,"_PCA_PC",factor.plot,"PC",factor2.plot,"_loadings.pdf"))
		par(mfrow=c(1,1))
		plot(loadings[,factor.plot],loadings[,factor2.plot],  type="p",xlim=c(min(loadings[,factor.plot]),max(loadings[,factor.plot])),ylim=c(min(loadings[,factor2.plot]),max(loadings[,factor2.plot])),
				xlab=paste0("PC ",factor.plot," - Compound Loadings"), 
				ylab=paste0("PC ",factor2.plot," - Compound Loadings"),
				pch=21,
				bg="grey",
				col=groupCol,
				cex=1,main=paste0("PCA - PC ",factor.plot," vs ",factor2.plot)
		)
		text(loadings[,factor.plot],loadings[,factor2.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#		if(length(legend.names)>0){legend(0.8,0.8,legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
		if(length(legend.names)>0){legend("topright",legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
		
		plot.out(plot.type)
	
		## PLOT: PC `factor.plot' vs PC 'factor2.plot' - Factors
		plot.in(plot.type,paste0(basefilename,"_PCA_PC",factor.plot,"PC",factor2.plot,"_factors.pdf"))
		plot(scores[,factor.plot],scores[,factor2.plot],  type="p",
				xlab=paste0("PC ",factor.plot,"- Gene Factor Scores"), 
				ylab=paste0("PC ",factor2.plot," - Gene Factor Scores"),
				pch=21,
				bg="grey",
				col="grey",
				cex=1,main=paste0("PCA  - PC ",factor.plot," vs ",factor2.plot)
		)
		text(scores[,factor.plot],scores[,factor2.plot], rownames(data),	pos=1,	cex=0.5,	col="grey")
		plot.out(plot.type)
	}
	
	if(4 %in% which){
		## PLOT: Genes Factors of 'factor.plot'
		bg.temp <-  col.temp <- rep("grey",length(scores[,factor.plot]))
		if(!is.null(gene.thresP)){
			temp.boolean <- (scores[,factor.plot]>=gene.thresP)
			bg.temp[temp.boolean] <- thresP.col
			col.temp[temp.boolean] <- thresP.col
		}
		if(!is.null(gene.thresN)){
			temp.boolean <- (scores[,factor.plot]<=gene.thresN)
			bg.temp[temp.boolean] <- thresN.col
			col.temp[temp.boolean] <- thresN.col
		}
		
		# highlighting genes
		if(!is.null(gene.highlight)){
			if(class(gene.highlight)=="numeric" | class(gene.highlight)=="integer"){
				col.temp[gene.highlight] <- "green"
			}
			if(class(gene.highlight)=="list"){
				if(length(gene.highlight)>5){stop("Too many different gene.highlight",call.=FALSE)}
				col.pool <- c("green","deeppink","darkorchid4","gold3","tan1")
				for(i.list in 1:length(gene.highlight)){
					col.temp[gene.highlight[[i.list]]] <- col.pool[i.list]
				}
			}
		}
		
		plot.in(plot.type,paste0(basefilename,"_PCAFactors.pdf"))
		plot(c(1:length(scores[,factor.plot])),scores[,factor.plot],  type="p",
				xlab="Gene Index", 
				ylab="Gene Factor Scores",
				pch=21,
				bg=bg.temp,
				col=col.temp,
				cex=1,main=paste0("PC ",factor.plot," - Gene Factor Scores")
		)
		text(c(1:length(scores[,factor.plot])),scores[,factor.plot], rownames(data),	pos=1,	cex=0.5,	col=col.temp)
		if(!is.null(gene.thresP)){abline(gene.thresP,0,lty=3)}
		if(!is.null(gene.thresN)){abline(gene.thresN,0,lty=3)}
		plot.out(plot.type)
	}
	
	## SELECTING THE GENES OF INTEREST (replot in device if necessary)
	if(6 %in% which & profile.type=="gene"){
		if(is.null(row.interest)){
			if(plot.type=="pdf" | !(4 %in% which)){
				dev.new()
				bg.temp <-  col.temp <- rep("grey",length(scores[,factor.plot]))
				if(!is.null(gene.thresP)){
					temp.boolean <- (scores[,factor.plot]>=gene.thresP)
					bg.temp[temp.boolean] <- thresP.col
					col.temp[temp.boolean] <- thresP.col
				}
				if(!is.null(gene.thresN)){
					temp.boolean <- (scores[,factor.plot]<=gene.thresN)
					bg.temp[temp.boolean] <- thresN.col
					col.temp[temp.boolean] <- thresN.col
				}
				
				# highlighting genes
				if(!is.null(gene.highlight)){
					if(class(gene.highlight)=="numeric" | class(gene.highlight)=="integer"){
						col.temp[gene.highlight] <- "green"
					}
					if(class(gene.highlight)=="list"){
						if(length(gene.highlight)>5){stop("Too many different gene.highlight",call.=FALSE)}
						col.pool <- c("green","deeppink","darkorchid4","gold3","tan1")
						for(i.list in 1:length(gene.highlight)){
							col.temp[gene.highlight[[i.list]]] <- col.pool[i.list]
						}
					}
				}
				plot(c(1:length(scores[,factor.plot])),scores[,factor.plot],  type="p",
						xlab="Gene Index", 
						ylab="Gene Factor Scores",
						pch=21,
						bg=bg.temp,
						col=col.temp,
						cex=1,main=paste0("PC ",factor.plot," - Gene Factor Scores")
				)
				text(c(1:length(scores[,factor.plot])),scores[,factor.plot], rownames(data),	pos=1,	cex=0.5,	col=col.temp)
				if(!is.null(gene.thresP)){abline(gene.thresP,0,lty=3)}
				if(!is.null(gene.thresN)){abline(gene.thresN,0,lty=3)}
			}
			cat("Select as many genes as desired with left mouse button. Right-click to end selection procedure.\n\n")
			genes_interest <- identify(c(1:length(scores[,factor.plot])),scores[,factor.plot],n=999,plot=TRUE,labels=rownames(data))
			
		}
		else{
			genes_interest <- row.interest
		}
	}
	
	if(5 %in% which){
	
		## PLOT: Compound Loadings of 'factor.plot'
		plot.in(plot.type,paste0(basefilename,"_PCALoadings.pdf"))
		par(mfrow=c(1,1))
		plot(c(1:length(loadings[,factor.plot])),loadings[,factor.plot],  type="p",
				xlab="Compound Index", 
				ylab="Compound Loadings",
				pch=21,
				bg="grey",
				col=groupCol,
				cex=1,main=paste0("PC ",factor.plot," - Compound Loadings")
		)
		text(c(1:length(loadings[,factor.plot])),loadings[,factor.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#		if(length(legend.names)>0){legend((length(loadings[,factor.plot])*(1-0.45)),max(loadings),legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
		if(length(legend.names)>0){legend("topright",legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
		
		plot.out(plot.type)
	}
	
	## SELECTING THE COMPOUNDS OF INTEREST (replot in device if necessary)
	if(6 %in% which){
		if(is.null(column.interest)){
			if(plot.type=="pdf" | !(5 %in% which)){
				## PLOT: Compound Loadings of 'factor.plot'
				dev.new()
				par(mfrow=c(1,1))
				plot(c(1:length(loadings[,factor.plot])),loadings[,factor.plot],  type="p",
						xlab="Compound Index", 
						ylab="Compound Loadings",
						pch=21,
						bg="grey",
						col=groupCol,
						cex=1,main=paste0("PC ",factor.plot," - Compound Loadings")
				)
				text(c(1:length(loadings[,factor.plot])),loadings[,factor.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#				if(length(legend.names)>0){legend((length(loadings[,factor.plot])*(1-0.45)),1,legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
				if(length(legend.names)>0){legend("topright",legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
				
			}
			cat("Select as many compounds as desired with left mouse button. Right-click to end selection procedure.\n\n")
			cmpds_interest <- identify(c(1:length(loadings[,factor.plot])),loadings[,factor.plot],n=999,plot=TRUE,labels=colnames(data))
		}
		else{
			cmpds_interest <- column.interest
		}
	}
	
	if( 6 %in% which){
		## PLOT: Profiles Plot
		CSprofiles(data=data,ref_index=ref.index,gene.select=genes_interest,cmpd.select=cmpds_interest,profile.type=profile.type,
				cmpd.loadings=loadings,gene.scores=scores,component.plot=factor.plot,gene.thresP=gene.thresP,gene.thresN=gene.thresN,
				basefilename=basefilename,plot.type=plot.type,thresP.col=thresP.col,thresN.col=thresN.col,main.base=paste0("PCA PC ",factor.plot))
		
		
	}
	
	if(7 %in% which){
		legend.bg <- "white"
		legend.names.csrank <- legend.names
		legend.cols.csrank <- legend.cols

		out_CS_rank <- list(CSrank2(loadings,ref.index,color.columns=groupCol,ref.plot=CSrank.refplot,loadings_names=colnames(data),component.plot=factor.plot,type.component="PC",plot=TRUE,plot.type=plot.type,basefilename=basefilename,legend.bg=legend.bg,legend.names=legend.names.csrank,legend.cols=legend.cols.csrank))
		names(out_CS_rank) <- paste0("Factor",factor.plot)
		
		
	}
	else{
		out_CS_rank <- list(CSrank2(loadings,ref.index,color.columns=groupCol,ref.plot=CSrank.refplot,loadings_names=colnames(data),component.plot=factor.plot,type.component="PC",plot=FALSE))
		names(out_CS_rank) <- paste0("Factor",factor.plot)
		
	}
	
	## Return the result of PCA
		
	out <- list(factor.select=factor.plot,result=resPCA,CSRank=out_CS_rank)	
	return(out)
	
	
}






analyse_Zhang_MFAPCA <- function(data,resZhang,resPCA=NULL,resMFA=NULL,ressparseMFA=NULL,resZhang2=NULL,
		ref.index,factor.plot,
		colour.columns=NULL,legend.names=NULL,legend.cols=unique(colour.columns),legend.pos="topright",
		plot.type="pdf",basefilename="Zhang_PCAMFA"){
	
	## Check if not both results
	#if(!is.null(resPCA) & !is.null(resMFA)){stop("Cannot PCA and MFA result at the same time.",call.=FALSE)}
	if(sum(c(!is.null(resZhang2),!is.null(resPCA),!is.null(resMFA),!is.null(ressparseMFA)))>1){stop("Cannot PCA, MFA or sMFA result at the same time.",call.=FALSE)}
	
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	
	## Adding Zhang scores for the references. They will all get 1! (Even if multiple references were used)
	zhang.temp <- rep(0,dim(data)[2])
	i.zhang <- 1
	for(ii in 1:dim(data)[2]){
		if(!(ii %in% ref.index)){
			zhang.temp[ii] <- resZhang$All[i.zhang,1]
			i.zhang <- i.zhang+1
		}
		else{
			zhang.temp[ii] <- 1
		}
	}
	
	##
	if(!is.null(colour.columns)){groupCol <- colour.columns} else { groupCol <- "black"}
	##	
	
	## PLOT: Zhang Scores VS PCA Loadings (Ref will get Zhang Score=1)
	if(!is.null(resPCA)){
		
		loadings <- resPCA$var$coord
		minX <- min(-1,zhang.temp,na.rm=TRUE)
		maxX <- max(1,zhang.temp,na.rm=TRUE)
		minY <- min(-1,loadings[,factor.plot])
		maxY <- max(1,loadings[,factor.plot])
		
		plot.in(plot.type,paste0(basefilename,"_loadingsPCA.pdf"))
		par(mfrow=c(1,1))
		plot(zhang.temp,loadings[,factor.plot],xlim=c(minX,maxX),ylim=c(minY,maxY),col=groupCol,bg="grey",main="Zhang VS PCA - 1 Ref Compound",xlab="Zhang Score",ylab=paste0("PC",factor.plot," Loadings"),pch=21)
		text(zhang.temp,loadings[,factor.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
		abline(a=0,b=1,lty=3)
#		if(is.null(legend.x)){legend.x <- max(zhang.temp)*(1-0.2)}
#		if(is.null(legend.y)){legend.y <- max(loadings[,factor.plot])*(1-0.2)}
		if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
		plot.out(plot.type)
		
	}
	
	
	if(!is.null(resMFA)){
		
		loadings <- resMFA$quanti.var$coord	# For compounds	
		minX <- min(-1,zhang.temp,na.rm=TRUE)
		maxX <- max(1,zhang.temp,na.rm=TRUE)
		minY <- min(-1,loadings[,factor.plot])
		maxY <- max(1,loadings[,factor.plot])
		
		plot.in(plot.type,paste0(basefilename,"_loadingsMFA.pdf"))
		par(mfrow=c(1,1))
		plot(zhang.temp,loadings[,factor.plot],xlim=c(minX,maxX),ylim=c(minY,maxY),col=groupCol,bg="grey",main=paste0("Zhang VS MFA - ",length(ref.index)," Ref Compound"),xlab="Zhang Score",ylab=paste0("Factor",factor.plot," Loadings"),pch=21)
		text(zhang.temp,loadings[,factor.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
		abline(a=0,b=1,lty=3)
#		if(is.null(legend.x)){legend.x <- max(zhang.temp)*(1-0.2)}
#		if(is.null(legend.y)){legend.y <- max(loadings[,factor.plot])*(1-0.2)}
		if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
				
		plot.out(plot.type)
		
	}
	
	if(!is.null(ressparseMFA)){
		loadings <- ressparseMFA$loadings	# For compounds	
		minX <- min(-1,zhang.temp,na.rm=TRUE)
		maxX <- max(1,zhang.temp,na.rm=TRUE)
		minY <- min(-1,loadings[,factor.plot])
		maxY <- max(1,loadings[,factor.plot])
		
		plot.in(plot.type,paste0(basefilename,"_loadingssMFA.pdf"))
		par(mfrow=c(1,1))
		plot(zhang.temp,loadings[,factor.plot],xlim=c(minX,maxX),ylim=c(minY,maxY),col=groupCol,bg="grey",main=paste0("Zhang VS sMFA - ",length(ref.index)," Ref Compound"),xlab="Zhang Score",ylab=paste0("Factor",factor.plot," Loadings"),pch=21)
		text(zhang.temp,loadings[,factor.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#		abline(a=0,b=1,lty=3)
#		if(is.null(legend.x)){legend.x <- max(zhang.temp)*(1-0.2)}
#		if(is.null(legend.y)){legend.y <- max(loadings[,factor.plot])*(1-0.2)}
		if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
		
		plot.out(plot.type)
		
	}
	
	if(!is.null(resZhang2)){
		loadings <- rep(0,dim(data)[2])
		i.zhang <- 1
		for(ii in 1:dim(data)[2]){
			if(!(ii %in% ref.index)){
				loadings[ii] <- resZhang2$All[i.zhang,1]
				i.zhang <- i.zhang+1
			}
			else{
				loadings[ii] <- 1
			}
		}
		
		minX <- min(-1,zhang.temp,na.rm=TRUE)
		maxX <- max(1,zhang.temp,na.rm=TRUE)
		minY <- min(-1,loadings,na.rm=TRUE)
		maxY <- max(1,loadings,na.rm=TRUE)
		
		plot.in(plot.type,paste0(basefilename,"_scoresZhang2.pdf"))
		par(mfrow=c(1,1))
		plot(zhang.temp,loadings,xlim=c(minX,maxX),ylim=c(minY,maxY),col=groupCol,bg="grey",main=paste0("Zhang VS Zhang - ",length(ref.index)," Ref Compound"),xlab="Zhang Score 1",ylab="Zhagn Score 2",pch=21)
		text(zhang.temp,loadings, colnames(data),	pos=1,	cex=0.5,	col=groupCol)
		if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
		
		plot.out(plot.type)
		
		
		
	}
	if(!is.null(resZhang2)){loadings.temp <- loadings}else{loadings.temp <- loadings[,factor.plot]}

	#print(cor(zhang.temp,loadings[,factor.plot]))
	cor.temp <- cor(zhang.temp,loadings.temp,use="complete.obs")
	names(cor.temp) <- "CS.cor"
	return(cor.temp)
}

#######################################################################################################################################
#######################################################################################################################################
###### FABIA Analysis  
## which:	-1: Information Content of Biclusters
##			-2: Loadings for reference compounds
##			-3: Bicluster 'BC1.plot' VS Factor 'BC2.plot' : Loadings & Genes
##			-4:  Genes for Factor 'BC.Plot'  (These can be multiple)
##			-5: Loadings for Factor 'BC.plot'
##			-6: Compound Profiles (Select if necessary)
##			-7: CS Rank Scores for 'BC.plot'



analyse_fabia <- function(data,p=13,alpha=0.01,cyc=500,spl=0,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0,
					basefilename="analyseFABIA",weighted.data=FALSE,
					ref.index=c(1),
					BC.plot=NULL,
					column.interest=NULL,row.interest=NULL,gene.thresP=NULL,gene.thresN=NULL,
					colour.columns=NULL,legend.names=NULL,legend.cols=unique(colour.columns),thresP.col="blue",thresN.col="red",
					result.available=NULL,plot.type="pdf",
					CSrank.refplot=FALSE,gene.highlight=NULL,profile.type="gene",
					which=c(1,2,3,4,5,6)){
	
			
		## Plot-in and -out functions
		plot.in <- function(plot.type,name){
			if(plot.type=="pdf"){pdf(name)}
			if(plot.type=="device"){dev.new()}
			if(plot.type=="sweave"){}
		}
		plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}			
						
		## Doing Fabia Analysis  (Note that the data will be transposed)
		if(weighted.data==TRUE){
			
			total.index <- c(1:dim(data)[2])
			Mat1 <- data[,ref.index,drop=FALSE]
			Mat2 <- data[,total.index[-ref.index]]
			
			data <- getWeightedDat(Mat1,Mat2,scale.unit=TRUE)[[1]]
				
			# FILL IN HERE THE PROCESS TO WEIGHT/NORMALIZE THE DATA
		}
		
		
		
		## Do the FABIA Analysis
		if(is.null(result.available)){
			resFAB <- fabia(X=t(data),p=p,alpha=alpha,cyc=cyc,spl=spl,spz=spz,non_negative=non_negative,random=random,center=center,norm=norm,scale=scale,lap=lap,nL=nL,lL=lL,bL=bL)
		}
		else{
			resFAB <- result.available@object
		}
		
		if(1 %in% which){
			## PLOT: Information Content 
			plot.in(plot.type,paste0(basefilename,"_FABIA_IC.pdf"))
			showSelected(resFAB,which=c(1))
			plot.out(plot.type)
		}
		
		##
		loadings <- resFAB@L
		scores <- t(resFAB@Z)
		if(!is.null(colour.columns)){groupCol <- colour.columns} else { groupCol <- "black"}
		##
		
		
		if(2 %in% which){
			## PLOT: Loadings for reference compounds
			plot.in(plot.type,paste0(basefilename,"_FABIA_RefLoadings"))
			plot(0,0,type="n",main=paste0("Loadings for Ref ",paste(ref.index,collapse=",")),xlab="BC Index",ylab="Loadings",ylim=c(min(loadings[ref.index,]),max(loadings[ref.index,])),xlim=c(1,dim(loadings)[2]))
			for(i.ref in ref.index){
				points(c(1:dim(loadings)[2]),loadings[i.ref,],col=i.ref)
			}
			abline(0,0,lty=3)
			legend(dim(loadings)[2]*(1-0.2),max(loadings[ref.index,]),colnames(data)[ref.index],col=ref.index,bty="n",pch=21)
			plot.out(plot.type)
		}
		
		 
		## Select Biclusters from 'loadings for reference compounds' (redraw if necessary)
		if(is.null(BC.plot)|(3 %in% which) | (4 %in% which)| (5 %in% which) | (is.null(column.interest)&(6 %in% which)) | ((6 %in% which) & ( (!is.null(gene.thresP))   | (!is.null(gene.thresN)) ) ) | (7 %in% which) ){
			
			if(is.null(BC.plot)){
				if(plot.type=="pdf" | !(2 %in% which)){
					dev.new()
					plot(0,0,type="n",main=paste0("Loadings for Ref ",paste(ref.index,collapse=",")),xlab="BC Index",ylab="Loadings",ylim=c(min(loadings[ref.index,]),max(loadings[ref.index,])),xlim=c(1,dim(loadings)[2]))
					for(i.ref in ref.index){
						points(c(1:dim(loadings)[2]),loadings[i.ref,],col=i.ref)
					}
					abline(0,0,lty=3)
					legend(dim(loadings)[2]*(1-0.2),max(loadings[ref.index,]),colnames(data)[ref.index],col=ref.index,bty="n",pch=21)
				}
				out.overlap <- fabia.overlap(resFAB,ref.index)
				cat("Bicluster Suggestions (BC's which overlap with Ref & Query with default thresholds):\n")
				cat("------------------------------------------------------------------------------------\n\n")
				print(out.overlap)
				cat("\n Please select with left mousebutton which Bicluster should be investigated. \n If multiple reference were used, click in the middle of the group of points.\nRight-click to end selection procedure.\n")
				y.mean <- apply(loadings[ref.index,,drop=FALSE],MARGIN=2,FUN=mean)
				BC.plot <- identify(x=c(1:dim(loadings)[2]),y=y.mean,plot=TRUE,n=dim(loadings)[2],tolerance=0.5)
			
			}
			if(length(BC.plot)>1){
				BC1.plot <- BC.plot[1]
				BC2.plot <- BC.plot[2]
			}
			else{
				BC1.plot <- BC.plot
				BC2.plot <- ifelse(BC1.plot==1,2,1)
			}
		}
			
		if(3 %in% which){
			## PLOT: PC 'BC1.plot'  vs PC 'BC2.plot' - Loadings
			plot.in(plot.type,paste0(basefilename,"_FABIA_BC",BC1.plot,"BC",BC2.plot,"_loadings.pdf"))
			par(mfrow=c(1,1))
			plot(loadings[,BC1.plot],loadings[,BC2.plot],  type="p",xlim=c(min(loadings[,BC1.plot]),max(loadings[,BC1.plot])),ylim=c(min(loadings[,BC2.plot]),max(loadings[,BC2.plot])),
					xlab=paste0("BC ",BC1.plot," - Compound Loadings"), 
					ylab=paste0("BC ",BC2.plot," - Compound Loadings"),
					pch=21,
					bg="grey",
					col=groupCol,
					cex=1,main=paste0("FABIA - BC ",BC1.plot," vs ",BC2.plot)
			)
			text(loadings[,BC1.plot],loadings[,BC2.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#			if(length(legend.names)>0){legend(0.8,0.8,legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
			if(length(legend.names)>0){legend("topright",legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
			plot.out(plot.type)
		

			## PLOT: PC `BC1.plot' vs PC 'BC2.plot' - Factors
			plot.in(plot.type,paste0(basefilename,"_FABIA_PC",BC1.plot,"PC",BC2.plot,"_factors.pdf"))
			plot(scores[,BC1.plot],scores[,BC2.plot],  type="p",
					xlab=paste0("BC ",BC1.plot,"- Gene Factor Scores"), 
					ylab=paste0("BC ",BC2.plot," - Gene Factor Scores"),
					pch=21,
					bg="grey",
					col="grey",
					cex=1,main=paste0("FABIA - BC ",BC1.plot," vs ",BC2.plot)
			)
			text(scores[,BC1.plot],scores[,BC2.plot], rownames(data),	pos=1,	cex=0.5,	col="grey")
			plot.out(plot.type)
		}

		
		##
#		row.interest.list <- list()
		row.interest.list <- vector("list",max(BC.plot))
		##
		
		for(i.BC in BC.plot){
			if(4 %in% which){
				
				bg.temp <-  col.temp <- rep("grey",length(scores[,i.BC]))
				if(!is.null(gene.thresP)){
					temp.boolean <- (scores[,i.BC]>=gene.thresP)
					bg.temp[temp.boolean] <- thresP.col
					col.temp[temp.boolean] <- thresP.col
				}
				if(!is.null(gene.thresN)){					
					temp.boolean <- (scores[,i.BC]<=gene.thresN)
					bg.temp[temp.boolean] <- thresN.col
					col.temp[temp.boolean] <- thresN.col
				}
				
				# highlighting genes
				if(!is.null(gene.highlight)){
					if(class(gene.highlight)=="numeric" | class(gene.highlight)=="integer"){
						col.temp[gene.highlight] <- "green"
					}
					if(class(gene.highlight)=="list"){
						if(length(gene.highlight)>5){stop("Too many different gene.highlight",call.=FALSE)}
						col.pool <- c("green","deeppink","darkorchid4","gold3","tan1")
						for(i.list in 1:length(gene.highlight)){
							col.temp[gene.highlight[[i.list]]] <- col.pool[i.list]
						}
					}
				}
				
				
				## PLOT: Genes Factors of 'BC.plot'
				plot.in(plot.type,paste0(basefilename,"_FABIA_BC",i.BC,"Factors.pdf"))
				plot(c(1:length(scores[,i.BC])),scores[,i.BC],  type="p",
						xlab="Gene Index", 
						ylab="Gene Factor Scores",
						pch=21,
						bg=bg.temp,
						col=col.temp,
						cex=1,main=paste0("BC ",i.BC," - Gene Factor Scores")
				)
				text(c(1:length(scores[,i.BC])),scores[,i.BC], rownames(data),	pos=1,	cex=0.5,	col=col.temp)
				if(!is.null(gene.thresP)){abline(gene.thresP,0,lty=3)}
				if(!is.null(gene.thresN)){abline(gene.thresN,0,lty=3)}
				plot.out(plot.type)
			
			}
			
			# SELECT GENES FOR PROFILES PLOT
			if(6 %in% which & profile.type=="gene"){
				if(is.null(row.interest)){
					if(plot.type=="pdf" | !(4 %in% which)){
						bg.temp <-  col.temp <- rep("grey",length(scores[,i.BC]))
						if(!is.null(gene.thresP)){
							temp.boolean <- (scores[,i.BC]>=gene.thresP)
							bg.temp[temp.boolean] <- thresP.col
							col.temp[temp.boolean] <- thresP.col
						}
						if(!is.null(gene.thresN)){					
							temp.boolean <- (scores[,i.BC]<=gene.thresN)
							bg.temp[temp.boolean] <- thresN.col
							col.temp[temp.boolean] <- thresN.col
						}
						
						# highlighting genes
						if(!is.null(gene.highlight)){
							if(class(gene.highlight)=="numeric" | class(gene.highlight)=="integer"){
								col.temp[gene.highlight] <- "green"
							}
							if(class(gene.highlight)=="list"){
								if(length(gene.highlight)>5){stop("Too many different gene.highlight",call.=FALSE)}
								col.pool <- c("green","deeppink","darkorchid4","gold3","tan1")
								for(i.list in 1:length(gene.highlight)){
									col.temp[gene.highlight[[i.list]]] <- col.pool[i.list]
								}
							}
						}
						dev.new()
						plot(c(1:length(scores[,i.BC])),scores[,i.BC],  type="p",
								xlab="Gene Index", 
								ylab="Gene Factor Scores",
								pch=21,
								bg=bg.temp,
								col=col.temp,
								cex=1,main=paste0("BC ",i.BC," - Gene Factor Scores")
						)
						text(c(1:length(scores[,i.BC])),scores[,i.BC], rownames(data),	pos=1,	cex=0.5,	col=col.temp)
						if(!is.null(gene.thresP)){abline(gene.thresP,0,lty=3)}
						if(!is.null(gene.thresN)){abline(gene.thresN,0,lty=3)}
					}
					cat("Select as many genes as desired with left mouse button. Right-click to end selection procedure.\n\n")
					id.temp <- identify(c(1:length(scores[,i.BC])),scores[,i.BC],n=999,labels=rownames(data))
					if(!(length(id.temp)==0)){
						row.interest.list[[i.BC]] <- id.temp
					}

				}
				else{
					row.interest.list[[i.BC]] <- row.interest
				}
			}
			
		}
		
		
		

		
		
		##
#		column.interest.list <- list()
		column.interest.list <- vector("list",max(BC.plot))
		##
			
		for(i.BC in c(1:length(BC.plot))){
			
			if(5 %in% which){
				## PLOT: Compound Loadings of 'BC.plot'
				plot.in(plot.type,paste0(basefilename,"_FABIA_BC",BC.plot[i.BC],"Loadings.pdf"))
				par(mfrow=c(1,1))
				plot(c(1:length(loadings[,BC.plot[i.BC]])),loadings[,BC.plot[i.BC]],  type="p",
						xlab="Compound Index", 
						ylab="Compound Loadings",
						pch=21,
						bg="grey",
						col=groupCol,
						cex=1,main=paste0("Bicluster ",BC.plot[i.BC]," - Compound Loadings")
				)
				text(c(1:length(loadings[,BC.plot[i.BC]])),loadings[,BC.plot[i.BC]], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#				if(length(legend.names)>0){legend((length(loadings[,BC.plot[i.BC]])*(1-0.45)),max(loadings),legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
				if(length(legend.names)>0){legend("topright",legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
				
				plot.out(plot.type)
			}
			
			## SELECTING THE COMPOUNDS OF INTEREST (replot in device if necessary)
			if(6 %in% which){
				if(is.null(column.interest)){
					if(plot.type=="pdf" | !(5 %in% which)){
						dev.new()
						par(mfrow=c(1,1))
						plot(c(1:length(loadings[,BC.plot[i.BC]])),loadings[,BC.plot[i.BC]],  type="p",
								xlab="Compound Index", 
								ylab="Compound Loadings",
								pch=21,
								bg="grey",
								col=groupCol,
								cex=1,main=paste0("Bicluster ",BC.plot[i.BC]," - Compound Loadings")
						)
						text(c(1:length(loadings[,BC.plot[i.BC]])),loadings[,BC.plot[i.BC]], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#						if(length(legend.names)>0){legend((length(loadings[,BC.plot[i.BC]])*(1-0.45)),max(loadings),legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
						if(length(legend.names)>0){legend("topright",max(loadings),legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
						
					}
				
					cat("Select as many compounds as desired with left mouse button. Right-click to end selection procedure.\n\n")
					id.temp <- identify(c(1:length(loadings[,BC.plot[i.BC]])),loadings[,BC.plot[i.BC]],n=999,labels=colnames(data))
					if(!(length(id.temp)==0)){
						column.interest.list[[BC.plot[i.BC]]] <- id.temp
					}

				
				}
				else{
					column.interest.list[[BC.plot[i.BC]]] <- column.interest
					
				}
			}
					
		}
		
		if( (6 %in% which) & (length(column.interest.list)>0) ){# Excluding the case in which you don't select anything + no pm
						 
			if(!(profile.type=="gene" & !length(row.interest.list)>0)){ # case of gene profiles but no selection not allowed
				
				for(i.BC in BC.plot){
					cmpds_interest <- column.interest.list[[i.BC]]
					
					if(length(row.interest.list)==0){
						genes_interest <- vector("list",length(cmpds_interest))
					}
					else{
						genes_interest <- row.interest.list[[i.BC]]
						
					}
					
					if(!is.null(cmpds_interest)){
						if(!(profile.type=="gene" & is.null(genes_interest))){
						
							
							base.name <- paste0(basefilename,"_BC",i.BC)
							main.base <- paste0("Fabia BC",i.BC)
							
							CSprofiles(data=data,ref_index=ref.index,gene.select=genes_interest,cmpd.select=cmpds_interest,profile.type=profile.type,cmpd.loadings=loadings,gene.scores=scores,component.plot=i.BC,gene.thresP=gene.thresP,gene.thresN=gene.thresN,basefilename=base.name,plot.type=plot.type,thresP.col=thresP.col,thresN.col=thresN.col,main.base=main.base)
								
						}
						
						
					}
				}
				
				
			}

			
		}
		
		
		
		if(7%in%which){
#			out_CS_rank <- replicate(length(BC.plot),list)
			out_CS_rank <- vector("list",length(BC.plot))
			
			legend.bg <- "white"
			legend.names.csrank <- legend.names
			legend.cols.csrank <- legend.cols
			
			for(i.BC in c(1:length(BC.plot))){
				out_CS_rank[[i.BC]] <- CSrank2(loadings,ref.index,color.columns=groupCol,ref.plot=CSrank.refplot,loadings_names=colnames(data),component.plot=BC.plot[i.BC],type.component="BC",plot=TRUE,plot.type=plot.type,basefilename=basefilename,legend.bg=legend.bg,legend.names=legend.names.csrank,legend.cols=legend.cols.csrank)
				
			}
			names(out_CS_rank) <- paste0("BC",BC.plot)
			
		}
		else{
			out_CS_rank <- replicate(length(BC.plot),list)
					
			for(i.BC in c(1:length(BC.plot))){
				out_CS_rank[[i.BC]] <- CSrank2(loadings,ref.index,color.columns=groupCol,ref.plot=CSrank.refplot,loadings_names=colnames(data),component.plot=BC.plot[i.BC],type.component="BC",plot=FALSE)
				
			}
			names(out_CS_rank) <- paste0("BC",BC.plot)
			
		}
	
		
	out <- list(BC.select=BC.plot,result=resFAB,CSRank=out_CS_rank)	
	return(out)
}

#"bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". 


#######################################################################################################################################
###### FABIA vs PCA/MFA
## which:	-1: Compound Loadings Plot
##			-2: Gene Scores Plot



analyse_fabia_MFAPCA <- function(data,resFAB,resPCA=NULL,resMFA=NULL,ressparseMFA=NULL,resFAB2=NULL,
		ref.index,BC.plot,factor.plot,BC2.plot,
		# Give 2 Pos and 2 Neg threshold (1: FAB, 2:PCA/MFA)
		# Give 3 Pos and 3 Neg Colors (1:FAB, 2:PCA/MFA)
		gene.thresP=NULL,gene.thresN=NULL,thresP.col=c("blue","light blue"),thresN.col=c("red","pink"),
		colour.columns=NULL,legend.names=NULL,legend.cols=unique(colour.columns),legend.pos='topright',
		plot.type="pdf",basefilename="FABIA_PCAMFA",which=c(1,2)){
	
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
	
	## Check if not both results
	#if(!is.null(resPCA) & !is.null(resMFA)){stop("Cannot PCA and MFA result at the same time.",call.=FALSE)}
	if(sum(c(!is.null(resPCA),!is.null(resMFA),!is.null(ressparseMFA),!is.null(resFAB2)))>1){stop("Cannot Fabia2, PCA, MFA or sMFA result at the same time.",call.=FALSE)}
	
	
	##
	if(!is.null(colour.columns)){groupCol <- colour.columns} else { groupCol <- "black"}
	##	
	
	
	## PLOT: FABIA VS PCA 
	if(!is.null(resPCA)){
		
		if(1 %in% which){
			## Plot Compound Loadings
				
			loadingsPCA <- resPCA$var$coord
			loadingsFAB <- resFAB@L
			minX <- min(-1,loadingsFAB[,BC.plot])
			maxX <- max(1,loadingsFAB[,BC.plot])
			minY <- min(-1,loadingsPCA[,factor.plot])
			maxY <- max(1,loadingsPCA[,factor.plot])
			
			plot.in(plot.type,paste0(basefilename,"_loadingsPCA.pdf"))
			par(mfrow=c(1,1))
			plot(loadingsFAB[,BC.plot],loadingsPCA[,factor.plot],xlim=c(minX,maxX),ylim=c(minY,maxY),col=groupCol,bg="grey",main="FABIA VS PCA Loadings - 1 Ref Compound",xlab=paste0("Fabia BC ",BC.plot,"Loadings"),ylab=paste0("PCA PC",factor.plot," Loadings"),pch=21)
			text(loadingsFAB[,BC.plot],loadingsPCA[,factor.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
			#abline(a=0,b=1,lty=3)
			if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
			plot.out(plot.type)
			#print(cor(loadingsFAB[,BC.plot],loadingsPCA[,factor.plot]))
			cor.CS <- cor(loadingsFAB[,BC.plot],loadingsPCA[,factor.plot])
		}			
		
		if(2 %in% which){
			## Plot Gene Scores
	
			scoresPCA <- resPCA$ind$coord
			scoresFAB <- t(resFAB@Z)
			minX <- min(-1,scoresFAB[,BC.plot])
			maxX <- max(1,scoresFAB[,BC.plot])
			minY <- min(-1,scoresPCA[,factor.plot])
			maxY <- max(1,scoresPCA[,factor.plot])
			
			# gene colors
			if(!is.null(gene.thresP) | !is.null(gene.thresN)){
				if(is.null(gene.thresP)){gene.thresP <- c(99999,99999)}
				if(is.null(gene.thresN)){gene.thresN <- c(-99999,-99999)}
								
				scores <- rbind(scoresFAB[,BC.plot],scoresPCA[,factor.plot])
				list.scores <- as.list(as.data.frame(scores))
				#(x,P1,N1,P2,N2,col.P1,col.N1,col.P2,col.N2)
				list.colors <- lapply(X=list.scores,FUN=.give.gene.color,P1=gene.thresP[1],P2=gene.thresP[2],N1=gene.thresN[1],N2=gene.thresN[2],col.P1=thresP.col[1],col.P2=thresP.col[2],col.N1=thresN.col[1],col.N2=thresN.col[2])
				list.colors <- t(as.data.frame(list.colors))
				colnames(list.colors) <- c("bg","col")
			}
			else{
				list.colors <- matrix("grey",ncol=2,nrow=length(scoresPCA[,1])) 
			}
						
			plot.in(plot.type,paste0(basefilename,"_scoresPCA.pdf"))
			par(mfrow=c(1,1))
			plot(scoresFAB[,BC.plot],scoresPCA[,factor.plot],xlim=c(minX,maxX),ylim=c(minY,maxY),col=list.colors[,2],bg=list.colors[,1],main="FABIA VS PCA Scores - 1 Ref Compound",xlab=paste0("Fabia BC ",BC.plot,"Scores"),ylab=paste0("PCA PC",factor.plot,"Scores"),pch=21)
			text(scoresFAB[,BC.plot],scoresPCA[,factor.plot], rownames(data),	pos=1,	cex=0.5,	col=list.colors[,2])
			#abline(a=0,b=1,lty=3)
			if(!is.null(gene.thresP)){
				abline(v=gene.thresP[1],lty=3)
				abline(h=gene.thresP[2],lty=3)
			}
			if(!is.null(gene.thresN)){
				abline(v=gene.thresN[1],lty=3)
				abline(h=gene.thresN[2],lty=3)
			}
			#if(length(legend.names)>0){legend(legend.pos,pch=21,col=legend.cols,pt.bg="grey",bty="n")}
			plot.out(plot.type)
#			print(cor(scoresFAB[,BC.plot],scoresPCA[,factor.plot]))
			cor.GS <- cor(scoresFAB[,BC.plot],scoresPCA[,factor.plot])
			
		}
	}
	
	if(!is.null(resMFA)){
		
		if(1 %in% which){
			## Plot Compound Loadings
			
			loadingsMFA <- resMFA$quanti.var$coord
			loadingsFAB <- resFAB@L
			minX <- min(-1,loadingsFAB[,BC.plot])
			maxX <- max(1,loadingsFAB[,BC.plot])
			minY <- min(-1,loadingsMFA[,factor.plot])
			maxY <- max(1,loadingsMFA[,factor.plot])
			
			plot.in(plot.type,paste0(basefilename,"_loadingsMFA.pdf"))
			par(mfrow=c(1,1))
			plot(loadingsFAB[,BC.plot],loadingsMFA[,factor.plot],xlim=c(minX,maxX),ylim=c(minY,maxY),col=groupCol,bg="grey",main=paste0("FABIA VS MFA Loadings - ",length(ref.index)," Ref Compound"),xlab=paste0("Fabia BC ",BC.plot,"Loadings"),ylab=paste0("MFA PC",factor.plot," Loadings"),pch=21)
			text(loadingsFAB[,BC.plot],loadingsMFA[,factor.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
			#abline(a=0,b=1,lty=3)
			if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
			plot.out(plot.type)
			#print(cor(loadingsFAB[,BC.plot],loadingsMFA[,factor.plot]))
			cor.CS <- cor(loadingsFAB[,BC.plot],loadingsMFA[,factor.plot])
			
		}
		
		if(2 %in% which){
			
			## Plot Gene Scores
			
			scoresMFA <- resMFA$ind$coord	
			scoresFAB <- t(resFAB@Z)
			minX <- min(-1,scoresFAB[,BC.plot])
			maxX <- max(1,scoresFAB[,BC.plot])
			minY <- min(-1,scoresMFA[,factor.plot])
			maxY <- max(1,scoresMFA[,factor.plot])
			
			# gene colors
			if(!is.null(gene.thresP) | !is.null(gene.thresN)){
				if(is.null(gene.thresP)){gene.thresP <- c(99999,99999)}
				if(is.null(gene.thresN)){gene.thresN <- c(-99999,-99999)}
								
				scores <- rbind(scoresFAB[,BC.plot],scoresMFA[,factor.plot])
				list.scores <- as.list(as.data.frame(scores))
				#(x,P1,N1,P2,N2,col.P1,col.N1,col.P2,col.N2)
				list.colors <- lapply(X=list.scores,FUN=.give.gene.color,P1=gene.thresP[1],P2=gene.thresP[2],N1=gene.thresN[1],N2=gene.thresN[2],col.P1=thresP.col[1],col.P2=thresP.col[2],col.N1=thresN.col[1],col.N2=thresN.col[2])
				list.colors <- t(as.data.frame(list.colors))
				colnames(list.colors) <- c("bg","col")
			}
			else{
				list.colors <- matrix("grey",ncol=2,nrow=length(scoresMFA[,1])) 
			}
			
			plot.in(plot.type,paste0(basefilename,"_scoresMFA.pdf"))
			par(mfrow=c(1,1))
			plot(scoresFAB[,BC.plot],scoresMFA[,factor.plot],xlim=c(minX,maxX),ylim=c(minY,maxY),col=list.colors[,2],bg=list.colors[,1],main=paste0("FABIA VS MFA Scores - ",length(ref.index)," Ref Compound"),xlab=paste0("Fabia BC ",BC.plot,"Scores"),ylab=paste0("MFA PC",factor.plot,"Scores"),pch=21)
			text(scoresFAB[,BC.plot],scoresMFA[,factor.plot], rownames(data),	pos=1,	cex=0.5,	col=list.colors[,2])
			#abline(a=0,b=1,lty=3)
			if(!is.null(gene.thresP)){
				abline(v=gene.thresP[1],lty=3)
				abline(h=gene.thresP[2],lty=3)
			}
			if(!is.null(gene.thresN)){
				abline(v=gene.thresN[1],lty=3)
				abline(h=gene.thresN[2],lty=3)
			}
			#if(length(legend.names)>0){legend(legend.pos,pch=21,col=legend.cols,pt.bg="grey",bty="n")}
			plot.out(plot.type)
			#print(cor(scoresFAB[,BC.plot],scoresMFA[,factor.plot]))
			cor.GS <- cor(scoresFAB[,BC.plot],scoresMFA[,factor.plot])
			
		}
		
	}
	
	if(!is.null(ressparseMFA)){
		if(1 %in% which){
			## Plot Compound Loadings
			
			loadingsMFA <- ressparseMFA$loadings
			loadingsFAB <- resFAB@L
			minX <- min(-1,loadingsFAB[,BC.plot])
			maxX <- max(1,loadingsFAB[,BC.plot])
			minY <- min(-1,loadingsMFA[,factor.plot])
			maxY <- max(1,loadingsMFA[,factor.plot])
			
			plot.in(plot.type,paste0(basefilename,"_loadingssMFA.pdf"))
			par(mfrow=c(1,1))
			plot(loadingsFAB[,BC.plot],loadingsMFA[,factor.plot],xlim=c(minX,maxX),ylim=c(minY,maxY),col=groupCol,bg="grey",main=paste0("FABIA VS sMFA Loadings - ",length(ref.index)," Ref Compound"),xlab=paste0("Fabia BC ",BC.plot,"Loadings"),ylab=paste0("sMFA PC",factor.plot," Loadings"),pch=21)
			text(loadingsFAB[,BC.plot],loadingsMFA[,factor.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
			#abline(a=0,b=1,lty=3)
			if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
			plot.out(plot.type)
			#print(cor(loadingsFAB[,BC.plot],loadingsMFA[,factor.plot]))
			cor.CS <- cor(loadingsFAB[,BC.plot],loadingsMFA[,factor.plot])
			
		}
		
		if(2 %in% which){
			## Weighted Data
			total.index <- c(1:dim(data)[2])
			Mat1 <- data[,ref.index,drop=FALSE]
			Mat2 <- data[,total.index[-ref.index]]
			
			data2 <- getWeightedDat(Mat1,Mat2,scale.unit=TRUE)[[1]]
			
			## Plot Gene Scores
			
			scoresMFA <- data2 %*% loadingsMFA
			scoresFAB <- t(resFAB@Z)
			minX <- min(-1,scoresFAB[,BC.plot])
			maxX <- max(1,scoresFAB[,BC.plot])
			minY <- min(-1,scoresMFA[,factor.plot])
			maxY <- max(1,scoresMFA[,factor.plot])
			
			# gene colors
			if(!is.null(gene.thresP) | !is.null(gene.thresN)){
				if(is.null(gene.thresP)){gene.thresP <- c(99999,99999)}
				if(is.null(gene.thresN)){gene.thresN <- c(-99999,-99999)}
				
				scores <- rbind(scoresFAB[,BC.plot],scoresMFA[,factor.plot])
				list.scores <- as.list(as.data.frame(scores))
				#(x,P1,N1,P2,N2,col.P1,col.N1,col.P2,col.N2)
				list.colors <- lapply(X=list.scores,FUN=.give.gene.color,P1=gene.thresP[1],P2=gene.thresP[2],N1=gene.thresN[1],N2=gene.thresN[2],col.P1=thresP.col[1],col.P2=thresP.col[2],col.N1=thresN.col[1],col.N2=thresN.col[2])
				list.colors <- t(as.data.frame(list.colors))
				colnames(list.colors) <- c("bg","col")
			}
			else{
				list.colors <- matrix("grey",ncol=2,nrow=length(scoresMFA[,1])) 
			}
			
			plot.in(plot.type,paste0(basefilename,"_scoressMFA.pdf"))
			par(mfrow=c(1,1))
			plot(scoresFAB[,BC.plot],scoresMFA[,factor.plot],xlim=c(minX,maxX),ylim=c(minY,maxY),col=list.colors[,2],bg=list.colors[,1],main=paste0("FABIA VS sMFA Scores - ",length(ref.index)," Ref Compound"),xlab=paste0("Fabia BC ",BC.plot,"Scores"),ylab=paste0("sMFA PC",factor.plot,"Scores"),pch=21)
			text(scoresFAB[,BC.plot],scoresMFA[,factor.plot], rownames(data),	pos=1,	cex=0.5,	col=list.colors[,2])
			#abline(a=0,b=1,lty=3)
			if(!is.null(gene.thresP)){
				abline(v=gene.thresP[1],lty=3)
				abline(h=gene.thresP[2],lty=3)
			}
			if(!is.null(gene.thresN)){
				abline(v=gene.thresN[1],lty=3)
				abline(h=gene.thresN[2],lty=3)
			}
			#if(length(legend.names)>0){legend(legend.pos,pch=21,col=legend.cols,pt.bg="grey",bty="n")}
			plot.out(plot.type)
			#print(cor(scoresFAB[,BC.plot],scoresMFA[,factor.plot]))
			cor.GS <- cor(scoresFAB[,BC.plot],scoresMFA[,factor.plot])
			
		}
		
	}
	
	if(!is.null(resFAB2)){
		if(1 %in% which){
			## Plot Compound Loadings
			
			loadingsFAB2 <- resFAB2@L
			loadingsFAB <- resFAB@L
			minX <- min(-1,loadingsFAB[,BC.plot])
			maxX <- max(1,loadingsFAB[,BC.plot])
			minY <- min(-1,loadingsFAB2[,BC2.plot])
			maxY <- max(1,loadingsFAB2[,BC2.plot])
			
			plot.in(plot.type,paste0(basefilename,"_loadingsFAB2.pdf"))
			par(mfrow=c(1,1))
			plot(loadingsFAB[,BC.plot],loadingsFAB2[,BC2.plot],xlim=c(minX,maxX),ylim=c(minY,maxY),col=groupCol,bg="grey",main=paste0("FABIA VS FABIA2 Loadings - ",length(ref.index)," Ref Compound"),xlab=paste0("Fabia BC ",BC.plot,"Loadings"),ylab=paste0("Fabia 2 BC",BC2.plot," Loadings"),pch=21)
			text(loadingsFAB[,BC.plot],loadingsFAB2[,BC2.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
			#abline(a=0,b=1,lty=3)
			if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
			plot.out(plot.type)
			#print(cor(loadingsFAB[,BC.plot],loadingsFAB2[,BC2.plot]))
			cor.CS <- cor(loadingsFAB[,BC.plot],loadingsFAB2[,BC2.plot])
		}
		
		if(2 %in% which){
			
			## Plot Gene Scores
			
			scoresFAB2 <- t(resFAB2@Z)
			scoresFAB <- t(resFAB@Z)
			minX <- min(-1,scoresFAB[,BC.plot])
			maxX <- max(1,scoresFAB[,BC.plot])
			minY <- min(-1,scoresFAB2[,BC2.plot])
			maxY <- max(1,scoresFAB2[,BC2.plot])
			
			# gene colors
			if(!is.null(gene.thresP) | !is.null(gene.thresN)){
				if(is.null(gene.thresP)){gene.thresP <- c(99999,99999)}
				if(is.null(gene.thresN)){gene.thresN <- c(-99999,-99999)}
				
				scores <- rbind(scoresFAB[,BC.plot],scoresFAB2[,BC2.plot])
				list.scores <- as.list(as.data.frame(scores))
				#(x,P1,N1,P2,N2,col.P1,col.N1,col.P2,col.N2)
				list.colors <- lapply(X=list.scores,FUN=.give.gene.color,P1=gene.thresP[1],P2=gene.thresP[2],N1=gene.thresN[1],N2=gene.thresN[2],col.P1=thresP.col[1],col.P2=thresP.col[2],col.N1=thresN.col[1],col.N2=thresN.col[2])
				list.colors <- t(as.data.frame(list.colors))
				colnames(list.colors) <- c("bg","col")
			}
			else{
				list.colors <- matrix("grey",ncol=2,nrow=length(scoresFAB2[,1])) 
			}
			
			plot.in(plot.type,paste0(basefilename,"_scoresFAB2.pdf"))
			par(mfrow=c(1,1))
			plot(scoresFAB[,BC.plot],scoresFAB2[,BC2.plot],xlim=c(minX,maxX),ylim=c(minY,maxY),col=list.colors[,2],bg=list.colors[,1],main=paste0("FABIA VS FABIA 2 Scores - ",length(ref.index)," Ref Compound"),xlab=paste0("Fabia BC ",BC.plot,"Scores"),ylab=paste0("Fabia BC ",BC2.plot,"Scores"),pch=21)
			text(scoresFAB[,BC.plot],scoresFAB2[,BC2.plot], rownames(data),	pos=1,	cex=0.5,	col=list.colors[,2])
			#abline(a=0,b=1,lty=3)
			if(!is.null(gene.thresP)){
				abline(v=gene.thresP[1],lty=3)
				abline(h=gene.thresP[2],lty=3)
			}
			if(!is.null(gene.thresN)){
				abline(v=gene.thresN[1],lty=3)
				abline(h=gene.thresN[2],lty=3)
			}
			#if(length(legend.names)>0){legend(legend.pos,pch=21,col=legend.cols,pt.bg="grey",bty="n")}
			plot.out(plot.type)
			#print(cor(scoresFAB[,BC.plot],scoresFAB2[,BC2.plot]))
			cor.GS <- cor(scoresFAB[,BC.plot],scoresFAB2[,BC2.plot])
		}
		
	}
	
	cor.temp <- c(cor.CS,cor.GS)
	names(cor.temp) <- c("CS.cor","GS.cor")
	return(cor.temp)
}

#analyse_fabia_MFAPCA <- function(data,resFAB,resPCA=NULL,resMFA=NULL,
#		ref.index,BC.plot,factor.plot,
#		# Give 2 Pos and 2 Neg threshold (1: FAB, 2:PCA/MFA)
#		# Give 3 Pos and 3 Neg Colors (1:FAB, 2:PCA/MFA)
#		gene.thresP=NULL,gene.thresN=NULL,thresP.col=c("blue","light blue"),thresN.col=c("red","pink"),
#		colour.columns=NULL,legend.names=NULL,legend.cols=unique(colour.columns),legend.pos='topright',
#		plot.type="pdf",basefilename="FABIA_PCAMFA",which=c(1,2)){

#analyse_fabia_MFAPCA(dataMFA,out2,resPCA=out,ref.index=2,BC.plot=1,factor.plot=1
#,gene.thresP=c(1.5,5)
#,gene.thresN=c(-1.5,-5)
#
#,colour.columns=groupCol,legend.names=c("ref","SP","WP","SN","c"),legend.cols=c("blue","green","red","cyan","black")
#
#
#,plot.type="device")


analyse_Zhang_fabia <- function(data,resZhang,resFAB=NULL,
		ref.index,BC.plot,
		colour.columns=NULL,legend.names=NULL,legend.cols=unique(colour.columns),legend.pos="topright",
		plot.type="pdf",basefilename="Zhang_fabia"){
	

	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	
	## Adding Zhang scores for the references. They will all get 1! (Even if multiple references were used)
	zhang.temp <- rep(0,dim(data)[2])
	i.zhang <- 1
	for(ii in 1:dim(data)[2]){
		if(!(ii %in% ref.index)){
			zhang.temp[ii] <- resZhang$All[i.zhang,1]
			i.zhang <- i.zhang+1
		}
		else{
			zhang.temp[ii] <- 1
		}
	}
	
	##
	if(!is.null(colour.columns)){groupCol <- colour.columns} else { groupCol <- "black"}
	##	
	
	## PLOT: Zhang Scores VS FABIA Loadings (Ref will get Zhang Score=1)
	loadings <- resFAB@L	# For compounds	

	minX <- min(-1,zhang.temp,na.rm=TRUE)
	maxX <- max(1,zhang.temp,na.rm=TRUE)
	minY <- min(-1,loadings[,BC.plot])
	maxY <- max(1,loadings[,BC.plot])
	
	plot.in(plot.type,paste0(basefilename,"_loadingsFABIA.pdf"))

	par(mfrow=c(1,1))
	plot(zhang.temp,loadings[,BC.plot],xlim=c(minX,maxX),ylim=c(minY,maxY),col=groupCol,bg="grey",main=paste0("Zhang VS FABIA - ",length(ref.index)," Ref Compound"),xlab="Zhang Score",ylab=paste0("BC",BC.plot," Loadings"),pch=21)
	text(zhang.temp,loadings[,BC.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#	abline(a=0,b=1,lty=3)
	if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
	
	plot.out(plot.type)
	
	#print(cor(zhang.temp,loadings[,BC.plot]))
	cor.temp <- cor(zhang.temp,loadings[,BC.plot],use="complete.obs")
	names(cor.temp) <- "CS.cor"
	return(cor.temp)
}


#######################################################################################################################################
#######################################################################################################################################
###### sMFA Analysis  (For multiple ref compounds)
## which:	-1: Percentage Variance Explained by factors
##			-2: Loadings for reference compounds
##			-3: Factor 'factor.plot' VS Factor '?' : Loadings & Genes
##			-4: Loadings & Genes for Factor 'factor.plot'
##			-5: Compound Profiles (Select if necessary)
##			-6: CS Rank Scores for 'factor.plot'


analyse_sMFA <- function(data,K=15,para,type=c("predictor","Gram"),sparse=c("penalty","varnum"),use.corr=FALSE,lambda=1e-6,max.iter=200,trace=FALSE,eps.conv=1e-3,
		basefilename="analyseMFA",ref.index=c(1),sparse.dim=2,
		factor.plot=1,column.interest=NULL,gene.thresP=NULL,gene.thresN=NULL,
		colour.columns=NULL,legend.names=NULL,legend.cols=unique(colour.columns),thresP.col="blue",thresN.col="red",
		CSrank.refplot=FALSE,gene.highlight=NULL,profile.type="gene",
		result.available=NULL,plot.type="pdf",row.interest=NULL,
		which=c(1,2,3,4,5)){
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	## Weighting the data
	total.index <- c(1:dim(data)[2])
	Mat1 <- data[,ref.index,drop=FALSE]
	Mat2 <- data[,total.index[-ref.index]]
	
	data <- getWeightedDat(Mat1,Mat2,scale.unit=TRUE)[[1]]
	
	## Doing sMFA Analysis
	if(is.null(result.available)){
		
		if(!(sparse.dim %in% c(1,2))){stop("Please use a correct sparse.dim")}
		
		if(sparse.dim==1){data.spca <- t(data)}else{data.spca <- data}
		
		if(lambda==Inf){
			if(dim(data.spca)[2]<dim(data.spca)[1]){warning("The to-be-reduced-with-sparsness dimension  is larger than the other dimension. Consider putting lambda at 0.")}
			
			
			resMFA <- arrayspc(x=data.spca,K=K,para=para,use.corr=use.corr,max.iter=max.iter,trace=trace,eps=eps.conv)
			
		}
		else{
			if(dim(data.spca)[2]>dim(data.spca)[1]){warning("The to-be-reduced-with-sparsness dimension  is larger than the other dimension. Consider putting lambda at Inf.")}
			
			resMFA <- spca(x=data.spca,K=K,para=para,type=type,sparse=sparse,use.corr=use.corr,lambda=lambda,max.iter=max.iter,trace=trace,eps.conv=eps.conv)
			
		}
		
	}
	else{
		resMFA <- result.available@object
	}
	

	if(1 %in% which){
		## PLOT: Amount of variance explained by PC's
		perc.var <- resMFA$pev				
		plot.in(plot.type,paste0(basefilename,"_sMFApercvar.pdf"))
		par(mfrow=c(1,1))
		plot(perc.var,main="Percentage of Variance Explained",xlab="Number of Components",ylab="Perc. Var. Explained")
		plot.out(plot.type)
	}
	
	##
	
	if(sparse.dim==2){
		loadings <- resMFA$loadings	# For compounds	
		scores <- data %*% loadings	# For genes
		
	}
	else if(sparse.dim==1){
		scores <- resMFA$loadings	# For genes	
		loadings <- data.spca %*% scores	# For compounds
	}
	
	##
	
	if(!is.null(colour.columns)){groupCol <- colour.columns} else { groupCol <- "black"}
	
	##
	
	if(2 %in% which){
		## PLOT: Loadings for reference compounds -> select which PC to use if factor.plot is NOT given
		plot.in(plot.type,paste0(basefilename,"_sMFA_RefLoadings"))
		plot(0,0,type="n",main=paste0("Loadings for Ref ",paste(ref.index,collapse=",")),xlab="Factor Index",ylab="Loadings",ylim=c(min(loadings[ref.index,]),max(loadings[ref.index,])),xlim=c(1,dim(loadings)[2]))
		for(i.ref in ref.index){
			points(c(1:dim(loadings)[2]),loadings[i.ref,],col=i.ref)
		}
		abline(0,0,lty=3)
		legend(dim(loadings)[2]*(1-0.4),max(loadings[ref.index,]),colnames(data)[ref.index],col=ref.index,bty="n",pch=21)
		plot.out(plot.type)
	}
	
	
	## Selecting the factor.plot
	if(is.null(factor.plot)|(3 %in% which) | (4 %in% which) | (5 %in% which) | (is.null(column.interest)&(6 %in% which)) | ((6 %in% which) & ( (!is.null(gene.thresP))   | (!is.null(gene.thresN)) ) ) | (7 %in% which) ){
		if(is.null(factor.plot)){
			if(plot.type=="pdf" | !(2 %in% which)){
				dev.new()
				plot(0,0,type="n",main=paste0("Loadings for Ref ",paste(ref.index,collapse=",")),xlab="Factor Index",ylab="Loadings",ylim=c(min(loadings[ref.index,]),max(loadings[ref.index,])),xlim=c(1,dim(loadings)[2]))
				for(i.ref in ref.index){
					points(c(1:dim(loadings)[2]),loadings[i.ref,],col=i.ref)
				}
				abline(0,0,lty=3)
				legend(dim(loadings)[2]*(1-0.4),max(loadings[ref.index,]),colnames(data)[ref.index],col=ref.index,bty="n",pch=21)
			}
			cat("Please select with left mousebutton which Factor should be investigated. (The Factor with highest loading for reference)\nIf multiple reference were used, click in the middle of the group of points.\n")
			y.mean <- apply(loadings[ref.index,,drop=FALSE],MARGIN=2,FUN=mean)
			factor.plot <- identify(x=c(1:dim(loadings)[2]),y=y.mean,plot=TRUE,n=1,tolerance=0.5)
		}
		
		##
		factor2.plot <- ifelse(factor.plot==1,2,1)
		##
	}
	
	
	if(3 %in% which){
		## PLOT: PC 'factor.plot'  vs PC 'factor2.plot' - Loadings
		plot.in(plot.type,paste0(basefilename,"_sMFA_PC",factor.plot,"PC",factor2.plot,"_loadings.pdf"))
		par(mfrow=c(1,1))
		plot(loadings[,factor.plot],loadings[,factor2.plot],  type="p",xlim=c(min(loadings[,factor.plot]),max(loadings[,factor.plot])),ylim=c(min(loadings[,factor2.plot]),max(loadings[,factor2.plot])),
				xlab=paste0("Factor ",factor.plot," - Compound Loadings"), 
				ylab=paste0("Factor ",factor2.plot," - Compound Loadings"),
				pch=21,
				bg="grey",
				col=groupCol,
				cex=1,main=paste0("sparse PCA weighted (sMFA) - Factor ",factor.plot," vs ",factor2.plot)
		)
		text(loadings[,factor.plot],loadings[,factor2.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#		if(length(legend.names)>0){legend(0.8,0.8,legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
		if(length(legend.names)>0){legend("topright",legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
		
		plot.out(plot.type)
		
		
		## PLOT: PC `factor.plot' vs PC 'factor2.plot' - Factors
		plot.in(plot.type,paste0(basefilename,"_sMFA_PC",factor.plot,"PC",factor2.plot,"_factors.pdf"))
		plot(scores[,factor.plot],scores[,factor2.plot],  type="p",
				xlab=paste0("Factor ",factor.plot,"- Gene Factor Scores"), 
				ylab=paste0("Factor ",factor2.plot," - Gene Factor Scores"),
				pch=21,
				bg="grey",
				col="grey",
				cex=1,main=paste0("sparse PCA weighted (sMFA) - Factor ",factor.plot," vs ",factor2.plot)
		)
		text(scores[,factor.plot],scores[,factor2.plot], rownames(data),	pos=1,	cex=0.5,	col="grey")
		plot.out(plot.type)
	}
	
	
	if(4 %in% which){
		## PLOT: Genes Factors of 'factor.plot'
		bg.temp <-  col.temp <- rep("grey",length(scores[,factor.plot]))
		if(!is.null(gene.thresP)){
			temp.boolean <- (scores[,factor.plot]>=gene.thresP)
			bg.temp[temp.boolean] <- thresP.col
			col.temp[temp.boolean] <- thresP.col
		}
		if(!is.null(gene.thresN)){
			temp.boolean <- (scores[,factor.plot]<=gene.thresN)
			bg.temp[temp.boolean] <- thresN.col
			col.temp[temp.boolean] <- thresN.col
		}
		
		# highlighting genes
		if(!is.null(gene.highlight)){
			if(class(gene.highlight)=="numeric" | class(gene.highlight)=="integer"){
				col.temp[gene.highlight] <- "green"
			}
			if(class(gene.highlight)=="list"){
				if(length(gene.highlight)>5){stop("Too many different gene.highlight",call.=FALSE)}
				col.pool <- c("green","deeppink","darkorchid4","gold3","tan1")
				for(i.list in 1:length(gene.highlight)){
					col.temp[gene.highlight[[i.list]]] <- col.pool[i.list]
				}
			}
		}
		
		plot.in(plot.type,paste0(basefilename,"_sMFAFactors.pdf"))
		plot(c(1:length(scores[,factor.plot])),scores[,factor.plot],  type="p",
				xlab="Gene Index", 
				ylab="Gene Factor Scores",
				pch=21,
				bg=bg.temp,
				col=col.temp,
				cex=1,main=paste0("Factor ",factor.plot," - Gene Factor Scores")
		)
		text(c(1:length(scores[,factor.plot])),scores[,factor.plot], rownames(data),	pos=1,	cex=0.5,	col=col.temp)
		if(!is.null(gene.thresP)){abline(gene.thresP,0,lty=3)}
		if(!is.null(gene.thresN)){abline(gene.thresN,0,lty=3)}
		plot.out(plot.type)
	}
	genes_interest <- row.interest #make the object in case of a cmpd profile
	
	## SELECTING THE GENES OF INTEREST (replot in device if necessary)
	if(6 %in% which & profile.type=="gene"){
		if(is.null(row.interest)){
			if(plot.type=="pdf" | !(4 %in% which)){
				dev.new()
				bg.temp <-  col.temp <- rep("grey",length(scores[,factor.plot]))
				if(!is.null(gene.thresP)){
					temp.boolean <- (scores[,factor.plot]>=gene.thresP)
					bg.temp[temp.boolean] <- thresP.col
					col.temp[temp.boolean] <- thresP.col
				}
				if(!is.null(gene.thresN)){
					temp.boolean <- (scores[,factor.plot]<=gene.thresN)
					bg.temp[temp.boolean] <- thresN.col
					col.temp[temp.boolean] <- thresN.col
				}
				
				# highlighting genes
				if(!is.null(gene.highlight)){
					if(class(gene.highlight)=="numeric" | class(gene.highlight)=="integer"){
						col.temp[gene.highlight] <- "green"
					}
					if(class(gene.highlight)=="list"){
						if(length(gene.highlight)>5){stop("Too many different gene.highlight",call.=FALSE)}
						col.pool <- c("green","deeppink","darkorchid4","gold3","tan1")
						for(i.list in 1:length(gene.highlight)){
							col.temp[gene.highlight[[i.list]]] <- col.pool[i.list]
						}
					}
				}
				plot(c(1:length(scores[,factor.plot])),scores[,factor.plot],  type="p",
						xlab="Gene Index", 
						ylab="Gene Factor Scores",
						pch=21,
						bg=bg.temp,
						col=col.temp,
						cex=1,main=paste0("Factor ",factor.plot," - Gene Factor Scores")
				)
				text(c(1:length(scores[,factor.plot])),scores[,factor.plot], rownames(data),	pos=1,	cex=0.5,	col=col.temp)
				if(!is.null(gene.thresP)){abline(gene.thresP,0,lty=3)}
				if(!is.null(gene.thresN)){abline(gene.thresN,0,lty=3)}
			}
			cat("Select as many genes as desired with left mouse button. Right-click to end selection procedure.\n\n")
			genes_interest <- identify(c(1:length(scores[,factor.plot])),scores[,factor.plot],n=999,plot=TRUE,labels=rownames(data))
			
		}
		else{
			genes_interest <- row.interest
		}
	}
	
	if(5 %in% which){
		
		## PLOT: Compound Loadings of 'factor.plot'
		plot.in(plot.type,paste0(basefilename,"_sMFALoadings.pdf"))
		par(mfrow=c(1,1))
		plot(c(1:length(loadings[,factor.plot])),loadings[,factor.plot],  type="p",
				xlab="Compound Index", 
				ylab="Compound Loadings",
				pch=21,
				bg="grey",
				col=groupCol,
				cex=1,main=paste0("Factor ",factor.plot," - Compound Loadings")
		)
		text(c(1:length(loadings[,factor.plot])),loadings[,factor.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#		if(length(legend.names)>0){legend((length(loadings[,factor.plot])*(1-0.45)),max(loadings),legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
		if(length(legend.names)>0){legend("topright",legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
		
		plot.out(plot.type)
	}
	
	## SELECTING THE COMPOUNDS OF INTEREST (replot in device if necessary)
	if(6 %in% which){
		if(is.null(column.interest)){
			if(plot.type=="pdf" | !(5 %in% which)){
				## PLOT: Compound Loadings of 'factor.plot'
				dev.new()
				par(mfrow=c(1,1))
				plot(c(1:length(loadings[,factor.plot])),loadings[,factor.plot],  type="p",
						xlab="Compound Index", 
						ylab="Compound Loadings",
						pch=21,
						bg="grey",
						col=groupCol,
						cex=1,main=paste0("Factor ",factor.plot," - Compound Loadings")
				)
				text(c(1:length(loadings[,factor.plot])),loadings[,factor.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
#				if(length(legend.names)>0){legend((length(loadings[,factor.plot])*(1-0.45)),max(loadings),legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
				if(length(legend.names)>0){legend("topright",legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
				
			}
			cat("Select as many compounds as desired with left mouse button. Right-click to end selection procedure.\n\n")
			cmpds_interest <- identify(c(1:length(loadings[,factor.plot])),loadings[,factor.plot],n=999,plot=TRUE,labels=colnames(data))
		}
		else{
			cmpds_interest <- column.interest
		}
	}
	
	if( 6 %in% which){
		## PLOT: Profiles Plot
		CSprofiles(data=data,ref_index=ref.index,gene.select=genes_interest,cmpd.select=cmpds_interest,profile.type=profile.type,
				cmpd.loadings=loadings,gene.scores=scores,component.plot=factor.plot,gene.thresP=gene.thresP,gene.thresN=gene.thresN,
				basefilename=basefilename,plot.type=plot.type,thresP.col=thresP.col,thresN.col=thresN.col,main.base=paste0("sMFA Factor ",factor.plot))
		
		
	}
	
	if(7 %in% which){
		legend.bg <- "white"
		legend.names.csrank <- legend.names
		legend.cols.csrank <- legend.cols

		
		out_CS_rank <- list(CSrank2(loadings,ref.index,color.columns=groupCol,ref.plot=CSrank.refplot,loadings_names=colnames(data),component.plot=factor.plot,type.component="Factor",plot=TRUE,plot.type=plot.type,basefilename=basefilename,legend.bg=legend.bg,legend.names=legend.names.csrank,legend.cols=legend.cols.csrank))
		names(out_CS_rank) <- paste0("Factor",factor.plot)
		
		
	}
	else{
		out_CS_rank <- list(CSrank2(loadings,ref.index,color.columns=groupCol,ref.plot=CSrank.refplot,loadings_names=colnames(data),component.plot=factor.plot,type.component="Factor",plot=FALSE))
		names(out_CS_rank) <- paste0("Factor",factor.plot)
		
	}
	
	
	
	## Returning the MFA result for further use
	
	resMFA$loadings <- loadings
	resMFA$scores <- scores
	out <- list(factor.select=factor.plot,result=resMFA,CSRank=out_CS_rank)	
	return(out)
	
	
}



CSrank <- function(loadings,ref.index,color.columns=NULL,ref.plot=FALSE,loadings_names=NULL,component.plot,type.component="Factor",plot=TRUE,plot.type="pdf",basefilename="base"){
	
	type <- "loadings"
#	type <- "loadings_abs1"
#	type <- "loadings_abs2"
#	type <- "loadings_0abs"
	
	
#	type <- "max"
#	type <- "median"
#	type <- "contributions"

	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	loadings <- loadings[,component.plot]
	names <- loadings_names
	
	
	SingleCS.matrix <- matrix(0,nrow=length(loadings[-ref.index]),ncol=length(ref.index))
	colnames(SingleCS.matrix) <- rep("",length(ref.index))
	
	if(type=="loadings"){
		avg_query2 <- mean(loadings[-ref.index]^2)
		
	}
	if(type=="loadings_abs"){
		avg_query2 <- mean(abs(loadings[-ref.index]))
	}
	
	if(type != "loadings_0abs"){
		if(type=="loadings"){refL <- loadings[ref.index]^2}
		if(type=="loadings_abs"){refL <- abs(loadings[ref.index])}
		if(sum(refL <= avg_query2)>=1){
			warning("No CS Rank Scores computed. One or more of the references is too low as a loading.")
#		return(NULL)
			return(NA)
		}
	}

	
	for(i.ref in 1:length(ref.index)){
		
		ref.loading <- loadings[ref.index[i.ref]]
		temp <- sapply(loadings[-ref.index],FUN=function(x){
											
					sign.temp <- ifelse(sign(ref.loading)==sign(x),1,-1)
					
					if(type=="loadings"){
						reflarger <- (ref.loading^2 >= x^2)
						
					}
					
					if(type=="loadings_abs" | type=="loadings_0abs"){
						reflarger <- (abs(ref.loading) >= abs(x))
					}
					
					if(reflarger){
#						return(
								if(type=="loadings"){return(((x^2 - avg_query2)^2/(ref.loading^2 - avg_query2)^2)*sign.temp)}
								if(type=="loadings_abs"){return(
											abs(abs(x)-avg_query2)/abs(abs(ref.loading)-avg_query2) * sign.temp
											)}
								if(type=="loadings_0abs"){return(
											abs(x)/abs(ref.loading) * sign.temp
											
											)}
								
								if(type=="loadings_abs"){return(((x^2 - avg_query2^2)^2/(ref.loading^2 - avg_query2^2)^2)*sign.temp)}
#								((x^2 - avg_query2)^2/(ref.loading^2 - avg_query2)^2)*sign.temp
#								((x^2 - avg_query2^2)^2/(ref.loading^2 - avg_query2^2)^2)*sign.temp
#								((x^2 - avg_query2^2)^2/(ref.loading^2 - avg_query2^2)^2)*sign.temp
								
#							)
					}
					else{
#						return(
								if(type=="loadings"){return(((ref.loading^2 - avg_query2)^2/(x^2 - avg_query2)^2)*sign.temp)}
								if(type=="loadings_abs"){return(
											abs(abs(ref.loading)-avg_query2)/abs(abs(x)-avg_query2) * sign.temp
											)}
								if(type=="loadings_0abs"){return(
											abs(ref.loading)/abs(x) * sign.temp
											)}
								
								if(type=="loadings_abs"){return(((ref.loading^2 - avg_query2^2)^2/(x^2 - avg_query2^2)^2)*sign.temp)}
#								((ref.loading^2 - avg_query2)^2/(x^2 - avg_query2)^2)*sign.temp
#								((ref.loading^2 - avg_query2^2)^2/(x^2 - avg_query2^2)^2)*sign.temp
#								((x^2 - avg_query2^2)^2/(ref.loading^2 - avg_query2^2)^2)*sign.temp
								
#							)
					}
					
					
			})
		SingleCS.matrix[,i.ref] <- temp
		colnames(SingleCS.matrix)[i.ref] <- paste0(names[ref.index[i.ref]]," (Ref ",i.ref,")")
	}
	
	if(type=="loadings"){weights <- loadings[ref.index]^2}
	if(type=="loadings_0abs"){weights <- abs(loadings[ref.index])}
	if(type=="loadings_abs"){weights <- abs(loadings[ref.index])}
	
	SingleCS.vector <- apply(SingleCS.matrix,MARGIN=1,FUN=weighted.mean,w=weights)
	
	if(plot){
		
		plot.in(plot.type,paste0(basefilename,"_CSRank.pdf"))
		if(ref.plot){
			col.plots <- ifelse((length(ref.index)+1)%%2 == 0, (length(ref.index)+1)%/%2, (length(ref.index)+1)%/%2 +1 )
			par(mfrow=c(2,col.plots))
			#dev.new()
			# Plot for each reference
			for(i.ref in 1:length(ref.index)){
				plot(1:dim(SingleCS.matrix)[1],SingleCS.matrix[,i.ref],xlab="Compound Index",ylab="CS Rankscore",col=color.columns[-ref.index],bg="grey",pch=21,main=paste0("Reference ",i.ref," (Cmpd ",names[ref.index[i.ref]],")"))
				text(1:dim(SingleCS.matrix)[1],SingleCS.matrix[,i.ref],names[-ref.index],col=color.columns[-ref.index],pos=2)
			}
			
			
		}
		else{
			par(mfrow=c(1,1))
			#dev.new()
		}
			
		## Plot weighted for references
		main.temp <- ifelse(ref.plot,"CS Rank Score (Weighted Mean)","CS Rank Score")
		plot(SingleCS.vector,col=color.columns[-ref.index],xlab="Compound Index",ylab="CS Rankscore",bg="grey",pch=21,main=paste0(type.component," ",component.plot," - ",main.temp))
		text(SingleCS.vector,labels=names[-ref.index],col=color.columns[-ref.index],pos=2)
		plot.out(plot.type)
	}
	
	
	out <- as.data.frame(cbind(SingleCS.vector,SingleCS.matrix))
	rownames(out) <- names[-ref.index]
	colnames(out)[1] <- "CSRankScores"
	
	
	return(out)
}


CSrank2 <- function(loadings,ref.index,color.columns=NULL,ref.plot=FALSE,loadings_names=NULL,component.plot,type.component="Factor",plot=TRUE,plot.type="pdf",basefilename="base",signCol="grey",legend.bg="grey",legend.names="",legend.cols="black"){
	
	term1.w <- 0.5
	term2.w <- 0.5
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	names <- loadings_names
	loadings <- loadings[,component.plot]
	
	RefCS.matrix <- matrix(0,nrow=length(loadings[-ref.index]),ncol=length(ref.index))
	colnames(RefCS.matrix) <- rep("",length(ref.index))
	
	query.loadings <- loadings[-ref.index]
	
	
	# |q|/|R| for each ref  
	for(i.ref in 1:length(ref.index)){
		ref.loading <- loadings[ref.index[i.ref]]	
		
		
		temp <- sapply(query.loadings,FUN=function(x){
					
					if(abs(ref.loading)>=abs(x)){
						return(abs(x)/abs(ref.loading))
					}
					else{
						return(abs(ref.loading)/abs(x))
					}
				})
		
		RefCS.matrix[,i.ref] <- temp
		colnames(RefCS.matrix)[i.ref] <- paste0(names(loadings)[ref.index[i.ref]]," (Ref ",i.ref,")")
	}
	
	weights <- abs(loadings[ref.index])
	term1 <- apply(RefCS.matrix,MARGIN=1,FUN=weighted.mean,w=weights)
	
	
	# Extra term for each q
	query.mean <- mean(abs(query.loadings))
	term2.temp <- abs(query.loadings)-query.mean
	
	term2 <- term2.temp/max(abs(term2.temp))
	
	# Combine both terms with (weighted) mean  +  change to 0 if result is negative
	combinetest1 <- sapply(1:length(term1),FUN=function(x){weighted.mean(c(term1[x],term2[x]),c(term1.w,term2.w))})
	combinetest2 <- sapply(combinetest1,FUN=function(x){if(x<0){return(0)}else{return(x)}})
	
	# Add Sign to final Score
	sign.temp <- sapply(query.loadings,FUN=function(x){
				ifelse(sign(mean(loadings[ref.index]))==sign(x),1,-1)
			})
	
	# Final
	combinetest3 <- combinetest2*sign.temp
	
	out <- as.data.frame(cbind(combinetest3,term1,RefCS.matrix,term2,combinetest1,combinetest2))
	colnames(out) <- c("CSRankScores","Term1.Weighted",paste0("Term1.Ref",c(1:length(ref.index))),"Term2","Combine1.wmean","Combine2.zero")
	rownames(out) <- names[-ref.index]
	
	
	## PLOTTING
	
	if(plot){
		
		if(ref.plot){
			
			nref <- dim(RefCS.matrix)[2]
			nplots <- dim(out)[2]-1
						
			col.plots <- ifelse(nplots<4,nplots,4)
			row.plots <- ifelse(nplots%%4==0,nplots%/%4,nplots%/%4 +1)
			par(mfrow=c(row.plots,col.plots))
						
			plot.in(plot.type,paste0(basefilename,"_CSRankExtra.pdf"))
			for(i.plot in 2:dim(out)[2]){
				plot(out[,i.plot],xlab="Query Compound Index",ylab="Score",col=color.columns[-ref.index],bg="grey",pch=21,main=paste0(type.component," ",component.plot," - ",colnames(out)[i.plot]))
				text(out[,i.plot],names[-ref.index],col=color.columns[-ref.index],pos=2)
			}
			plot.out(plot.type)
		
		}
		
		par(mfrow=c(1,1))
		
		## Plot Final
		main.temp <- ifelse(ref.plot,"CS Rank Score (Final)","CS Rank Score")
		plot.in(plot.type,paste0(basefilename,"_CSRank.pdf"))
		plot(combinetest3,col=color.columns[-ref.index],xlab="Query Compound Index",ylab="CS Rankscore",bg=signCol,pch=21,main=paste0(type.component," ",component.plot," - ",main.temp))
		text(combinetest3,labels=names[-ref.index],col=color.columns[-ref.index],pos=2)
		if(length(legend.names)>0){legend("topright",legend.names,pch=21,col=legend.cols,pt.bg=legend.bg,bty="n")}
		
		plot.out(plot.type)
	}
	
	

	return(out)

}
