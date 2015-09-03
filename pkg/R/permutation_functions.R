# Project: CSFA
# 
# Author: lucp8394
###############################################################################


#' Permute CS results
#' 
#' Apply permutation on MFA or Zhang results to obtain p-values. 
#' The function asks for a CSresult object which is returned by CSanalysis. The CSpermute function will return the same CSresult object with added information such as p-values.
#' If asked, the CSpermute function will also draw a volcanoplot and/or histograms of the p-values. If you simply want to redraw these plots, simply use the returned CSresult object by CSpermute again in the CSpermute function.
#' If the number of permutations was not changed, this will prevent the entire permutation analysis from being redone.
#' 
#' @export
#' @param refMat Reference matrix (Rows = genes and columns = compounds).
#' @param querMat Query matrix
#' @param CSresult A CSresult class object.
#' @param B Number of permutations.
#' @param mfa.factor If permuting a CSmfa result, mfa.factor will decide of which factor the p-values should be computed. If \code{NULL}, the factor chosen in CSanalysis will be chosen (the factor chosen in the CS slot of the CSresult). NOTE: If the mfa.factor is different from the factor in the CS slot, the CS slot will be overwritten with this new factor.
#' @param method.adjust Correction method of multiplicity adjusted p-values: "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" or "fdr". (Raw p-values are also always provided)
#' @param verbose If \code{TRUE}, progression dots of the permutation analysis will be printed.
#' @param querMat.Perm Possible to provide user-created permuted querMat data. Should be a list object of B times permuted querMat matrices.
#' @param save.querMat.Perm If \code{TRUE}, the list of permuted querMat matrices will be saved in the permutation.object slot of the CSresult.
#' @param which Choose which plot to draw. 1: A volcano plot of the -log(p-values) versus the observed connection scores. 2: A histogram of the permuted connection scores under the null hypothesis for a specific compound. A vertical line(s) is added for the observed CS and its p-value. The \code{cmpd.hist} parameter determines which compounds are drawn like this.
#' @param cmpd.hist Decides which compounds are plotted for the histogram distribution under null hypothesis (\code{which=2}). If \code{NULL}, you can select which compounds you want interactively on the volcano plot.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param color.columns Option to color the compounds on the volcano plot (\code{which=1}). Should be a vector of colors with the length of number of references and queries together.
#' @param basefilename Basename of the graphs if saved in pdf files
#' @return Returns the same CSresult object with added p-values to the CS slot and added information to the permutation.object slot. This CSresult can be reused in CSpermute to redraw the plots without calculation.
#' @examples
#' \dontrun{
#' data("dataSIM",package="CSFA")
#' Mat1 <- dataSIM[,c(1:6)]
#' Mat2 <- dataSIM[,-c(1:6)]
#' 
#' MFA_analysis <- CSanalysis(Mat1,Mat2,"CSmfa")
#' MFA_analysis <- CSpermute(Mat1,Mat2,MFA_analysis,B=200)
#' }
CSpermute <- function(refMat,querMat,CSresult,B=500,mfa.factor=NULL,method.adjust="none",verbose=TRUE,querMat.Perm=NULL,save.querMat.Perm=FALSE,
		which=c(1,2),cmpd.hist=NULL,color.columns=NULL,plot.type="device",basefilename="CSpermute"){
	
	if(class(CSresult)!="CSresult"){stop("CSresult is not a class object of CSresult")}
	
	type <- CSresult@type
	ref.index <- c(1:dim(refMat)[2])
	
	if(!(type %in% c("CSmfa","CSzhang"))){stop("Permutation is only available for MFA and Zhang results")}
	
	
	# redo permutation and analysis if object not available of different number of permutations is asked
	if(is.null(CSresult@permutation.object) | (length(CSresult@permutation.object$CS.Perm)!=B)){ 
		
		permutation.object <- list()
		
		if(is.null(querMat.Perm)){
			###### MAKING PERMUTATION DATA - unless perm data is provided ######
			
			if(verbose){
				cat("Permuting Data:\n")
			}
			
			querMat.Perm <- list()
			
			for(i in 1:B){
				perm <- sample(c(1:dim(querMat)[1]),dim(querMat)[1],replace=FALSE)
				querMat.Perm[[i]] <- querMat[perm,]
				
				
				if(verbose){
					cat(".")
					if(i%%100==0){cat(" ",i,"\n")}
				}
			}
			if(verbose){cat("\nDONE\n")}
		}
		
		if(save.querMat.Perm){permutation.object$querMat.Perm <- querMat.Perm}
	
	
		###### APPLYING ANALYSIS TO DATA ########
		
		
		## MFA Analysis ##
		if(type=="CSmfa"){
		
			if(verbose){cat("Analysing Permuted Data with MFA\n(Factor is chosen based on highest average reference loadings):\n")}
		
			MFA.Perm <- list()
#			MFA.bc.Perm <- c(1:B)
		
			for(i in 1:B){
				data_comb <- cbind(refMat,querMat.Perm[[i]])
				rownames(data_comb) <- rownames(refMat)
				colnames(data_comb) <- c(colnames(refMat),colnames(querMat))
			
				out_MFA <- MFA(base=data_comb,group=c(ncol(refMat),ncol(querMat)),type=c("s","s"),ind.sup=NULL,row.w=CSresult@call$analysis.pm$row.w,weight.col.mfa=CSresult@call$analysis.pm$weight.col.mfa,ncp=CSresult@call$analysis.pm$ncp,name.group=c("Reference","Query"),graph=FALSE)
				ref.loadings <- apply(out_MFA$quanti.var$coord,MARGIN=2,FUN=function(x){mean(x[ref.index])})
				factor.choose <- which(abs(ref.loadings)==max(abs(ref.loadings)))
#				MFA.bc.Perm[i] <- factor.choose
				MFA.Perm[[i]] <- out_MFA$quanti.var$coord[,factor.choose]
		
				if(verbose){
					cat(".")
					if(i%%100==0){cat(" ",i,"\n")}
				}
			}
			if(verbose){cat("\nDONE\n")}
			
			permutation.object$CS.Perm <- MFA.Perm
			
		}
	
		## Zhang Analysis ##
		if(type=="CSzhang"){
			if(verbose){cat("Analysing Permuted Data with Zhang and Gant:\n")}
			
			Zhang.Perm <- list()
		
			for(i in 1:B){
				rownames(querMat.Perm[[i]]) <- rownames(refMat)
				colnames(querMat.Perm[[i]]) <- colnames(querMat)
			
				out_zhang <- analyse_zhang(refMat,querMat.Perm[[i]],
						nref=CSresult@call$analysis.pm$nref,nquery=CSresult@call$analysis.pm$nquery,ord.query=CSresult@call$analysis.pm$ord.query,ntop.scores=CSresult@call$analysis.pm$ntop.scores,
						basefilename="analyseZhang",which=c(),plot.type="device",print.top=FALSE)
				
				dim2.data <- (dim(refMat)[2]+dim(querMat)[2])
				zhang.temp <- rep(0,dim2.data)
				i.zhang <- 1
				for(ii in 1:dim2.data){
					if(!(ii %in% ref.index)){
						zhang.temp[ii] <- out_zhang$All[i.zhang,1]
						i.zhang <- i.zhang+1
					}
					else{
						zhang.temp[ii] <- 1
					}
				}
				Zhang.Perm[[i]] <- zhang.temp

				if(verbose){
					cat(".")
					if(i%%100==0){cat(" ",i,"\n")}
				}
			}
			
			if(verbose){cat("\nDONE\n")}
			
			permutation.object$CS.Perm <- Zhang.Perm
			
			mfa.factor <- NULL
		}
	}
	else{
		permutation.object <- CSresult@permutation.object # If the previous part was not necessary to do, extract the existing permutation.objects
	}
	
	
	##### COMPUTING THE P-VALUES  + FINISHING PERMUTATION.OBJECT #####
	
	if(type=="CSmfa"){
		# Choosing the mfa.factor and giving the proper warning later (namely that CS and GS will be overwritten)
		
		if(is.null(mfa.factor)){ # If no factor given, take the one chosen in CSanalysis
			mfa.factor <- CSresult@call$factor.select
		}
		
		if(!is.null(mfa.factor) & (mfa.factor != CSresult@call$factor.select)){
			CS.warning <- TRUE # Prepare for warning if mfa.factor is given and different from the one in CSanalysis. CSresult will need to be updated
			CSresult@call$factor.select <- mfa.factor
		}
		else{
			CS.warning <- FALSE
		}
	}
			
	# Note: pvalues will always be re-computed to allow for different mfa.factor
	pval.dataframe <- pvalue_compute(obs.result=CSresult,list.h0.result=permutation.object$CS.Perm,mfa.factor=mfa.factor,method.adjust=method.adjust)
	
	permutation.object$pval.dataframe <- pval.dataframe
	
	CSresult@permutation.object <- permutation.object # Saving the updated permutation.object in CSresult
		
	
	#### UPDATE THE CSresult object which will be returned in the end
	if(type=="CSzhang"){
		
		## Adding p-values to CS slot
		
		CS <- CSresult@CS
		
		CS$CS.ref$pvalues <- pval.dataframe$pvalues[ref.index]
		CS$CS.query$pvalues <- pval.dataframe$pvalues[-ref.index]
		
		if(method.adjust!="none"){
			CS$CS.ref$pvalues.adjusted <- pval.dataframe$pvalues.adjusted[ref.index]
			CS$CS.query$pvalues.adjusted <- pval.dataframe$pvalues.adjusted[-ref.index]
		}
		
		CSresult@CS <- CS
		
		## Adding p-values to CS.top.query slot
		
		oldtop <- CSresult@CS$CS.top.query
		posname <- as.character(oldtop$posname)
		negname <- as.character(oldtop$negname)
		
		index.pos <- sapply(posname,FUN=function(x){which(x==pval.dataframe$Cmpd)})
		index.neg <- sapply(negname,FUN=function(x){which(x==pval.dataframe$Cmpd)})
		p.pos <- pval.dataframe$pvalues[index.pos]
		p.neg <- pval.dataframe$pvalues[index.neg]
		if(method.adjust!="none"){
			p.pos.adj <- pval.dataframe$pvalues.adjusted[index.pos]
			p.neg.adj <- pval.dataframe$pvalues.adjusted[index.neg]
			newtop <- data.frame(posname=oldtop$posname,posscore=oldtop$posscore,pospvalue=p.pos,pospvalue.adj=p.pos.adj,negname=oldtop$negname,negscore=oldtop$negscore,negpvalue=p.neg,negvalue.adj=p.neg.adj)
		}
		else{
			newtop <- data.frame(posname=oldtop$posname,posscore=oldtop$posscore,pospvalue=p.pos,negname=oldtop$negname,negscore=oldtop$negscore,negpvalue=p.neg)
			
		}
		CSresult@CS$CS.top.query <- newtop
		
	}
	if(type=="CSmfa"){
		
		
		CS <- CSresult@CS
		
		if(CS.warning){
			# If warning, re-extract CS and GS 
			warning("Due to choice of mfa.factor, the CS and GS slot in CSresult will be overwritten with another factor. ")
						
			loadings <- CSresult@object$quanti.var$coord	
			scores <- CSresult@object$ind$coord			
			
			CS.loadings <- data.frame(loadings[,mfa.factor])
			CS.query <- data.frame(CS.loadings[-c(1:dim(refMat)[2]),])
			CS.ref <- data.frame(CS.loadings[c(1:dim(refMat)[2]),])
			
			colnames(CS.query) <- colnames(CS.ref) <- sapply(mfa.factor,FUN=function(x){paste0("Factor",x)})
			rownames(CS.query) <- rownames(loadings)[-c(1:dim(refMat)[2])]
			rownames(CS.ref) <- rownames(loadings)[c(1:dim(refMat)[2])]	
			
			CS <- list(CS.query=CS.query,CS.ref=CS.ref)
			
			GS <- data.frame(scores[,mfa.factor])
			rownames(GS) <- rownames(scores)
			colnames(GS) <- sapply(mfa.factor,FUN=function(x){paste0("Factor",x)})
					
			CSresult@GS <- GS
		}
		# Add p-values to CS slot			
		CS$CS.ref$pvalues <- pval.dataframe$pvalues[ref.index]
		CS$CS.query$pvalues <- pval.dataframe$pvalues[-ref.index]
			
		if(method.adjust!="none"){
			CS$CS.ref$pvalues.adjusted <- pval.dataframe$pvalues.adjusted[ref.index]
			CS$CS.query$pvalues.adjusted <- pval.dataframe$pvalues.adjusted[-ref.index]
		}
		
		CSresult@CS <- CS
	}
			
		
	
	##### POSSIBLE PLOTS #####
	
	# NOTE: Still need to check all cases with 2 plot.types 
	
	hist.drawn <- FALSE
	
	# Volcano plot
	if(1 %in% which){
		
		if((2 %in% which) & (is.null(cmpd.hist))){
			# Volcano plot with selecting compounds
			pvalue_volc(pval.dataframe=permutation.object$pval.dataframe,type=type,color.columns=color.columns,list.h0.result=permutation.object$CS.Perm,make.hist=TRUE,plot.type=plot.type,plot.type.hist=plot.type,basefilename=paste0(basefilename,"_volcanoplot"))
						
			hist.drawn <- TRUE
		}else{
			# Volcano plots without selecting compounds
			pvalue_volc(pval.dataframe=permutation.object$pval.dataframe,type=type,color.columns=color.columns,list.h0.result=NULL,make.hist=FALSE,plot.type=plot.type,plot.type.hist=plot.type,basefilename=paste0(basefilename,"_volcanoplot"))
		
		}

	}
	
	
	# Histogram plot
	if((2 %in% which) & !hist.drawn){
		if(is.null(cmpd.hist)){
			#volc with plot device
			pvalue_volc(pval.dataframe=permutation.object$pval.dataframe,type=type,color.columns=color.columns,list.h0.result=permutation.object$CS.Perm,make.hist=TRUE,plot.type="device",plot.type.hist=plot.type,basefilename=paste0(basefilename))
		
		}else{
			# just hist
			pvalue_hist(pval.dataframe=permutation.object$pval.dataframe[cmpd.hist,],list.h0.result=permutation.object$CS.Perm,type=type,plot.type=plot.type,basefilename=paste0(basefilename,"_histogram"))
			
		}

		hist.drawn <- TRUE
	}
	
	
	
	##### RETURN OBJECT ######
	
#	# Also add the method parameters used for this permutation...
#	CSresult@permutation.object$analysis.perm.pm <- CSresult@call$analysis.pm
	
	return(CSresult)
	# Contains:
	# - updated CS with pvalues 
	# - permutation.object slot with permuted CS, pval.dataframe and if asked the permuted data.
	
	
}



# add adjusted pvalues
pvalue_compute <- function(obs.result,list.h0.result,cmpd.index=NULL,mfa.factor=1,method.adjust=c("none","holm","hochberg","hommel","bonferroni","BH","BY","fdr")){
	if(class(obs.result)!="CSresult"){stop("obs.result is not a \"CSresult\" class object")}
	
	type <- obs.result@type
	
	if(!(type %in% c("CSzhang","CSmfa"))){stop("Only Zhang and MFA results can be used")}
	
	if(type=="CSzhang"){
		
		obs.scores <- rbind(obs.result@CS$CS.ref,obs.result@CS$CS.query)
		
		if(is.null(cmpd.index)){cmpd.index <- c(1:dim(obs.scores)[1])}
		
		pval.dataframe <- data.frame(Cmpd=rownames(obs.scores)[cmpd.index],cmpd.index=cmpd.index)
		temp.pval <- c()
		temp.obs <- c()
		
		for(i.cmpd in cmpd.index){
			h0.data <- unlist(lapply(list.h0.result,FUN=function(x){return(x[i.cmpd])}))
			obs <- obs.scores[i.cmpd,1]
			pvalue <- ifelse(obs>=0 ,sum(h0.data>obs)/(length(h0.data)+1), sum(h0.data<obs)/(length(h0.data)+1))
			temp.pval <- c(temp.pval,pvalue)
			temp.obs <- c(temp.obs,obs)
		}
		
		pval.dataframe$pvalues <- temp.pval
		if(method.adjust!="none"){pval.dataframe$pvalues.adjusted <- p.adjust(temp.pval,method=method.adjust)}
		pval.dataframe$observed <- temp.obs
		
		
		return(pval.dataframe)
	}
	
	if(type=="CSmfa"){
		
		obs.scores <- obs.result@object$quanti.var$coord[,mfa.factor]
		
		if(is.null(cmpd.index)){cmpd.index <- c(1:length(obs.scores))}
		
		pval.dataframe <- data.frame(Cmpd=names(obs.scores)[cmpd.index],cmpd.index=cmpd.index)
		temp.pval <- c()
		temp.pval1 <- c()
		temp.pval2 <- c()
		temp.obs <- c()
		
		for(i.cmpd in cmpd.index){
			
			h0.data <- unlist(lapply(list.h0.result,FUN=function(x){return(x[i.cmpd])}))
			obs <- obs.scores[i.cmpd]
			
			# PROBLEM WITH MFA: Solution 2-sided p-value? (FOR MFA SIGN ALSO NOT ALWAYS POS OR NEG CONN, DEPENDS ON DATA!!!!!)
			pvalue1 <- ifelse(obs>=0 ,sum(h0.data>obs)/(length(h0.data)+1), sum(h0.data<obs)/(length(h0.data)+1))
			pvalue2 <- ifelse(obs>=0 ,sum(h0.data<(-obs))/(length(h0.data)+1), sum(h0.data>(-obs))/(length(h0.data)+1))
			
			temp.pval1 <- c(temp.pval1,pvalue1)
			temp.pval2 <- c(temp.pval2,pvalue2)
			temp.pval <- c(temp.pval,pvalue1+pvalue2)
			temp.obs <- c(temp.obs,obs)
			
		}
		pval.dataframe$pvalues <- temp.pval
		pval.dataframe$pvalue1 <- temp.pval1
		pval.dataframe$pvalue2 <- temp.pval2
		if(method.adjust!="none"){pval.dataframe$pvalues.adjusted <- p.adjust(temp.pval,method=method.adjust)}
		pval.dataframe$observed <- temp.obs
		
		
		return(pval.dataframe)
		
	}
}

# function for all pvalues and h0 given
pvalue_hist <- function(pval.dataframe,list.h0.result,type=c("CSmfa","CSzhang"),plot.type="device",basefilename="pvalue_hist"){
	
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(paste0(name,".pdf"))}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	
	#if(dim(pval.dataframe)[1]!=length(list.h0.result[[1]])){stop("Number of pvalues is not same as number h0 results")}
	cmpd.names <- pval.dataframe$Cmpd
	cmpd.index <- pval.dataframe$cmpd.index
	obs.scores <- pval.dataframe$observed
	
	
	if(type=="CSzhang"){
		
		for(i in 1:dim(pval.dataframe)[1]){
			
			h0.data <- unlist(lapply(list.h0.result,FUN=function(x){return(x[cmpd.index[i]])}))
			
			obs <- obs.scores[i]
			pvalue <- pval.dataframe$pvalues[i]
			
			plot.in(plot.type,paste0(basefilename,"_cmpd",cmpd.index[i]))
			hist(h0.data,xlim=c(min(c(h0.data,obs)),max(c(h0.data,obs))),nclass=15,col='lightblue',main=paste0('Distribution of Cmpd ',cmpd.index[i]," (",cmpd.names[i],") under H0"),xlab="Zhang Score")
			abline(v=obs,lw=2,col="red")
			if(obs>=0){text(obs,50, paste0("Observed Value (p-value=",pvalue,")"),pos=2,col="red",offset=1)} else{text(obs,50, paste0("Observed Value (p-value=",pvalue,")"),pos=4,col="red",offset=1)}
			plot.out(plot.type)
		}
		
	}
	
	if(type=="CSmfa"){
		
		
		for(i in 1:dim(pval.dataframe)[1]){
			h0.data <- unlist(lapply(list.h0.result,FUN=function(x){return(x[cmpd.index[i]])}))
			obs <- obs.scores[i]
			pvalue1 <- pval.dataframe$pvalue1[i]
			pvalue2 <- pval.dataframe$pvalue2[i]
			
			
			# HERE
			plot.in(plot.type,paste0(basefilename,"_cmpd",cmpd.index[i]))
			hist(h0.data,xlim=c(min(c(h0.data,obs,-obs)),max(c(h0.data,obs,-obs))),nclass=15,col='lightblue',main=paste0('Distribution of Cmpd ',cmpd.index[i]," (",cmpd.names[i],") under H0"),xlab="MFA CScore")
			abline(v=obs,lw=2,col="red")
			if(obs>=0){text(obs,50, paste0("Observed Value (",pvalue1,")"),pos=2,col="red",offset=1)} else{text(obs,50, paste0("Observed Value (p-value=",pvalue1,")"),pos=4,col="red",offset=1)}
			abline(v=-obs,lw=2,col="red",lty=2)
			if(obs>=0){text(-obs,30, paste0("- Observed Value (",pvalue2,")"),pos=4,col="red",offset=1)} else{text(-obs,50, paste0("- Observed Value (",pvalue2,")"),pos=2,col="red",offset=1)}
			text(obs,60,paste0("P-Value = ",pvalue1+pvalue2),col="red",pos=ifelse(obs>=0,2,4))
			plot.out(plot.type)
		}
		
	}
	
	
}


# Use compute pvalue plot
pvalue_volc <- function(pval.dataframe,type=c("CSmfa","CSzhang"),color.columns=NULL,list.h0.result=NULL,make.hist=FALSE,plot.type="device",plot.type.hist="device",basefilename="pvalue_volc"){
	
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(paste0(name,".pdf"))}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	pvalues <- pval.dataframe$pvalues
	obs.scores <- pval.dataframe$observed
	
	plot.pvalues <- -log(pvalues)
	inf.index <- which(plot.pvalues == Inf)
	if(length(inf.index)>0){inf.present <- TRUE}else{inf.present <- FALSE}
	
	max.real <- max(plot.pvalues[-inf.index])
	plot.pvalues[inf.index] <- max.real+1
	
	
	if(is.null(color.columns)){color.columns <- "black"}
	
	plot.in(plot.type,basefilename)
	plot(obs.scores,plot.pvalues,col=color.columns,bg="grey",pch=21,xlab="Observed CS for cmpds",ylab="-log(pvalues)",main=paste0("Volcano Plot for ",type," result"))
	text(obs.scores,plot.pvalues,pval.dataframe$Cmpd,pos=1,col=color.columns)
	if(inf.present){axis(4,at=max.real+1,"Inf Value",col="red")}
	plot.out(plot.type)
	
	if(make.hist){
		if(!is.null(list.h0.result)){
			#if(dim(pval.dataframe)[1]!=length(list.h0.result[[1]])){stop("Number of pvalues is not same as number h0 results")}
			
			if(plot.type!="device"){
				dev.new()
				plot(obs.scores,plot.pvalues,col=color.columns,bg="grey",pch=21,xlab="Observed CS for cmpds",ylab="-log(pvalues)",main=paste0("Volcano Plot for ",type," result"))
				text(obs.scores,plot.pvalues,pval.dataframe$Cmpd,pos=1,col=color.columns)
				if(inf.present){axis(4,at=max.real+1,"Inf Value",col="red")}
			}
						
			cat("Please choose one or more compounds with left-click.\n To end the selection, press right-click.")
			choose.cmpd <- identify(obs.scores,plot.pvalues,n=9999,labels=pval.dataframe$Cmpd,col="slateblue3")
			
			pvalue_hist(pval.dataframe[choose.cmpd,],list.h0.result,type,plot.type=plot.type.hist,basefilename=paste(basefilename,"_histogram"))
			
		}
	}
}

# add this to the CScompare with a pvalcompare option
# special case for only 2 pval.dataframe and when transforming already happened before
pvalue2_compare <- function(list.pval.dataframe,threshold=0.05){
	if(length(list.pval.dataframe)!=2){stop("Need exactly 2 pval.dataframe results")}
	
	pvalues <- lapply(list.pval.dataframe,FUN=function(x){x$pvalues})
	m.pvalues <- matrix(unlist(pvalues),ncol=length(pvalues))
	compare.val <- apply(m.pvalues,MARGIN=1,FUN=function(x){paste0("S",c(1:length(pvalues))[x<=threshold],collapse="")})
	compare.val <- sapply(compare.val,FUN=function(x){if(x=="S"){return("NS")}else{return(x)}})
	
	table.compare.val <- table(compare.val)
	
	
	m11 <- ifelse(!is.na(table.compare.val["S1S2"]),table.compare.val["S1S2"],0)
	m21 <- ifelse(!is.na(table.compare.val["S1"]),table.compare.val["S1"],0)
	m12 <- ifelse(!is.na(table.compare.val["S2"]),table.compare.val["S2"],0)
	m22 <- ifelse(!is.na(table.compare.val["NS"]),table.compare.val["NS"],0)
	
	pval.matrix <- matrix(c(m11,m21,m12,m22),nrow=2)
	colnames(pval.matrix) <- c("Result1.Sign","Result1.NotSign")
	rownames(pval.matrix) <- c("Result2.Sign","Result2.NotSign")
		
	# adjusted pvalues
	adjusted.available <- unlist(lapply(list.pval.dataframe,FUN=function(x){"pvalues.adjusted"%in%colnames(x)}))
	if(sum(adjusted.available)==1){warning("Only 1 of the results has adjusted p-values.")}
	if(sum(adjusted.available)==2){
		
		pvalues.adj <- lapply(list.pval.dataframe,FUN=function(x){x$pvalues.adjusted})
		m.pvalues.adj <- matrix(unlist(pvalues.adj),ncol=length(pvalues.adj))
		compare.val.adj <- apply(m.pvalues.adj,MARGIN=1,FUN=function(x){paste0("S",c(1:length(pvalues.adj))[x<=threshold],collapse="")})
		compare.val.adj <- sapply(compare.val.adj,FUN=function(x){if(x=="S"){return("NS")}else{return(x)}})
		
		table.compare.val.adj <- table(compare.val.adj)
		
		
		m11 <- ifelse(!is.na(table.compare.val.adj["S1S2"]),table.compare.val.adj["S1S2"],0)
		m21 <- ifelse(!is.na(table.compare.val.adj["S1"]),table.compare.val.adj["S1"],0)
		m12 <- ifelse(!is.na(table.compare.val.adj["S2"]),table.compare.val.adj["S2"],0)
		m22 <- ifelse(!is.na(table.compare.val.adj["NS"]),table.compare.val.adj["NS"],0)
		
		pval.adj.matrix <- matrix(c(m11,m21,m12,m22),nrow=2)
		colnames(pval.adj.matrix) <- c("Result1.Sign","Result1.NotSign")
		rownames(pval.adj.matrix) <- c("Result2.Sign","Result2.NotSign")
		
	}
	else{pval.adj.matrix <- NULL}
		
	# add to pval.dataframes
	for(i in 1:length(list.pval.dataframe)){
		list.pval.dataframe[[i]]$pval.sign <- compare.val
		if(sum(adjusted.available)==2){list.pval.dataframe[[i]]$pval.adj.sign <- compare.val.adj}
	}
		
	# return 1 (or 2) matrices and the new list.pval.dataframes
	return(list(list.pval.dataframe=list.pval.dataframe,compare.pvalues=pval.matrix,compare.pvalues.adjusted=pval.adj.matrix))
	
	
}


