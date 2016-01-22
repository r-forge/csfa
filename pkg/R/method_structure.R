# Project: Connectivity
# 
# Author: lucp8394
###############################################################################

## EXAMPLE DATA ##
#' Simulated Microarray Data
#'
#' A matrix containing some simulated example microarray data. 
#' The first 6 columns of this matrix make up the reference matrix part.
#'
#' @format A matrix with 1000 rows and 341 columns.
#' @name dataSIM
NULL


## IMPORTS ##

#' @import methods
#' @importFrom fabia fabia extractBic showSelected
#' @importFrom FactoMineR PCA MFA
#' @importFrom pls stdize
#' @importFrom elasticnet spca arrayspc



## CLASSES ##
setClass("CSfabia",slots=c(call="call"))
setClass("CSmfa",slots=c(call="call"))
setClass("CSpca",slots=c(call="call"))
setClass("CSsmfa",slots=c(call="call"))
setClass("CSzhang",slots=c(call="call"))

#' An S4 class in which the results of the Connectivity Scores by Factor Analysis are stored.
#' 
#' @export
#' @slot type A character string containing the analysis type.
#' @slot CS List containing the connectivity scores for the query (and reference loadings). If permutation was applied, will also contain p-values.
#' @slot GS Dataframe containing the gene scores.
#' @slot CSRankScores Dataframe containing the CS Rank Scores (Not for CSzhang). First column is the weighted mean and the subsequent ones per reference.
#' @slot object Object containing the FA or Zhang result.
#' @slot permutation.object Contains permuted data, CS for permuted data and a dataframe with the p-values (only for MFA and Zhang).
#' @slot call List object containing the original call of \code{CSanalysis} as well as the parameters for the chosen method.
setClass("CSresult",slots=list(type="character",CS="list",GS="data.frame",CSRankScores="list",object="ANY",permutation.object="ANY",call="ANY"))

setClass("CSzhangCompare",slots=c(CSresult1="CSresult",CSresult2="CSresult"))
setClass("CSfabiaCompare",slots=c(CSresult1="CSresult",CSresult2="CSresult"))


## METHODS  - CSANALYSIS ##
#' Connectivity Score Analysis.
#' 
#' Doing a CS analysis, interactively generating graphs. See specific type for additional parameteres.
#' 
#' @export
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type Type of Factor Analysis or Zhang & Gant ( \code{"CSfabia"}, \code{"CSmfa"}, \code{"CSpca"}, \code{"CSsmfa"} or \code{"CSzhang"})
#' @param ... Additional parameters for analysis
#' @return An object of the S4 Class \code{CSresult}.
#' @examples
#' \dontrun{
#' data("dataSIM",package="CSFA")
#' Mat1 <- dataSIM[,c(1:6)]
#' Mat2 <- dataSIM[,-c(1:6)]
#' 
#' MFA_analysis <- CSanalysis(Mat1,Mat2,"CSmfa")
#' FABIA_analysis <- CSanalysis(Mat1,Mat2,"CSfabia")
#' ZHANG_analysis <- CSanalysis(Mat1,Mat2,"CSzhang")
#' }
setGeneric('CSanalysis', function(refMat,querMat,type, ...){standardGeneric('CSanalysis')})


#' Connectivity Score Analysis.
#' 
#' Doing a CS analysis, interactively generating graphs. See specific type for additional parameteres.
#' 
#' @export
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type Type of Factor Analysis or Zhang & Gant ( \code{"CSfabia"}, \code{"CSmfa"}, \code{"CSpca"}, \code{"CSsmfa"} or \code{"CSzhang"})
#' @param ... Additional parameters for analysis
#' @return An object of the S4 Class \code{CSresult}.
#' @examples
#' \dontrun{
#' data("dataSIM",package="CSFA")
#' Mat1 <- dataSIM[,c(1:6)]
#' Mat2 <- dataSIM[,-c(1:6)]
#' 
#' MFA_analysis <- CSanalysis(Mat1,Mat2,"CSmfa")
#' FABIA_analysis <- CSanalysis(Mat1,Mat2,"CSfabia")
#' ZHANG_analysis <- CSanalysis(Mat1,Mat2,"CSzhang")
#' }
setMethod('CSanalysis', c('matrix','matrix','character'),
		function(refMat,querMat,type, ...) {
			if(type %in% c("CSfabia","CSmfa","CSpca","CSsmfa","CSzhang")){
				type <- new(type,call=match.call())
				CSanalysis(refMat,querMat,type,...)
			}
			else{
				stop("This method type is not available.")
			}
		})


# FABIA

#' "CSfabia"
#' 
#' Doing interactive CS analysis with FABIA (Factor Analysis for Bicluster Acquisition). One or multiple reference compounds are possible in this analysis.
#' 
#' @export 
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type \code{"CSfabia"}
#' @param p \emph{Fabia Parameter:} number of hidden factors = number of biclusters; default = 13
#' @param alpha \emph{Fabia Parameter:} sparseness loadings (0 - 1.0); default = 0.01
#' @param cyc \emph{Fabia Parameter:} number of iterations; default = 500 
#' @param spl \emph{Fabia Parameter:} sparseness prior loadings (0 - 2.0); default = 0 (Laplace)
#' @param spz \emph{Fabia Parameter:} sparseness factors (0.5 - 2.0); default = 0.5 (Laplace)
#' @param non_negative \emph{Fabia Parameter:} Non-negative factors and loadings if non_negative > 0; default = 0
#' @param random \emph{Fabia Parameter:} <=0: by SVD, >0: random initialization of loadings in [-random,random]; default = 1.0
#' @param center \emph{Fabia Parameter:} data centering: 1 (mean), 2 (median), > 2 (mode), 0 (no); default = 2
#' @param norm \emph{Fabia Parameter:} data normalization: 1 (0.75-0.25 quantile), >1 (var=1), 0 (no); default = 1
#' @param scale \emph{Fabia Parameter:} loading vectors are scaled in each iteration to the given variance. 0.0 indicates non scaling; default = 0.0
#' @param lap \emph{Fabia Parameter:} minimal value of the variational parameter; default = 1.0
#' @param nL \emph{Fabia Parameter:} maximal number of biclusters at which a row element can participate; default = 0 (no limit)
#' @param lL \emph{Fabia Parameter:} maximal number of row elements per bicluster; default = 0 (no limit)
#' @param bL \emph{Fabia Parameter:} cycle at which the nL or lL maximum starts; default = 0 (start at the beginning)
#' @param which Choose one or more plots to draw:\cr 1. Information Content of Biclusters \cr 2. Loadings for reference compounds\cr 3. Bicluster \code{BC1.plot} VS Bicluster \code{BC2.plot} : Loadings & Genes\cr 4. Gene Scores for Bicluster \code{BC.Plot}\cr 5. Loadings for Factor \code{BC.plot}\cr 6. Profile plot (see \code{profile.type})\cr 7. CS Rank Scores for Factor \code{BC.plot}
#' @param BC.plot Which biclusters should be investigated? Can be a vector of multiple (e.g. \code{c(1,3,5)}). If \code{NULL}, you can choose biclusters of interest interactively from reference loadings plot.
#' @param column.interest Numeric vector of indices of query columns which should be in the profiles plots (\code{which=6}). If \code{NULL}, you can interactively select genes on the Compound Loadings plot (\code{which=5}).
#' @param row.interest Numeric vector of gene indices to be plotted in gene profiles plot (\code{which=6}, \code{profile.type="gene"}). If \code{NULL}, you can interactively select them in the gene scores plot (\code{which=4}).
#' @param profile.type Type of \code{which=6} plot:\cr - \code{"gene"}: Gene profiles plot of selected genes in \code{row.interest} with the reference compounds and those selected in \code{column.interest} ordered first on the x axis. The other compounds are ordered in decreasing CScore. \cr - \code{"cmpd"}: Compound profiles plot of reference and selected compounds (\code{column.interest}) and only those genes on the x-axis which beat the thresholds (\code{gene.thresP}, \code{gene.thresN})
#' @param color.columns Vector of colors for the reference and query columns (compounds). If \code{NULL}, blue will be used for reference and black for query. Use this option to highlight reference columns and query columns of interest.
#' @param gene.highlight Single numeric vector or list of maximum 5 numeric vectors. This highlights gene of interest in gene scores plot (\code{which=4}) up to 5 different colors.
#' @param gene.thresP Threshold for genes with a high score (\code{which=4}).
#' @param gene.thresN Threshold for genes with a low score (\code{which=4}).
#' @param thresP.col Color of genes above \code{gene.thresP}.
#' @param thresN.col Color of genes below \code{gene.thresN}.
#' @param legend.names Option to draw a legend of colored columns in Compound Loadings plot (\code{which=5}). If \code{NULL}, only "References" will be in the legend.
#' @param legend.cols Colors to be used in legend in \code{which=5}. If \code{NULL}, only blue for "References is used".
#' @param result.available You can a previously returned object by \code{CSanalysis} in order to only draw graphs, not recompute the scores.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param basefilename Basename of the graphs if saved in pdf files
#' @param CSrank.refplot Logical value deciding if the CS Rank Scores (\code{which=6}) should also be plotted per reference (instead of only the weighted mean).
#' @return An object of the S4 Class \code{CSresult}.
setMethod('CSanalysis',c('matrix','matrix','CSfabia'),
		function(refMat,querMat,type="CSfabia",p=13,alpha=0.01,cyc=500,spl=0,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0
				,which=c(2,4,5,7)
				,BC.plot=NULL
				,column.interest=NULL,row.interest=NULL,profile.type="gene"
				,color.columns=NULL,gene.highlight=NULL
				,gene.thresP=1,gene.thresN=-1,thresP.col="blue",thresN.col="red"
				,legend.names=NULL,legend.cols=NULL
				,CSrank.refplot=FALSE
				,result.available=NULL
				,basefilename="analyseFABIA",plot.type="device"
		) {

			data <- cbind(as.matrix(refMat),as.matrix(querMat))
			
			if(is.null(color.columns)){colour.columns <- c(rep("blue",dim(refMat)[2]),rep("black",dim(querMat)[2]))} else {colour.columns <- color.columns}
			
			if(is.null(legend.names)){legend.names <- c("References")}
			if(is.null(legend.cols)){legend.cols <- c("blue")}

			analysis.pm.temp <- list(p=p,alpha=alpha,cyc=cyc,spl=spl,spz=spz,non_negative=non_negative,random=random,center=center,norm=norm,scale=scale,lap=lap,nL=nL,lL=lL,bL=bL)
			
			if(!is.null(result.available)){
				if(class(result.available) != "CSresult"){
					stop("result.available is not of the correct class") 
				}
				else if(result.available@type != class(type)){
					stop("result.available is from a different analysis")
				}
				else if(!identical(result.available@call$analysis.pm,analysis.pm.temp)){
					stop("result.available used other analysis parameters")
				}
				else{result.available2 <- result.available}
			}
			else{
				result.available2 <- NULL
			}
			
			if(!is.null(column.interest)){column.interest <- column.interest + dim(refMat)[2]}

			out <- analyse_fabia(data=data,p=p,alpha=alpha,cyc=cyc,spl=spl,spz=spz,non_negative=non_negative,random=random,center=center,norm=norm,scale=scale,lap=lap,nL=nL,lL=lL,bL=bL,
					basefilename=basefilename,weighted.data=TRUE,
					ref.index=c(1:dim(refMat)[2]),
					BC.plot=BC.plot,
					column.interest=column.interest,row.interest=row.interest,gene.thresP=gene.thresP,gene.thresN=gene.thresN,
					colour.columns=colour.columns,legend.names=legend.names,legend.cols=legend.cols,thresP.col=thresP.col,thresN.col=thresN.col,
					result.available=result.available2,plot.type=plot.type,CSrank.refplot=CSrank.refplot,
					which=which,
					gene.highlight=gene.highlight,profile.type=profile.type)

			BC.select <- out[[1]]
			result <- out[[2]]
			
			loadings <- result@L
			scores <- t(result@Z)
			
			CS.loadings <- data.frame(loadings[,BC.select])
			CS.query <- data.frame(CS.loadings[-c(1:dim(refMat)[2]),])
			CS.ref <- data.frame(CS.loadings[c(1:dim(refMat)[2]),])
			
			colnames(CS.query) <- colnames(CS.ref) <- sapply(BC.select,FUN=function(x){paste0("BC",x)})
			rownames(CS.query) <- rownames(loadings)[-c(1:dim(refMat)[2])]
			rownames(CS.ref) <- rownames(loadings)[c(1:dim(refMat)[2])]	
			
			CS <- list(CS.query=CS.query,CS.ref=CS.ref)
			
			
			
			GS <- data.frame(scores[,BC.select])
			rownames(GS) <- rownames(scores)
			colnames(GS) <- sapply(BC.select,FUN=function(x){paste0("BC",x)})
			
			call.object <- list(match.call=type@call,analysis.pm=list(p=p,alpha=alpha,cyc=cyc,spl=spl,spz=spz,non_negative=non_negative,random=random,center=center,norm=norm,scale=scale,lap=lap,nL=nL,lL=lL,bL=bL),BC.select=BC.select)

			if(!is.null(result.available)){
				if(BC.select != result.available@call$BC.select){
					if(!is.null(result.available@permutation.object)){
						warning("CS, GS and CSRankScores Slot in CSresult will be overwritten due to different factor choice. To add p-values again, run CSpermute again on this object!\nThe permutation.object was also deleted due to a different choice of BC.")
						
					}else{
						warning("CS, GS and CSRankScores Slot in CSresult will be overwritten due to different BC choice.")
						
					}
					
					return(new("CSresult",type="CSfabia",CS=CS,CSRankScores=out[[3]],GS=GS,object=result,permutation.object=NULL,call=call.object))
					
				}
				else{
					return(new("CSresult",type="CSfabia",CS=result.available@CS,CSRankScores=result.available@CSRankScores,GS=result.available@GS,object=result.available@object,permutation.object=result.available@permutation.object,call=call.object))
					
				}
				
			}
			else{
				return(new("CSresult",type="CSfabia",CS=CS,CSRankScores=out[[3]],GS=GS,object=result,permutation.object=NULL,call=call.object))
				
			}
		})




# MFA
#' "CSmfa"
#' 
#' Doing interactive CS analysis with MFA (Multiple Factor Analysis). Should use multiple references for this analysis.
#' 
#' @export 
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type \code{"CSmfa"}
#' @param ncp \emph{MFA Parameter:} Number of dimensions kept in the results (by default 5).
#' @param weight.col.mfa \emph{MFA Parameter:} Vector of weights, useful for HMFA method (by default, \code{NULL} and an MFA is performed).
#' @param row.w \emph{MFA Parameter:} An optional row weights (by default, a vector of 1 for uniform row weights).
#' @param which Choose one or more plots to draw:\cr 1. Percentage of variance explained by factors\cr 2. Loadings for reference compounds\cr 3. Factor 1 VS Factor 2 : Loadings & Genes\cr 4. Genes Scores for Factor \code{factor.plot}\cr 5. Compound Loadings for Factor \code{factor.plot} \cr 6. Profile plot (see \code{profile.type})\cr 7. CS Rank Scores for Factor \code{factor.plot}
#' @param factor.plot Which factor (only one) should be investigated? If \code{NULL}, you can choose a factor of interest interactively from reference loadings plot.
#' @param column.interest Numeric vector of indices of query columns which should be in the profiles plot (\code{which=5}). If \code{NULL}, you can interactively select genes in the Compound Loadings plot (\code{which=4}).
#' @param row.interest Numeric vector of gene indices to be plotted in gene profiles plot (\code{which=6}, \code{profile.type="gene"}). If \code{NULL}, you can interactively select them in the gene scores plot (\code{which=4}).
#' @param profile.type Type of \code{which=6} plot:\cr - \code{"gene"}: Gene profiles plot of selected genes in \code{row.interest} with the reference compounds and those selected in \code{column.interest} ordered first on the x axis. The other compounds are ordered in decreasing CScore. \cr - \code{"cmpd"}: Compound profiles plot of reference and selected compounds (\code{column.interest}) and only those genes on the x-axis which beat the thresholds (\code{gene.thresP}, \code{gene.thresN})
#' @param color.columns Vector of colors for the reference and query columns (compounds). If \code{NULL}, blue will be used for reference and black for query. Use this option to highlight reference columns and query columns of interest.
#' @param gene.highlight Single numeric vector or list of maximum 5 numeric vectors. This highlights gene of interest in gene scores plot (\code{which=4}) up to 5 different colors.
#' @param gene.thresP Threshold for genes with a high score (\code{which=4}).
#' @param gene.thresN Threshold for genes with a low score (\code{which=4}).
#' @param thresP.col Color of genes above \code{gene.thresP}.
#' @param thresN.col Color of genes below \code{gene.thresN}.
#' @param legend.names Option to draw a legend of colored columns in Compound Loadings plot (\code{which=5}). If \code{NULL}, only "References" will be in the legend.
#' @param legend.cols Colors to be used in legend in \code{which=5}. If \code{NULL}, only blue for "References is used".
#' @param result.available You can a previously returned object by \code{CSanalysis} in order to only draw graphs, not recompute the scores. If this object also contains the permutation object, in the loadings and CSRank plot (\code{which=c(5,7)}) the values with a (adjusted) pvalue smaller than 0.05 will be colored purple.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param basefilename Basename of the graphs if saved in pdf files
#' @param CSrank.refplot Logical value deciding if the CS Rank Scores (\code{which=7}) should also be plotted per reference (instead of only the weighted mean).
#' @return An object of the S4 Class \code{CSresult}.
setMethod("CSanalysis",c("matrix","matrix","CSmfa"),function(
				refMat,querMat,type="CSmfa",ncp=5,weight.col.mfa=NULL,row.w=NULL,
#				which=c(1,2,3,4,5,6),
				which=c(2,4,5,7),
				factor.plot=NULL,
				column.interest=NULL,row.interest=NULL,profile.type="gene",
				color.columns=NULL,gene.highlight=NULL,				
				gene.thresP=1,gene.thresN=-1,thresP.col="blue",thresN.col="red",
				legend.names=NULL,legend.cols=NULL,
				CSrank.refplot=FALSE,
				result.available=NULL,
				basefilename="analyseMFA",plot.type="device"
				){
				
			if(!(dim(refMat)[2]>1)){
				stop("Reference matrix should have more than 1 reference")
			}
			
			data <- cbind(as.matrix(refMat),as.matrix(querMat))
			
			
			if(is.null(color.columns)){colour.columns <- c(rep("blue",dim(refMat)[2]),rep("black",dim(querMat)[2]))} else {colour.columns <- color.columns}
			if(is.null(legend.names)){legend.names <- c("References")}
			if(is.null(legend.cols)){legend.cols <- c("blue")}
			
			analysis.pm.temp <- list(ncp=ncp,weight.col.mfa=weight.col.mfa,row.w=row.w)
			
			if(!is.null(result.available)){
				if(class(result.available) != "CSresult"){
					stop("result.available is not of the correct class") 
				}
				else if(result.available@type != class(type)){
					stop("result.available is from a different analysis")
				}
				else if(!identical(result.available@call$analysis.pm,analysis.pm.temp)){
					stop("result.available used other analysis parameters")
				}
				else{result.available2 <- result.available}
			}
			else{
				result.available2 <- NULL
			}
			
			if(!is.null(column.interest)){column.interest <- column.interest + dim(refMat)[2]}
			
			out <- analyse_MFA(data=data,group=c(ncol(refMat),ncol(querMat)),type=rep("s",2),ind.sup=NULL,ncp=ncp,name.group=c("Reference","Query"),num.group.sup=NULL,graph=FALSE,weight.col.mfa=weight.col.mfa,row.w=row.w,axes=c(1,2),tab.comp=NULL,
					basefilename=basefilename,
					factor.plot=factor.plot,column.interest=column.interest,row.interest=row.interest,gene.thresP=gene.thresP,gene.thresN=gene.thresN,
					colour.columns=colour.columns,legend.names=legend.names,legend.cols=legend.cols,thresP.col=thresP.col,thresN.col=thresN.col,
					result.available=result.available2,plot.type=plot.type,CSrank.refplot=CSrank.refplot,
					which=which,
					gene.highlight=gene.highlight,profile.type=profile.type)
			
			factor.select <- out[[1]]
			result <- out[[2]]
			
			loadings <- result$quanti.var$coord	
			scores <- result$ind$coord			
			
			CS.loadings <- data.frame(loadings[,factor.select])
			CS.query <- data.frame(CS.loadings[-c(1:dim(refMat)[2]),])
			CS.ref <- data.frame(CS.loadings[c(1:dim(refMat)[2]),])
			
			colnames(CS.query) <- colnames(CS.ref) <- sapply(factor.select,FUN=function(x){paste0("Factor",x)})
			rownames(CS.query) <- rownames(loadings)[-c(1:dim(refMat)[2])]
			rownames(CS.ref) <- rownames(loadings)[c(1:dim(refMat)[2])]	
				
			CS <- list(CS.query=CS.query,CS.ref=CS.ref)
			
			GS <- data.frame(scores[,factor.select])
			rownames(GS) <- rownames(scores)
			colnames(GS) <- sapply(factor.select,FUN=function(x){paste0("Factor",x)})
			
			call.object <- list(match.call=type@call,analysis.pm=list(ncp=ncp,weight.col.mfa=weight.col.mfa,row.w=row.w),factor.select=factor.select)
			

			if(!is.null(result.available)){
				if(factor.select != result.available@call$factor.select){
					if(!is.null(result.available@permutation.object)){
						warning("CS, GS and CSRankScores Slot in CSresult will be overwritten due to different factor choice. To add p-values again, run CSpermute again on this object!\nThe permutation.object was also deleted due to a different choice of factor.")
						
					}else{
						warning("CS, GS and CSRankScores Slot in CSresult will be overwritten due to different factor choice.")
						
					}
					
					return(new("CSresult",type="CSmfa",CS=CS,CSRankScores=out[[3]],GS=GS,object=result,permutation.object=NULL,call=call.object))
					
				}
				else{
					return(new("CSresult",type="CSmfa",CS=result.available@CS,CSRankScores=result.available@CSRankScores,GS=result.available@GS,object=result.available@object,permutation.object=result.available@permutation.object,call=call.object))
					
				}
			
			}
			else{
				return(new("CSresult",type="CSmfa",CS=CS,CSRankScores=out[[3]],GS=GS,object=result,permutation.object=NULL,call=call.object))
				
			}
			
			

					
		})



# PCA
#' "CSpca"
#' 
#' Doing interactive CS analysis with PCA (Principal Component Analysis). This analysis is meant for 1 reference signature.
#' 
#' @export 
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type \code{"CSpca"}
#' @param ncp \emph{PCA Parameter:} Number of dimensions kept in the results (by default 5).
#' @param scale.unit \emph{PCA Parameter:} A boolean, if TRUE (value set by default) then data are scaled to unit variance.
#' @param row.w \emph{PCA Parameter:} An optional row weights (by default, a vector of 1 for uniform row weights).
#' @param col.w \emph{PCA Parameter:} An optional column weights (by default, uniform column weights).
#' @param which Choose one or more plots to draw\cr 1. Percentage of variance explained by PC's\cr 2. Loadings for reference compounds\cr 3. PC 1 VS PC 2 : Loadings & Genes\cr 4. Genes for PC \code{factor.plot}\cr 5. Loadings for PC \code{factor.plot}\cr 6. Profile plot (see \code{profile.type})\cr 7. CS Rank Scores for Factor \code{factor.plot}
#' @param factor.plot Which PC (only one) should be investigated? If \code{NULL}, you can choose a PC of interest interactively from reference loadings plot.
#' @param column.interest Numeric vector of indices of query columns which should be in the profiles plot (\code{which=5}). If \code{NULL}, you can interactively select genes on the Compound Loadings plot (\code{which=4}).
#' @param row.interest Numeric vector of gene indices to be plotted in gene profiles plot (\code{which=6}, \code{profile.type="gene"}). If \code{NULL}, you can interactively select them in the gene scores plot (\code{which=4}).
#' @param profile.type Type of \code{which=6} plot:\cr - \code{"gene"}: Gene profiles plot of selected genes in \code{row.interest} with the reference compounds and those selected in \code{column.interest} ordered first on the x axis. The other compounds are ordered in decreasing CScore. \cr - \code{"cmpd"}: Compound profiles plot of reference and selected compounds (\code{column.interest}) and only those genes on the x-axis which beat the thresholds (\code{gene.thresP}, \code{gene.thresN})
#' @param color.columns Vector of colors for the reference and query columns (compounds). If \code{NULL}, blue will be used for reference and black for query. Use this option to highlight reference columns and query columns of interest.
#' @param gene.highlight Single numeric vector or list of maximum 5 numeric vectors. This highlights gene of interest in gene scores plot (\code{which=4}) up to 5 different colors.
#' @param gene.thresP Threshold for genes with a high score (\code{which=4}).
#' @param gene.thresN Threshold for genes with a low score (\code{which=4}).
#' @param thresP.col Color of genes above \code{gene.thresP}.
#' @param thresN.col Color of genes below \code{gene.thresN}.
#' @param legend.names Option to draw a legend of colored columns in Compound Loadings plot (\code{which=5}). If \code{NULL}, only "References" will be in the legend.
#' @param legend.cols Colors to be used in legend in \code{which=5}. If \code{NULL}, only blue for "References is used".
#' @param result.available You can a previously returned object by \code{CSanalysis} in order to only draw graphs, not recompute the scores.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param basefilename Basename of the graphs if saved in pdf files
#' @param CSrank.refplot Logical value deciding if the CS Rank Scores (\code{which=7}) should also be plotted per reference (instead of only the weighted mean).
#' @return An object of the S4 Class \code{CSresult}.
setMethod("CSanalysis",c("matrix","matrix","CSpca"),function(
				refMat,querMat,type="CSpca",ncp=5,scale.unit=TRUE,row.w=NULL,col.w=NULL,
				which=c(2,4,5,7),
				factor.plot=NULL,
				column.interest=NULL,row.interest=NULL,profile.type="gene",
				color.columns=NULL,gene.highlight=NULL,				
				gene.thresP=1,gene.thresN=-1,thresP.col="blue",thresN.col="red",
				legend.names=NULL,legend.cols=NULL,
				CSrank.refplot=FALSE,
				result.available=NULL,
				basefilename="analysePCA",plot.type="device"
				){
			
			if((dim(refMat)[2]!=1)){
				stop("Reference matrix should have only 1 reference")
			}
					
			data <- cbind(as.matrix(refMat),as.matrix(querMat))
									
			if(is.null(color.columns)){colour.columns <- c(rep("blue",dim(refMat)[2]),rep("black",dim(querMat)[2]))} else {colour.columns <- color.columns}
			if(is.null(legend.names)){legend.names <- c("References")}
			if(is.null(legend.cols)){legend.cols <- c("blue")}
					
			analysis.pm.temp <- list(ncp=ncp,scale.unit=scale.unit,row.w=row.w,col.w=col.w)
			
			if(!is.null(result.available)){
				if(class(result.available) != "CSresult"){
					stop("result.available is not of the correct class") 
				}
				else if(result.available@type != class(type)){
					stop("result.available is from a different analysis")
				}
				else if(!identical(result.available@call$analysis.pm,analysis.pm.temp)){
					stop("result.available used other analysis parameters")
				}
				else{result.available2 <- result.available}
			}
			else{
				result.available2 <- NULL
			}
			
			if(!is.null(column.interest)){column.interest <- column.interest + dim(refMat)[2]}
				
			out <- analyse_PCA(data=data, scale.unit = scale.unit, ncp = ncp, ind.sup = NULL,
					quanti.sup = NULL, quali.sup = NULL, row.w = row.w,
					col.w = col.w, graph = FALSE, axes = c(1,2),
					basefilename=basefilename,
					ref.index=1,factor.plot=factor.plot,column.interest=column.interest,gene.thresP=gene.thresP,gene.thresN=gene.thresN,
					colour.columns=colour.columns,legend.names=legend.names,legend.cols=legend.cols,thresP.col=thresP.col,thresN.col=thresN.col,
					result.available=result.available2,plot.type=plot.type,which=which,CSrank.refplot=CSrank.refplot,
					profile.type=profile.type,gene.highlight=gene.highlight)
			
			factor.select <- out[[1]]
			result <- out[[2]]
			
			loadings <- result$var$coord
			scores <- result$ind$coord		
			
			CS.loadings <- data.frame(loadings[,factor.select])
			CS.query <- data.frame(CS.loadings[-c(1:dim(refMat)[2]),])
			CS.ref <- data.frame(CS.loadings[c(1:dim(refMat)[2]),])
			
			colnames(CS.query) <- colnames(CS.ref) <- sapply(factor.select,FUN=function(x){paste0("PC",x)})
			rownames(CS.query) <- rownames(loadings)[-c(1:dim(refMat)[2])]
			rownames(CS.ref) <- rownames(loadings)[c(1:dim(refMat)[2])]	
			
			CS <- list(CS.query=CS.query,CS.ref=CS.ref)
			
			GS <- data.frame(scores[,factor.select])
			rownames(GS) <- rownames(scores)
			colnames(GS) <- sapply(factor.select,FUN=function(x){paste0("PC",x)})
			
			call.object <- list(match.call=type@call,analysis.pm=list(ncp=ncp,scale.unit=scale.unit,row.w=row.w,col.w=col.w),factor.select=factor.select)
			
			if(!is.null(result.available)){
				if(factor.select != result.available@call$factor.select){
					if(!is.null(result.available@permutation.object)){
						warning("CS, GS and CSRankScores Slot in CSresult will be overwritten due to different factor choice. To add p-values again, run CSpermute again on this object!\nThe permutation.object was also deleted due to a different choice of factor.")
						
					}else{
						warning("CS, GS and CSRankScores Slot in CSresult will be overwritten due to different factor choice.")
						
					}
					
					return(new("CSresult",type="CSpca",CS=CS,CSRankScores=out[[3]],GS=GS,object=result,permutation.object=NULL,call=call.object))
					
				}
				else{
					return(new("CSresult",type="CSpca",CS=result.available@CS,CSRankScores=result.available@CSRankScores,GS=result.available@GS,object=result.available@object,permutation.object=result.available@permutation.object,call=call.object))
					
				}
				
			}
			else{
				return(new("CSresult",type="CSpca",CS=CS,CSRankScores=out[[3]],GS=GS,object=result,permutation.object=NULL,call=call.object))
				
			}
			
					
		})

		
# sMFA
#' "CSsmfa"
#' 
#' Doing interactive CS analysis with sMFA (Sparse Multiple Factor Analysis). Should use multiple references for this analysis.
#' 
#' @export 
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type \code{"CSsmfa"}
#' @param K \emph{sMFA Parameters:} Number of components.
#' @param para \emph{sMFA Parameters:} A vector of length K. All elements should be positive. If \code{sparse="varnum"}, the elements integers. 
#' @param lambda \emph{sMFA Parameters:} Quadratic penalty parameter. Default value is 1e-6. If the target dimension of the sparsness is higher than the other dimension (p > n), it is advised to put \code{lambda} to \code{Inf} which uses the \code{arrayspc} algorithm optimized for this case. For the other case, p < n, a zero or positive \code{lambda} is sufficient and will utilize the normal \code{spca} algorithm.
#' @param sparse.dim \emph{sMFA Parameters:} Which dimension should be sparse? 1: Rows, 2: Columns (default) (Note: For Connectivity Scores it is advised to apply sparsity on the compounds/columns)
#' @param sparse \emph{sMFA Parameters (\code{lambda < Inf} only):} If \code{sparse="penalty"}, \code{para} is a vector of 1-norm penalty parameters. If \code{sparse="varnum"}, \code{para} defines the number of sparse loadings to be obtained.
#' @param max.iter \emph{sMFA Parameters:} Maximum number of iterations.
#' @param eps.conv \emph{sMFA Parameters:} Convergence criterion.
#' @param which Choose one or more plots to draw.\cr 1. Percentage of variance explained by PC's\cr 2. Loadings for reference compounds\cr 3. Factor 1 VS PC Factor 2 : Loadings & Genes\cr 4. Genes for Factor \code{factor.plot}\cr 5. Loadings for Factor \code{factor.plot}\cr 6. Profile plot (see \code{profile.type})\cr 7. CS Rank Scores for Factor \code{factor.plot}
#' @param factor.plot Which factor (only one) should be investigated? If \code{NULL}, you can choose a factor of interest interactively from reference loadings plot.
#' @param column.interest Numeric vector of indices of query columns which should be in the profiles plot (\code{which=5}). If \code{NULL}, you can interactively select genes on the Compound Loadings plot (\code{which=5}).
#' @param row.interest Numeric vector of gene indices to be plotted in gene profiles plot (\code{which=6}, \code{profile.type="gene"}). If \code{NULL}, you can interactively select them in the gene scores plot (\code{which=4}).
#' @param profile.type Type of \code{which=6} plot:\cr - \code{"gene"}: Gene profiles plot of selected genes in \code{row.interest} with the reference compounds and those selected in \code{column.interest} ordered first on the x axis. The other compounds are ordered in decreasing CScore. \cr - \code{"cmpd"}: Compound profiles plot of reference and selected compounds (\code{column.interest}) and only those genes on the x-axis which beat the thresholds (\code{gene.thresP}, \code{gene.thresN})
#' @param color.columns Vector of colors for the reference and query columns (compounds). If \code{NULL}, blue will be used for reference and black for query. Use this option to highlight reference columns and query columns of interest.
#' @param gene.highlight Single numeric vector or list of maximum 5 numeric vectors. This highlights gene of interest in gene scores plot (\code{which=4}) up to 5 different colors.
#' @param gene.thresP Threshold for genes with a high score (\code{which=4}).
#' @param gene.thresN Threshold for genes with a low score (\code{which=4}).
#' @param thresP.col Color of genes above \code{gene.thresP}.
#' @param thresN.col Color of genes below \code{gene.thresN}.
#' @param legend.names Option to draw a legend of colored columns in Compound Loadings plot (\code{which=5}). If \code{NULL}, only "References" will be in the legend.
#' @param legend.cols Colors to be used in legend in \code{which=5}. If \code{NULL}, only blue for "References is used".
#' @param result.available You can a previously returned object by \code{CSanalysis} in order to only draw graphs, not recompute the scores.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param basefilename Basename of the graphs if saved in pdf files
#' @param CSrank.refplot Logical value deciding if the CS Rank Scores (\code{which=7}) should also be plotted per reference (instead of only the weighted mean).
#' @return An object of the S4 Class \code{CSresult}.
setMethod("CSanalysis",c("matrix","matrix","CSsmfa"),function(
				refMat,querMat,type="Csmfa",K=15,para,lambda=1e-6,sparse.dim=2,sparse="penalty",max.iter=200,eps.conv=1e-3,
				which=c(2,4,5,7),
				factor.plot=NULL,
				column.interest=NULL,row.interest=NULL,profile.type="gene",
				color.columns=NULL,gene.highlight=NULL,
				gene.thresP=1,gene.thresN=-1,thresP.col="blue",thresN.col="red",
				legend.names=NULL,legend.cols=NULL,
				CSrank.refplot=FALSE,
				result.available=NULL,
				basefilename="analysesMFA",plot.type="device"
				){
			
			if(!(dim(refMat)[2]>1)){
				stop("Reference matrix should have more than 1 reference")
			}
					
			data <- cbind(as.matrix(refMat),as.matrix(querMat))
					
					
			if(is.null(color.columns)){colour.columns <- c(rep("blue",dim(refMat)[2]),rep("black",dim(querMat)[2]))} else {colour.columns <- color.columns}
			if(is.null(legend.names)){legend.names <- c("References")}
			if(is.null(legend.cols)){legend.cols <- c("blue")}
					
			analysis.pm.temp <- list(K=K,lambda=lambda,sparse=sparse,max.iter=max.iter,eps.conv=eps.conv,para=para,sparse.dim=sparse.dim)
			
			if(!is.null(result.available)){
				if(class(result.available) != "CSresult"){
					stop("result.available is not of the correct class") 
				}
				else if(result.available@type != class(type)){
					stop("result.available is from a different analysis")
				}
				else if(!identical(result.available@call$analysis.pm,analysis.pm.temp)){
					stop("result.available used other analysis parameters")
				}
				else{result.available2 <- result.available}
			}
			else{
				result.available2 <- NULL
			}
			
			if(!is.null(column.interest)){column.interest <- column.interest + dim(refMat)[2]}
			
			out <- analyse_sMFA(data=data,K=K,para=para,type="predictor",sparse=sparse,use.corr=FALSE,lambda=lambda,max.iter=max.iter,trace=FALSE,eps.conv=eps.conv,sparse.dim=sparse.dim,
					basefilename=basefilename,ref.index=c(1:dim(refMat)[2]),
					factor.plot=factor.plot,column.interest=column.interest,gene.thresP=gene.thresP,gene.thresN=gene.thresN,
					colour.columns=colour.columns,legend.names=legend.names,legend.cols=legend.cols,thresP.col=thresP.col,thresN.col=thresN.col,
					result.available=result.available2,plot.type=plot.type,CSrank.refplot=CSrank.refplot,
					which=which,profile.type=profile.type,gene.highlight=gene.highlight,row.interest=row.interest)
			
			factor.select <- out[[1]]
			result <- out[[2]]
			
			loadings <- result$loadings	
			scores <- result$scores	
			
			CS.loadings <- data.frame(loadings[,factor.select])
			CS.query <- data.frame(CS.loadings[-c(1:dim(refMat)[2]),])
			CS.ref <- data.frame(CS.loadings[c(1:dim(refMat)[2]),])
			
			colnames(CS.query) <- colnames(CS.ref) <- sapply(factor.select,FUN=function(x){paste0("Factor",x)})
			rownames(CS.query) <- rownames(loadings)[-c(1:dim(refMat)[2])]
			rownames(CS.ref) <- rownames(loadings)[c(1:dim(refMat)[2])]	
			
			CS <- list(CS.query=CS.query,CS.ref=CS.ref)
			
			GS <- data.frame(scores[,factor.select])
			rownames(GS) <- rownames(scores)
			colnames(GS) <- sapply(factor.select,FUN=function(x){paste0("Factor",x)})
			
			call.object <- list(match.call=type@call,analysis.pm=list(K=K,lambda=lambda,sparse=sparse,max.iter=max.iter,eps.conv=eps.conv,para=para,sparse.dim=sparse.dim),factor.select=factor.select)
			
			if(!is.null(result.available)){
				if(factor.select != result.available@call$factor.select){
					if(!is.null(result.available@permutation.object)){
						warning("CS, GS and CSRankScores Slot in CSresult will be overwritten due to different factor choice. To add p-values again, run CSpermute again on this object!\nThe permutation.object was also deleted due to a different choice of factor.")
						
					}else{
						warning("CS, GS and CSRankScores Slot in CSresult will be overwritten due to different factor choice.")
						
					}
					
					return(new("CSresult",type="CSsmfa",CS=CS,CSRankScores=out[[3]],GS=GS,object=result,permutation.object=NULL,call=call.object))
					
				}
				else{
					return(new("CSresult",type="CSsmfa",CS=result.available@CS,CSRankScores=result.available@CSRankScores,GS=result.available@GS,object=result.available@object,permutation.object=result.available@permutation.object,call=call.object))
					
				}
				
			}
			else{
				return(new("CSresult",type="CSsmfa",CS=CS,CSRankScores=out[[3]],GS=GS,object=result,permutation.object=NULL,call=call.object))
				
			}
			
			
			
		})


# Zhang and Gant
#' "CSzhang"
#' 
#' Compute the Connectivity Scores by Zhang and Gant (2008). One or multiple reference compounds are possible in this analysis.
#' 
#' @export 
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type \code{"CSzhang"}
#' @param nref \emph{Zhang Parameter:} Number of top up- and downregulated genes in reference signature. If \code{NULL}, all rows (genes) are used.
#' @param nquery \emph{Zhang Parameter:} Number of top up- and downregulated genes in query signature. If \code{NULL}, all rows (genes) are used. (Note that \eqn{nref >= nquery})
#' @param ord.query \emph{Zhang Parameter:} Logical value. Should the query signature be treated as ordered?
# @param permute \emph{Zhang Parameter:} Logical value. Should p-values be computed through permutation?
# @param B \emph{Zhang Parameter:} Number of permutations for p-value calculation.
# @param ntop.pvalues \emph{Zhang Parameter:} Number of top p-values to be reported first. 
#' @param ntop.scores \emph{Zhang Parameter:} Number of top positive and negative CS to be reported first.
#' @param color.query Vector of colors for the query columns. You can use this option to highlight columns(compounds) of interest in the CS plot. (This does not include the reference columns since they are not included in the CS plot.)
#' @param legend.names Option to draw a legend (about the highlights in \code{color.query}) in the CS plot. If \code{NULL}, no legend will be drawn.
#' @param legend.cols Colors to be used for the \code{legend.names}.
#' @param result.available You can a previously returned object by \code{CSanalysis} in order to only draw graphs, not recompute the scores. If this object also contains the permutation object, in the score plot the values with a (adjusted) pvalue smaller than 0.05 will be colored purple.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param basefilename Basename of the graphs if saved in pdf files
#' @return An object of the S4 Class \code{CSresult}. The CS slot will also contain the top positive and negative scores as well as the top p-values. The GS slot will be empty for Zhang and Gant.
setMethod("CSanalysis",c("matrix","matrix","CSzhang"),function(refMat,querMat,type="CSzhang",
				nref=NULL,nquery=NULL,ord.query=TRUE,ntop.scores=20,
#				B=100000,ntop.pvalues=20,permute=FALSE,
				color.query=NULL,
				legend.names=NULL,legend.cols=NULL,
				result.available=NULL,plot.type="device",basefilename="analyseZhang"
				){
			
			
			colour.query <- color.query
			if(is.null(legend.cols)){legend.cols <- "black"}
			
			if(!is.null(result.available)){
				if(class(result.available) != "CSresult"){
					stop("result.available is not of the correct class") 
				}
				else if(result.available@type != class(type)){
					stop("result.available is from a different analysis")
				}
				else{result.available2 <- result.available}
			}
			else{
				result.available2 <- NULL
			}
			
			out <- analyse_zhang(dataref=refMat,dataquery=querMat,nref=nref,nquery=nquery,ord.query=ord.query,ntop.pvalues=20,ntop.scores=ntop.scores,
#					permute=permute,B=B,ntop.pvalues=ntop.values,
					basefilename=basefilename,
					colour.query=colour.query,legend.names=legend.names,legend.cols=legend.cols,legend.x=NULL,legend.y=NULL,
					result.available=result.available2,plot.type=plot.type,
					which=c(1))
			
			CS.ref <-  data.frame(ZhangScore=rep(1,dim(refMat)[2]))
			rownames(CS.ref) <- colnames(refMat)
					
#			CS <- list(CS.query=out$All,CS.ref=CS.ref,CS.top.query=list(TopCS=out$Top,TopPvalues=out$Toppvalues))
			CS <- list(CS.query=out$All,CS.ref=CS.ref,CS.top.query=out$Top)
			
			
			call.object <- list(match.call=type@call,analysis.pm=list(nref=nref,nquery=nquery,ord.query=ord.query,ntop.scores=ntop.scores))
			
			if(!is.null(result.available)){
				if(!is.null(result.available@permutation.object)){
					
#					if(identical(result.available@permutation.object$analysis.perm.pm , call.object$analysis.pm   )){
						warning("CS Slot in CSresult will be overwritten. To add p-values again, run CSpermute again on this object!")
						return(new("CSresult",type="CSzhang",CS=CS,GS=data.frame(),object=out,permutation.object=result.available@permutation.object,call=call.object))
						
#					}
#					else{
#						warning("CS Slot in CSresult will be overwritten. To add p-values again, run CSpermute again on this object!\nThe permutation.object was also deleted due to a different choice of method parameters.")
#						return(new("CSresult",type="CSzhang",CS=CS,GS=data.frame(),object=out,permutation.object=NULL,call=call.object))
						
#					}
					
				}
				
				
			}
			else{
				return(new("CSresult",type="CSzhang",CS=CS,GS=data.frame(),object=out,permutation.object=NULL,call=call.object))
			}
			
		})


## METHODS - CSCOMPARE ##
##' Compare CS Results.
##' 
##' After applying different CSanalysis on the same data, you can compare the connectivity scores (and gene scores), based on a certain component (factor, PC, bicluster).
##' Right now it is possible to compare Zhang with any other result and Fabia with any other result.
##' 
##' @export
##' @param CSresult1 First result.
##' @param CSresult2 Second result.
##' @param ... Additional parameters for analysis, depending on which results you are comparing.
##' @return Vector with correlation between Connectivity Scores (and Gene Scores if available).
#setGeneric('CScompare', function(CSresult1,CSresult2, ...){standardGeneric('CScompare')})
#
##' Compare CS Results.
##' 
##' After applying different CSanalysis on the same data, you can compare the connectivity scores (and gene scores), based on a certain component (factor, PC, bicluster).
##' Right now it is possible to compare Zhang with any other result and Fabia with any other result.
##' 
##' @export
##' @param CSresult1 First result.
##' @param CSresult2 Second result.
##' @param ... Additional parameters for analysis, depending on which results you are comparing.
##' @return Vector with correlation between Connectivity Scores (and Gene Scores if available).
#setMethod("CScompare",c("CSresult","CSresult"),function(CSresult1,CSresult2,...){
#			
#			type1 <- CSresult1@type
#			type2 <- CSresult2@type
#			
#			if("CSzhang" %in% c(type1,type2)){
#				if(type1=="CSzhang"){
#					CSresultSet <- new("CSzhangCompare",CSresult1=CSresult1,CSresult2=CSresult2)
#				}
#				else{
#					CSresultSet <- new("CSzhangCompare",CSresult1=CSresult2,CSresult2=CSresult1)
#					
#				}
#				
#				CScompare(CSresultSet,CSresult1,...)
#				
#			}
#			
#			else if("CSfabia" %in% c(type1,type2)){
#				
#				if(type1=="CSfabia"){
#					CSresultSet <- new("CSfabiaCompare",CSresult1=CSresult1,CSresult2=CSresult2)
#				}
#				else{
#					CSresultSet <- new("CSfabiaCompare",CSresult1=CSresult2,CSresult2=CSresult1)
#				}
#				
#			}
#			
#			else{
#				stop("This comparison is not yet possible")
#			}
#			
#			
#		})
#
## Zhang with...
##' Compare Zhang and Gant CS with...
##' 
##' This compare method will be used if you are comparing a Zhang and Gant result with any of the other results (including Zhang and Gant). 
##' The connectivity scores will be drawn on a scatter plot and the pearson correlation will be returned.
##' @export
##' @param CSresult1 First result.
##' @param CSresult2 Second result.
##' @param component.plot If you are using a non-Zhang&Gant result, specify the bicluster, factor or principal component which should be used to derive connectivity scores from.
##' @param color.columns Vector of colors for the reference and query columns (compounds). If \code{NULL}, blue will be used for reference and black for query. Use this option to highlight reference columns and query columns of interest.
##' @param legend.names Option to draw a legend (about the highlights in \code{color.columns}) in the CS plot. If \code{NULL}, only references are in the legend.
##' @param legend.cols Colors to be used for the \code{legend.names}.
##' @param legend.pos The location of the legend: \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"} and \code{"center"}.
##' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
##' @param basefilename Basename of the graphs if saved in pdf files
##' @return Numeric with correlation between Connectivity Scores.
#setMethod("CScompare",c("CSzhangCompare","CSresult"),function(CSresult1,CSresult2,
#				component.plot,
#				color.columns=NULL,
#				legend.names=NULL,legend.cols=NULL,legend.pos="topright",
#				plot.type="device",basefilename="ZhangCompare"){
#			
#			
#			result1 <- CSresult1@CSresult1
#			result2 <- CSresult1@CSresult2
#			
#			if(missing(component.plot) & !(result1@type=="CSzhang" & result2@type=="CSzhang")){stop("argument 'component.plot' is missing, with no default")}
#			
#			
#			refdim1 <- dim(result1@CS$CS.ref)[1]
#			refdim2 <- dim(result2@CS$CS.ref)[1]
#			querdim1 <- dim(result1@CS$CS.query)[1]
#			querdim2 <- dim(result2@CS$CS.query)[1]
#			
#			if(refdim1!=refdim2){stop("Using 2 different reference matrices")}
#			if(querdim1!=querdim2){stop("Using 2 different query matrices")}
#			
#			ref.index <- c(1:refdim1)
#			
#			if(is.null(color.columns)){colour.columns <- c(rep("blue",refdim1),rep("black",querdim1))} else {colour.columns <- color.columns}
#			if(is.null(legend.names)){legend.names <- c("References")}
#			if(is.null(legend.cols)){legend.cols <- c("blue")}
#			
#			# Make dummy data (for dimensions and colnames for older functions
#			if(result2@type=="CSzhang"){
#				data <- matrix(0,nrow=1,ncol=refdim1+querdim1)
#				colnames(data) <- c(rownames(result2@CS$CS.ref),rownames(result2@CS$CS.query))
#			}
#			else{
#				data <- matrix(0,nrow=dim(result2@GS)[1],ncol=refdim1+querdim1)
#				rownames(data) <- rownames(result2@GS)
#				colnames(data) <- c(rownames(result2@CS$CS.ref),rownames(result2@CS$CS.query))
#			}
#			
#			if(result2@type=="CSfabia"){
#				out <- analyse_Zhang_fabia(data=data,resZhang=result1@object,resFAB=result2@object,
#						ref.index=ref.index,BC.plot=component.plot,
#						colour.columns=colour.columns,legend.names=legend.names,legend.cols=legend.cols,legend.pos=legend.pos,
#						plot.type=plot.type,basefilename=basefilename)
#				return(out)
#			}
#			if(result2@type=="CSmfa"){
#				
#				out <- analyse_Zhang_MFAPCA(data=data,resZhang=result1@object,resPCA=NULL,resMFA=result2@object,ressparseMFA=NULL,
#						ref.index=ref.index,factor.plot=component.plot,
#						colour.columns=colour.columns,legend.names=legend.names,legend.cols=legend.cols,legend.pos=legend.pos,
#						plot.type=plot.type,basefilename=basefilename)
#				return(out)
#				
#			}
#			if(result2@type=="CSpca"){
#				out <- analyse_Zhang_MFAPCA(data=data,resZhang=result1@object,resPCA=result2@object,resMFA=NULL,ressparseMFA=NULL,
#						ref.index=ref.index,factor.plot=component.plot,
#						colour.columns=colour.columns,legend.names=legend.names,legend.cols=legend.cols,legend.pos=legend.pos,
#						plot.type=plot.type,basefilename=basefilename)
#				return(out)
#			}
#			if(result2@type=="CSsmfa"){
#				out <- analyse_Zhang_MFAPCA(data=data,resZhang=result1@object,resPCA=NULL,resMFA=NULL,ressparseMFA=result2@object,
#						ref.index=ref.index,factor.plot=component.plot,
#						colour.columns=colour.columns,legend.names=legend.names,legend.cols=legend.cols,legend.pos=legend.pos,
#						plot.type=plot.type,basefilename=basefilename)
#				return(out)
#			}
#			
#			if(result2@type=="CSzhang"){
#				out <- analyse_Zhang_MFAPCA(data=data,resZhang=result1@object,resPCA=NULL,resMFA=NULL,resZhang2=result2@object,
#						ref.index=ref.index,factor.plot=1,
#						colour.columns=colour.columns,legend.names=legend.names,legend.cols=legend.cols,legend.pos=legend.pos,
#						plot.type=plot.type,basefilename=basefilename)
#				return(out)
#			}
#			
#		})
#
#
#
## Fabia with...
##' Compare Fabia CS and GS with...
##' 
##' This compare method will be used if you are comparing a Fabia result with MFA, PCA, sMFA or another Fabia result. 
##' The connectivity scores and gene scores will be drawn on a scatter plot and the pearson correlation of both will be returned.
##' @export
##' @param CSresult1 First result.
##' @param CSresult2 Second result.
##' @param component1.plot Component (Bicluster, factor or principal component) of the first result of which CS and GS are derived.
##' @param component2.plot Component (Bicluster, factor or principal component) of the second result of which CS and GS are derived.
#setMethod("CScompare",c("CSfabiaCompare","CSresult"),function(CSresult1,CSresult2,
#				component1.plot,component2.plot,
#				gene.thresP=c(1,1),gene.thresN=c(-1,-1),thresP.col=c("blue","light blue"),thresN.col=c("red","pink"),
#				
#				){
#			#if(com)
#			
#		})
#
######## FABIA vs PCA/MFA
#### which:	-1: Compound Loadings Plot
####			-2: Gene Scores Plot
##analyse_fabia_MFAPCA <- function(data,resFAB,resPCA=NULL,resMFA=NULL,ressparseMFA=NULL,resFAB2=NULL,
##		ref.index,BC.plot,factor.plot,BC2.plot,
##		# Give 2 Pos and 2 Neg threshold (1: FAB, 2:PCA/MFA)
##		# Give 3 Pos and 3 Neg Colors (1:FAB, 2:PCA/MFA)
#		gene.thresP=NULL,gene.thresN=NULL,thresP.col=c("blue","light blue"),thresN.col=c("red","pink"),
##		colour.columns=NULL,legend.names=NULL,legend.cols=unique(colour.columns),legend.pos='topright',
##		plot.type="pdf",basefilename="FABIA_PCAMFA",which=c(1,2)){
#
#
## Make classes like in biclust package: CSfabia, CSmfa, CSzhang -> these have a @function inside


# HOW TO DO COMPARE: All results same class "CSresult"
# compare(CSresult,CSresult,...) -> combine both results in a set and give class FABandMFA..
# compare(FABandMFA,function(FABandMFA,ref.index))