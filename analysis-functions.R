# ICW Index function

# Function to standardize columns of a matrix
# where you designate a standardization group
# (e.g., the control group in an experiment)
# with "sgroup", a logical vector.

matStand <- function(x, sgroup = rep(TRUE, nrow(x))){
  for(j in 1:ncol(x)){
    x[,j] <- (x[,j] - mean(x[sgroup,j]))/sd(x[sgroup,j])
  }
  return(x)
}

# Function that takes in data in matrix format and returns
# (i) IC weights and (ii) ICW index.
# Weights can be incorporated using the "wgts" argument.
# The "revcols" argument takes a vector indicating which columns,
# if any, should have values reversed (that is, standardized 
# values are multiplied by -1) prior to construction of the index. 

icwIndex <- function(	xmat,
                      wgts=rep(1, nrow(xmat)),
                      revcols = NULL,
                      sgroup = rep(TRUE, nrow(xmat))){
  X <- matStand(xmat, sgroup)
  if(length(revcols)>0){
    X[,revcols] <-  -1*X[,revcols]
  }
  i.vec <- as.matrix(rep(1,ncol(xmat)))
  Sx <- cov.wt(X, wt=wgts)[[1]]
  weights <- solve(t(i.vec)%*%solve(Sx)%*%i.vec)%*%t(i.vec)%*%solve(Sx)
  index <- t(solve(t(i.vec)%*%solve(Sx)%*%i.vec)%*%t(i.vec)%*%solve(Sx)%*%t(X))
  return(list(weights = weights, index = index))
}

# For printing R screen output in HTML

print_output <- function(output, cex = 0.45, yrange = c(0,1)) {
  tmp <- capture.output(output)
  plot.new()
  plot.window(xlim=c(0,1), ylim=yrange)
  text(0, 1, paste(tmp, collapse='\n'), adj = c(0,1), family = 'mono', cex = cex)
  box()
}

# Properly formatted cross-tab

crossTab <- function(rowvar, colvar, rowlab, collab){
  crosstab <- table(rowvar, colvar)
  colproptab <- crosstab/as.matrix(rep(1,
                             nrow(crosstab)))%x%t(as.matrix(apply(crosstab, 
                                                                  2, sum)))
  forPrint.sc <- rbind(c("",colnames(crosstab),""),
                       cbind(rbind(cbind(c(rbind(rownames(crosstab), 
                                                 rep("", nrow(crosstab)))),
                                         interleave(matrix(as.character(crosstab), 
                                                           nrow=nrow(crosstab)),
                                                    matrix(as.character(round(colproptab, 2)), 
                                                           nrow=nrow(colproptab)))),
                                   c("",apply(crosstab, 2, sum))),
                             c(rbind(apply(crosstab, 1, sum), 
                                     rep("",nrow(crosstab))),
                               sum(crosstab))))
  forPrint <- rbind(c("","",collab,rep("", ncol(forPrint.sc)-2)),
                    cbind(c("",rowlab,rep("", nrow(forPrint.sc)-2)), 
                          forPrint.sc))
  colnames(forPrint) <- rep("", ncol(forPrint))
  crossTabresult <- list(crosstab, forPrint)
  return(crossTabresult)
}


# Estimating and Reporting Treatment Effects

fitR <- function(vUp, label_name=""){
  listOut <- list(NA)
  formUp <- as.formula(paste0(vUp, "~treated+as.factor(blocknum)"))
  fitHH <- lm_robust(formUp,cluster=distrito, data=subset(HH, leader==0))
  fitLE <- lm_robust(formUp,cluster=distrito, data=subset(HH, leader==1))
  listOut[[1]] <- tidy(fitHH)[2,]
  listOut[[1]]$controlmean <- mean(HH[HH$treated==0,vUp])
  listOut[[1]]$controlsd <- sd(HH[HH$treated==0,vUp])
  listOut[[1]]$N <- fitHH$N
  listOut[[1]]$label <- label_name
  listOut[[2]] <- tidy(fitLE)[2,]
  listOut[[2]]$controlmean <- mean(subset(HH, leader==1)[subset(HH, leader==1)$treated==0,vUp])
  listOut[[2]]$controlsd <- sd(subset(HH, leader==1)[subset(HH, leader==1)$treated==0,vUp])  
  listOut[[2]]$N <- fitLE$N
  listOut[[2]]$label <- label_name
    return(listOut)
}


fitR_2 <- function(vUp, label_name="", HHdata=NULL){
  listOut <- list(NA)
  formUp <- as.formula(paste0(vUp, "~treated+as.factor(blocknum)"))
  fitHH <- lm_robust(formUp,cluster=distrito, data=subset(HHdata, leader==0))
  fitLE <- lm_robust(formUp,cluster=distrito, data=subset(HHdata, leader==1))
  listOut[[1]] <- tidy(fitHH)[2,]
  listOut[[1]]$controlmean <- mean(HHdata[HHdata$treated==0,vUp])
  listOut[[1]]$controlsd <- sd(HHdata[HHdata$treated==0,vUp])
  listOut[[1]]$N <- fitHH$N
  listOut[[1]]$label <- label_name
  listOut[[2]] <- tidy(fitLE)[2,]
  listOut[[2]]$controlmean <- mean(subset(HHdata, leader==1)[subset(HHdata, leader==1)$treated==0,vUp])
  listOut[[2]]$controlsd <- sd(subset(HHdata, leader==1)[subset(HHdata, leader==1)$treated==0,vUp])  
  listOut[[2]]$N <- fitLE$N
  listOut[[2]]$label <- label_name
    return(listOut)
}

fitR_dist <- function(vUp,
                      label_name="",
                      ddata=NULL){
  formUp <- as.formula(paste0(vUp, "~treated+as.factor(blocknum)"))
  fit <- lm_robust(formUp,cluster=distrito, data=ddata)
  listOut <- tidy(fit)[2,]
  listOut$controlmean <- mean(ddata[ddata$treated==0,vUp])
  listOut$controlsd <- sd(ddata[ddata$treated==0,vUp])
  listOut$N <- fit$N
  listOut$label <- label_name
  return(listOut)
}


resVec_dist <- function(x){
  if("coefficients"%in%names(x)){
    resVecOut <- c(as.character(round(x["coefficients"], 2)),
                     paste0("(",round(x["se"],2),")"),
                     paste0("[",round(x["p"],2),"]"),
                     as.character(round(x["controlmean"], 2)),
                     paste0("(",round(x["controlsd"],2),")"),
                     as.character(round(x["N"], 0)),
                     x["label"])
  }
  if("estimate"%in%names(x)){
    resVecOut <- c(as.character(round(x["estimate"], 2)),
                     paste0("(",round(x["std.error"],2),")"),
                     paste0("[",round(x["p.value"],2),"]"),
                     as.character(round(x["controlmean"], 2)),
                     paste0("(",round(x["controlsd"],2),")"),
                     as.character(round(x["N"], 0)),
                     x["label"])
  }
  resVecOut <- as.matrix(resVecOut)
  rownames(resVecOut) <- c("Program effect","(SE)","[p]","Control mean","(Control SD)","N","Outcome")
  return(resVecOut)  
}



resVec <- function(x){
  if("coefficients"%in%names(x[[1]])){
  resVecOutHH <- c(as.character(round(x[[1]]["coefficients"], 2)),
                   paste0("(",round(x[[1]]["se"],2),")"),
                   paste0("[",round(x[[1]]["p"],2),"]"),
                   as.character(round(x[[1]]["controlmean"], 2)),
                   paste0("(",round(x[[1]]["controlsd"],2),")"),
                   as.character(round(x[[1]]["N"], 0)),
                   x[[1]]["label"])
  resVecOutLE <- c(as.character(round(x[[2]]["coefficients"], 2)),
                   paste0("(",round(x[[2]]["se"],2),")"),
                   paste0("[",round(x[[2]]["p"],2),"]"),
                   as.character(round(x[[2]]["controlmean"], 2)),
                   paste0("(",round(x[[2]]["controlsd"],2),")"),
                   as.character(round(x[[2]]["N"], 0)),
                   x[[2]]["label"])
  }
  if("estimate"%in%names(x[[1]])){
    resVecOutHH <- c(as.character(round(x[[1]]["estimate"], 2)),
                     paste0("(",round(x[[1]]["std.error"],2),")"),
                     paste0("[",round(x[[1]]["p.value"],2),"]"),
                     as.character(round(x[[1]]["controlmean"], 2)),
                     paste0("(",round(x[[1]]["controlsd"],2),")"),
                     as.character(round(x[[1]]["N"], 0)),
                     x[[1]]["label"])
    resVecOutLE <- c(as.character(round(x[[2]]["estimate"], 2)),
                     paste0("(",round(x[[2]]["std.error"],2),")"),
                     paste0("[",round(x[[2]]["p.value"],2),"]"),
                     as.character(round(x[[2]]["controlmean"], 2)),
                     paste0("(",round(x[[2]]["controlsd"],2),")"),
                     as.character(round(x[[2]]["N"], 0)),
                     x[[2]]["label"])
  }
  resVecOut <- cbind(resVecOutHH, resVecOutLE)
  colnames(resVecOut) <- c("Non-ldr.","Ldr.")
  rownames(resVecOut) <- c("Program effect","(SE)","[p]","Control mean","(Control SD)","N","Outcome")
  return(resVecOut)
}

# Create a binary variable with no missingness

makeBinary <- function(varIn, data, yesValue=1){
  varOut <- as.numeric(data[varIn] == yesValue)
  varOut[is.na(varOut)] <- 0
  return(varOut)
}


# Reverse code and clean Likert scales

revCleanLik <- function(xIn, orignegend=4, origposend=1){
  xIn[xIn==99] <- (orignegend+origposend)/2
  xIn[xIn==100] <- (orignegend+origposend)/2
  xIn[is.na(xIn)] <- (orignegend+origposend)/2
  xOut <- orignegend - xIn
}


# For cleaning numeric variables with -99 and NAs and imputing mean values

cleanNeg99mean <- function(varUp){
  varOut <- HH[,varUp]
  varOut[HH[,varUp]==-99] <- mean(varOut[HH[,varUp]!=-99], na.rm=TRUE)
  varOut[is.na(HH[,varUp])] <- mean(varOut[!is.na(HH[,varUp])], na.rm=TRUE)
  return(varOut)
}  

# For cleaning numeric variables with 99 and NAs and imputing median values

clean99median <- function(varUp){
  varOut <- HH[,varUp]
  varOut[HH[,varUp]==99] <- median(varOut[HH[,varUp]!=99], na.rm=TRUE)
  varOut[is.na(HH[,varUp])] <- median(varOut[!is.na(HH[,varUp])], na.rm=TRUE)
  return(varOut)
}  

# For cleaning numeric variables with 98, -99, and NAs and imputing modal values

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}


cleanMissMode <- function(varUp, missVals){
  varOut <- HH[,varUp]
  varOut[HH[,varUp]%in%missVals] <- NA
  varOut[is.na(varOut)] <- getmode(varOut[!is.na(varOut)])
  return(varOut)
}  


# For cleaning numeric variables with -99 and NAs and imputing 0

cleanNeg99_to_0 <- function(varUp){
  varOut <- HH[,varUp]
  varOut[HH[,varUp]==-99] <- 0
  varOut[is.na(HH[,varUp])] <- 0
  return(varOut)
}  

cleanNeg99_to_0_gen <- function(varUp, dataUp){
  varOut <- dataUp[,varUp]
  varOut[dataUp[,varUp]==-99] <- 0
  varOut[is.na(dataUp[,varUp])] <- 0
  return(varOut)
}  


# For computing roster ag income data

rostIncome <- function(	rosterIn=NULL,
						quantVar=NULL,
						priceVar=NULL,
						incomeOut=NULL){
	quantClean <- paste0("use",quantVar)
	priceClean <- paste0("use",priceVar)
	rosterIn[,quantClean] <- cleanNeg99_to_0_gen(quantVar, rosterIn)
	rosterIn[,priceClean] <- cleanNeg99_to_0_gen(priceVar, rosterIn)
	rosterIn[,incomeOut] <- rosterIn[,quantClean]*rosterIn[,priceClean]
	if(length(incomeOut)==1){
		rosterAgg <- aggregate(as.formula(paste0(incomeOut,"~parent_key")),
							data= rosterIn[,c("parent_key",incomeOut)],
							sum)
	}
	if(length(incomeOut)>1){
		rosterAgg <- aggregate(as.formula(paste0("cbind(",
												paste0(incomeOut,
														collapse=","),
												")",
												"~parent_key")),
							data= rosterIn[,c("parent_key",incomeOut)],
							sum)
	}
	HH_incomeOut <- rep(0, nrow(HH))
	HH_incomeOut[match(rosterAgg$parent_key, 
									HH$key)] <- rosterAgg[,
													incomeOut]
	return(HH_incomeOut)
}


# For computing roster ag contract data

rostContract <- function(	rosterIn=NULL,
							contrVar=NULL,
							anyOut=NULL){
	indOut <- paste0("use", contrVar, "_ind")
	colOut <- paste0("use", contrVar, "_col")
	rosterIn[, anyOut] <- 1
	rosterIn[, indOut] <- makeBinary(contrVar, rosterIn, 1)
	rosterIn[, colOut] <- makeBinary(contrVar, rosterIn, 2)
	rosterIn_agg <- aggregate(as.formula(
								paste0("cbind(",
										paste0(
											c(anyOut,
											indOut,
											colOut),
											collapse=","),
										")",
										"~parent_key")
									), 
                  data= rosterIn[,c("parent_key",
                                     anyOut,
                                     indOut,
                                     colOut)], 
                  sum)
	HH_add <- matrix(0, nrow=nrow(HH),ncol=3)
	colnames(HH_add) <- c(anyOut, indOut, colOut)
	HH_add <- as.data.frame(HH_add)
	HH_add[,anyOut][match(rosterIn_agg$parent_key, 
                        HH$key)] <- rosterIn_agg[,anyOut]
	HH_add[,indOut][match(rosterIn_agg$parent_key, 
                        HH$key)] <- rosterIn_agg[,indOut]
	HH_add[,colOut][match(rosterIn_agg$parent_key, 
                        HH$key)] <- rosterIn_agg[,colOut]
	return(HH_add)
}

# For finding variables that start with a certain pattern

findNames <- function(x, dataIn){
	names(dataIn)[apply(as.matrix(names(dataIn)), 
			1, 
			function(a){substr(a, 1, nchar(x))})  == x]
}


# TE graphs

ci_plot <- function(b_up=NULL,
                    ci_low=NULL,
                    ci_high=NULL,
                    yArg=NULL,
                    lwdArg=2,
                    alphaArg=.25){
  segments(ci_low,
           yArg,
           ci_high,
           yArg,
           lwd=lwdArg*.8,
           col=gray(0, alpha=alphaArg))
  segments(b_up-0.84*abs(b_up-ci_low),
           yArg,
           b_up+0.84*abs(b_up-ci_high),
           yArg,
           lwd=lwdArg,
           col=gray(0, alpha=alphaArg))
}



te_graph <- function(in1=NULL, in2=NULL){
  c_nl <- c(in1[[1]]$controlmean,
               in2[[1]]$controlmean)
  t_nl <- c(in1[[1]]$controlmean+in1[[1]]$estimate,
            in2[[1]]$controlmean+in2[[1]]$estimate)
  i_nl <- cbind(c(in1[[1]]$conf.low+in1[[1]]$controlmean,
                  in1[[1]]$conf.high+in1[[1]]$controlmean),
                c(in2[[1]]$conf.low+in2[[1]]$controlmean,
                  in2[[1]]$conf.high+in2[[1]]$controlmean))
  c_le <- c(in1[[2]]$controlmean,
            in2[[2]]$controlmean)
  t_le <- c(in1[[2]]$controlmean+in1[[2]]$estimate,
            in2[[2]]$controlmean+in2[[2]]$estimate)
  i_le <- cbind(c(in1[[2]]$conf.low+in1[[2]]$controlmean,
                  in1[[2]]$conf.high+in1[[2]]$controlmean),
                c(in2[[2]]$conf.low+in2[[2]]$controlmean,
                  in2[[2]]$conf.high+in2[[2]]$controlmean))
  
  tVec <- c(1,2)
  yrange <- range(c(c_nl,c_le,t_nl,t_le,i_nl,i_le))
  xlim_arg <- c(.75,2.25)
  c_shift <- -.05
  cex_arg <- 1.5
  par(mfrow=c(1,2))
  plot(tVec+c_shift,
       c_nl,
       ylim=yrange,
       xlim= xlim_arg,
       type="b",
       main="Non-leaders",
       ylab=in1[[1]]$label,
       axes=F,
       xlab="Outcome wave",
       cex=cex_arg)
  points(tVec,
         t_nl,
         ylim=yrange,
         type="b",
         pch=19,
       	 cex=cex_arg)
  segments(tVec,
           i_nl[1,],
           tVec,
           i_nl[2,],
           lwd=1.5,
           col=gray(0, .25))
  segments(tVec,
           t_nl-.84*abs(t_nl-i_nl[1,]),
           tVec,
           t_nl+.84*abs(t_nl-i_nl[2,]),
           lwd=2.5,
           col=gray(0, .35))

  
#    segments(tVec, i_nl[1,],tVec, i_nl[2,])
  axis(1, c(1,2))
  axis(2, seq(from=round(yrange[1], 1),
              to=round(yrange[2], 1),
              length=5),
       las=1)
  box()
  plot(tVec+c_shift, 
  		c_le,
         ylim=yrange,
       xlim= xlim_arg,
         type="b",
          main="Leaders",
          ylab="",
       axes=F,
       xlab="Outcome wave",
       cex=cex_arg)
  points(tVec,
         t_le,
         ylim=yrange,
         type="b",
         pch=19,
       cex=cex_arg)

  segments(tVec,
           i_le[1,],
           tVec,
           i_le[2,],
           lwd=1.5,
           col=gray(0, .25))
  segments(tVec,
           t_le-.84*abs(t_le-i_le[1,]),
           tVec,
           t_le+.84*abs(t_le-i_le[2,]),
           lwd=2.5,
           col=gray(0, .35))
  

  axis(1, c(1,2))
  axis(2, seq(from=round(yrange[1], 1),
              to=round(yrange[2], 1),
              length=5),
       las=1)
  box()
}


te_est_graph <- function(fileName=NULL,
                         fileFolder="tabs-figs/endline/te/",
                         heightArg=3,
                         widthArg=6,
                         yUp_1 = NULL,
                         yUp_2=yUp_1,
                         yLab="",
                         data1=end1,
                         data2=end2){
  
  te_1 <- fitR_2(yUp_1,
                 yLab,
                 HHdata=data1)
  te_2 <- fitR_2(yUp_2,
                 yLab,
                 HHdata=data2)
  
  pdf(file=paste0(fileFolder, fileName, ".pdf"),
      height=heightArg,
      width=widthArg)  
  te_graph(te_1, te_2)
  dev.off()
}

