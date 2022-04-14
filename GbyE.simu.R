`GbyE.simu` <- 
function(GD,h2,covg=1,ysdmean=0,qtn=FALSE,nqtn=NULL){
	D=GD[,2:ncol(GD)]
	D1=as.data.frame(as.matrix(D))
	D2=t(D1)
	#QTL <- 100*(1:20) #pick 20 QTL
	#构建真实QTL
	if(qtn==FALSE){
		if(is.null(nqtn)) stop("Please enter the number of QTN")
		QTL <- sample(1:ncol(D),nqtn,replace = FALSE)
		}else{
			QTL=qtn
		}
	u <-rep(0,ncol(D)) #marker effects
	u[QTL] <- covg
	#求SNP效应
	g <-as.vector(crossprod(D2,u))
	h2 <- h2      #heritability 
	y <- g +rnorm(nrow(D),mean=ysdmean,sd=sqrt((1-h2)/h2*var(g)))
	y = data.frame(GD[,1],y)
	colnames(y)=c("Taxa","trait")
	#Saving simulated file
	return(list(Y=y,QTL.position=QTL))
	#write.table(y,"H60Q20.txt", sep="\t")
}

