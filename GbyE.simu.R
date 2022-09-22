`GbyE.simu` <- 
function(GD,h2,covg=1,ysdmean=0,qtn=NULL,nqtn=NULL,qtn.cor=1){
	D=GD[,2:ncol(GD)]
	D1=as.data.frame(as.matrix(D))
	D2=t(D1)
	#构建真实QTL
	if(is.null(qtn)==TRUE){
		if(is.null(nqtn)) stop("Please enter the number of QTN")
		QTL <- sample(1:ncol(D),nqtn,replace = FALSE) #pick up QTL
		}else{
			if(qtn.cor==1){
				QTL=qtn
				}else{
					if(qtn.cor==0){
						qtn=sample(setdiff(1:ncol(GD),qtn),length(qtn),replace = FALSE)
					}else{
						setpoint=round(qtn.cor*length(qtn))
						cor=sample(qtn,setpoint,replace = FALSE)
						sqtn=sample(setdiff(1:ncol(GD),cor),(length(qtn)-setpoint),replace = FALSE)
						QTL=c(cor,sqtn)
					}
				}
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
}

