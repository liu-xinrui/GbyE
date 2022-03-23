GbyE <- function(
	GD=NULL,
	GM=NULL,
	GbyE.GD=NULL,
	GbyE.GM=NULL,
	Y=NULL,
	ha2=NULL,
	NQTN,
	NE=NULL,
	nrep=NULL,
	nfold=NULL,
	cov_g=NULL,
	cov_e=NULL,
	SimuY=NULL,
	gwas=NULL,
	gs=NULL,
	nIter=NULL,
	burnIn=NULL,
	plot=NULL,
	file.output=NULL){

if(is.null(GD)){stop("The value of GD is null")}
if(is.null(GM)){stop("The value of GM is null")}
if(is.null(NE)){NE=2}
if(is.null(nrep)){stop("nrep is NULL")}
if(is.null(nfold)){nfold=5
					print("nfold is NULL, default 5")}
if(is.null(SimuY)){SimuY=FALSE}
if(is.null(Y) && SimuY==TRUE){
if(is.null(ha2))
	{stop("The value of ha2 is null")}
if(is.null(NQTN))
	{stop("The value of NQTN is null")}
if(is.null(cov_g))
	{stop("The value of cov_g is null")}
if(is.null(cov_e)){cov_e=1
	print("cov_e is NULL, default 1")}
		}else{
			if(is.null(Y)){stop("The value of Y is null")}
			if(is.null(ha2)){ha2="real"}
			if(is.null(NQTN)){ha2="real"}
			if(is.null(cov_g)){ha2="real"}
			if(is.null(cov_e)){cov_e="real"}
		}
if(is.null(gwas))
	{gwas=TRUE
		print("gwas test is NULL, default TRUE")}
if(is.null(gs))
	{gs=TRUE
		print("gs test is NULL, default TRUE")}
if(is.null(nIter)){nIter=12000}
if(is.null(burnIn)){burnIn=2000}
if(is.null(file.output)){file.output=FALSE}
if(is.null(plot))
	{plot=FALSE
		print("plot default FALSE")}

#样本数据
n=nrow(GD)
nq=c(1:n) #Mean样本总体数量
nq2=c(1:(n*2)) #GbyE样本总体数量

#计算GbyE文件
if((is.null(GbyE.GM) && is.null(GbyE.GD=NULL))
{
print("Calculating the GbyE file")
mydata=GbyE.file.Calculate(GM=GM,GD=GD,file.output=file.output)
GbyE.GD=mydata$GbyE.GD
GbyE.GM=mydata$GbyE.GM
print("The GbyE file Calculate done")
}
#绘图参数
if(plot==TRUE){
pdf("Power FDR and Heritability.pdf",width=6*length(cov_g),height=6*length(ha2))
par(mfrow=c(length(ha2),length(cov_g)))
}
#遗传力遗传相关循环
for(g in cov_g){
	for(h in ha2){
		for(NNQTN in NQTN){
		kh=h
		print(paste(kh,"Heritability and",g,"Genetic Correlation",sep=" "))
		h=c(h,h)
	#模拟函数的设计矩阵
	rg=matrix(g,NE,NE)
	    diag(rg)=1
	re=matrix(cov_e,NE,NE)
	    diag(re)=1

mean.rep=NULL
gbye.rep=NULL
gexd.rep=NULL

StatRep=replicate(nrep,{
	krep=1
	#报告循环到第几次了
	print(paste("The GbyE program is running its",krep,"cycle",sep=" "))
	krep=krep+1

	if(SimuY==TRUE){
	#模拟表型
	mysimu=G_E.Simulation(GD=GD,ha2=h,NQTN=NNQTN,NE=NE,rg=rg,re=re,file.output=file.output)
	#GbyE.Y=mysimu$GbyE.Y #gbyey
	#Mean.Y=mysimu$Mean.Y
	myY=mysimu$Y
	QTN.Position=mysimu$QTN.Position
	#Mean.QTN.Position=mysimu$QTN.Position #QTN position
	#GbyE.QTN.Position=mysimu$GbyE.QTN.Position
	}else{
		myY=data.frame(Y,Y[,2])
		}
	Mean.Y=myY[,-3]
	Mean.Y[,2]=(myY[,2]+myY[,3])/2
	colnames(Mean.Y)=c("taxa","trait")
	myY.1=cbind(paste(myY[,1],"-1",sep=""),myY[,2])
	myY.2=cbind(paste(myY[,1],"-2",sep=""),myY[,3])
	GbyE.Y=rbind(myY.1,myY.2)
	colnames(GbyE.Y)=c("taxa","trait")
	GbyE.Y=as.data.frame(as.matrix(data.table(GbyE.Y)))
	GbyE.Y[,2]=as.numeric(GbyE.Y[,2])


	#A Part of GWAS	
	if(gwas==TRUE){
			print("###########################  Running GWAS  ###########################")
			#Mean				
			myGAPIT<-GAPIT(
					Y=Mean.Y, #fist column is individual ID, the third columns is days to pollination
					GD=GD,
					GM=GM,
					PCA.total=3,
					#QTN.position=QTN.position,
					model="MLM",
					file.output=FALSE
			)

			print("Mean GWAS Power VS FDR Calculate")
						
			myGWAS=myGAPIT$GWAS[,1:4]
			myGWAS=myGWAS[order(myGWAS[,4]),]
			cc=as.vector(GM[QTN.Position,][,1])
			Mean.data=Power.FDR.Calculate(cc=cc,myGWAS=myGWAS,myGM=GM,NQTN=NNQTN)

			#GbyE
			myGAPIT<-GAPIT(
					Y=GbyE.Y, #fist column is individual ID, the third columns is days to pollination
					GD=GbyE.GD,
					GM=GbyE.GM,
					PCA.total=3,
					#QTN.position=QTN.position,
					model="MLM",
					file.output=FALSE
			)

			print("###########################  GbyE GWAS Power VS FDR Calculate  ###########################")
					
			#GbyE
			#GbyE.GWAS=Comparison.GWAS.Reslut.PValue(GD=GbyE.GD,GM=GbyE.GM,GWASresult=myGAPIT)
			GbyE.GWAS=Comparison.GWAS.Reslut.PValue(GM=GbyE.GM,GWASresult=myGAPIT)
			rownames(GbyE.GWAS)=c(1:(ncol(GD)-1))
			cc=as.vector(GbyE.GWAS[QTN.Position,][,1])
			GbyE.GWAS=GbyE.GWAS[order(GbyE.GWAS[,4]),]
			GbyE.data=Power.FDR.Calculate(cc=cc,myGWAS=GbyE.GWAS,myGM=GbyE.GM[1:(nrow(GbyE.GM)/2),],NQTN=NNQTN)

			print("###########################  AddE GWAS Power VS FDR Calculate  ###########################")
					
			#Additive Eﬀects
			AddE.GWAS=myGAPIT$GWAS[1:(nrow(GbyE.GM)/2),1:4]
			cc=as.vector(AddE.GWAS[QTN.Position,][,1])
			AddE.GWAS=AddE.GWAS[order(AddE.GWAS[,4]),]
			AddE.data=Power.FDR.Calculate(cc=cc,myGWAS=AddE.GWAS,myGM=GbyE.GM[1:(nrow(GbyE.GM)/2),],NQTN=NNQTN)
								
			print("###########################  GGEE GWAS Power VS FDR Calculate  ###########################")
					
			#GXE
			GGEE.GWAS=myGAPIT$GWAS[-(1:(nrow(GbyE.GM)/2)),1:4]
			rownames(GGEE.GWAS)=c(1:(ncol(GD)-1))
			cc=as.vector(GGEE.GWAS[QTN.Position,][,1])
			GGEE.GWAS=GGEE.GWAS[order(GGEE.GWAS[,4]),]
			GGEE.data=Power.FDR.Calculate(cc=cc,myGWAS=GGEE.GWAS,myGM=GbyE.GM[1:(nrow(GbyE.GM)/2),],NQTN=NNQTN)
		}else{
			print("###########################  Don't calculate GWAS  ###########################")
		}

		if(gs==TRUE){
			print("###########################  Running GS  ###########################")
		#A part of GS
		sets=sample(cut(1:n,nfold,labels=FALSE),n)
		#Set up database (load cyclic data)
		mean=NULL
		gbye=NULL
		gexd=NULL

		for(xfold in 1:nfold)
		{
			#Phenotypic backup
			Mean.Y.raw=Mean.Y[,c(1,2)] #choos a trait
			GbyE.Y.raw=GbyE.Y[,c(1,2)]
			#K折提取
			Mean.Y.X=sets==xfold
			GbyE.Y.X=c(Mean.Y.X,Mean.Y.X)
			Mean.Y.raw[Mean.Y.X,2]=NA
			GbyE.Y.raw[GbyE.Y.X,2]=NA
			Mean.Y.raw=Mean.Y.raw[!is.na(Mean.Y.raw[,2]),]
			GbyE.Y.raw=GbyE.Y.raw[!is.na(GbyE.Y.raw[,2]),]
			#other methods
			#Mean.Y.raw=Mean.Y.raw[sets!=xfold,]
			#GbyE.Y.raw=GbyE.Y.raw[c(sets!=xfold,sets!=xfold),]
			
			#保存提取的原始行
			raw.M=c(as.numeric(rownames(Mean.Y.raw)),nq)
			raw.M=raw.M[duplicated(raw.M)]
			raw.M <- nq[!(nq %in% raw.M)]
			raw.G=c(as.numeric(rownames(GbyE.Y.raw)),nq2)
			raw.G=raw.G[duplicated(raw.G)]
			raw.G <- nq2[!(nq2 %in% raw.G)]

			#Mean
			Mean.method=mixed.solve(y=Mean.Y.raw[,2],Z=GD[-raw.M,-1],method="ML")
     		mean.cor=as.matrix(GD[raw.M,-1]) %*% Mean.method$u
      		mean=append(mean,cor(Mean.Y[raw.M,2],mean.cor))

      		#GbyE
      		GbyE.method=mixed.solve(y=GbyE.Y.raw[,2],Z=GbyE.GD[-raw.G,-1],method="ML")
     		gbye.cor=as.matrix(GbyE.GD[raw.G,-1]) %*% GbyE.method$u
      		gbye=append(gbye,cor(GbyE.Y[raw.G,2],gbye.cor))

      		#GExd
      		GExd.GD<-list(list(X=GbyE.GD[,2:ncol(GD)],model='BRR'),list(X=GbyE.GD[,(ncol(GD)+1):(ncol(GD)*2-1)],model='BRR')) #list两部分
      		#Additive effects and interaction effects are calculated separately
      		brrY=GbyE.Y 
      		brrY[raw.G,2]=NA #GExd的表型缺失，BGLR包要求为缺失对应
      		GExd.method<-BGLR(y=as.numeric(brrY[,2]),ETA=GExd.GD, nIter=nIter, burnIn=burnIn)
      		gexd=append(gexd,cor(GbyE.Y[raw.G,2],GExd.method$yHat[raw.G]))
		}

		#the Mean ValueGbyE
		mean.rep=append(mean.rep,mean(mean)) #mean 
		gbye.rep=append(gbye.rep,mean(gbye)) #GbyE
		gexd.rep=append(gexd.rep,mean(gexd)) 
		}else{
			print("###########################  Don't calculate GS  ###########################")
		 }
		#if(gwas==TRUE&gs=TRUE){myStat=list(Mean.data,GbyE.data,AddE.data,GGEE.data,mean.rep,gbye.rep,gexd.rep)}
		if(gwas==FALSE && gs==TRUE){Mean.data=NULL
									GbyE.data=NULL
									AddE.data=NULL
									GGEE.data=NULL}
		if(gwas==TRUE && gs==FALSE){mean.rep=NULL
									gbye.rep=NULL
									gexd.rep=NULL}

		myStat=list(Mean.data,GbyE.data,AddE.data,GGEE.data,mean.rep,gbye.rep,gexd.rep)
		}) #This is the end of the loop

#Computed loop value

#A part of GWAS
if(gwas==TRUE){
	Mean_power_FDR=StatRep[1,]
	GbyE_power_FDR=StatRep[2,]
	AddE_power_FDR=StatRep[3,]
	GGEE_power_FDR=StatRep[4,]
	Mean_power=NULL
	Mean_FDR=NULL
	GbyE_power=NULL
	GbyE_FDR=NULL
	AddE_power=NULL
	AddE_FDR=NULL
	GGEE_power=NULL
	GGEE_FDR=NULL

	for (i in 1:nrep)
		{
			#Meandata
			Mean_power_A=as.matrix(as.data.frame(Mean_power_FDR[i])$Power)
			Mean_FDR_A=as.matrix(as.data.frame(Mean_power_FDR[i])$FDR)
			#GbyEdata
			GbyE_power_A=as.matrix(as.data.frame(GbyE_power_FDR[i])$Power)
			GbyE_FDR_A=as.matrix(as.data.frame(GbyE_power_FDR[i])$FDR)
			#AddEdata
			AddE_power_A=as.matrix(as.data.frame(AddE_power_FDR[i])$Power)
			AddE_FDR_A=as.matrix(as.data.frame(AddE_power_FDR[i])$FDR)
			#GGEEdata
			GGEE_power_A=as.matrix(as.data.frame(GGEE_power_FDR[i])$Power)
			GGEE_FDR_A=as.matrix(as.data.frame(GGEE_power_FDR[i])$FDR)
			
			if(i==1)
			{
				#Meandata
				Mean_power=Mean_power_A
				Mean_FDR=Mean_FDR_A
				#GbyEdata
				GbyE_power=GbyE_power_A
				GbyE_FDR=GbyE_FDR_A
				#AddEdata
				AddE_power=AddE_power_A
				AddE_FDR=AddE_FDR_A
				#GGEEdata
				GGEE_power=GGEE_power_A
				GGEE_FDR=GGEE_FDR_A
				
				}else{
						#Meandata
						Mean_power=cbind(Mean_power,Mean_power_A)
						Mean_FDR=cbind(Mean_FDR,Mean_FDR_A)
						#GbyEdata
						GbyE_power=cbind(GbyE_power,GbyE_power_A)
						GbyE_FDR=cbind(GbyE_FDR,GbyE_FDR_A)
						#AddEdata
						AddE_power=cbind(AddE_power,AddE_power_A)
						AddE_FDR=cbind(AddE_FDR,AddE_FDR_A)
						#GGEEdata
						GGEE_power=cbind(GGEE_power,GGEE_power_A)
						GGEE_FDR=cbind(GGEE_FDR,GGEE_FDR_A)
					} #end for i=1
		} #end for i in nrep
		rm(Mean_power_A)
		rm(GbyE_power_A)
		rm(AddE_power_A)
		rm(GGEE_power_A)
		rm(Mean_FDR_A)
		rm(GbyE_FDR_A)
		rm(AddE_FDR_A)
		rm(GGEE_FDR_A)
		
		#GbyE
		GbyE_power=apply(GbyE_power,2,as.numeric)
		GbyE_power=as.data.frame(rowMeans(GbyE_power))
		colnames(GbyE_power)="GbyE_power"
		
		GbyE_FDR=apply(GbyE_FDR,2,as.numeric)
		GbyE_FDR=as.data.frame(rowMeans(GbyE_FDR))
		colnames(GbyE_FDR)="GbyE_FDR"
		
		#Mean
		Mean_power=apply(Mean_power,2,as.numeric)
		Mean_power=as.data.frame(rowMeans(Mean_power))
		colnames(Mean_power)="Mean_power"
		
		Mean_FDR=apply(Mean_FDR,2,as.numeric)
		Mean_FDR=as.data.frame(rowMeans(Mean_FDR))
		colnames(Mean_FDR)="Mean_FDR"
		
		#AddE
		AddE_power=apply(AddE_power,2,as.numeric)
		AddE_power=as.data.frame(rowMeans(AddE_power))
		colnames(AddE_power)="AddE_power"
		
		AddE_FDR=apply(AddE_FDR,2,as.numeric)
		AddE_FDR=as.data.frame(rowMeans(AddE_FDR))
		colnames(AddE_FDR)="AddE_FDR"
		
		#GGEE
		GGEE_power=apply(GGEE_power,2,as.numeric)
		GGEE_power=as.data.frame(rowMeans(GGEE_power))
		colnames(GGEE_power)="GGEE_power"
		
		GGEE_FDR=apply(GGEE_FDR,2,as.numeric)
		GGEE_FDR=as.data.frame(rowMeans(GGEE_FDR))
		colnames(GGEE_FDR)="GGEE_FDR"
		
		
	Power.VS.FDR=cbind(GbyE_power,GbyE_FDR,Mean_power,Mean_FDR,AddE_power,AddE_FDR,GGEE_power,GGEE_FDR)
	#write.table(Power.VS.FDR,file="Power.VS.FDR.txt",col.names=F,row.names=T,quote=F)
	
	#Power.VS.FDR.Plot.Read=as.data.frame(rbind(
	#											as.matrix(cbind(as.data.frame(matrix("GbyE",nrow(Power.VS.FDR[1:2]),1)),Power.VS.FDR[1:2])),
	#											as.matrix(cbind(as.data.frame(matrix("Mean",nrow(Power.VS.FDR[3:4]),1)),Power.VS.FDR[3:4])),
	#											as.matrix(cbind(as.data.frame(matrix("AddE",nrow(Power.VS.FDR[5:6]),1)),Power.VS.FDR[5:6])),
	#											as.matrix(cbind(as.data.frame(matrix("GGEE",nrow(Power.VS.FDR[7:8]),1)),Power.VS.FDR[7:8]))
	#											))
	#colnames(Power.VS.FDR.Plot.Read)=c("Method","Power","FDR")

	#GbyE.data=cbind(as.data.frame(matrix("GbyE",nrow(Power.VS.FDR[1:2]),1)),Power.VS.FDR[1:2])
	#Mean.data=cbind(as.data.frame(matrix("Mean",nrow(Power.VS.FDR[3:4]),1)),Power.VS.FDR[3:4])
	#AddE.data=cbind(as.data.frame(matrix("AddE",nrow(Power.VS.FDR[5:6]),1)),Power.VS.FDR[5:6])
	#GGEE.data=cbind(as.data.frame(matrix("GGEE",nrow(Power.VS.FDR[7:8]),1)),Power.VS.FDR[7:8])
		rm(Mean_power_FDR)
		rm(GbyE_power_FDR)
		rm(AddE_power_FDR)
		rm(GGEE_power_FDR)
		rm(Mean_power)
		rm(Mean_FDR)
		rm(GbyE_power)
		rm(GbyE_FDR)
		rm(AddE_power)
		rm(AddE_FDR)
		rm(GGEE_power)
		rm(GGEE_FDR)
if(file.output==TRUE){
write.table(Power.VS.FDR,file=paste(kh,".Heritability.and.",g,".Genetic.Correlation.","NQTN=",NNQTN,".Power.VS.FDR.",nrep,".times.txt",sep=""),
			col.names=T,row.names=F,quote=F,sep="\t")
} #output file
	}#end of the gwas result wirte


if(gs==TRUE){
#A part of GS
mean.rep=as.numeric(StatRep[5,])
gbye.rep=as.numeric(StatRep[6,])
gexd.rep=as.numeric(StatRep[7,])

GS.results=cbind(as.data.frame(mean.rep),as.data.frame(gbye.rep),as.data.frame(gexd.rep))
colnames(GS.results)=c("Mean","GbyE","GExd")
if(file.output==TRUE){
write.table(GS.results,file=paste(kh,".Heritability.and.",g,".Genetic.Correlation.","NQTN=",NNQTN,".GS.results.",nrep,".times.txt"),
			quote = FALSE, sep = "\t",row.names=F)
} #output file
#write.csv(GS.results,file="GS.results.csv",quote = FALSE,row.names=F)
print(paste(kh," Heritability and ",g," Genetic Correlation ","NQTN=",NNQTN,sep=""))
print(paste("Mean:",sprintf("%0.5f",mean(mean.rep)),"GbyE:",
					sprintf("%0.5f",mean(gbye.rep)),"GExd:",
					sprintf("%0.5f",mean(gexd.rep)),sep=" "))

		} #end of the gs result wirte


#part of Power FDR plot
if(plot==TRUE & gwas==TRUE){
#pdf("Power FDR and Heritability.pdf",width=4*length(g),height=4*length(kh))
#par(mfrow=c(length(kh),length(g)))

	#pdf(file=paste(i,".Genetic.Correlation.pdf",sep=""))
	#Power.VS.FDR=read.table(paste(i,".Genetic.Correlation.Power.VS.FDR.txt",sep=""),head=T)
	#Power.VS.FDR=apply(Power.VS.FDR,2,as.numeric)
	plot(x=Power.VS.FDR[,2],y=Power.VS.FDR[,1],type="l",  #GbyE
						xlab="FDR",ylab="Power",
						xlim = c(0,1),ylim = c(0,1),
						col=rgb(253,10,7,255,maxColorValue=255),lwd=2,
						main=paste("ha2=",kh," cov_e=",g," NQTN=",NNQTN,"  Power VS FDR",sep=""))
	lines(x=Power.VS.FDR[,8],y=Power.VS.FDR[,7],type="l",col=rgb(141,172,208,255,maxColorValue=255),lwd=2)  #GGEE
	lines(x=Power.VS.FDR[,6],y=Power.VS.FDR[,5],type="l",col=rgb(11,11,11,255,maxColorValue=255),lwd=2)  #AddE
	lines(x=Power.VS.FDR[,4],y=Power.VS.FDR[,3],type="l",col=rgb(3,175,78,255,maxColorValue=255),lwd=2) #Mean
	legend("bottomright", inset=.05, c("GbyE","GXE","AddE","Mean"),lty=1,lwd=3,
						col=c(rgb(253,10,7,255,maxColorValue=255),
						rgb(141,172,208,255,maxColorValue=255),
						rgb(11,11,11,255,maxColorValue=255),
						rgb(3,175,78,255,maxColorValue=255)
						))
			} #end of NQTN
		} #end of plot 
	} #end of ha2
} #End the cyclic calculation of heritability and heredity (ha2 & cov_g)
	#end of cov_g
if(plot==TRUE){dev.off()}
if(gwas==FALSE){Power.VS.FDR=NULL}
if(gs==FALSE){GS.results=NULL}
file.remove(dir(".", pattern="(.dat)$")) #remove Cache data
#part of GS boxs plot
return(list(GWAS=Power.VS.FDR,GS=GS.results))
} #end of all

