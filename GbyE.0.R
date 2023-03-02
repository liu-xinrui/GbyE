`GbyE` <- function(
		GD=NULL, #Genotype data 0 1 2
		GM=NULL, #SNP location information
		GbyE.GD=NULL, #If you calculate the GbyE file, then you can add between
		GbyE.GM=NULL,
		Y=NULL,  #The phenotype file must contain two environments or two traits
		gwas=FALSE,
		PCA.total=3, #PCA default 3
		CV=NULL, #gwas 协变量
		gs=FALSE,
		gwas.model=NULL, #GWAS MLM-CMLM-BLINK etc; GS
		gs.model=NULL, #GS model: rrblup-default-"ML", gapit "gblup" "mablup" "cblup" "sblup", bglr "RRB" "BA" "BC" "BB" "BL"
		method="gblup",  #GS method: "gapit", "rrblup", "bglr" can select
		cutOff=0.05,  #manhattan cutoff
		nIter=12000,  #BGLR
		burnIn=10000, #BGLR-
		plot=FALSE, #GWAS plot manhattan
		file.output=FALSE, 
		file.type="pdf", #输出图片的格式
		dpi=600){
#=============================================================================================
 #Object: Construction of interactive genotype matrices using multi-trait or multi-environment information,
 #			complete with improved statistical power of GWAS and GS
 #Output: Computational genotype and phenotype files for GbyE, GWAS and GS results using GbyE analysis
 #intput: Input genotype and phenotype files
 #Authors: Xinrui Liu, Jiabo Wang
 #Designed by Jiabo Wang
 #Last update: March 2, 2023
##############################################################################################

	#Data Judgment
	if(is.null(GD)){stop("The value of GD is null")}
	if(is.null(GM)){stop("The value of GM is null")}
	if(gwas==FALSE & gs==FALSE){print("No arithmetic analysis, only GbyE files are calculated")}
	if(is.null(Y)){stop("The value of Y is null")}

	#Conditional Judgment
	if(!require(data.table)) install.packages("data.table")
	if(!require(rrBLUP)) install.packages("rrBLUP")
	if(!require(BGLR)) install.packages("BGLR")
	source("https://raw.githubusercontent.com/liu-xinrui/GbyE/main/GbyE.function.R")
	#Handling data types
	GD[,1]=as.factor(GD[,1])
	Y[,1]=as.factor(Y[,1])
	#Calculating GbyE files
	if(is.null(GbyE.GM) & is.null(GbyE.GD)){
		print("Calculating the GbyE file")
		mydata=GbyE.file.Calculate(GM=GM,GD=GD,file.output=file.output)
		GbyE.GD=mydata$GbyE.GD
		GbyE.GM=mydata$GbyE.GM
		print("The GbyE file Calculate done")
	}else{
		print("The GbyE file has been entered")
		GbyE.GD[,1]=as.factor(GbyE.GD[,1])
	}
	myY=Y[,1:3]
	print("Only the first two phenotypes of the phenotype data are read")
	print("The first is assumed to be additive effect and the second is assumed to be interactive effect")
	myY.1=cbind(paste(myY[,1],"-1",sep=""),myY[,2])
	myY.2=cbind(paste(myY[,1],"-2",sep=""),myY[,3])
	GbyE.Y=data.frame(rbind(myY.1,myY.2))
	colnames(GbyE.Y)=c("taxa","trait")
	GbyE.Y[,2]=as.numeric(as.character(GbyE.Y[,2]))
	AddE=NULL
	GXE=NULL
	gbye.gwas=NULL
	GMM=NULL
	if(gwas==TRUE){
		#running GWAS
		if(is.null(gwas.model)==TRUE){
			gwas.model="MLM"
			print("GbyE selects 'MLM' model for GWAS analysis by default")
		}
		gbye.gwas<-GAPIT(
						Y=GbyE.Y, #fist column is individual ID, the third columns is ..
						GD=GbyE.GD,
						GM=GbyE.GM,
						PCA.total=PCA.total,
						CV=CV,
						model=gwas.model,
						file.output=F
				)
		GbyE=gbye.gwas$GWAS
		raw=unlist(strsplit(GbyE[,1],"-"))
		GbyE[,1]=raw[seq(1,length(raw),2)]
		GbyE=data.frame(GbyE,raw[seq(2,length(raw),2)])
		colnames(GbyE)[10]="part"
		AddE=GbyE[which(GbyE$part==1),-11]
		GXE=GbyE[which(GbyE$part==2),-11]
		GXE[,2]=GXE[,2]-10
		if(file.output==TRUE){
			write.csv(AddE,"GAPIT.Association.GWAS_Results.AddE.Traits.csv",row.names=F)
			write.csv(GXE,"GAPIT.Association.GWAS_Results.GbyE.Traits.csv",row.names=F)
		} #end of outpout gwas
		if(plot==TRUE){
			#Ploting Manhattan
			write.csv(AddE,"GAPIT.Association.GWAS_Results.AddE.Traits.csv",row.names=F)
			write.csv(GXE,"GAPIT.Association.GWAS_Results.GbyE.Traits.csv",row.names=F)
			GMM=GAPIT.Multiple.Manhattan(model_store=c("AddE","GbyE"),
				Y.names="Traits",
				cutOff=cutOff,
				GM=GM)
			GAPIT.Circle.Manhattan.Plot(GMM$multip_mapP,
			            plot.type=c("m","q"),
			            signal.line=1,
			            xz=GMM$xz,
			            threshold=0.01,
			            file=file.type,
			            ylim=c(0,15),
			            memo="",
			            dpi=dpi)
		} #end of plot manhattan
	} #end of GWAS

	pred.y=NULL
	gbye.rrbup=NULL
	gbye.gapit=NULL
	gbye.bglr=NULL
	if(gs==TRUE){
		#GS
		print("Run genome-wide association selection analysis")
		if(method=="rrblup"){
			#rrblup
			print("GS method selects 'rrBLUP'")
			if(is.null(gs.model)==TRUE){
				gs.model="ML"
				print("GbyE selects 'ML' model for rrblup analysis by default")
			}
			gbye.rrbup=mixed.solve(
				y=GbyE.Y[!is.na(GbyE.Y[,2]),2],
				Z=GbyE.GD[!is.na(GbyE.Y[,2]),-1],
				method=gs.model)
			pred=gbye.rrbup$beta+as.vector(as.matrix(GbyE.GD[,-1]) %*% as.numeric(gbye.rrbup$u))
			#pred.y=data.frame(cbind(as.character(Y[,1]),pred[1:nrow(GD)],pred[-(1:nrow(GD))]))
			#colnames(pred.y)=c("taxa",colnames(Y)[2:3])
			pred.y=data.frame(cbind(as.character(GbyE.Y[,1]),rep(1,1,nrow(GbyE.GD)),pred))
			pred.y[is.na(GbyE.Y[,2]),2]=2
			colnames(pred.y)=c("Taxa","RefInf","Prediction")
		}

		if(method=="gapit"){
			#GAPIT
			print("GS method selects 'GAPIT'")
			#include gblup mablup
			if(is.null(gs.model)==TRUE){
				gs.model="gblup"
				print("GbyE selects 'gblup' model for gapit-gs analysis by default")
			}

			gbye.gapit <-GAPIT(
								Y=GbyE.Y, #fist column is individual ID, the third columns is ..
								GD=GbyE.GD,
								GM=GbyE.GM,
								PCA.total=PCA.total,
								CV=CV,
								model=gs.model,
								file.output=F,
								SNP.test=F
						)
			gapit=gbye.gapit$Pred[gbye.gapit$Pred[,3]==2,]
			pred.y=gapit[,c(1,3,8)]
			#pred=gbye.gapit$Pred[]
			#inf.y=merge(Y[!training,],inf.y,by.x=colnames(Y)[1],by.y=colnames(inf.y)[1])
			#pred.Y=cor(Y[!training,2],inf.y[,9])
			#pred.Y.rep=append(pred.Y.rep,pred.Y)
		}

		if(method=="bglr"){
			#BGLR
			print("GS method selects 'BGLR'")
			if(is.null(gs.model)==TRUE){
				gs.model="BRR"
				print("GbyE selects 'BRR' model for BGLR analysis by default")
			}
			#Additive genetic and reciprocal effects are predicted in two parts in BGLR prediction
			bglr.gd=GbyE.GD[,-1]
			gbye.bglr <- BGLR(
				y=GbyE.Y[!is.na(GbyE.Y[,2]),2],
				ETA=list(
					list(X=bglr.gd[!is.na(GbyE.Y[,2]),1:nrow(GD)],model=gs.model), #Additive effect
					list(X=bglr.gd[!is.na(GbyE.Y[,2]),-(1:nrow(GD))],model=gs.model)), #Interaction effect
				nIter=nIter,
				burnIn=burnIn)

			#Calculate predicted phenotype values
			pred=gbye.bglr$mu+as.vector(as.matrix(bglr.gd[is.na(GbyE.Y[,2]),])%*%as.numeric(c(gbye.bglr$ETA[[1]]$b,gbye.bglr$ETA[[2]]$b))) #gbye
			bglr.env1=gbye.bglr$mu+as.vector(as.matrix(bglr.gd[is.na(GbyE.Y[,2]),1:nrow(GD)])%*%as.numeric(gbye.bglr$ETA[[1]]$b)) #adde
			bglr.env2=gbye.bglr$mu+as.vector(as.matrix(bglr.gd[is.na(GbyE.Y[,2]),-(1:nrow(GD))])%*%as.numeric(gbye.bglr$ETA[[2]]$b)) #inter
			#all prediciton Y traits value
			bglr.pred.y=gbye.bglr$mu+as.vector(as.matrix(bglr.gd)%*%as.numeric(c(gbye.bglr$ETA[[1]]$b,gbye.bglr$ETA[[2]]$b))) #gbye
			pred.y=data.frame(cbind(as.character(GbyE.Y[,1]),rep(1,1,nrow(GbyE.GD)),bglr.pred.y))
			pred.y[is.na(GbyE.Y[,2]),2]=2
			colnames(pred.y)=c("Taxa","RefInf","Prediction")
		}

	} #end of GS
	if(file.output==TRUE){write.csv(pred.y,file="pred.y.csv",row.names=F)}

file.remove(dir(".", pattern="(.dat)$")) #remove Cache data
return(list(GbyE=mydata,GbyE.Y=GbyE.Y,
			GWAS=list(adde=AddE,GbyE=GXE,gwas.summary=gbye.gwas,GMM=GMM),
			GS=list(pred=pred.y,rrblup=gbye.rrbup,gapit=gbye.gapit,bglr=gbye.bglr)))
}