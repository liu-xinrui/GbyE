GbyE.file.Calculate <- function(GD,GM,file.output){
#library(data.table)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GbyE.GM的处理
	if(is.numeric(GM[,2])==TRUE){
	print("The Chr is numeric")
	GM.L1=GM
	GM.L2=GM
	GM.L2[,2]=GM.L2[,2]+max(GM.L2[,2])
	GM.L1[,1]=paste(GM.L1[,1],"-1",sep="")
	GM.L2[,1]=paste(GM.L2[,1],"-2",sep="")
	GbyE.GM=rbind(GM.L1,GM.L2)
	#rm(GM.L1)
	#rm(GM.L2)
		}else{
		print("The Chr is not numeric")
		GM.L1=GM
		GM.L2=GM
		#GM.L2[,2]=GM.L2[,2]+10*(length(unlist(strsplit(as.character(max(GM[,3])),"")))+1)
		GM.L1[,1]=paste(GM.L1[,1],"-1",sep="")
		GM.L2[,1]=paste(GM.L2[,1],"-2",sep="")
		GM.L2[,2]=paste(GM.L2[,2],"_2",sep="")
		GbyE.GM=rbind(GM.L1,GM.L2)
		}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GbyE.GD的处理
	GD.d=apply(GD[,-1],2,as.numeric)
	Matrix.r=matrix(1,2,2)
	Matrix.r[1,2]=0
	GbyE.GD=kronecker(Matrix.r,as.matrix(GD.d))
	#rm(GD.d)
	#rm(Matrix.r)
	Taxa="taxa"
	#row.name=rbind(as.data.frame(c(as.character(Taxa),paste(GD[,1],"-1",sep=""),paste(GD[,1],"-2",sep=""))))
	row.name=rbind(as.data.frame(c(paste(GD[,1],"-1",sep=""),paste(GD[,1],"-2",sep=""))))
	#row.name=rbind(as.matrix(Taxa),as.matrix(paste(GD[,1],"-1",sep="")),as.matrix(paste(GD[,1],"-2",sep="")))
	#col.name=t(as.matrix(GbyE.GM[,1]))
	#GbyE.GD=rbind(col.name,apply(GbyE.GD,2,as.numeric))
	#GbyE.GD=data.frame(col.name,GbyE.GD)
	#col.name=t(as.matrix(as.character(GbyE.GM[,1])))
	#lm.GD=data.frame(matrix(NA,nrow(GbyE.GD)+1,ncol(GbyE.GD)))
	#lm.GD[1,]=col.name
	#lm.GD[-1,]=GbyE.GD
	#GbyE.GD=as.martix(lm.GD)
	GbyE.GD=cbind(row.name,GbyE.GD)
	colnames(GbyE.GD)=c(paste(as.character(Taxa)),as.character(GbyE.GM[,1]))
	#GbyE.GD=rbind(col.name,apply(GbyE.GD,2,as.numeric))
	#GbyE.GD=cbind(row.name,GbyE.GD)
	#GbyE.GD=data.frame(GbyE.GD)
	#col.names=as.character(GbyE.GD[1,])
	#GbyE.GD=GbyE.GD[-1,]
	#colnames(GbyE.GD)=col.names
	#rm(col.name)
	#rm(row.name)
	#转化字符串
	GbyE.GD[,1]=as.character(GbyE.GD[,1])
	#转化数值型
	#if(is.numeric(GbyE.GD[,-1])==FALSE){
	#GbyE.GD[,-1]=apply(GbyE.GD[,-1],1,as.numeric)
	#}
	if(file.output==TRUE){
						fwrite(GbyE.GD,paste(getwd(),"/","GbyE.GD.txt",sep=""),row.names=F,col.names=T,sep="\t",quote=FALSE)
						fwrite(GbyE.GM,paste(getwd(),"/","GbyE.GM.txt",sep=""),row.names=F,col.names=T,sep="\t",quote=FALSE)
						}else{
								print("~~~~~~~~~The GbyE file is not output~~~~~~~~~~")
								print("~~~~~~~~~If you want the file output~~~~~~~~~~")
								print("~~~~~~~~~Please set file.output=TRUE~~~~~~~~~~")
								}
						#data=list(GbyE.GD,GbyE.GM)
						#names(data)=c("GbyE.GD","GbyE.GM")
						return(list(GbyE.GD=GbyE.GD,GbyE.GM=GbyE.GM))
						}
