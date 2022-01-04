GbyE.file.Calculate <- function(GD,GM,file.output){
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GbyE.GM的处理
	GM.L1=GM
	GM.L2=GM
	GM.L2[,2]=GM.L2[,2]+max(GM.L2[,2])
	GM.L1[,1]=paste(GM.L1[,1],"-1",sep="")
	GM.L2[,1]=paste(GM.L2[,1],"-2",sep="")
	GbyE.GM=rbind(GM.L1,GM.L2)
	rm(GM.L1)
	rm(GM.L2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GbyE.GD的处理
	GD.d=apply(GD[-1,-1],2,as.numeric)
	Matrix.r=matrix(1,2,2)
	Matrix.r[1,2]=0
	GbyE.GD=kronecker(Matrix.r,as.matrix(GD.d))
	rm(GD.d)
	rm(Matrix.r)
	Taxa="taxa"
	row.name=rbind(as.matrix(Taxa),as.matrix(paste(GD[-1,1],"-1",sep="")),as.matrix(paste(GD[-1,1],"-2",sep="")))
	col.name=t(as.matrix(GbyE.GM[,1]))
	GbyE.GD=rbind(col.name,GbyE.GD)
	GbyE.GD=cbind(row.name,GbyE.GD)
	rm(col.name)
	rm(row.name)
	if(file.output==TRUE){
						fwrite(GbyE.GD,paste(O_Path,"/","GbyE.GD.txt",sep=""),row.names=F,sep="\t",quote=FALSE,col.names=F)
						fwrite(GbyE.GM,paste(O_Path,"/","GbyE.GM.txt",sep=""),row.names=F,sep="\t",quote=FALSE,col.names=T)
						}else{
								print("~~~~~~~~~~There GbyE file not output~~~~~~~~~~")
								print("~~~~~~~~~If you want the file output~~~~~~~~~~")
								print("~~~~~~~~~~Please the file.output=1~~~~~~~~~~~~")
								}
						data=list(GbyE.GD,GbyE.GM)
						names(data)=c("GbyE.GD","GbyE.GM")
						return(data)
						}