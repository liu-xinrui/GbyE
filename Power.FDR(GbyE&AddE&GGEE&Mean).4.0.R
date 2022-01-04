rm(list=ls())

setwd("C:/Users/VEN MIKIR/Desktop/GbyE.test")
setwd("/home/student/workshop/lxr/GbyE")
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("G_E.Duplicated.R")
source("G&E_Simulation.R")
#source("F:/swu-青研院/实验室/code/Power.FDR.Calculate.R")
#source("GWASbyCor.R")
library(data.table)
library(ggplot2)
library('MASS')
O_Path=getwd()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Power.FDR.Calculate
Power.FDR.Calculate <- function(cc,myGWAS,NQTN,myGM){
result_transcripts=myGWAS
data=NULL
				#for (q in c(1:2)){
				power=NULL
				fdr=NULL
				s=0
				S=0	
				m=0
				n=0
				q=nrow(myGM)
				j=NQTN
				
					for (r in 1:q){
						if(result_transcripts[r,1] %in% cc[1:j]){
						   power[[r]]=(m+1)/j
						   m=m+1
						}else{
						   if(r == 1){		       
						   power[[r]]=n/j
						   }else{
						   power[[r]]=power[[r-1]]		   
						   }
						power=c(power,power[[r]])
						}
						if(result_transcripts[r,1] %in% cc[1:j]){
						   if(r == 1){		        
						   fdr[[r]]=s/(q-j)
						   }else{
						   fdr[[r]]=fdr[[r-1]]	
						   }
						}else{
						   fdr[[r]]=(S+1)/(q-j)
						   S=S+1
						}
						fdr=c(fdr,fdr[[r]])     
					}
				data=data.frame(gene_id=result_transcripts$SNP,Power=as.matrix(power[-(q+1)]),FDR=as.matrix(fdr[-(q+1)]))
				
				return(data)
				}
				

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
	if(file.output==1){
						fwrite(GbyE.GD,paste(O_Path,"/","GbyE.GD.txt",sep=""),row.names=F,sep=" ",quote=FALSE,col.names=F)
						fwrite(GbyE.GM,paste(O_Path,"/","GbyE.GM.txt",sep=""),row.names=F,sep=" ",quote=FALSE,col.names=T)
						}else{
								print("~~~~~~~~~~There GbyE file not output~~~~~~~~~~")
								print("~~~~~~~~~If you want the file output~~~~~~~~~~")
								print("~~~~~~~~~~Please the file.output=1~~~~~~~~~~~~")
								}
						data=list(GbyE.GD,GbyE.GM)
						names(data)=c("GbyE.GD","GbyE.GM")
						return(data)
						}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~参数设置


nrep=100
#set.seed(99164)
ha2=c(0.5,0.5)
NQTN=20
NE=2 #NE*NE的矩阵
cov_g=c(0.001,0.999)
cov_e=0.6
file.output=1

GD=read.table("mdp_numeric.txt",head=F)
GM=read.table("mdp_SNP_information.txt",head=T)


mydata=GbyE.file.Calculate(GM=GM,GD=GD,file.output=file.output)
	if(file.output==1){
			myGD=read.table(paste(O_Path,"/","GbyE.GD.txt",sep=""),head=T)
			myGM=read.table(paste(O_Path,"/","GbyE.GM.txt",sep=""),head=T)
			}else{
						myGD=mydata[[1]]
						myGM=mydata[[2]]
						colnames(myGD)=myGD[1,]
						myGD=myGD[-1,]
						myGD=as.data.frame(myGD)
						myGM=as.data.frame(myGM)
						}


for (g in cov_g){
									print(paste("##############",g," & ",g," & ",g," & ",g," & ",g," & ",g,"##############",sep=""))
									print("######################################################################")
									print("######################################################################")
									print(paste("##############",g," Genetic.Correlation.Power.VS.FDR##############",sep=""))
									print("######################################################################")
									print("######################################################################")
									print(paste("##############",g," & ",g," & ",g," & ",g," & ",g," & ",g,"##############",sep=""))
									
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~replicate重复POWER运算
	StatRep=replicate(nrep,{
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GbyE.Y的模拟表型###############
		#GD=read.table(paste(O_Path,"/","mdp_numeric.txt",sep=""),head=F)
			fwrite(GD[-1,],"GD.temp.txt",row.names=F,sep=" ",quote=FALSE,col.names=F)
			GD.n=read.table("GD.temp.txt",head=F)
			file.remove("GD.temp.txt")

				rg=matrix(g,NE,NE)
				diag(rg)=1
				re=matrix(cov_e,NE,NE)
				diag(re)=1

			mysimu=G_E.Simulation(GD=GD.n,ha2=ha2,NQTN=NQTN,NE=NE,rg=rg,re=re)
			rm(rg)
			rm(re)
			rm(GD.n)
			myY=mysimu$Y
			Mean.Y=myY[,-3]
			Mean.Y[,2]=(myY[,2]+myY[,3])/2
			myY.1=cbind(paste(myY[,1],"-1",sep=""),myY[,2])
			myY.2=cbind(paste(myY[,1],"-2",sep=""),myY[,3])
			GbyE.Y=rbind(myY.1,myY.2)
			colnames(GbyE.Y)=c("taxa","trait")
			rm(myY.1)
			rm(myY.2)
			rm(myY)
									print("Simulation Y & Gbye.Simulation Y output")
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~创建处理结果文件夹及路径设置#####

			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~常规GWAS分析。
			GD=read.table(paste(O_Path,"/","mdp_numeric.txt",sep=""),head=T)
			GM=read.table(paste(O_Path,"/","mdp_SNP_information.txt",sep=""),head=T)
			Y=Mean.Y
			
									cat("Input GD file ing\nInput GM file ing\nInput Mean.Y file ing")
									cat("Running GAPIT to Mean GWAS")
									
			myGAPIT<-GAPIT(
					Y=Y, #fist column is individual ID, the third columns is days to pollination
					GD=GD,
					GM=GM,
					PCA.total=3,
					#QTN.position=mysimu$QTN.position,
					model="MLM",
					file.output=FALSE
			)
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~power&fdr结果运算
								cat("Mean GWAS Running\nMean GWAS Running\n
										\nGWAS Result order with P-Value\n
										\nUse the GWAS result about $SNP $Chrosome $Position $P-Value\n
										\nRunning FDR&Power")
										
					print("###########################Mean GWAS Power VS FDR Calculate###########################")
										
			myGWAS=myGAPIT$GWAS[,1:4]
			myGWAS=myGWAS[order(myGWAS[,4]),]
			cc=as.vector(GM[mysimu$QTN.position,][,1])
			Mean.data=Power.FDR.Calculate(cc=cc,myGWAS=myGWAS,myGM=GM,NQTN=NQTN)

			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GbyE
								cat("GbyE GWAS Running\nGbyE GWAS Running\n 
									\nGWAS Result order with P-Value\n
									\nUse the GWAS result about $SNP $Chrosome $Position $P-Value\n
									\nRunning FDR&Power")

			fwrite(GbyE.Y,"temp.txt",row.names=F,sep=" ",quote=FALSE,col.names=F)
			myY=read.table("temp.txt",head=F)
			file.remove("temp.txt")
			#myY=as.data.frame(GbyE.Y)
			QTN.position=mysimu$QTN.position+nrow(myGM)/2
			QTN.position=as.integer(c(mysimu$QTN.position,QTN.position))


			myGAPIT<-GAPIT(
					Y=myY, #fist column is individual ID, the third columns is days to pollination
					GD=myGD,
					GM=myGM,
					PCA.total=3,
					QTN.position=QTN.position,
					model="MLM",
					file.output=FALSE
			)
			
			#GbyE
			myGWAS=myGAPIT$GWAS[,1:4]
			GbyE.GWAS=as.data.frame(matrix(NA,nrow(myGM)/2,4))
			colnames(GbyE.GWAS)=c("SNP","Chrosome","Position","P.value")
			repp=nrow(myGM)/2
			for (Chr.pos in 1:repp){
					if(myGWAS[Chr.pos,4]<=myGWAS[Chr.pos+repp,4])
						{
						GbyE.GWAS[Chr.pos,]=myGWAS[Chr.pos,]
							}else{
									GbyE.GWAS[Chr.pos,]=myGWAS[Chr.pos+repp,]
									}
				}
					
					print("###########################GbyE GWAS Power VS FDR Calculate###########################")
					
			#result_transcripts=myGAPIT$GWAS[,1:4]
			rownames(GbyE.GWAS)=c(1:repp)
			cc=as.vector(GbyE.GWAS[mysimu$QTN.position,][,1])
			GbyE.GWAS=GbyE.GWAS[order(GbyE.GWAS[,4]),]
			GbyE.data=Power.FDR.Calculate(cc=cc,myGWAS=GbyE.GWAS,myGM=myGM[1:(nrow(myGM)/2),],NQTN=NQTN)
							
					print("###########################AddE GWAS Power VS FDR Calculate###########################")
					
			#Additive Eﬀects
			AddE.GWAS=myGAPIT$GWAS[1:(nrow(myGM)/2),1:4]
			cc=as.vector(AddE.GWAS[mysimu$QTN.position,][,1])
			AddE.GWAS=AddE.GWAS[order(AddE.GWAS[,4]),]
			AddE.data=Power.FDR.Calculate(cc=cc,myGWAS=AddE.GWAS,myGM=myGM[1:(nrow(myGM)/2),],NQTN=NQTN)
								
					print("###########################GGEE GWAS Power VS FDR Calculate###########################")
					
			#GXE
			GGEE.GWAS=myGAPIT$GWAS[-(1:(nrow(myGM)/2)),1:4]
			rownames(GGEE.GWAS)=c(1:repp)
			cc=as.vector(GGEE.GWAS[mysimu$QTN.position,][,1])
			GGEE.GWAS=GGEE.GWAS[order(GGEE.GWAS[,4]),]
			GGEE.data=Power.FDR.Calculate(cc=cc,myGWAS=GGEE.GWAS,myGM=myGM[1:(nrow(myGM)/2),],NQTN=NQTN)
			
			
			myStat=list(Mean.data,GbyE.data,AddE.data,GGEE.data)
		})


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
						}
		}
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
	
	Power.VS.FDR.Plot.Read=as.data.frame(rbind(
												as.matrix(cbind(as.data.frame(matrix("GbyE",nrow(Power.VS.FDR[1:2]),1)),Power.VS.FDR[1:2])),
												as.matrix(cbind(as.data.frame(matrix("Mean",nrow(Power.VS.FDR[3:4]),1)),Power.VS.FDR[3:4])),
												as.matrix(cbind(as.data.frame(matrix("AddE",nrow(Power.VS.FDR[5:6]),1)),Power.VS.FDR[5:6])),
												as.matrix(cbind(as.data.frame(matrix("GGEE",nrow(Power.VS.FDR[7:8]),1)),Power.VS.FDR[7:8]))
												))
	colnames(Power.VS.FDR.Plot.Read)=c("Method","Power","FDR")

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


			if(g>0.5){
				HH2="GHigh"
				}else{
						if(g==0.5){
								HH2="Mibble"
								}else{
										HH2="GLow"}
			}
			write.table(Power.VS.FDR,file=paste(HH2,".Genetic.Correlation.Power.VS.FDR.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}


pdf("1.pdf",width=8,height=8)
pdf("1.pdf",width=12,height=4)

par(mfrow=c(1,3))
	for (i in c("High","Middle","Low"))
			{
			if(i=="High"){
						cov_e=0.8
						}else{
							if(i=="Middle"){
									cov_e=0.5}else{
													cov_e=0.2}
						}
				
				#pdf(file=paste(i,".Genetic.Correlation.pdf",sep=""))
				Power.VS.FDR=read.table(paste(i,".Genetic.Correlation.Power.VS.FDR.txt",sep=""),head=T)
				#Power.VS.FDR=apply(Power.VS.FDR,2,as.numeric)
				plot(x=Power.VS.FDR[,2],y=Power.VS.FDR[,1],type="l",  #GbyE
						xlab="FDR",ylab="Power",
						xlim = c(0,1),ylim = c(0,1),
						col=rgb(253,10,7,255,maxColorValue=255),lwd=2,
						main=paste(i," Genetic Correlation, Value=",cov_e,sep=""))
				lines(x=Power.VS.FDR[,8],y=Power.VS.FDR[,7],type="l",col=rgb(141,172,208,255,maxColorValue=255),lwd=2)  #GGEE
				lines(x=Power.VS.FDR[,6],y=Power.VS.FDR[,5],type="l",col=rgb(11,11,11,255,maxColorValue=255),lwd=2)  #AddE
				lines(x=Power.VS.FDR[,4],y=Power.VS.FDR[,3],type="l",col=rgb(3,175,78,255,maxColorValue=255),lwd=2) #Mean
				legend("bottomright", inset=.05, c("GbyE","GXE","AddE","Mean"),lty=1,lwd=3,
						col=c(rgb(253,10,7,255,maxColorValue=255),
						rgb(141,172,208,255,maxColorValue=255),
						rgb(11,11,11,255,maxColorValue=255),
						rgb(3,175,78,255,maxColorValue=255)
						))
			}
			dev.off()
			
			
pdf("Heritability.pdf",width=12,height=4)

par(mfrow=c(1,3))
	for (i in c("High","Middle","Low"))
			{
			if(i=="High"){
						cov_e=0.8
						}else{
							if(i=="Middle"){
									cov_e=0.5}else{
													cov_e=0.2}
						}
				#pdf(file=paste(i,".Genetic.Correlation.pdf",sep=""))
				Power.VS.FDR=read.table(paste(i,".ha2.Power.VS.FDR.txt",sep=""),head=T)
				#Power.VS.FDR=apply(Power.VS.FDR,2,as.numeric)
				plot(x=Power.VS.FDR[,2],y=Power.VS.FDR[,1],type="l",  #GbyE
						xlab="FDR",ylab="Power",
						xlim = c(0,1),ylim = c(0,1),
						col=rgb(253,10,7,255,maxColorValue=255),lwd=2,
						main=paste(i," Heritability, Value=",cov_e,sep=""))
				lines(x=Power.VS.FDR[,8],y=Power.VS.FDR[,7],type="l",col=rgb(141,172,208,255,maxColorValue=255),lwd=2)  #GGEE
				lines(x=Power.VS.FDR[,6],y=Power.VS.FDR[,5],type="l",col=rgb(11,11,11,255,maxColorValue=255),lwd=2)  #AddE
				lines(x=Power.VS.FDR[,4],y=Power.VS.FDR[,3],type="l",col=rgb(3,175,78,255,maxColorValue=255),lwd=2) #Mean
				legend("bottomright", inset=.05, c("GbyE","GXE","AddE","Mean"),lty=1,lwd=3,
						col=c(rgb(253,10,7,255,maxColorValue=255),
						rgb(141,172,208,255,maxColorValue=255),
						rgb(11,11,11,255,maxColorValue=255),
						rgb(3,175,78,255,maxColorValue=255)
						))
			}
			dev.off()
		
		
		
		




################################################################################GGPLOT2##################
p1=ggplot(Power.VS.FDR) +
			geom_line(aes(GbyE_FDR,GbyE_power,color='red'),size=1)+
			geom_line(aes(Mean_FDR,Mean_power,colour='blue'),size=1)+
			geom_line(aes(AddE_FDR,AddE_power,colour="green"),size=1)+
			geom_line(aes(GGEE_FDR,GGEE_power,colour="black"),size=1)+
			#geom_line(size=1)+
			#scale_fill_manual(c("red","blue","green","black"),size=2)+
			#geom_smooth(aes(GbyE_FDR,GbyE_power),color="red")+
			#geom_point()+
			#geom_smooth(aes(Mean_FDR,Mean_power),color="blue")+
			#geom_point()+
			#geom_point(aes(GbyE_FDR,GbyE_power,colour="red"),shape=16,size=0.5)+
			#geom_point(aes(Mean_FDR,Mean_power,color="blue"),shape=16,size=0.5)+
			#geom_point(aes(AddE_FDR,AddE_power,color="green"),shape=16,size=0.5)+
			#geom_point(aes(GGEE_FDR,GGEE_power,color="black"),shape=16,size=0.5)+
            labs(x = 'FDR',y = 'Power')+
			xlim(0,1)+ylim(0,1)+
			scale_colour_discrete(name="Method",labels=c("GbyE","Mean","AddE","GGEE"))+
			#guides(shape = guide_legend(override.aes = list(size = 15)))+
			#scale_color_manual(c("red","blue","green","black"))+
			theme(
					axis.title=element_text(size = 14), #坐标轴的标题
					axis.text =element_text(size = 12, color = 'black'), #坐标轴的大小标题
					axis.text.x = element_text(angle = 65, hjust = 1), #X轴字体的大小偏转角度等
					panel.background = element_blank(), #去掉绘图背景
					axis.line = element_line(colour = "black"), #线条标题
					legend.key.size = unit(30, "pt"), #图例大小
					legend.key = element_blank(), #图例背景去掉
					legend.text = element_text(face = "bold",size = 12,margin = margin(r=20)), #图例文本大小
					legend.title=element_blank(),
					legend.position=c(0.9, 0.25), #图例位置
					text=element_text(size=16,  family="sans") #设置字体
					)
p1
ggsave(p1, file="ratings.pdf", width=16, height=12)
#plot(P.O.mean[,4], pogwas, type="b", col="blue",xlim=c(0,1))
#lines(P.G.mean[,4], pggwas, type="b", col= "red")

