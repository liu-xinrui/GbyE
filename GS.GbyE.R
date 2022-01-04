#Import files
#######################################################################################
#myY <- read.table("mdp_traits.txt", head = TRUE)
#myKI <- read.table("KSN.txt", head = FALSE)
#myCV <- read.table("Copy of Q_First_Three_Principal_Components.txt", head = TRUE)
#Initial
#GD、GM为Mean
#myGD、myGM为GbyE的
#######################################################################################
rm(list=ls())

setwd("C:/Users/VEN MIKIR/Desktop/GbyE.test")
#setwd("/home/student/workshop/lxr/GbyE")
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("G&E_Simulation.R")
source("F:/swu-青研院/实验室/code/Power.FDR.Calculate.R")
source("F:/swu-青研院/实验室/code/GbyE.file.Calculate.R")
#source("/home/student/workshop/lxr/software/code/Power.FDR.Calculate.R")
#source("/home/student/workshop/lxr/software/code/GbyE.file.Calculate.R")
library('MASS')
library(data.table)

CV.function <- function(GD,GM,NQTN,GWAS)
				{
													ntop=NQTN
													index.raw=GWAS
													index=index.raw[order(index.raw[,4]),]
													index=subset(index,index[,4]<=0.01/nrow(GM))
													if(nrow(index)>=ntop)
													{
																index=index[1:ntop,]
																indexs <- GM[,1] %in% index[,1]
																CV=cbind(data.frame(GD[,1]),data.frame(GD[,indexs]))
																colnames(CV)=c("taxa",index[,1])
													}else  {
																		index=order(index.raw[,4])
																		top=index[1:ntop]
																		CV=GD[,c(1,top+1)]
															}
					return(CV)
					}


Comparison.GWAS.Reslut.PValue <- function(GD,GM,GWASresult){
													index.raw=GWASresult$GWAS[1:4]
													index.tem=as.data.frame(matrix(NA,nrow(GM)/2,4))
													colnames(index.tem)=c("SNP","Chrosome","Position","P.value")
													repp=nrow(GM)/2
													for (Chr.pos in 1:repp){
															if(index.raw[Chr.pos,4]<=index.raw[Chr.pos+repp,4])
																{
																index.tem[Chr.pos,]=index.raw[Chr.pos,]
																	}else{
																			index.tem[Chr.pos,]=index.raw[Chr.pos+repp,]
																			}
														}
														return(index.tem)
														}
		
					
					
O_Path=getwd()


t=3 #total replicates
nfold=3 #缺失比例
#ntop=10
#nrep=100
#set.seed(99164)
ha2=c(0.8,0.8)
NQTN=20
NE=2 #NE*NE的矩阵
g=0.2
cov_g=0.5
cov_e=0.6
file.output=T

set.seed(15555)

GD=read.table(file="http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=F)#导入基因型数据
GM=read.table(file="http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T) #导入SNP信息（SNP的染色体编号及位置）
#myCV=read.table(file="http://zzlab.net/GAPIT/data/mdp_env.txt",head=T)
#GD=read.table("mdp_numeric.txt",head=F)
#GM=read.table("mdp_SNP_information.txt",head=T)

mydata=GbyE.file.Calculate(GM=GM,GD=GD,file.output=file.output)
	if(file.output==TRUE)
			{
			myGD=read.table(paste(O_Path,"/","GbyE.GD.txt",sep=""),head=T,check.names=FALSE)
			myGM=read.table(paste(O_Path,"/","GbyE.GM.txt",sep=""),head=T)
			}else{
					myGD=mydata[[1]]
					myGM=mydata[[2]]
					colnames(myGD)=myGD[1,]
					myGD=myGD[-1,]
					myGD=as.data.frame(myGD)
					myGM=as.data.frame(myGM)
					}

			fwrite(GD[-1,],"GD.temp.txt",row.names=F,sep=" ",quote=FALSE,col.names=F)
			GD.n=read.table("GD.temp.txt",head=F)
			file.remove("GD.temp.txt")

				rg=matrix(g,NE,NE)
				diag(rg)=1
				re=matrix(cov_e,NE,NE)
				diag(re)=1
			
#rm(rg)
#rm(re)
#rm(GD.n)
#rm(myY.1)
#rm(myY.2)
#rm(myY)
#rm(mydata)

GD=read.table("http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)

n=nrow(GD)
nq=c(1:n) #Mean样本总体数量
nq2=c(1:(n*2)) #GbyE样本总体数量
#s=1/5 #sample of inference, e.g. set it to 1/5 for five fold cross validation
#n.missing=round(n*s) #四舍五入
result=matrix(NA,t,6)
result.AVO=matrix(NA,t,6)
result.rep=matrix(NA,6,nfold)
Mean.real.Matrix=as.data.frame(matrix(NA,n,9))
Mean.Pred.Matrix=as.data.frame(matrix(NA,n,9))
Mean.noCV.Matrix=as.data.frame(matrix(NA,n,9))
GbyE.real.Matrix=as.data.frame(matrix(NA,2*n,9))
GbyE.Pred.Matrix=as.data.frame(matrix(NA,2*n,9))
GbyE.noCV.Matrix=as.data.frame(matrix(NA,2*n,9))

#Loop on replicates
for(rep in 1:t)
	{
		mysimu=G_E.Simulation(GD=GD.n,ha2=ha2,NQTN=NQTN,NE=NE,rg=rg,re=re)
			myY=mysimu$Y
			Mean.Y=myY[,-3]
			Mean.Y[,2]=(myY[,2]+myY[,3])/2
			myY.1=cbind(paste(myY[,1],"-1",sep=""),myY[,2])
			myY.2=cbind(paste(myY[,1],"-2",sep=""),myY[,3])
			GbyE.Y=rbind(myY.1,myY.2)
			colnames(GbyE.Y)=c("taxa","trait")
			GbyE.Y=as.data.frame(as.matrix(data.table(GbyE.Y)))
			GbyE.Y[,2]=as.numeric(GbyE.Y[,2])
			#GbyE.Y[,2]=apply(GbyE.Y[,2],1,as.numeric)
			if(file.output==TRUE)
				{
					fwrite(myY,"myY.txt",row.names=F,sep=" ",quote=FALSE,col.names=T)
					print("Output the original simulated phenotype")
					fwrite(Mean.Y,"Mean.Y.txt",row.names=F,sep=" ",quote=FALSE,col.names=T)
					print("Output the Mean simulated phenotype")
					fwrite(GbyE.Y,"GbyE.Y.txt",row.names=F,sep=" ",quote=FALSE,col.names=T)
					print("Output the GbyE simulated phenotype")
					
				}else{
						print("No simulated phenotype output")
						}
			
			QTN.Position=mysimu$QTN.position

		#Mean.Y.raw=Mean.Y.raw[!is.na(Mean.Y.raw[,2]),]
		#GbyE.Y.raw=GbyE.Y.raw[!is.na(GbyE.Y.raw[,2]),] #Remove missing data
		#n=nrow(Mean.Y.raw)
		sets=sample(cut(1:n,nfold,labels=FALSE),n)
	
		for(xfold in 1:nfold)
			{
				Mean.Y.raw=Mean.Y[,c(1,2)] #choos a trait
				GbyE.Y.raw=GbyE.Y[,c(1,2)]
				
				Mean.Y.X=sets==xfold
				GbyE.Y.X=c(Mean.Y.X,Mean.Y.X)
				Mean.Y.raw[Mean.Y.X,2]=NA
				GbyE.Y.raw[GbyE.Y.X,2]=NA
				GbyE.Y.111=GbyE.Y.raw
				Mean.Y.111=Mean.Y.raw
				Mean.Y.raw=Mean.Y.raw[!is.na(Mean.Y.raw[,2]),]
				GbyE.Y.raw=GbyE.Y.raw[!is.na(GbyE.Y.raw[,2]),]
		
				raw.M=c(as.numeric(rownames(Mean.Y.raw)),nq)
				raw.M=raw.M[duplicated(raw.M)]
				raw.M <- nq[!(nq %in% raw.M)]
				raw.G=c(as.numeric(rownames(GbyE.Y.raw)),nq2)
				raw.G=raw.G[duplicated(raw.G)]
				raw.G <- nq2[!(nq2 %in% raw.G)]
		#testing1=sample(n,round(n/nfold),replace=F)#设置20%的个体为试验群体
		#testing2=sample(n-testing1)
		#folds <- createFolds(Mean.Y.raw,k=3)	#
		#Mean.Y.raw=Mean.Y.raw[-testing,]
		#GbyE.Y.raw=GbyE.Y.raw[-c(testing,(testing+n)),]
#a=as.data.frame(cbind(as.data.frame(colnames(GD))[-1,],GM[,1]))
#a=cbind(as.data.frame(data.frame(colnames(GD))[QTN.Position+1,]),GM[QTN.Position,1])

		#Mean
				#GWAS	  
				Mean.GWAS <- GAPIT(
									Y=Mean.Y.raw,
									GD=GD,
									GM=GM,
									#QTN.position=QTN.Position,
									PCA.total=3,
									model="MLM",
									file.output=FALSE
								)
				#GWAS.Result do CV SNP			
				#提取我们需要的P值最小的20个值，作为协变量CV，加入GAPIT运算，GS
				PCA=Mean.GWAS$PCA
																		
				#Real SNP CV
				#CV1=cbind(GD[,c(1,QTN.Position+1)])
				CV1=cbind(PCA,GD[,QTN.Position+1])

				CV2=CV.function(GD=GD,GM=GM,GWAS=Mean.GWAS$GWAS,NQTN=NQTN)
				CV2=cbind(PCA,CV2[,-1])

				#Prediction
				#1.Real CV
				Mean.real <- GAPIT(
									Y=Mean.Y.raw, #fist column is individual ID, the third columns is days to pollination
									GD=GD,
									GM=GM,
									CV=CV1,
									#PCA.total=3,
									model="gBLUP",
									SNP.test=FALSE,
									file.output=FALSE
								)
								
								
					Mean.real.prediction=Mean.real$Pred
					#Separate reference (with phenotype) and inference (without phenotype)
					Mean.real=Mean.real.prediction[Mean.real.prediction[,3]==2,]
					#Merge prediction with original Y
					Mean.real <- merge(Mean.Y, Mean.real, by.x = "taxa", by.y = "Taxa")
					Mean.real.Matrix[raw.M,]=Mean.real
				
					
					#Mean.real=rbind
					
					#Calculate correlation and store them
					result.rep[1,xfold]=cor(as.numeric(as.vector(Mean.real[,2])),as.numeric(as.vector(Mean.real[,9])))
					if(xfold==nfold)
							{
								Mean.real.raw=cor(as.numeric(as.vector(Mean.real.Matrix[,2])),as.numeric(as.vector(Mean.real.Matrix[,9])))
								result[rep,1]=Mean.real.raw
								result.AVO[rep,1]=rowMeans(result.rep)[1]
								}else{
										print(paste("No",nfold,"Times",sep=" "))
										}

				#2.Pred CV	
				Mean.Pred <- GAPIT(
									Y=Mean.Y.raw, #fist column is individual ID, the third columns is days to pollination
									GD=GD,
									GM=GM,
									CV=CV2,
									#PCA.total=3,
									model=c("gBLUP"),
									SNP.test=FALSE,
									file.output=FALSE
								)
								
								
					Mean.Pred.prediction=Mean.Pred$Pred
					#Separate reference (with phenotype) and inference (without phenotype)
					Mean.Pred=Mean.Pred.prediction[Mean.Pred.prediction[,3]==2,]
					#Merge prediction with original Y
					Mean.Pred <- merge(Mean.Y, Mean.Pred, by.x = "taxa", by.y = "Taxa")
					#Calculate correlation and store them
					Mean.Pred.Matrix[raw.M,]=Mean.Pred
				
					
					#Mean.real=rbind
					
					#Calculate correlation and store them
					
					result.rep[2,xfold]=cor(as.numeric(as.vector(Mean.Pred[,2])),as.numeric(as.vector(Mean.Pred[,9])))
					if(xfold==nfold)
							{
								Mean.Pred.raw=cor(as.numeric(as.vector(Mean.Pred.Matrix[,2])),as.numeric(as.vector(Mean.Pred.Matrix[,9])))
								result[rep,2]=Mean.Pred.raw
								result.AVO[rep,2]=rowMeans(result.rep)[2]
								}else{
										print(paste("No",nfold,"Times",sep=" "))
										}
										
				#3.NO CV	
				Mean.noCV <- GAPIT(
									Y=Mean.Y.raw, #fist column is individual ID, the third columns is days to pollination
									GD=GD,
									GM=GM,
									#CV=CV2,
									#PCA.total=3,
									model=c("gBLUP"),
									SNP.test=FALSE,
									file.output=FALSE
								)	
					Mean.noCV.prediction=Mean.noCV$Pred
					#Separate reference (with phenotype) and inference (without phenotype)
					Mean.noCV=Mean.noCV.prediction[Mean.noCV.prediction[,3]==2,]
					#Merge prediction with original Y
					Mean.noCV <- merge(Mean.Y, Mean.noCV, by.x = "taxa", by.y = "Taxa")
					#Calculate correlation and store them
					Mean.noCV.Matrix[raw.M,]=Mean.noCV

					#Calculate correlation and store them
					
					result.rep[3,xfold]=cor(as.numeric(as.vector(Mean.noCV[,2])),as.numeric(as.vector(Mean.noCV[,9])))
					if(xfold==nfold)
							{
								Mean.noCV.raw=cor(as.numeric(as.vector(Mean.noCV.Matrix[,2])),as.numeric(as.vector(Mean.noCV.Matrix[,9])))
								result[rep,3]=Mean.noCV.raw
								result.AVO[rep,3]=rowMeans(result.rep)[3]
								}else{
										print(paste("No",nfold,"Times",sep=" "))
										}

		#GbyE
				#GWAS	  
				GbyE.GWAS <- GAPIT(
									Y=GbyE.Y.raw,
									GD=myGD,
									GM=myGM,
									PCA.total=3,
									model="MLM",
									file.output=FALSE
								)
			#Function do Comparison Pvalue by GbyE.Pvalue, the smart we used 
				#Comparison.GWAS.Reslut.PValue
		
				#Running the function,
				#GWAS.Pvalue=Comparison.GWAS.Reslut.PValue(GD=myGD,GM=myGM,GWASresult=GbyE.GWAS)
					
				PCA=GbyE.GWAS$PCA
				#Real SNP CV
				CV1=myGD[,c(1,QTN.Position+1,(QTN.Position+(nrow(myGM)/2))+1)]
				CV1=cbind(PCA,CV1[,-1])
				
										
				#提取我们需要的P值最小的10个值，作为协变量CV，加入GAPIT运算，GS											
				CV2=CV.function(GD=myGD,GM=myGM,GWAS=GbyE.GWAS$GWAS,NQTN=NQTN*2)
				CV2=cbind(PCA,CV2[,-1])
										
				#Prediction
				#1.Real CV
				GbyE.real <- GAPIT(
									Y=GbyE.Y.raw, 
									GD=myGD,
									GM=myGM,
									CV=CV1,
									#PCA.total=3,
									model=c("gBLUP"),
									SNP.test=FALSE,
									file.output=FALSE
								)
								
				GbyE.real.prediction=GbyE.real$Pred
				#Separate reference (with phenotype) and inference (without phenotype)
				GbyE.real=GbyE.real.prediction[GbyE.real.prediction[,3]==2,]
				#Merge prediction with original Y
				GbyE.real <- merge(GbyE.Y, GbyE.real, by.x = "taxa", by.y = "Taxa")
				#Calculate correlation and store them
				GbyE.real.Matrix[raw.G,]=GbyE.real

				#Calculate correlation and store them
				
				result.rep[4,xfold]=cor(as.numeric(as.vector(GbyE.real[,2])),as.numeric(as.vector(GbyE.real[,9])))
				if(xfold==nfold)
						{	
							GbyE.real.raw=cor(as.numeric(as.vector(GbyE.real.Matrix[,2])),as.numeric(as.vector(GbyE.real.Matrix[,9])))
							result[rep,4]=GbyE.real.raw
							result.AVO[rep,4]=rowMeans(result.rep)[4]
							}else{
									print(paste("No",nfold,"Times",sep=" "))
									}

				#2.Pred CV
				GbyE.Pred <- GAPIT(
									Y=GbyE.Y.raw, 
									GD=myGD,
									GM=myGM,
									CV=CV2,
									#PCA.total=3,
									model=c("gBLUP"),
									SNP.test=FALSE,
									file.output=FALSE
								)

				GbyE.Pred.prediction=GbyE.Pred$Pred
				#Separate reference (with phenotype) and inference (without phenotype)
				GbyE.Pred=GbyE.Pred.prediction[GbyE.Pred.prediction[,3]==2,]
				#Merge prediction with original Y
				GbyE.Pred <- merge(GbyE.Y, GbyE.Pred, by.x = "taxa", by.y = "Taxa")
				#Calculate correlation and store them
				GbyE.Pred.Matrix[raw.G,]=GbyE.Pred

				#Calculate correlation and store them
				result.rep[5,xfold]=cor(as.numeric(as.vector(GbyE.Pred[,2])),as.numeric(as.vector(GbyE.Pred[,9])))
				if(xfold==nfold)
						{	
							GbyE.Pred.raw=cor(as.numeric(as.vector(GbyE.Pred.Matrix[,2])),as.numeric(as.vector(GbyE.Pred.Matrix[,9])))
							result[rep,5]=GbyE.Pred.raw
							result.AVO[rep,5]=rowMeans(result.rep)[5]
							}else{
									print(paste("No",nfold,"Times",sep=" "))
									}

				#3.no CV
				GbyE.noCV <- GAPIT(
									Y=GbyE.Y.raw, 
									GD=myGD,
									GM=myGM,
									#CV=CV1,
									#PCA.total=3,
									model=c("gBLUP"),
									SNP.test=FALSE,
									file.output=FALSE
								)

				GbyE.noCV.prediction=GbyE.noCV$Pred
				#Separate reference (with phenotype) and inference (without phenotype)
				GbyE.noCV=GbyE.noCV.prediction[GbyE.noCV.prediction[,3]==2,]
				#Merge prediction with original Y
				GbyE.noCV <- merge(GbyE.Y, GbyE.noCV, by.x = "taxa", by.y = "Taxa")
				#Calculate correlation and store them
				GbyE.noCV.Matrix[raw.G,]=GbyE.noCV

				#Calculate correlation and store them
				result.rep[6,xfold]=cor(as.numeric(as.vector(GbyE.noCV[,2])),as.numeric(as.vector(GbyE.noCV[,9])))
				if(xfold==nfold)
						{	
							GbyE.noCV.raw=cor(as.numeric(as.vector(GbyE.noCV.Matrix[,2])),as.numeric(as.vector(GbyE.noCV.Matrix[,9])))
							result[rep,6]=GbyE.noCV.raw
							result.AVO[rep,6]=rowMeans(result.rep)[6]
							}else{
									print(paste("No",nfold,"Times",sep=" "))
									}
				}
} #End of for (rep in 1:t)

colnames(result)=c("Mean-realCV","Mean-predCV","Mean-noCV","GbyE-realCV","GbyE-predCV","GbyE-noCV")
colnames(result.AVO)=c("Mean-realCV","Mean-predCV","Mean-noCV","GbyE-realCV","GbyE-predCV","GbyE-noCV")
write.table(result.AVO, paste("Average.",t,".Times.",nfold,".Group",".txt",sep=""), quote = FALSE, sep = "\t")
write.table(result, paste("Overall.",t,".Times.",nfold,".Group",".txt",sep=""), quote = FALSE, sep = "\t")


#ggplot(alpha2, aes(x = group2, y = value, fill = group1)) + 
#  geom_boxplot(outlier.size = 1) +
#  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank()) +
#  labs(x = '', y = 'Chao1') 


				YY=GbyE.Y.111[-1]
				rownames(YY)=GbyE.Y.111[,1]
				GG=myGD[,-1]
				colnames(GG)=c(1:ncol(GG))
				GG <- mapvalues(apply(GG,2,as.numeric),from=c(0,2,1),to=c(1,-1,0))
				markers_impute=GG
				yield_answer <- mixed.solve(as.matrix(YY), Z=markers_impute)
				YLD<-yield_answer$u
				e=as.matrix(YLD)
				pred_yield_valid=markers_impute %*% e
				pred_yield=pred_yield_valid[,1]+c(yield_answer$beta)
				YLD_accuracy<-cor(pred_yield_valid[raw.G,],GbyE.Y[raw.G,2],use="complete")
				YLD_accuracy
				
				#EBV0=as.matrix(GG)%*%GbyE.noCV$u
				#EBV0=cbind(GbyE.Y.111,EBV0)
				#EBV0=cbind(GbyE.Y[raw.G,],EBV0[raw.G,3])
				#cor(EBV0[,2],EBV0[,3])
				

				YY=Mean.Y.111[-1]
				rownames(YY)=Mean.Y.111[,1]
				GG=GD[,-1]
				colnames(GG)=c(1:ncol(GG))
				GG <- mapvalues(apply(GG,2,as.numeric),from=c(0,2,1),to=c(1,-1,0))
				#impute<-A.mat(GG,max.missing = 0.5,impute.method = "mean",return.imputed = T)
				#markers_impute=impute$imputed
				markers_impute=GG
				yield_answer <- mixed.solve(as.matrix(YY), Z=markers_impute)
				
				
				YLD<-yield_answer$u
				e=as.matrix(YLD)
				pred_yield_valid=markers_impute %*% e
				pred_yield=pred_yield_valid[,1]+c(yield_answer$beta)
				YLD_accuracy<-cor(pred_yield_valid[raw.M,],Mean.Y[raw.M,2],use="complete")
				YLD_accuracy
				#EBV0=markers_impute %*% Mean.noCV$u
				#EBV0=cbind(Mean.Y.111,EBV0)
				#EBV0=cbind(Mean.Y[raw.M,],EBV0[raw.M,3])
				#cor(EBV0[,2],EBV0[,3])





