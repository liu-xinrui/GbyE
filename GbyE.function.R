GbyE.file.Calculate <- function(GD,GM,file.output){
if(!require(data.table)) install.packages("data.table")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GbyE.GM的处理
	GM.L1=GM
	GM.L2=GM
	GM.L2[,2]=GM.L2[,2]+max(GM.L2[,2])
	GM.L1[,1]=paste(GM.L1[,1],"-1",sep="")
	GM.L2[,1]=paste(GM.L2[,1],"-2",sep="")
	GbyE.GM=rbind(GM.L1,GM.L2)
	#rm(GM.L1)
	#rm(GM.L2)
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


`GAPIT.Multiple.Manhattan` <-
function(model_store,DPP=50000,chor_taxa=NULL,cutOff=0.01,band=5,seqQTN=NULL,byTraits=FALSE,
    Y.names=NULL,GM=NULL,interQTN=NULL,WS=10e5,outpch=NULL,inpch=NULL,
    plot.style="Oceanic",plot.line=TRUE,plot.type=c("h","s","w")){
    #Object: Make a Manhattan Plot
    #Output: pdfs of the Multiple Manhattan Plot
    #Authors: Zhiwu Zhang and Jiabo Wang
    # Last update: AUG 24, 2022
    ##############################################################################################
  Nenviron=length(model_store)*length(Y.names)
  environ_name=NULL
  new_xz=NULL
  if(byTraits)
  {
    for(i in 1:length(Y.names))
    {
       for(j in 1:length(model_store))
       {
      # environ_name=c(environ_name,paste(model_store[i],".",Y.names[j],sep=""))
          environ_name=c(environ_name,paste(model_store[j],".",Y.names[i],sep=""))
       }
    }
  }else{
    for(i in 1:length(model_store))
    {
       for(j in 1:length(Y.names))
       {
      # environ_name=c(environ_name,paste(model_store[i],".",Y.names[j],sep=""))
          environ_name=c(environ_name,paste(model_store[i],".",Y.names[j],sep=""))
       }
    }
  }
sig_pos=NULL
simulation=FALSE
    if(!is.null(seqQTN)){    
        #seqQTN=-seqQTN
        simulation=TRUE    
    }
themax.y0=NULL
store.x=NULL
y_filter0=NULL
for(i in 1:length(environ_name))
{
  print(paste("Reading GWAS result with ",environ_name[i],sep=""))
  environ_result=read.csv(paste("GAPIT.Association.GWAS_Results.",environ_name[i],".csv",sep=""),head=T)
  environ_result=environ_result[order(environ_result[,3]),]
  environ_result=environ_result[order(environ_result[,2]),]
  environ_filter=environ_result[!is.na(environ_result[,4]),]
  themax.y=round(max(-log10(environ_filter[,4])),0)
  themax.y0=round(max(c(themax.y,themax.y0)),0)
  chm.to.analyze <- unique(environ_result[,2])
  nchr=length(chm.to.analyze)
  # xx=environ_result[,-1]
  # ticks=NULL
  # lastbase=0
  # max.x=NULL

  # for(ii in chm.to.analyze)
  #   {
  #   index=(xx[,1]==ii)
  #   ticks <- c(ticks, lastbase+mean(xx[index,2]))
  #   xx[index,2]=xx[index,2]+lastbase
  #   lastbase=max(xx[index,2])
  #   max.x=c(max.x,max(as.numeric(xx[index,2])))
  #   }
  # max.x=c(min(as.numeric(xx[,2])),max.x)
  # store.x=c(store.x,as.numeric(xx[,2]))
  y_filter=environ_filter[environ_filter[,4]<(cutOff/(nrow(environ_filter))),,drop=FALSE]
  traits=environ_name[i]
  # print(head(y_filter))
  # print(traits)
  if(nrow(y_filter)>0)y_filter=cbind(y_filter[,1:5],traits)
  y_filter0=rbind(y_filter0,y_filter)
  # write.table(y_filter,paste("GAPIT.Filter_",environ_name[i],"_GWAS_result.txt",sep=""))

  result=environ_result[,1:4]
  result=result[match(as.character(GM[,1]),as.character(result[,1])),]
  rownames(result)=1:nrow(result)
  #print(i)
  if(i==1){
    result0=result
    colnames(result0)[4]=environ_name[i]
    }
  if(i!=1){
    result0=merge(result0,result[,c(1,4)],by.x=colnames(result0)[1],by.y=colnames(result)[1])
    colnames(result0)[i+3]=environ_name[i]
    }
  rownames(result)=1:nrow(result)
  result[is.na(result[,4]),4]=1
  # map_store=max.x
  sig_pos=append(sig_pos,as.numeric(rownames(result[result[!is.na(result[,4]),4]<(cutOff/nrow(result)),,drop=FALSE])))
}
  write.csv(y_filter0,paste("GAPIT.Association.Filter_GWAS_results.csv",sep=""),quote=FALSE)

# print(sig_pos)
#if(length(sig_pos)!=0)sig_pos=sig_pos[!duplicated(sig_pos)]
 if(length(sig_pos[!is.na(sig_pos)])>1)
 {
 # {     x_matrix=as.matrix(table(sig_pos))
 #       x_matrix=cbind(as.data.frame(rownames(x_matrix)),x_matrix)
       #print(x_matrix)
        lastbase=0
        map_store=cbind(as.data.frame(GM[,2]),as.numeric(GM[,3]))
        ticks=NULL
        # print(head(map_store))
        max.x=NULL
        for (j in unique(map_store[,1]))
        {
            index=map_store[,1]==j
            ticks <- c(ticks, lastbase+mean(map_store[index,2]))
            map_store[index,2]=as.numeric(map_store[index,2])+lastbase
            lastbase=max(as.numeric(map_store[index,2]))
            max.x=c(max.x,max(as.numeric(map_store[index,2])))
        }
        # print(ticks)
       max.x=c(min(as.numeric(map_store[,2])),max.x)
       store.x=c(store.x,as.numeric(map_store[,2]))
       # colnames(x_matrix)=c("pos","times")
       new_xz0=cbind(sig_pos,map_store[as.numeric(as.character(sig_pos)),,drop=FALSE])
       common=as.numeric(new_xz0[,3])
       scom=sort(common)
       de.sc=scom[-1]-scom[-length(scom)]
       dayu1.index=duplicated(scom)|c(abs(de.sc)<WS,FALSE)

       print(table(dayu1.index))
       if(sum(dayu1.index)>0)
       {
       scom2=scom[dayu1.index]
       # scom2=scom2[!duplicated(scom2)]
       # print(new_xz0)
       # print(scom2)
       sc.index=as.character(new_xz0[,3])%in%scom2
       # print(table(sc.index))
       new_xz=new_xz0[sc.index,,drop=FALSE]
       # print(new_xz)
       new_xz=cbind(new_xz[,1],2,new_xz[,-1])
       new_xz[duplicated(new_xz[,4]),2]=1
       colnames(new_xz)=c("pos","times","chro","xlab")
       new_xz=new_xz[!duplicated(new_xz),]
       new_xz=as.matrix(new_xz)
       new_xz=new_xz[new_xz[,2]!="0",]
       new_xz=matrix(as.numeric(new_xz),length(as.vector(new_xz))/4,4)
       }
       # print(head(new_xz))
}else{
        lastbase=0
        map_store=cbind(as.data.frame(GM[,2]),as.numeric(GM[,3]))
        ticks=NULL
        max.x=NULL
        # print(head(map_store))
        for (j in unique(map_store[,1]))
        {
            index=map_store[,1]==j
            ticks <- c(ticks, lastbase+mean(map_store[index,2]))
            map_store[index,2]=as.numeric(map_store[index,2])+lastbase
            lastbase=max(as.numeric(map_store[index,2]))
            max.x=c(max.x,max(as.numeric(map_store[index,2])))
        }
       max.x=c(min(as.numeric(map_store[,2])),max.x)
       store.x=c(store.x,as.numeric(map_store[,2]))
}

# print(new_xz)
# setup colors
# print(head(result))
# chm.to.analyze <- unique(result[,2])
nchr=length(chm.to.analyze)
size=1 #1
ratio=10 #5
base=1 #1
numCHR=nchr
numMarker=nrow(GM)
bonferroniCutOff=-log10(cutOff/numMarker)

ncycle=ceiling(nchr/5)
ncolor=band*ncycle
thecolor=seq(1,nchr,by= ncycle)
col.Rainbow=rainbow(ncolor+1)     
col.FarmCPU=rep(c("#CC6600","deepskyblue","orange","forestgreen","indianred3"),ceiling(numCHR/5))
col.Rushville=rep(c("orangered","navyblue"),ceiling(numCHR/2))    
col.Congress=rep(c("deepskyblue3","firebrick"),ceiling(numCHR/2))
col.Ocean=rep(c("steelblue4","cyan3"),ceiling(numCHR/2))    
col.PLINK=rep(c("gray10","gray70"),ceiling(numCHR/2))     
col.Beach=rep(c("turquoise4","indianred3","darkolivegreen3","red","aquamarine3","darkgoldenrod"),ceiling(numCHR/5))
col.Oceanic=rep(c(  '#EC5f67',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5'),ceiling(numCHR/5))
col.cougars=rep(c(  '#990000',    'dimgray'),ceiling(numCHR/2))  
if(plot.style=="Rainbow")plot.color= col.Rainbow
if(plot.style =="FarmCPU")plot.color= col.Rainbow
if(plot.style =="Rushville")plot.color= col.Rushville
if(plot.style =="Congress")plot.color= col.Congress
if(plot.style =="Ocean")plot.color= col.Ocean
if(plot.style =="PLINK")plot.color= col.PLINK
if(plot.style =="Beach")plot.color= col.Beach
if(plot.style =="Oceanic")plot.color= col.Oceanic
if(plot.style =="cougars")plot.color= col.cougars  

if("h"%in%plot.type)
{
    Max.high=6*Nenviron
    if(Max.high>8)Max.high=40
    pdf(paste("GAPIT.Association.Manhattans_High",".pdf" ,sep = ""), width = 20,height=6*Nenviron)
    par(mfrow=c(Nenviron,1))
    mypch=1
    for(k in 1:Nenviron)
    { 
       if(k==Nenviron)
        {
        par(mar = c(3.5,8,0,8))
         # par(pin=c(10,((8-mtext.h)/Nenviron)+mtext.h))

        }else{
            #par(mfrow=c(Nenviron,1))
        par(mar = c(1.5,8,0.5,8))    
        }
       environ_result=read.csv(paste("GAPIT.Association.GWAS_Results.",environ_name[k],".csv",sep=""),head=T)
       result=environ_result[,1:4]
       result=result[order(result[,3]),]
       result=result[order(result[,2]),]
       result=result[match(as.character(GM[,1]),as.character(result[,1])),]
       rownames(result)=1:nrow(result)
       GI.MP=result[,c(2:4)]
       borrowSlot=4
       GI.MP[,borrowSlot]=0 #Inicial as 0
       GI.MP[,5]=1:(nrow(GI.MP))
       GI.MP[,6]=1:(nrow(GI.MP)) 
       GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
       GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
       GI.MP[is.na(GI.MP[,3]),3]=1
    #Retain SNPs that have P values between 0 and 1 (not na etc)
       GI.MP <- GI.MP[GI.MP[,3]>0,]
       GI.MP <- GI.MP[GI.MP[,3]<=1,]
    #Remove chr 0 and 99
       GI.MP <- GI.MP[GI.MP[,1]!=0,]
       total_chromo=length(unique(GI.MP[,1]))
    # print(dim(GI.MP))
       if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1
       numMarker=nrow(GI.MP)
       GI.MP[,3] <-  -log10(GI.MP[,3])
       GI.MP[,5]=1:numMarker
       y.lim <- ceiling(max(GI.MP[,3]))  
       chm.to.analyze <- unique(GI.MP[,1])
       nchr=length(chm.to.analyze)
       GI.MP[,6]=1:(nrow(GI.MP))
       MP_store=GI.MP
       index_GI=MP_store[,3]>=0
       MP_store <- MP_store[index_GI,]
       ticks=NULL
       lastbase=0
       # print(head(MP_store))
       for(i in chm.to.analyze)
          {
           index=(MP_store[,1]==i)
           ticks <- c(ticks, lastbase+mean(MP_store[index,2]))
           MP_store[index,2]=MP_store[index,2]+lastbase
           lastbase=max(MP_store[index,2])
          }
       x0 <- as.numeric(MP_store[,2])
       y0 <- as.numeric(MP_store[,3])
       z0 <- as.character(MP_store[,1])
       # convert chromosome character to number
       chor_taxa=as.character(unique(MP_store[,1]))
       chor_taxa=chor_taxa[order(as.numeric(as.character(chor_taxa)))]
       chr_letter=grep("[A-Z]|[a-z]",chor_taxa)
       if(!setequal(integer(0),chr_letter))
         {     
           z0=as.character(MP_store[,1])
           for(i in 1:(length(chor_taxa)))
              {
                index=z0==chor_taxa[i]
                z0[index]=i    
              }
          }
       z0=as.numeric(z0)
       # print(ticks)
       x1=sort(x0)
       position=order(y0,decreasing = TRUE)
       values=y0[position]
       if(length(values)<=DPP)
         {
         index=position[c(1:length(values))]
         }else{       
          # values=sqrt(values)  #This shift the weight a little bit to the low building.
        #Handler of bias plot
        cut0=ceiling(-log10(cutOff/length(values))/2)
        rv=runif(length(values))
        values=values+rv*(values+cut0)
        index=position[which(values>cut0)]
         }        
       x=x0[index]
       y=y0[index]
       z=z0[index]
        #Extract QTN
       QTN=MP_store[which(MP_store[,borrowSlot]==1),]
        #Draw circles with same size and different thikness
       
       themax=ceiling(max(y))
       themax2=ceiling((ceiling(themax/4))*4)

       themin=floor(min(y))
       wd=((y-themin+base)/(themax-themin+base))*size*ratio
       s=size-wd/ratio/2
       plot(y~x,xlab="",ylab="" ,ylim=c(0,themax2),xlim=c(min(x),max(x)),
           cex.axis=4, cex.lab=4, ,col=plot.color[z],axes=FALSE,type = "p",
           pch=mypch,lwd=wd,cex=s+2.5,cex.main=2)
       mtext(side=2,expression(-log[10](italic(p))),line=3.5, cex=2.5)
       if(!simulation)
         {
          abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")}else{
          points(QTN[,2], QTN[,3], pch=21, cex=2.8,lwd=1.5,col="dimgrey")
          points(QTN[,2], QTN[,3], pch=20, cex=1.8,lwd=1.5,col="dimgrey")
         }        
       if(plot.line)
         {
          if(!is.null(nrow(new_xz)))  
            {
             abline(v=as.numeric(new_xz[,4]),col="grey",lty=as.numeric(new_xz[,2]),untf=T,lwd=3)
            }else{
             abline(v=as.numeric(new_xz[1]),col=plot.color[as.numeric(new_xz[3])],lty=as.numeric(new_xz[2]),untf=T,lwd=3)
            }
         }
        #Add a horizontal line for bonferroniCutOff
       # print(ticks)
       abline(h=bonferroniCutOff,lty=1,untf=T,lwd=3,col="forestgreen")
       axis(2, yaxp=c(0,themax2,4),cex.axis=2.3,tick=T,las=1,lwd=2.5)
       if(k==Nenviron)axis(1, at=max.x,cex.axis=2.5,labels=rep("",length(max.x)),tick=T,lwd=2.5)
       if(k==Nenviron)axis(1, at=ticks,cex.axis=2.5,labels=chm.to.analyze,tick=F,line=1)
       mtext(side=4,paste(environ_name[k],sep=""),line=3.2,cex=2)
    }#end of environ_name
       dev.off()
}#end of plot.type

if("w"%in%plot.type)
{
 pdf(paste("GAPIT.Association.Manhattans_Wide",".pdf" ,sep = ""), width = 16,height=8.5)
 par(mfrow=c(Nenviron,1))
 mtext.h=0.5
 size=2
 ratio=5
 for(k in 1:Nenviron)
 { 
  if(k==Nenviron)
        {#par(mfrow=c(Nenviron,1))
          # print(par())
        par(mar = c(2.5,8,0,8))
         # par(pin=c(10,((8-mtext.h)/Nenviron)+mtext.h))

        }else{
            #par(mfrow=c(Nenviron,1))
        par(mar = c(2,8,0.5,8))    
         # par(pin=c(10,(8-mtext.h)/Nenviron))

        }
  environ_result=read.csv(paste("GAPIT.Association.GWAS_Results.",environ_name[k],".csv",sep=""),head=T)
  #print(environ_result[as.numeric(new_xz[,1]),])
  result=environ_result[,1:4]
    result=result[order(result[,3]),]
    result=result[order(result[,2]),]
    result=result[match(as.character(GM[,1]),as.character(result[,1])),]
    rownames(result)=1:nrow(result)
    GI.MP=result[,c(2:4)]
    borrowSlot=4
    GI.MP[,borrowSlot]=0 #Inicial as 0
    GI.MP[,5]=1:(nrow(GI.MP))
    GI.MP[,6]=1:(nrow(GI.MP))
    
    
    GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
    GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
    GI.MP[is.na(GI.MP[,3]),3]=1
    
    #Retain SNPs that have P values between 0 and 1 (not na etc)
    GI.MP <- GI.MP[GI.MP[,3]>0,]
    GI.MP <- GI.MP[GI.MP[,3]<=1,]
    #Remove chr 0 and 99
    GI.MP <- GI.MP[GI.MP[,1]!=0,]
    total_chromo=length(unique(GI.MP[,1]))
    # print(dim(GI.MP))
    if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1
    numMarker=nrow(GI.MP)
    bonferroniCutOff=-log10(cutOff/numMarker)
    GI.MP[,3] <-  -log10(GI.MP[,3])
    GI.MP[,5]=1:numMarker
    y.lim <- ceiling(max(GI.MP[,3]))
    
    chm.to.analyze <- unique(GI.MP[,1])
    # chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
    nchr=length(chm.to.analyze)
    GI.MP[,6]=1:(nrow(GI.MP))
    MP_store=GI.MP
        index_GI=MP_store[,3]>=0
        MP_store <- MP_store[index_GI,]
        ticks=NULL
        lastbase=0
        for (i in chm.to.analyze)
        {
            index=(MP_store[,1]==i)
            ticks <- c(ticks, lastbase+mean(MP_store[index,2]))
            MP_store[index,2]=MP_store[index,2]+lastbase
            lastbase=max(MP_store[index,2])
        }
        
        x0 <- as.numeric(MP_store[,2])
        y0 <- as.numeric(MP_store[,3])
        z0 <- as.character(MP_store[,1])
       # convert chromosome character to number
       chor_taxa=as.character(unique(MP_store[,1]))
       chor_taxa=chor_taxa[order(as.numeric(as.character(chor_taxa)))]
       chr_letter=grep("[A-Z]|[a-z]",chor_taxa)
       if(!setequal(integer(0),chr_letter))
         {     
           z0=as.character(MP_store[,1])
           for(i in 1:(length(chor_taxa)))
              {
                index=z0==chor_taxa[i]
                z0[index]=i    
              }
          }
        z0=as.numeric(z0)
        x1=sort(x0)

        position=order(y0,decreasing = TRUE)
        values=y0[position]
        if(length(values)<=DPP)
        {
         index=position[c(1:length(values))]
            }else{      
        # values=sqrt(values)  #This shift the weight a little bit to the low building.
        #Handler of bias plot
        cut0=ceiling(-log10(cutOff/length(values))/2)
        rv=runif(length(values))
        values=values+rv*(values+cut0)
      
        index=position[which(values>cut0)]
        }     
        x=x0[index]
        y=y0[index]
        z=z0[index]
        # print(length(x))
        #Extract QTN
        #if(!is.null(seqQTN))MP_store[seqQTN,borrowSlot]=1
        #if(!is.null(interQTN))MP_store[interQTN,borrowSlot]=2
        QTN=MP_store[which(MP_store[,borrowSlot]==1),]
        #Draw circles with same size and different thikness
        themax=ceiling(max(y))
        themax2=ceiling((ceiling(themax/4))*4)
        themin=floor(min(y))
        # size=5
        wd=((y-themin+base)/(themax-themin+base))*size*ratio
        # wd=0.5
        s=size-wd/ratio/2
        mypch=1
        bamboo=4
        # plot(y~x,xlab="",ylab="" ,ylim=c(0,themax),
        #     cex.axis=4, cex.lab=4, ,col=plot.color[z],axes=FALSE,type = "p",pch=mypch,lwd=0.5,cex=0.7,cex.main=2)
        plot(y~x,xlab="",ylab="" ,ylim=c(0,themax2),xlim=c(min(x),max(x)),
           cex.axis=4, cex.lab=4, ,col=plot.color[z],axes=FALSE,type = "p",
           pch=mypch,lwd=wd,cex=s,cex.main=2)
        mtext(side=2,expression(-log[10](italic(p))),line=3, cex=1)
        if(plot.line)
        {
          if(!is.null(nrow(new_xz)))  {abline(v=as.numeric(new_xz[,4]),col="grey",lty=as.numeric(new_xz[,2]),untf=T,lwd=2)
             }else{abline(v=as.numeric(new_xz[1]),col=plot.color[as.numeric(new_xz[3])],lty=as.numeric(new_xz[2]),untf=T,lwd=2)
             }
        }
        if(!simulation){abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")}else{
          # print("$$$")
          points(QTN[,2], QTN[,3], type="p",pch=21, cex=2.8,lwd=1.5,col="dimgrey")
          points(QTN[,2], QTN[,3], type="p",pch=20, cex=1.5,lwd=1.5,col="dimgrey")
          }
        #Add a horizontal line for bonferroniCutOff
        abline(h=bonferroniCutOff,lty=1,untf=T,lwd=1,col="forestgreen")
        axis(2, yaxp=c(0,themax2,bamboo),cex.axis=1.5,las=1,tick=F)
        if(k==Nenviron)axis(1, at=ticks,cex.axis=1.5,line=0.001,labels=chm.to.analyze,tick=F)
        mtext(side=4,paste(environ_name[k],sep=""),line=2,cex=1,base_family="Arial")
 box()
 }#end of environ_name
 dev.off()
}#end of plot.type

if("s"%in%plot.type)
{
    # wd=((y-themin+base)/(themax-themin+base))*size*ratio
 wd=2
 if(is.null(outpch))
 {
   allpch0=c(0,1,2,5,6)
 }else{
   allpch0=outpch
 }
 if(is.null(inpch))
 {
   add.pch=c("+","*","-","#","<",">","^","$","=","|","?",as.character(1:9),letters[1:26],LETTERS[1:26]) 
 }else{
   add.pch=inpch
 }
 n.vals=ceiling(Nenviron/length(allpch0))-1
 s=size-wd/ratio/2
 DPP=500
 


 pdf(paste("GAPIT.Association.Manhattans_Symphysic",".pdf" ,sep = ""), width = 30,height=18)
 par(mfrow=c(1,1))
 par(mar = c(5,8,5,1))
 themax.y02=ceiling((ceiling(themax.y0/4))*4)

 plot(1~1,col="white",xlab="",ylab="" ,ylim=c(0,themax.y02),xlim=c(min(store.x,na.rm=TRUE),max(store.x,na.rm=TRUE)),yaxp=c(0,themax.y02,4),
    cex.axis=4, cex.lab=4,axes=FALSE,
    pch=1,cex.main=4)
 # print(ticks)   
        #Add a horizontal line for bonferroniCutOff
 axis(1, at=max.x,cex.axis=2,labels=rep("",length(max.x)),tick=T,lwd=2.5)
 axis(1, at=ticks,cex.axis=2,labels=chm.to.analyze,tick=F,line=1)
 axis(2, yaxp=c(0,themax.y02,4),cex.axis=2,tick=T,las=1,lwd=2.5)
 if(!is.null(cutOff))abline(h=bonferroniCutOff,lty=1,untf=T,lwd=3,col="forestgreen")
 if(plot.line)
    {
    if(!is.null(nrow(new_xz)))  
        {
        abline(v=as.numeric(new_xz[,4]),col="grey",lty=as.numeric(new_xz[,2]),untf=T,lwd=3)
        }else{
        abline(v=as.numeric(new_xz[1]),col=plot.color[as.numeric(new_xz[3])],lty=as.numeric(new_xz[2]),untf=T,lwd=3)
        }
    }
 mtext(side=2,expression(-log[10](italic(p))),line=4, cex=2.5)
 # legend("top",legend=paste(environ_name,sep=""),ncol=length(environ_name),
 #       col="black",pch=allpch[1:Nenviron],lty=0,lwd=1,cex=2,
 #       bty = "o", bg = "white",box.col="white")
 # step.vals=0
 # if(1>2){
 for(k in 1:Nenviron)
  { 
    step.vals=ceiling(k/length(allpch0))-1

    environ_result=read.csv(paste("GAPIT.Association.GWAS_Results.",environ_name[k],".csv",sep=""),head=T)
    result=environ_result[,1:4]
    result=result[order(result[,3]),]
    result=result[order(result[,2]),]
    result=result[match(as.character(GM[,1]),as.character(result[,1])),]
    rownames(result)=1:nrow(result)
    GI.MP=result[,c(2:4)]
    borrowSlot=4
    GI.MP[,borrowSlot]=0 #Inicial as 0
    GI.MP[,5]=1:(nrow(GI.MP))
    GI.MP[,6]=1:(nrow(GI.MP)) 
    GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
    GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
    GI.MP[is.na(GI.MP[,3]),3]=1
    
    #Retain SNPs that have P values between 0 and 1 (not na etc)
    GI.MP <- GI.MP[GI.MP[,3]>0,]
    GI.MP <- GI.MP[GI.MP[,3]<=1,]
    #Remove chr 0 and 99
    GI.MP <- GI.MP[GI.MP[,1]!=0,]
    total_chromo=length(unique(GI.MP[,1]))
    # print(dim(GI.MP))
    if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1
    numMarker=nrow(GI.MP)
    bonferroniCutOff=-log10(cutOff/numMarker)
    GI.MP[,3] <-  -log10(GI.MP[,3])
    GI.MP[,5]=1:numMarker
    y.lim <- ceiling(max(GI.MP[,3]))
    
    chm.to.analyze <- unique(GI.MP[,1])
    # print(chm.to.analyze)
    # chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
    nchr=length(chm.to.analyze)
    # print(chm.to.analyze)
    GI.MP[,6]=1:(nrow(GI.MP))
    MP_store=GI.MP
    index_GI=MP_store[,3]>=0
    MP_store <- MP_store[index_GI,]
    ticks=NULL
    lastbase=0
    for (i in chm.to.analyze)
        {
            index=(MP_store[,1]==i)
            ticks <- c(ticks, lastbase+mean(MP_store[index,2]))
            MP_store[index,2]=MP_store[index,2]+lastbase
            lastbase=max(MP_store[index,2])
        }
        
    x0 <- as.numeric(MP_store[,2])
    y0 <- as.numeric(MP_store[,3])
    z0 <- as.character(MP_store[,1])
       # convert chromosome character to number
       chor_taxa=as.character(unique(MP_store[,1]))
       chor_taxa=chor_taxa[order(as.numeric(as.character(chor_taxa)))]
       chr_letter=grep("[A-Z]|[a-z]",chor_taxa)
       if(!setequal(integer(0),chr_letter))
         {     
           z0=as.character(MP_store[,1])
           for(i in 1:(length(chor_taxa)))
              {
                index=z0==chor_taxa[i]
                z0[index]=i    
              }
          }
       z0=as.numeric(z0)
       max.x=NULL
    for (i in chm.to.analyze)
        {
            index=(MP_store[,1]==i)
            max.x=c(max.x,max(x0[index]))
        }
    max.x=c(min(x0),max.x)
    x1=sort(x0)

    position=order(y0,decreasing = TRUE)
    values=y0[position]
    if(length(values)<=DPP)
        {
         index=position[c(1:length(values))]
        }else{       
          # values=sqrt(values)  #This shift the weight a little bit to the low building.
        #Handler of bias plot
        cut0=ceiling(-log10(cutOff/length(values))/2)
        rv=runif(length(values))
        values=values+rv*(values+cut0)

        index=position[which(values>cut0)]
        }        
    x=x0[index]
    y=y0[index]
    z=z0[index]
        # print(length(x))

        #Extract QTN
        #if(!is.null(seqQTN))MP_store[seqQTN,borrowSlot]=1
        #if(!is.null(interQTN))MP_store[interQTN,borrowSlot]=2
    QTN=MP_store[which(MP_store[,borrowSlot]==1),]
        #Draw circles with same size and different thikness
    themax=ceiling(max(y))
    # themax.y02=ceiling((ceiling(themax.y0/4)+1)*4)
    # print(themax.y02)
    themin=floor(min(y))
    mypch=allpch0[k-step.vals*length(allpch0)]
    
   # if(k!=1) par(new=T)
    
    par(new=T)
    plot(y~x,xlab="",ylab="" ,ylim=c(0,themax.y02),xlim=c(min(x),max(x)),yaxp=c(0,themax.y02,4),
    cex.axis=4, cex.lab=4,col=plot.color[z],axes=FALSE,
    pch=mypch,lwd=1,cex=s+2.5,cex.main=4)
    if(step.vals!=0)
    {
      points(y~x,pch=add.pch[step.vals],col=plot.color[z],cex=s+0.5,cex.main=4)
    }
    if(!simulation)
       {
        abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")
        }else{
        points(QTN[,2], QTN[,3], pch=20, cex=2,lwd=2.5,col="dimgrey")
       }        
    
 }#end of environ_name
 dev.off()
 # }
 ## Plot legend
 nchar.traits=1.5
 # environ_name=paste(environ_name,"1234",sep="")
 nchar0=max(nchar(environ_name))
 print(nchar0)
 if(Nenviron>5)
 {
  yourpch=c(rep(allpch0,n.vals),allpch0[1:(Nenviron-length(allpch0)*n.vals)])
 }else{
  yourpch=allpch0[1:Nenviron]
 }
  yourpch2=NULL
  for(pp in 1:n.vals)
  {
    yourpch2=c(yourpch2,rep(add.pch[pp],length(allpch0)))
  }
  # yourpch2=c(rep(allpch0,n.vals),allpch0[Nenviron-length(allpch0)*n.vals])
  if(Nenviron>5){
  yourpch2=yourpch2[1:(Nenviron-length(allpch0))]
  yourpch2=c(rep(NA,length(allpch0)),yourpch2)
  }
 
 max.row=25
 max.pch=ifelse(Nenviron<max.row,Nenviron,max.row)
 n.col.pch=ceiling(Nenviron/max.row)
 ratio.cex=ceiling(Nenviron/5)
  if(ratio.cex>5)ratio.cex=5
  cex.Ne=ratio.cex*(0.05*ratio.cex+0.35) #the size of cex
  if(ratio.cex<3)cex.Ne=1
  c.t.d=c(0.5,1,1,1,1.5)[ratio.cex] # the different between size of cex and text
  cex.di=0.3*ratio.cex # the different size between cex and signal
  text.di=c(.02,0.01,0.01,0.02,0.02)[ratio.cex] #the different distance between cex and text

  high.Ne=2*ratio.cex # the total highth of figure
  cex.betw=c(.38,.5,.7,0.9,1)[ratio.cex] # distance between cexes
  x.di=c(1.12,1.09,0.5,0.8,1.3)[ratio.cex]/2 # distance between markers in x axis
  # print(nchar0)
  # x.di0=c(0)
  if(n.col.pch>1){
    text.di=(0.01*nchar0)/3+0.02
    # x.di=0.52*n.col.pch
    x.di=(0.1*(ceiling(nchar0/5)-1)+(n.col.pch-1)*0.1)*ceiling(nchar0/5)#*n.col.pch
  }
  # print(x.di)
 # if(Nenviron>5){
 #  cex.Ne=3
 #  cex.di=1.5
 #  text.di=.02
 #  high.Ne=10
 #  cex.betw=0.9
 #  }else{
 #  cex.Ne=1
 #  cex.di=0.3
 #  high.Ne=Nenviron/2 
 #  text.di=.02
 #  cex.betw=0.9
 #  }


 write.csv(environ_name,"GAPIT.Association.Manhattans_Symphysic_Traitsnames.csv",quote=FALSE)
 pdf(paste("GAPIT.Association.Manhattans_Symphysic_Legend",".pdf" ,sep = ""), width = 4+(x.di*(n.col.pch+1)),height=high.Ne)
 par(mfrow=c(1,1))
 par(mar = c(cex.Ne+1,2,cex.Ne+1,2))
 # print(length(yourpch))
 # print(length(yourpch2))
 plot(0,0,xlab="",ylab="" ,axes=FALSE,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),col="white")
 for(kk in 1:n.col.pch)
 {
 par(new=T)

 if(kk==n.col.pch)
 {
  if(n.col.pch==1)
  {
  # print(kk)
  max.pch2=Nenviron-(n.col.pch-1)*max.row
  
  plot(rep(0,max.pch2),(max.pch:(max.pch-max.pch2+1))*cex.betw,xlab="",ylab="" ,axes=FALSE,col="black",
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne,
  pch=yourpch[((kk-1)*max.row+1):Nenviron])
  if(Nenviron>5) points(rep(0,max.pch2),((max.pch):(max.pch-max.pch2+1))*cex.betw,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne-cex.di,
  pch=yourpch2[((kk-1)*max.row+1):Nenviron])
  text(rep((0+text.di),max.pch2),(max.pch:(max.pch-max.pch2+1))*cex.betw,labels=environ_name[((kk-1)*max.row+1):Nenviron],pos=4,cex=cex.Ne-c.t.d)
  }else{
  # print(kk)
  max.pch2=Nenviron-(n.col.pch-1)*max.row
  
  plot(rep((kk-1)*x.di,max.pch2),(max.pch:(max.pch-max.pch2+1))*cex.betw,xlab="",ylab="" ,axes=FALSE,col="black",
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne,
  pch=yourpch[((kk-1)*max.row+1):Nenviron])
  if(Nenviron>5) points(rep((kk-1)*x.di,max.pch2),((max.pch):(max.pch-max.pch2+1))*cex.betw,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne-cex.di,
  pch=yourpch2[((kk-1)*max.row+1):Nenviron])
  text(rep(((kk-1)*x.di+text.di),max.pch2),(max.pch:(max.pch-max.pch2+1))*cex.betw,labels=environ_name[((kk-1)*max.row+1):Nenviron],pos=4,cex=cex.Ne-c.t.d)
    
  }
 }else{
  if(kk==1)
  {
  print(kk)
  plot(rep(0,max.pch),(max.pch:1)*cex.betw,xlab="",ylab="" ,axes=FALSE,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne,
  pch=yourpch[((kk-1)*max.row+1):((kk-1)*max.row+max.row)])
  if(Nenviron>5) points(rep(0,max.pch),((max.pch):1)*cex.betw,
  xlim=c(1,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne-cex.di,
  pch=yourpch2[((kk-1)*max.row+1):((kk-1)*max.row+max.row)])
  text(rep((0+text.di),max.pch),(max.pch:1)*cex.betw,labels=environ_name[((kk-1)*max.row+1):((kk-1)*max.row+max.row)],pos=4,cex=cex.Ne-c.t.d)
  }else{
  print(kk)
  plot(rep((kk-1)*x.di,max.pch),(max.pch:1)*cex.betw,xlab="",ylab="" ,axes=FALSE,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne,
  pch=yourpch[((kk-1)*max.row+1):((kk-1)*max.row+max.row)])
  if(Nenviron>5) points(rep((kk-1)*x.di,max.pch),((max.pch):1)*cex.betw,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne-cex.di,
  pch=yourpch2[((kk-1)*max.row+1):((kk-1)*max.row+max.row)])
  text(rep(((kk-1)*x.di+text.di),max.pch),(max.pch:1)*cex.betw,labels=environ_name[((kk-1)*max.row+1):((kk-1)*max.row+max.row)],pos=4,cex=cex.Ne-c.t.d)
   
  }

 }
}#end of plot.type
 dev.off()

}
print("GAPIT.Association.Manhattans has done !!!")
return(list(multip_mapP=result0,xz=new_xz))
} #end of GAPIT.Manhattan
#=============================================================================================

########## These three functions come from MVP package, Jiabo did some modifications
########## Following Apache License, we thank MVP developper to build these functions.
    ########## 1 creat P value scale in addtitional chromsome
    ########## 2 set col  same as GAPIT
    ########## 3 
circle.plot <- function(myr,type="l",x=NULL,lty=1,lwd=1,col="black",add=TRUE,n.point=1000)
    {
        graphics::curve(sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=add)
        graphics::curve(-sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=TRUE)
    }
Densitplot <- function(
        map,
        col=c("darkblue", "white", "red"),
        main="SNP Density",
        bin=1e6,
        band=3,
        width=5,
        legend.len=10,
        legend.max=NULL,
        legend.pt.cex=3,
        legend.cex=1,
        legend.y.intersp=1,
        legend.x.intersp=1,
        plot=TRUE
    )
    {   #print(head(map))
        map <- as.matrix(map)
        map <- map[!is.na(map[, 2]), ]
        map <- map[!is.na(map[, 3]), ]
        map <- map[map[, 2] != 0, ]
        #map <- map[map[, 3] != 0, ]
        options(warn = -1)
        max.chr <- max(as.numeric(map[, 2]), na.rm=TRUE)
        if(is.infinite(max.chr))    max.chr <- 0
        map.xy.index <- which(!as.numeric(map[, 2]) %in% c(0 : max.chr))
        if(length(map.xy.index) != 0){
            chr.xy <- unique(map[map.xy.index, 2])
            for(i in 1:length(chr.xy)){
                map[map[, 2] == chr.xy[i], 2] <- max.chr + i
            }
        }
        map <- map[order(as.numeric(map[, 2]), as.numeric(map[, 3])), ]
        chr <- as.numeric(map[, 2])
        pos <- as.numeric(map[, 3])
        chr.num <- unique(chr)
        #print(chr.num)
        chorm.maxlen <- max(pos)
        if(plot)    plot(NULL, xlim=c(0, chorm.maxlen + chorm.maxlen/10), ylim=c(0, length(chr.num) * band + band), main=main,axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
        pos.x <- list()
        col.index <- list()
        maxbin.num <- NULL
        #print(chr.num)
        for(i in 1 : length(chr.num)){
            pos.x[[i]] <- pos[which(chr == chr.num[i])]
            cut.len <- ceiling((max(pos.x[[i]]) - min(pos.x[[i]])) / bin)
            if(cut.len <= 1){
                col.index[[i]] = 1
            }else{
                cut.r <- cut(pos.x[[i]], cut.len, labels=FALSE)
                eachbin.num <- table(cut.r)
                #print(eachbin.num)

                maxbin.num <- c(maxbin.num, max(eachbin.num))
                col.index[[i]] <- rep(eachbin.num, eachbin.num)
            }
        }

        Maxbin.num <- max(maxbin.num)
        maxbin.num <- Maxbin.num
        if(!is.null(legend.max)){
            maxbin.num <- legend.max
        }
        #print(col)
        #print(maxbin.num)
        col = grDevices::colorRampPalette(col)(maxbin.num)
        col.seg=NULL
        for(i in 1 : length(chr.num)){
            if(plot)    graphics::polygon(c(0, 0, max(pos.x[[i]]), max(pos.x[[i]])), 
                c(-width/5 - band * (i - length(chr.num) - 1), width/5 - band * (i - length(chr.num) - 1), 
                width/5 - band * (i - length(chr.num) - 1), -width/5 - band * (i - length(chr.num) - 1)), col="grey", border="grey")
            if(!is.null(legend.max)){
                if(legend.max < Maxbin.num){
                    col.index[[i]][col.index[[i]] > legend.max] <- legend.max
                }
            }
            col.seg <- c(col.seg, col[round(col.index[[i]] * length(col) / maxbin.num)])
            if(plot)    graphics::segments(pos.x[[i]], -width/5 - band * (i - length(chr.num) - 1), pos.x[[i]], width/5 - band * (i - length(chr.num) - 1), 
            col=col[round(col.index[[i]] * length(col) / maxbin.num)], lwd=1)
        }
        if(length(map.xy.index) != 0){
            for(i in 1:length(chr.xy)){
                chr.num[chr.num == max.chr + i] <- chr.xy[i]
            }
        }
        chr.num <- rev(chr.num)
        if(plot)    graphics::mtext(at=seq(band, length(chr.num) * band, band),text=paste("Chr", chr.num, sep=""), side=2, las=2, font=1, cex=0.6, line=0.2)
        if(plot)    graphics::axis(3, at=seq(0, chorm.maxlen, length=10), labels=c(NA, paste(round((seq(0, chorm.maxlen, length=10))[-1] / 1e6, 0), "Mb", sep="")),
            font=1, cex.axis=0.8, tck=0.01, lwd=2, padj=1.2)
        # image(c(chorm.maxlen-chorm.maxlen * legend.width / 20 , chorm.maxlen), 
        # round(seq(band - width/5, (length(chr.num) * band + band) * legend.height / 2 , length=maxbin.num+1), 2), 
        # t(matrix(0 : maxbin.num)), col=c("white", rev(heat.colors(maxbin.num))), add=TRUE)
        legend.y <- round(seq(0, maxbin.num, length=legend.len))
        len <- legend.y[2]
        legend.y <- seq(0, maxbin.num, len)
        if(!is.null(legend.max)){
            if(legend.max < Maxbin.num){
                if(!maxbin.num %in% legend.y){
                    legend.y <- c(legend.y, paste(">=", maxbin.num, sep=""))
                    legend.y.col <- c(legend.y[c(-1, -length(legend.y))], maxbin.num)
                }else{
                    legend.y[length(legend.y)] <- paste(">=", maxbin.num, sep="")
                    legend.y.col <- c(legend.y[c(-1, -length(legend.y))], maxbin.num)
                }
            }else{
                if(!maxbin.num %in% legend.y){
                    legend.y <- c(legend.y, maxbin.num)
                }
                legend.y.col <- c(legend.y[-1])
            }
        }else{
            if(!maxbin.num %in% legend.y){
                legend.y <- c(legend.y, paste(">", max(legend.y), sep=""))
                legend.y.col <- c(legend.y[c(-1, -length(legend.y))], maxbin.num)
            }else{
                legend.y.col <- c(legend.y[-1])
            }
        }
        legend.y.col <- as.numeric(legend.y.col)
        legend.col <- c("grey", col[round(legend.y.col * length(col) / maxbin.num)])
        if(plot)    graphics::legend(x=(chorm.maxlen + chorm.maxlen/100), y=( -width/2.5 - band * (length(chr.num) - length(chr.num) - 1)), title="", legend=legend.y, pch=15, pt.cex = legend.pt.cex, col=legend.col,
            cex=legend.cex, bty="n", y.intersp=legend.y.intersp, x.intersp=legend.x.intersp, yjust=0, xjust=0, xpd=TRUE)
        if(!plot)   return(list(den.col=col.seg, legend.col=legend.col, legend.y=legend.y))
    }

GAPIT.Circle.Manhattan.Plot <- function(
    Pmap,
    col=c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
    #col=c("darkgreen", "darkblue", "darkyellow", "darkred"),
    
    bin.size=1e6,
    bin.max=NULL,
    pch=19,
    band=1,
    cir.band=0.5,
    H=1.5,
    ylim=NULL,
    cex.axis=1,
    plot.type="c",
    multracks=TRUE,
    cex=c(0.5,0.8,1),
    r=0.3,
    xlab="Chromosome",
    ylab=expression(-log[10](italic(p))),
    xaxs="i",
    yaxs="r",
    outward=TRUE,
    threshold = 0.01, 
    threshold.col="red",
    threshold.lwd=1,
    threshold.lty=2,
    amplify= TRUE,     # is that available for remark signal pch col
    chr.labels=NULL,
    signal.cex = 2,
    signal.pch = 8,
    signal.col="red",
    signal.line=NULL,
    cir.chr=TRUE,
    cir.chr.h=1.3,
    chr.den.col=c("darkgreen", "yellow", "red"),
    #chr.den.col=c(126,177,153),
    cir.legend=TRUE,
    cir.legend.cex=0.8,
    cir.legend.col="grey45",
    LOG10=TRUE,
    box=FALSE,
    conf.int.col="grey",
    file.output=TRUE,
    file="pdf",
    dpi=300,
    xz=NULL,
    memo=""
)
{       #print("Starting Circular-Manhattan plot!",quote=F)
    taxa=colnames(Pmap)[-c(1:3)]
    if(!is.null(memo) && memo != "")    memo <- paste("_", memo, sep="")
    if(length(taxa) == 0)   taxa <- "Index"
    taxa <- paste(taxa, memo, sep="")
    col=rep(c( '#FF6A6A',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5'),ceiling(length(taxa)/5))
    legend.bit=round(nrow(Pmap)/30)

    numeric.chr <- as.numeric(Pmap[, 1])
    options(warn = 0)
    max.chr <- max(numeric.chr, na.rm=TRUE)
    aa=Pmap[1:legend.bit,]
    aa[,2]=max.chr+1
    #print(aa[,3])
    aa[,3]=sample(1:10^7.5,legend.bit)
    aa[,-c(1:3)]=0
    Pmap=rbind(Pmap,aa)
    #print(unique(Pmap[,2]))
    #SNP-Density plot
    if("d" %in% plot.type){
        print("SNP_Density Plotting...")
        if(file.output){
            if(file=="jpg") grDevices::jpeg(paste("SNP_Density.",paste(taxa,collapse="."),".jpg",sep=""), width = 9*dpi,height=7*dpi,res=dpi,quality = 100)
            if(file=="pdf") grDevices::pdf(paste("GAPIT.Association.SNP_Density", taxa,".pdf" ,sep=""), width = 9,height=7)
            if(file=="tiff")    grDevices::tiff(paste("SNP_Density.",paste(taxa,collapse="."),".tiff",sep=""), width = 9*dpi,height=7*dpi,res=dpi)
            graphics::par(xpd=TRUE)
        }else{
            if(is.null(grDevices::dev.list()))  grDevices::dev.new(width = 9,height=7)
            graphics::par(xpd=TRUE)
        }

        Densitplot(map=Pmap[,c(1:3)], col=col, bin=bin.size, legend.max=bin.max, main=paste("The number of SNPs within ", bin.size/1e6, "Mb window size", sep=""))
        if(file.output) grDevices::dev.off()
    }

    if(length(plot.type) !=1 | (!"d" %in% plot.type)){
    
        #order Pmap by the name of SNP
        #Pmap=Pmap[order(Pmap[,1]),]
        Pmap <- as.matrix(Pmap)

        #delete the column of SNPs names
        Pmap <- Pmap[,-1]
        Pmap[is.na(Pmap)]=1
        #print(dim(Pmap))

        #scale and adjust the parameters
        cir.chr.h <- cir.chr.h/5
        cir.band <- cir.band/5
        threshold=threshold/nrow(Pmap)
        if(!is.null(threshold)){
            threshold.col <- rep(threshold.col,length(threshold))
            threshold.lwd <- rep(threshold.lwd,length(threshold))
            threshold.lty <- rep(threshold.lty,length(threshold))
            signal.col <- rep(signal.col,length(threshold))
            signal.pch <- rep(signal.pch,length(threshold))
            signal.cex <- rep(signal.cex,length(threshold))
        }
        if(length(cex)!=3) cex <- rep(cex,3)
        if(!is.null(ylim)){
            if(length(ylim)==1) ylim <- c(0,ylim)
        }
        
        if(is.null(conf.int.col))   conf.int.col <- NA
        if(is.na(conf.int.col)){
            conf.int=FALSE
        }else{
            conf.int=TRUE
        }

        #get the number of traits
        R=ncol(Pmap)-2

        #replace the non-euchromosome
        options(warn = -1)
        numeric.chr <- as.numeric(Pmap[, 1])
        options(warn = 0)
        max.chr <- max(numeric.chr, na.rm=TRUE)
        if(is.infinite(max.chr))    max.chr <- 0
        map.xy.index <- which(!numeric.chr %in% c(0:max.chr))
        if(length(map.xy.index) != 0){
            chr.xy <- unique(Pmap[map.xy.index, 1])
            for(i in 1:length(chr.xy)){
                Pmap[Pmap[, 1] == chr.xy[i], 1] <- max.chr + i
            }
        }

        Pmap <- matrix(as.numeric(Pmap), nrow(Pmap))

        #order the GWAS results by chromosome and position
        Pmap <- Pmap[order(Pmap[, 1], Pmap[,2]), ]

        #get the index of chromosome
        chr <- unique(Pmap[,1])
        chr.ori <- chr
        if(length(map.xy.index) != 0){
            for(i in 1:length(chr.xy)){
                chr.ori[chr.ori == max.chr + i] <- chr.xy[i]
            }
        }

        pvalueT <- as.matrix(Pmap[,-c(1:2)])
        #print(dim(pvalueT))
        pvalue.pos <- Pmap[, 2]
        p0.index <- Pmap[, 1] == 0
        if(sum(p0.index) != 0){
            pvalue.pos[p0.index] <- 1:sum(p0.index)
        }
        pvalue.pos.list <- tapply(pvalue.pos, Pmap[, 1], list)
        
        #scale the space parameter between chromosomes
        if(!missing(band)){
            band <- floor(band*(sum(sapply(pvalue.pos.list, max))/100))
        }else{
            band <- floor((sum(sapply(pvalue.pos.list, max))/100))
        }
        if(band==0) band=1
        
        if(LOG10){
            pvalueT[pvalueT <= 0] <- 1
            pvalueT[pvalueT > 1] <- 1
        }

        #set the colors for the plot
        #palette(heat.colors(1024)) #(heatmap)
        #T=floor(1024/max(pvalue))
        #plot(pvalue,pch=19,cex=0.6,col=(1024-floor(pvalue*T)))
        
        #print(col)
        if(is.vector(col)){
            col <- matrix(col,R,length(col),byrow=TRUE)
        }
        if(is.matrix(col)){
            #try to transform the colors into matrix for all traits
            col <- matrix(as.vector(t(col)),R,dim(col)[2],byrow=TRUE)
        }

        Num <- as.numeric(table(Pmap[,1]))
        Nchr <- length(Num)
        N <- NULL
        #print(Nchr)
        #set the colors for each traits
        for(i in 1:R){
            colx <- col[i,]
            colx <- colx[!is.na(colx)]
            N[i] <- ceiling(Nchr/length(colx))
        }
        
        #insert the space into chromosomes and return the midpoint of each chromosome
        ticks <- NULL
        pvalue.posN <- NULL
        #pvalue <- pvalueT[,j]
        for(i in 0:(Nchr-1)){
            if (i==0){
                #pvalue <- append(pvalue,rep(Inf,band),after=0)
                pvalue.posN <- pvalue.pos.list[[i+1]] + band
                ticks[i+1] <- max(pvalue.posN)-floor(max(pvalue.pos.list[[i+1]])/2)
            }else{
                #pvalue <- append(pvalue,rep(Inf,band),after=sum(Num[1:i])+i*band)
                pvalue.posN <- c(pvalue.posN, max(pvalue.posN) + band + pvalue.pos.list[[i+1]])
                ticks[i+1] <- max(pvalue.posN)-floor(max(pvalue.pos.list[[i+1]])/2)
            }
        }
        pvalue.posN.list <- tapply(pvalue.posN, Pmap[, 1], list)
        #NewP[[j]] <- pvalue
        
        #merge the pvalues of traits by column
        if(LOG10){
            logpvalueT <- -log10(pvalueT)
        }else{
            pvalueT <- abs(pvalueT)
            logpvalueT <- pvalueT
        }

        add <- list()
        for(i in 1:R){
            colx <- col[i,]
            colx <- colx[!is.na(colx)]
            add[[i]] <- c(Num,rep(0,N[i]*length(colx)-Nchr))
        }

        TotalN <- max(pvalue.posN)

        if(length(chr.den.col) > 1){
            cir.density=TRUE
            den.fold <- 20
            density.list <- Densitplot(map=Pmap[,c(1,1,2)], col=chr.den.col, plot=FALSE, bin=bin.size, legend.max=bin.max)
            #list(den.col=col.seg, legend.col=legend.col, legend.y=legend.y)
        }else{
            cir.density=FALSE
        }


        #print(dim(pvalueT))

        if(is.null(xz)){
        signal.line.index <- NULL
        if(!is.null(threshold)){
            if(!is.null(signal.line)){
                for(l in 1:R){
                    if(LOG10){
                        signal.line.index <- c(signal.line.index,which(pvalueT[,l] < min(threshold)))
                    }else{
                        signal.line.index <- c(signal.line.index,which(pvalueT[,l] > max(threshold)))
                    }
                }
                signal.line.index <- unique(signal.line.index)
            }
        }
        signal.lty=rep(2,length(signal.line.index))
        }else{
        signal.line.index=as.numeric(as.vector(xz[,1]))
        signal.lty=as.numeric(as.vector(xz[,2]))
        }#end is.null(xz)
        
        signal.line.index <- pvalue.posN[signal.line.index]
    }
        


    if("c" %in% plot.type)
    {
        if(file.output){
            if(file=="jpg") grDevices::jpeg(paste("GAPIT.Manhattan.Multiple.Plot.circular.jpg",sep=""), width = 8*dpi,height=8*dpi,res=dpi,quality = 100)
            if(file=="pdf") grDevices::pdf(paste("GAPIT.Association.Manhattans_Circular.pdf" ,sep=""), width = 10,height=10)
            if(file=="tiff")    grDevices::tiff(paste("GAPIT.Manhattan.Multiple.Plot.circular.tiff",sep=""), width = 8*dpi,height=8*dpi,res=dpi)
        }
        if(!file.output){
            if(!is.null(grDevices::dev.list())) grDevices::dev.new(width=8, height=8)
            graphics::par(pty="s", xpd=TRUE, mar=c(1,1,1,1))
        }
        graphics::par(pty="s", xpd=TRUE, mar=c(1,1,1,1))
        RR <- r+H*R+cir.band*R
        if(cir.density){
            plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.1*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.1*(RR+4*cir.chr.h)),axes=FALSE,xlab="",ylab="")
        }else{
            plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),axes=FALSE,xlab="",ylab="")
        }
        if(!is.null(signal.line)){
            if(!is.null(signal.line.index)){
                X1chr <- (RR)*sin(2*pi*(signal.line.index-round(band/2))/TotalN)
                Y1chr <- (RR)*cos(2*pi*(signal.line.index-round(band/2))/TotalN)
                X2chr <- (r)*sin(2*pi*(signal.line.index-round(band/2))/TotalN)
                Y2chr <- (r)*cos(2*pi*(signal.line.index-round(band/2))/TotalN)
                #print(signal.line)

                #print(dim(pvalueT))
                #print(head(pvalueT))
                #print(dim(xz))
                #print(xz)
                #print(head(pvalue.posN))
                graphics::segments(X1chr,Y1chr,X2chr,Y2chr,lty=signal.lty,lwd=signal.line,col="grey")
            }
        }
        for(i in 1:R){
        
            #get the colors for each trait
            colx <- col[i,]
            colx <- colx[!is.na(colx)]
            
            #debug
            #print(colx)
            
            #print(paste("Circular_Manhattan Plotting ",taxa[i],"...",sep=""))
            pvalue <- pvalueT[,i]
            logpvalue <- logpvalueT[,i]
            if(is.null(ylim)){
                if(LOG10){
                    Max <- ceiling(-log10(min(pvalue[pvalue!=0])))
                }else{
                    Max <- ceiling(max(pvalue[pvalue!=Inf]))
                    if(Max<=1)
                    Max <- max(pvalue[pvalue!=Inf])
                }
            }else{
                Max <- ylim[2]
            }
            Cpvalue <- (H*logpvalue/Max)

            if(outward==TRUE){
                if(cir.chr==TRUE){
                    
                    #plot the boundary which represents the chromosomes
                    polygon.num <- 1000
                    #print(length(chr))
                    for(k in 1:length(chr)){
                        if(k==1){
                            polygon.index <- seq(round(band/2)+1,-round(band/2)+max(pvalue.posN.list[[1]]), length=polygon.num)
                            #change the axis from right angle into circle format
                            X1chr=(RR)*sin(2*pi*(polygon.index)/TotalN)
                            Y1chr=(RR)*cos(2*pi*(polygon.index)/TotalN)
                            X2chr=(RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
                            Y2chr=(RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)

                            #print(length(X1chr))
                            if(is.null(chr.den.col)){
                                graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k]) 
                            }else{
                                if(cir.density){
                                        graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
                                }else{
                                        graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
                                }
                            }
                        }else{
                            polygon.index <- seq(1+round(band/2)+max(pvalue.posN.list[[k-1]]),-round(band/2)+max(pvalue.posN.list[[k]]), length=polygon.num)
                            X1chr=(RR)*sin(2*pi*(polygon.index)/TotalN)
                            Y1chr=(RR)*cos(2*pi*(polygon.index)/TotalN)
                            X2chr=(RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
                            Y2chr=(RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
                            if(is.null(chr.den.col)){
                                graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])
                            }else{
                                if(cir.density){
                                        graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
                                }else{
                                        graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
                                }
                            }       
                        }
                    }
                    
                    if(cir.density){

                        graphics::segments(
                            (RR)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
                            (RR)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
                            (RR+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
                            (RR+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
                            col=density.list$den.col, lwd=0.1
                        )
                        graphics::legend(
                            x=RR+4*cir.chr.h,
                            y=(RR+4*cir.chr.h)/2,
                            horiz=F,
                            title="Density", legend=density.list$legend.y, pch=15, pt.cex = 3, col=density.list$legend.col,
                            cex=1, bty="n",
                            y.intersp=1,
                            x.intersp=1,
                            yjust=0.5, xjust=0, xpd=TRUE
                        )
                        
                    }
                    
                    # XLine=(RR+cir.chr.h)*sin(2*pi*(1:TotalN)/TotalN)
                    # YLine=(RR+cir.chr.h)*cos(2*pi*(1:TotalN)/TotalN)
                    # lines(XLine,YLine,lwd=1.5)
                    if(cir.density){
                        circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE,col='grey')
                        circle.plot(myr=RR,lwd=1.5,add=TRUE,col='grey')
                    }else{
                        circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE)
                        circle.plot(myr=RR,lwd=1.5,add=TRUE)
                    }

                }
                
                #plot the y axis of legend for each trait
                if(cir.legend==TRUE){
                    #try to get the number after radix point
                    if(Max<=1) {
                        round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
                    }else{
                        round.n=1
                    }
                    graphics::segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    graphics::segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    circle.plot(myr=r+H*(i-1)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    graphics::segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    circle.plot(myr=r+H*(i-0.75)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    graphics::segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    circle.plot(myr=r+H*(i-0.5)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    graphics::segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    circle.plot(myr=r+H*(i-0.25)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    graphics::segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    circle.plot(myr=r+H*(i-0)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    #text(-r/15,r+H*(i-0.75)+cir.band*(i-1),round(Max*0.25,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    graphics::text(-r/15,r+H*(i-0.5)+cir.band*(i-1),round(Max*0.5,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    graphics::text(-r/15,r+H*(i-0.25)+cir.band*(i-1),round(Max*0.75,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    #text(-r/15,r+H*(i-0)+cir.band*(i-1),round(Max*1,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    #text(r/5,0.4*(i-1),taxa[i],adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    
                }
                X=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN-round(band/2))/TotalN)
                Y=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN-round(band/2))/TotalN)
                # plot point in figure
                graphics::points(X[1:(length(X)-legend.bit)],Y[1:(length(Y)-legend.bit)],pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
                
                # plot significant line
                if(!is.null(threshold)){
                    if(sum(threshold!=0)==length(threshold)){
                        for(thr in 1:length(threshold)){
                            significantline1=ifelse(LOG10, H*(-log10(threshold[thr]))/Max, H*(threshold[thr])/Max)
                            #s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
                            #s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
                            # plot significant line
                            if(significantline1<H){
                                #lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
                                #if(thr==length(threshold))circle.plot(myr=(significantline1+r+H*(i-1)+cir.band*(i-1)),col="black",lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                                #print("!!!!!")
                                circle.plot(myr=(significantline1+r+H*(i-1)+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                                #circle.plot(myr=(significantline1+r+H*(i-1)+cir.band*(i-1)),col="black",lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                            }else{
                                warning(paste("No significant points for ",taxa[i]," pass the threshold level using threshold=",threshold[thr],"!",sep=""))
                            }
                        }
                    }
                }
                
                if(!is.null(threshold)){
                    if(sum(threshold!=0)==length(threshold)){
                        if(amplify==TRUE){
                            if(LOG10){
                                threshold <- sort(threshold)
                                significantline1=H*(-log10(max(threshold)))/Max
                            }else{
                                threshold <- sort(threshold, decreasing=TRUE)
                                significantline1=H*(min(threshold))/Max
                            }
                            
                            p_amp.index <- which(Cpvalue>=significantline1)
                            HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
                            HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
                            
                            #cover the points that exceed the threshold with the color "white"
                            graphics::points(HX1,HY1,pch=19,cex=cex[1],col="white")
                            
                                for(ll in 1:length(threshold)){
                                    if(ll == 1){
                                        if(LOG10){
                                            significantline1=H*(-log10(threshold[ll]))/Max
                                        }else{
                                            significantline1=H*(threshold[ll])/Max
                                        }
                                        p_amp.index <- which(Cpvalue>=significantline1)
                                        HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
                                        HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
                                    }else{
                                        if(LOG10){
                                            significantline0=H*(-log10(threshold[ll-1]))/Max
                                            significantline1=H*(-log10(threshold[ll]))/Max
                                        }else{
                                            significantline0=H*(threshold[ll-1])/Max
                                            significantline1=H*(threshold[ll])/Max
                                        }
                                        p_amp.index <- which(Cpvalue>=significantline1 & Cpvalue < significantline0)
                                        HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
                                        HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
                                    }
                                
                                    if(is.null(signal.col)){
                                        # print(signal.pch)
                                        graphics::points(HX1,HY1,pch=signal.pch,cex=signal.cex[ll]*cex[1],col=rep(rep(colx,N[i]),add[[i]])[p_amp.index])
                                    }else{
                                        # print(signal.pch)
                                        graphics::points(HX1,HY1,pch=signal.pch,cex=signal.cex[ll]*cex[1],col=signal.col[ll])
                                    }
                                }
                        }
                    }
                }
                if(cir.chr==TRUE){
                    ticks1=1.07*(RR+cir.chr.h)*sin(2*pi*(ticks-round(band/2))/TotalN)
                    ticks2=1.07*(RR+cir.chr.h)*cos(2*pi*(ticks-round(band/2))/TotalN)
                    if(is.null(chr.labels)){
                        #print(length(ticks))
                        for(i in 1:(length(ticks)-1)){
                            angle=360*(1-(ticks-round(band/2))[i]/TotalN)
                            graphics::text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
                        }
                    }else{
                        for(i in 1:length(ticks)){
                            angle=360*(1-(ticks-round(band/2))[i]/TotalN)
                            graphics::text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
                        }
                    }
                }else{
                    ticks1=(0.9*r)*sin(2*pi*(ticks-round(band/2))/TotalN)
                    ticks2=(0.9*r)*cos(2*pi*(ticks-round(band/2))/TotalN)
                    if(is.null(chr.labels)){
                        for(i in 1:length(ticks)){
                        angle=360*(1-(ticks-round(band/2))[i]/TotalN)
                        graphics::text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
                        }
                    }else{
                        for(i in 1:length(ticks)){
                            angle=360*(1-(ticks-round(band/2))[i]/TotalN)
                            graphics::text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
                        }
                    }
                }
            }
            if(outward==FALSE){
                if(cir.chr==TRUE){
                    # XLine=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(1:TotalN)/TotalN)
                    # YLine=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(1:TotalN)/TotalN)
                    # lines(XLine,YLine,lwd=1.5)

                    polygon.num <- 1000
                    for(k in 1:length(chr)){
                        if(k==1){
                            polygon.index <- seq(round(band/2)+1,-round(band/2)+max(pvalue.posN.list[[1]]), length=polygon.num)
                            X1chr=(2*cir.band+RR)*sin(2*pi*(polygon.index)/TotalN)
                            Y1chr=(2*cir.band+RR)*cos(2*pi*(polygon.index)/TotalN)
                            X2chr=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
                            Y2chr=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
                                if(is.null(chr.den.col)){
                                    graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k]) 
                                }else{
                                    if(cir.density){
                                        graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
                                    }else{
                                        graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
                                    }
                                }
                        }else{
                            polygon.index <- seq(1+round(band/2)+max(pvalue.posN.list[[k-1]]),-round(band/2)+max(pvalue.posN.list[[k]]), length=polygon.num)
                            X1chr=(2*cir.band+RR)*sin(2*pi*(polygon.index)/TotalN)
                            Y1chr=(2*cir.band+RR)*cos(2*pi*(polygon.index)/TotalN)
                            X2chr=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
                            Y2chr=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
                            if(is.null(chr.den.col)){
                                graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k]) 
                            }else{
                                    if(cir.density){
                                        graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
                                    }else{
                                        graphics::polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
                                    }
                            }   
                        }
                    }
                    if(cir.density){

                        graphics::segments(
                            (2*cir.band+RR)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
                            (2*cir.band+RR)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
                            (2*cir.band+RR+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
                            (2*cir.band+RR+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
                            col=density.list$den.col, lwd=0.1
                        )
                        graphics::legend(
                            x=RR+4*cir.chr.h,
                            y=(RR+4*cir.chr.h)/2,
                            title="Density", legend=density.list$legend.y, pch=15, pt.cex = 3, col=density.list$legend.col,
                            cex=1, bty="n",
                            y.intersp=1,
                            x.intersp=1,
                            yjust=0.5, xjust=0, xpd=TRUE
                        )
                        
                    }
                    
                    if(cir.density){
                        circle.plot(myr=2*cir.band+RR+cir.chr.h,lwd=1.5,add=TRUE,col='grey')
                        circle.plot(myr=2*cir.band+RR,lwd=1.5,add=TRUE,col='grey')
                    }else{
                        circle.plot(myr=2*cir.band+RR+cir.chr.h,lwd=1.5,add=TRUE)
                        circle.plot(myr=2*cir.band+RR,lwd=1.5,add=TRUE)
                    }

                }


                if(cir.legend==TRUE){
                    
                    #try to get the number after radix point
                    if(Max<=1) {
                        round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
                    }else{
                        round.n=2
                    }
                    graphics::segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    graphics::segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    circle.plot(myr=r+H*(i-1)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    graphics::segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    circle.plot(myr=r+H*(i-0.75)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    graphics::segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    circle.plot(myr=r+H*(i-0.5)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    graphics::segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    circle.plot(myr=r+H*(i-0.25)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    graphics::segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    circle.plot(myr=r+H*(i-0)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    graphics::text(-r/15,r+H*(i-0.25)+cir.band*(i-1),round(Max*0.25,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    #text(-r/15,r+H*(i-0.5)+cir.band*(i-1),round(Max*0.5,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    graphics::text(-r/15,r+H*(i-0.75)+cir.band*(i-1),round(Max*0.75,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    #text(-r/15,r+H*(i-1)+cir.band*(i-1),round(Max*1,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    #text(r,0.4*(i-1),taxa[i],adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)

                }
                
                X=(-Cpvalue+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN-round(band/2))/TotalN)
                Y=(-Cpvalue+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN-round(band/2))/TotalN)
                #points(X,Y,pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
                graphics::points(X[1:(length(X)-legend.bit)],Y[1:(length(Y)-legend.bit)],pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
                
                if(!is.null(threshold)){
                    if(sum(threshold!=0)==length(threshold)){
                    
                        for(thr in 1:length(threshold)){
                            significantline1=ifelse(LOG10, H*(-log10(threshold[thr]))/Max, H*(threshold[thr])/Max)
                            #s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
                            #s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
                            if(significantline1<H){
                                #lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
                                circle.plot(myr=(-significantline1+r+H*i+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                            }else{
                                warning(paste("No significant points for ",taxa[i]," pass the threshold level using threshold=",threshold[thr],"!",sep=""))
                            }
                        }
                        if(amplify==TRUE){
                            if(LOG10){
                                threshold <- sort(threshold)
                                significantline1=H*(-log10(max(threshold)))/Max
                            }else{
                                threshold <- sort(threshold, decreasing=TRUE)
                                significantline1=H*(min(threshold))/Max
                            }
                            p_amp.index <- which(Cpvalue>=significantline1)
                            HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
                            HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
                            
                            #cover the points that exceed the threshold with the color "white"
                            graphics::points(HX1,HY1,pch=19,cex=cex[1],col="white")
                            
                                for(ll in 1:length(threshold)){
                                    if(ll == 1){
                                        if(LOG10){
                                            significantline1=H*(-log10(threshold[ll]))/Max
                                        }else{
                                            significantline1=H*(threshold[ll])/Max
                                        }
                                        p_amp.index <- which(Cpvalue>=significantline1)
                                        HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
                                        HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
                                    }else{
                                        if(LOG10){
                                            significantline0=H*(-log10(threshold[ll-1]))/Max
                                            significantline1=H*(-log10(threshold[ll]))/Max
                                        }else{
                                            significantline0=H*(threshold[ll-1])/Max
                                            significantline1=H*(threshold[ll])/Max
                                        }
                                        p_amp.index <- which(Cpvalue>=significantline1 & Cpvalue < significantline0)
                                        HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
                                        HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2))/TotalN)
                                    
                                    }
                                
                                    if(is.null(signal.col)){
                                        graphics::points(HX1,HY1,pch=signal.pch,cex=signal.cex[ll]*cex[1],col=rep(rep(colx,N[i]),add[[i]])[p_amp.index])
                                    }else{
                                        graphics::points(HX1,HY1,pch=signal.pch,cex=signal.cex[ll]*cex[1],col=signal.col[ll])
                                    }
                                }
                        }
                    }
                }
                
                if(cir.chr==TRUE){
                    ticks1=1.1*(2*cir.band+RR)*sin(2*pi*(ticks-round(band/2))/TotalN)
                    ticks2=1.1*(2*cir.band+RR)*cos(2*pi*(ticks-round(band/2))/TotalN)
                    if(is.null(chr.labels)){
                        for(i in 1:(length(ticks)-1)){
                          angle=360*(1-(ticks-round(band/2))[i]/TotalN)
                          graphics::text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
                        }
                    }else{
                        for(i in 1:length(ticks)){
                            angle=360*(1-(ticks-round(band/2))[i]/TotalN)
                            graphics::text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
                        }
                    }
                }else{
                    ticks1=1.0*(RR+cir.band)*sin(2*pi*(ticks-round(band/2))/TotalN)
                    ticks2=1.0*(RR+cir.band)*cos(2*pi*(ticks-round(band/2))/TotalN)
                    if(is.null(chr.labels)){
                        for(i in 1:length(ticks)){
                        
                            #adjust the angle of labels of circle plot
                            angle=360*(1-(ticks-round(band/2))[i]/TotalN)
                            graphics::text(ticks1[i],ticks2[i],chr.ori[i],srt=angle,font=2,cex=cex.axis)
                        }
                    }else{
                        for(i in 1:length(ticks)){
                            angle=360*(1-(ticks-round(band/2))[i]/TotalN)
                            graphics::text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2,cex=cex.axis)
                        }
                    }   
                }
            }
        }
        taxa=append("Centre",taxa,)
        taxa_col=rep("black",R)
        taxa_col=append("red",taxa_col)
        for(j in 1:(R+1)){
            graphics::text(r/5,0.4*(j-1),taxa[j],adj=1,col=taxa_col[j],cex=cir.legend.cex,font=2)
                    
        }
        taxa=taxa[-1]
        if(file.output) grDevices::dev.off()
    }

    if("q" %in% plot.type){
        #print("Starting QQ-plot!",quote=F)
        amplify=FALSE
        if(multracks){
            if(file.output){
                if(file=="jpg") grDevices::jpeg(paste("GAPIT.Multracks.QQ.plot.jpg",sep=""), width = R*2.5*dpi,height=5.5*dpi,res=dpi,quality = 100)
                if(file=="pdf") grDevices::pdf(paste("GAPIT.Association.QQs_Tracks.pdf",sep=""), width = R*2.5,height=5.5)
                if(file=="tiff")    grDevices::tiff(paste("GAPIT.Multracks.QQ.plot.tiff",sep=""), width = R*2.5*dpi,height=5.5*dpi,res=dpi)
                graphics::par(mfcol=c(1,R),mar = c(0,1,4,1.5),oma=c(3,5,0,0),xpd=TRUE)
            }else{
                if(is.null(grDevices::dev.list()))  grDevices::dev.new(width = 2.5*R, height = 5.5)
                graphics::par(xpd=TRUE)
            }
            for(i in 1:R){
                print(paste("Multracks_QQ Plotting ",taxa[i],"...",sep=""))     
                P.values=as.numeric(Pmap[,i+2])
                P.values=P.values[!is.na(P.values)]
                if(LOG10){
                    P.values=P.values[P.values>0]
                    P.values=P.values[P.values<=1]
                    N=length(P.values)
                    P.values=P.values[order(P.values)]
                }else{
                    N=length(P.values)
                    P.values=P.values[order(P.values,decreasing=TRUE)]
                }
                p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
                log.Quantiles <- -log10(p_value_quantiles)
                if(LOG10){
                    log.P.values <- -log10(P.values)
                }else{
                    log.P.values <- P.values
                }
                
                #calculate the confidence interval of QQ-plot
                if(conf.int){
                    N1=length(log.Quantiles)
                    c95 <- rep(NA,N1)
                    c05 <- rep(NA,N1)
                    for(j in 1:N1){
                        xi=ceiling((10^-log.Quantiles[j])*N)
                        if(xi==0)xi=1
                        c95[j] <- stats::qbeta(0.95,xi,N-xi+1)
                        c05[j] <- stats::qbeta(0.05,xi,N-xi+1)
                    }
                    index=length(c95):1
                }else{
                    c05 <- 1
                    c95 <- 1
                }
                
                YlimMax <- max(floor(max(max(-log10(c05)), max(-log10(c95)))+1), floor(max(log.P.values)+1))
                plot(NULL, xlim = c(0,floor(max(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=1.2,ylim=c(0,YlimMax),xlab ="", ylab="", main = taxa[i])
                graphics::axis(1, at=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), labels=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), cex.axis=cex.axis)
                graphics::axis(2, at=seq(0,YlimMax,ceiling(YlimMax/10)), labels=seq(0,YlimMax,ceiling(YlimMax/10)), cex.axis=cex.axis)
                
                #plot the confidence interval of QQ-plot
                
                if(conf.int)    graphics::polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=conf.int.col,border=conf.int.col)
                
                if(!is.null(threshold.col)){
                    graphics::par(xpd=FALSE);
                    graphics::abline(a = 0, b = 1, col = threshold.col[1],lwd=2);
                    graphics::par(xpd=TRUE)
                }
                graphics::points(log.Quantiles, log.P.values, col = col[1],pch=1,cex=cex[3])
                #print(max(log.Quantiles))
                #   print(length(log.Quantiles))
                #   print(length(log.P.values))
                if(!is.null(threshold)){
                    if(sum(threshold!=0)==length(threshold)){
                        thre.line=-log10(min(threshold))
                        if(amplify==TRUE){
                            thre.index=which(log.P.values>=thre.line)
                            if(length(thre.index)!=0){
                            
                                #cover the points that exceed the threshold with the color "white"
                                graphics::points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
                                if(is.null(signal.col)){
                                    graphics::points(log.Quantiles[thre.index],log.P.values[thre.index],col = col[1],pch=signal.pch[1],cex=signal.cex[1])
                                }else{
                                    graphics::points(log.Quantiles[thre.index],log.P.values[thre.index],col = signal.col[1],pch=signal.pch[1],cex=signal.cex[1])
                                }
                            }
                        }
                    }
                }
            }
            if(box) box()
            if(file.output) grDevices::dev.off()
            if(R > 1){
                #qq_col=rainbow(R)
                qq_col=rep(c( '#FF6A6A',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5'),ceiling(R/5))

                signal.col <- NULL
                if(file.output){
                    if(file=="jpg") grDevices::jpeg(paste("GAPIT.Multiple.QQ.plot.symphysic.jpg",sep=""), width = 5.5*dpi,height=5.5*dpi,res=dpi,quality = 100)
                    if(file=="pdf") grDevices::pdf(paste("GAPIT.Association.QQs_Symphysic.pdf",sep=""), width = 5.5,height=5.5)
                    if(file=="tiff")    grDevices::tiff(paste("GAPIT.Multiple.QQ.plot.symphysic.tiff",sep=""), width = 5.5*dpi,height=5.5*dpi,res=dpi)
                    graphics::par(mar = c(5,5,4,2),xpd=TRUE)
                }else{
                    grDevices::dev.new(width = 5.5, height = 5.5)
                    graphics::par(xpd=TRUE)
                }
                P.values=as.numeric(Pmap[,i+2])
                P.values=P.values[!is.na(P.values)]
                if(LOG10){
                    P.values=P.values[P.values>0]
                    P.values=P.values[P.values<=1]
                    N=length(P.values)
                    P.values=P.values[order(P.values)]
                }else{
                    N=length(P.values)
                    P.values=P.values[order(P.values,decreasing=TRUE)]
                }
                p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
                log.Quantiles <- -log10(p_value_quantiles)
                                            
                # calculate the confidence interval of QQ-plot
                if(conf.int){
                    N1=length(log.Quantiles)
                    c95 <- rep(NA,N1)
                    c05 <- rep(NA,N1)
                    for(j in 1:N1){
                        xi=ceiling((10^-log.Quantiles[j])*N)
                        if(xi==0)xi=1
                        c95[j] <- stats::qbeta(0.95,xi,N-xi+1)
                        c05[j] <- stats::qbeta(0.05,xi,N-xi+1)
                    }
                    index=length(c95):1
                }
                
                if(!conf.int){c05 <- 1; c95 <- 1}
                
                Pmap.min <- Pmap[,3:(R+2)]

                YlimMax <- max(floor(max(max(-log10(c05)), max(-log10(c95)))+1), -log10(min(Pmap.min[Pmap.min > 0])))
                plot(NULL, xlim = c(0,floor(max(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=1.2,ylim=c(0, floor(YlimMax+1)),xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = "QQ plot")
                #legend("topleft",taxa,col=t(col)[1:R],pch=1,pt.lwd=2,text.font=6,box.col=NA)           
                graphics::legend("topleft",taxa,col=qq_col[1:R],pch=1,pt.lwd=3,text.font=6,box.col=NA)
                graphics::axis(1, at=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), labels=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), cex.axis=cex.axis)
                graphics::axis(2, at=seq(0,floor(YlimMax+1),ceiling((YlimMax+1)/10)), labels=seq(0,floor((YlimMax+1)),ceiling((YlimMax+1)/10)), cex.axis=cex.axis)
                #print(log.Quantiles[index])
                #print(index)
                #print(length(log.Quantiles))

                # plot the confidence interval of QQ-plot
                if(conf.int)    graphics::polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=conf.int.col,border=conf.int.col)
                
                for(i in 1:R){
                    #print(paste("Multraits_QQ Plotting ",taxa[i],"...",sep=""))
                    P.values=as.numeric(Pmap[,i+2])
                    P.values=P.values[!is.na(P.values)]
                    if(LOG10){
                    P.values=P.values[P.values>0]
                    P.values=P.values[P.values<=1]
                    N=length(P.values)
                    P.values=P.values[order(P.values)]
                    }else{
                    N=length(P.values)
                    P.values=P.values[order(P.values,decreasing=TRUE)]
                    }
                    p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
                    log.Quantiles <- -log10(p_value_quantiles)
                    if(LOG10){
                    log.P.values <- -log10(P.values)
                    }else{
                    log.P.values <- P.values
                    }
                
                        
                    if((i == 1) & !is.null(threshold.col)){
                        graphics::par(xpd=FALSE);
                        graphics::abline(a = 0, b = 1, col = threshold.col[1],lwd=2);
                        graphics::par(xpd=TRUE)}
                    #print(length(log.Quantiles))
                    #print("!!!!!") 
                    #points(log.Quantiles, log.P.values, col = t(col)[i],pch=1,lwd=3,cex=cex[3])
                    graphics::points(log.Quantiles, log.P.values, col = qq_col[i],pch=1,lwd=3,cex=cex[3])
                    
                    #print(max(log.Quantiles))
                    #
    
                    if(!is.null(threshold)){
                        if(sum(threshold!=0)==length(threshold)){
                            thre.line=-log10(min(threshold))
                            if(amplify==TRUE){
                                thre.index=which(log.P.values>=thre.line)
                                if(length(thre.index)!=0){
                                
                                    # cover the points that exceed the threshold with the color "white"
                                    graphics::points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,lwd=3,cex=cex[3])
                                    if(is.null(signal.col)){
                                        graphics::points(log.Quantiles[thre.index],log.P.values[thre.index],col = t(col)[i],pch=signal.pch[1],cex=signal.cex[1])
                                    }else{
                                        graphics::points(log.Quantiles[thre.index],log.P.values[thre.index],col = signal.col[1],pch=signal.pch[1],cex=signal.cex[1])
                                    }
                                }
                            }
                        }
                    }
                }
                    box()
                if(file.output) grDevices::dev.off()
            }
        }else{
            for(i in 1:R){
                print(paste("Q_Q Plotting ",taxa[i],"...",sep=""))
                if(file.output){
                    if(file=="jpg") grDevices::jpeg(paste("QQplot.",taxa[i],".jpg",sep=""), width = 5.5*dpi,height=5.5*dpi,res=dpi,quality = 100)
                    if(file=="pdf") grDevices::pdf(paste("GAPIT.Association.Q_Q.",taxa[i],".pdf",sep=""), width = 5.5,height=5.5)
                    if(file=="tiff")    grDevices::tiff(paste("QQplot.",taxa[i],".tiff",sep=""), width = 5.5*dpi,height=5.5*dpi,res=dpi)
                    graphics::par(mar = c(5,5,4,2),xpd=TRUE)
                }else{
                    if(is.null(grDevices::dev.list()))  grDevices::dev.new(width = 5.5, height = 5.5)
                    graphics::par(xpd=TRUE)
                }
                P.values=as.numeric(Pmap[,i+2])
                P.values=P.values[!is.na(P.values)]
                if(LOG10){
                    P.values=P.values[P.values>0]
                    P.values=P.values[P.values<=1]
                    N=length(P.values)
                    P.values=P.values[order(P.values)]
                }else{
                    N=length(P.values)
                    P.values=P.values[order(P.values,decreasing=TRUE)]
                }
                p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
                log.Quantiles <- -log10(p_value_quantiles)
                if(LOG10){
                    log.P.values <- -log10(P.values)
                }else{
                    log.P.values <- P.values
                }
                
                #calculate the confidence interval of QQ-plot
                if(conf.int){
                    N1=length(log.Quantiles)
                    c95 <- rep(NA,N1)
                    c05 <- rep(NA,N1)
                    for(j in 1:N1){
                        xi=ceiling((10^-log.Quantiles[j])*N)
                        if(xi==0)xi=1
                        c95[j] <- stats::qbeta(0.95,xi,N-xi+1)
                        c05[j] <- stats::qbeta(0.05,xi,N-xi+1)
                    }
                    index=length(c95):1
                }else{
                    c05 <- 1
                    c95 <- 1
                }
                #print(max(log.Quantiles))
                #print("@@@@@")
                YlimMax <- max(floor(max(max(-log10(c05)), max(-log10(c95)))+1), floor(max(log.P.values)+1))
                plot(NULL, xlim = c(0,floor(max(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=1.2,ylim=c(0,YlimMax),xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = paste("QQplot of",taxa[i]))
                graphics::axis(1, at=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), labels=seq(0,floor(max(log.Quantiles)+1),ceiling((max(log.Quantiles)+1)/10)), cex.axis=cex.axis)
                graphics::axis(2, at=seq(0,YlimMax,ceiling(YlimMax/10)), labels=seq(0,YlimMax,ceiling(YlimMax/10)), cex.axis=cex.axis)
                
                #plot the confidence interval of QQ-plot
                #print(log.Quantiles[index])
                qq_col = grDevices::rainbow(R)
                #if(conf.int)   polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=conf.int.col,border=conf.int.col)
                if(conf.int)    graphics::polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=qq_col[i],border=conf.int.col)
                
                if(!is.null(threshold.col)){
                    graphics::par(xpd=FALSE);
                    graphics::abline(a = 0, b = 1, col = threshold.col[1],lwd=2);
                    graphics::par(xpd=TRUE)
                    }
                 
                graphics::points(log.Quantiles, log.P.values, col = col[1],pch=19,cex=2)
                
                if(!is.null(threshold)){
                    if(sum(threshold!=0)==length(threshold)){
                        thre.line=-log10(min(threshold))
                        if(amplify==TRUE){
                            thre.index=which(log.P.values>=thre.line)
                            if(length(thre.index)!=0){
                                #print("!!!!")
                                #cover the points that exceed the threshold with the color "white"
                                graphics::points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,lwd=3,cex=cex[3])
                                if(is.null(signal.col)){
                                    graphics::points(log.Quantiles[thre.index],log.P.values[thre.index],col = col[1],pch=signal.pch[1],cex=signal.cex[1])
                                }else{
                                    graphics::points(log.Quantiles[thre.index],log.P.values[thre.index],col = signal.col[1],pch=signal.pch[1],cex=signal.cex[1])
                                }
                            }
                        }
                    }
                }
                box()
                if(file.output) grDevices::dev.off()
            }
        }
        print("Multiple QQ plot has been finished!",quote=F)
    }

        





    }#End of Whole function

#}



#===========BGLR====================================================================================
#This function creates an incidence matrix that will be included in the 
#linear term of the model
#Arguments: LT, Linear term, an object of the class "formula" that also includes
#optionally a data.frame to obtain the information
#It returns the incidence matrix
set.X=function(LT)
{

    flag=TRUE
    n_elements=length(LT)
    i=0
    while(i<=n_elements & flag)
    {
       i=i+1;
       if(is(LT[[i]],"formula"))
       {
            flag=FALSE
            rhs=LT[[i]]
            if(is.null(LT$data))
            {
                mf = model.frame(formula=rhs)
            }else{
                mf = model.frame(formula=rhs,data=LT$data)
            }
            X = model.matrix(attr(mf, "terms"), data=mf)
            Xint = match("(Intercept)", colnames(X), nomatch=0L)
            if(Xint > 0L) X = X[, -Xint, drop=FALSE]
       }
    }
    if(flag) stop("Unable to build incidence matrix, wrong formula or data")
    return(X)
}

## Fixed Effects ##################################################################
#Function for initializing regression coefficients for Fixed effects.
#All the arguments are defined in the function BGLR
setLT.Fixed=function(LT,n,j,y,weights,nLT,saveAt,rmExistingFiles,groups,nGroups)
{

    if(is.null(LT$X)) LT$X=set.X(LT)

    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
    LT$colNames=colnames(LT$X)
    
    if(any(is.na(LT$X)))
    { 
    stop("LP ",j," has NAs in X")
    }

    if(nrow(LT$X)!=n)
    {
        stop("Number of rows of LP ",j,"  not equal to the number of phenotypes")
    }

    #weight inputs if necessary
        
    LT$X=sweep(LT$X,1L,weights,FUN="*")        #weights

    if(!is.null(groups))
    {
        x2=matrix(NA,nrow=nGroups,ncol=ncol(LT$X))
        for(g in 1:nGroups)
        {
                x2[g,]=apply(LT$X[groups==g,,drop=FALSE],2L,function(x) sum(x^2)) #the sum of the square of each of the columns for each group
        }
        LT$x2=x2;
    }else{
        LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    }   
    
    #Objects for saving posterior means from MCMC
    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)   

    fname=paste(saveAt,LT$Name,"_b.dat",sep="")

    LT$NamefileOut=fname; 

    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    LT$fileOut=file(description=fname,open="w")
    tmp=LT$colNames
    write(tmp, ncolumns = LT$p, file = LT$fileOut, append = TRUE)
    LT$X=as.vector(LT$X)
    LT$x2=as.vector(LT$x2)
    LT$varB=1e10
    return(LT)
}

## Gaussian Regression ############################################################
#Function for initializing regression coefficients for Ridge Regression.
#All the arguments are defined in the function BGLR
setLT.BRR=function(LT,y,n,j,weights,nLT,R2,saveAt,rmExistingFiles,groups,nGroups,verbose,thin,nIter,burnIn,lower_tri){ #*#

    #Check inputs

    if(is.null(LT$X)) LT$X=set.X(LT)
       
    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
    LT$colNames=colnames(LT$X)
    
    if(any(is.na(LT$X)))
    { 
      stop("LP ",j," has NAs in X")
    }
    
    if(nrow(LT$X)!=n)
    {
      stop("Number of rows of LP ",j,"  not equal to the number of phenotypes")
    }   

    #Weight inputs if necessary
    LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights

    if(!is.null(groups))
    {
    
    x2=matrix(NA,nrow=nGroups,ncol=ncol(LT$X))
    for(g in 1:nGroups)
    {
        x2[g,]=apply(LT$X[groups==g,,drop=FALSE],2L,function(x) sum(x^2)) #the sum of the square of each of the columns for each group
    }
        LT$x2=x2;
    }else{
    LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    }  

    sumMeanXSq = sum((apply(LT$X,2L,mean))^2)

    #Default df for the prior assigned to the variance of marker effects
    if(is.null(LT$df0))
    {
    LT$df0=5

    if(verbose)
    {
        #message("Degree of freedom of LP ",j,"  set to default value (",LT$df0,")")
    }
    }

    if(is.null(LT$R2))
    { 
        LT$R2=R2/nLT
    }

 
    #Default scale parameter for the prior assigned to the variance of marker effects
    if(is.null(LT$S0))
    {
        if(LT$df0<=0) stop("df0>0 in BRR in order to set S0")

    LT$MSx=sum(LT$x2)/n-sumMeanXSq       
    LT$S0=((var(y,na.rm=TRUE)*LT$R2)/(LT$MSx))*(LT$df0+2)  
    
    if(verbose)
    {
        #message("Scale parameter of LP ",j,"  set to default value (",LT$S0,")")
    }
    }

    
    #Objects for saving posterior means from MCMC
    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    LT$varB=LT$S0/(LT$df0+2)
    LT$post_varB=0                 
    LT$post_varB2=0
    
    fname=paste(saveAt,LT$Name,"_varB.dat",sep=""); 
    
    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    LT$NamefileOut=fname
    LT$fileOut=file(description=fname,open="w")

    if(is.null(LT$lower_tri)) LT$lower_tri=FALSE;

    if(LT$lower_tri)
    {
    #message("You have provided a lower triangular matrix for LP ", j)
    #message("Checking dimmensions...")
    if(ncol(LT$X)==nrow(LT$X))
    {
        #message("Ok.")
        LT$X=LT$X[lower.tri(LT$X,diag=TRUE)]
        }
    }else{
     LT$X=as.vector(LT$X)
    }

    LT$x2=as.vector(LT$x2)
    
    #*#
    if(is.null(LT$saveEffects)){LT$saveEffects=FALSE}
    if(LT$saveEffects){
        if(is.null(LT$storageMode)){LT$storageMode="double"}
        if(!LT$storageMode%in%c("single","double")) {
            stop("storageMode of LP ",j," can either be 'single' or 'double' (default)")
        }
        if(is.null(LT$thin)){ LT$thin=thin }
        fname=paste(saveAt,LT$Name,"_b.bin",sep="")
        if(rmExistingFiles){ unlink(fname) }
        LT$fileEffects=file(fname,open='wb')
        nRow=floor((nIter-burnIn)/LT$thin)
        writeBin(object=c(nRow,LT$p),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
    }#*#

    return(LT)
}

#Ridge regression using sets of markers
#This is just a Ridge Regression set-specific variances, 
#LT has an extra attribute: sets

setLT.BRR_sets=function(LT,y,n,j,weights,nLT,R2,saveAt,rmExistingFiles,verbose,thin,nIter,burnIn)
{

    
    #Check the inputs
    if(is.null(LT$X)) LT$X=set.X(LT)
   
    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
    LT$colNames=colnames(LT$X)

    if(is.null(LT$sets)) stop("Argument sets (a vector grouping effects into sets) is required in BRR_sets");
    if(length(LT$sets)!=LT$p){ stop("The length of sets must be equal to the number of predictors") }
    
    LT$sets<-as.integer(factor(LT$sets,ordered=TRUE,levels=unique(LT$sets)))
        
    LT$n_sets=length(unique(LT$sets))   
    
    if(LT$n_sets>=LT$p){ stop("The number of sets is greater or equal than the number of effects!") }
    
    if(any(is.na(LT$X))){ 
      stop("LP ",j," has NAs in X")
    }
    
    if(nrow(LT$X)!=n){
      stop("Number of rows of LP ",j,"  not equal to the number of phenotypes")
    }   

    #Weight inputs if necessary
    LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
    LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    sumMeanXSq = sum((apply(LT$X,2L,mean))^2)
    

    if(is.null(LT$df0)){
        LT$df0=5
        if(verbose){
            #message("Degree of freedom of LP ",j,"  set to default value (",LT$df0,")")
        }
    }

    if(is.null(LT$R2)){ 
        LT$R2=R2/nLT
    }

    if(is.null(LT$S0))    {
        if(LT$df0<=0) stop("df0 must be greater than 0")

        LT$MSx=sum(LT$x2)/n-sumMeanXSq       
        LT$S0=((var(y,na.rm=TRUE)*LT$R2)/(LT$MSx))*(LT$df0+2) 
        if(verbose){ 
            #message("Scale parameter of LP ",j,"  set to default value (",LT$S0,")")
        }
    }
    
    LT$DF1=table(LT$sets)+LT$df0
    
    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    LT$varB=rep(LT$S0/(LT$df0+2),LT$p)
    LT$varSets=rep(0,LT$n_sets)
    LT$post_varSets=rep(0,LT$n_sets)
    LT$post_varSets2<-rep(0,LT$n_sets)
    LT$post_varB=rep(0 ,LT$p)                
    LT$post_varB2=rep(0,LT$p)

    fname=paste(saveAt,LT$Name,"_varB.dat",sep=""); 
    
    if(rmExistingFiles){ 
       unlink(fname) 
    }

    LT$NamefileOut=fname
    LT$fileOut=file(description=fname,open="w")
    LT$X=as.vector(LT$X)

    if(is.null(LT$saveEffects)){LT$saveEffects=FALSE}
    if(LT$saveEffects){
        if(is.null(LT$storageMode)){LT$storageMode="double"}
        if(!LT$storageMode%in%c("single","double")) {
            stop("storageMode of LP ",j," can either be 'single' or 'double' (default)")
        }
        if(is.null(LT$thin)){ LT$thin=thin }
        fname=paste(saveAt,LT$Name,"_b.bin",sep="")
        if(rmExistingFiles){ unlink(fname) }
        LT$fileEffects=file(fname,open='wb')
        nRow=floor((nIter-burnIn)/LT$thin)
        writeBin(object=c(nRow,LT$p),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
    }
    return(LT)
}

## Bayesian LASSO ############################################################
## The well known Bayesian LASSO (Park and Casella, 2008) and 
## de los Campos et al (2009). 
#  This functions simply sets hyper-parameters for quantities involved in BL regression

setLT.BL=function(LT,y,n,j,weights,nLT,R2,saveAt,rmExistingFiles,verbose,thin,nIter,burnIn)
{
    #Check the inputs
    if(is.null(LT$minAbsBeta)) LT$minAbsBeta=1e-9

    if(is.null(LT$X)) LT$X=set.X(LT)

    LT$X=as.matrix(LT$X)
    LT$p=ncol(LT$X)
    LT$colNames=colnames(LT$X)
        
    if(any(is.na(LT$X)))
    {
       stop("LP ",j," has NAs in X")
    }

    if(nrow(LT$X)!=n)
    {
        stop("Number of rows of LP ",j,"  not equal to the number of phenotypes")
    } 

    #Wheight inputs if necessary
    LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
    LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
    sumMeanXSq = sum((apply(LT$X,2L,mean))^2)   

    LT$MSx=sum(LT$x2)/n-sumMeanXSq

    # Prior
    if(is.null(LT$R2))
    { 
       LT$R2=R2/nLT
    }
    
    # Setting default value of lambda
    if(!is.null(LT$lambda))
    { 
        if(LT$lambda<0)
        {
        stop("lambda should be positive")
    }  
    }
    if(is.null(LT$lambda))
    {
        LT$lambda2=2*(1-R2)/(LT$R2)*LT$MSx
        LT$lambda=sqrt(LT$lambda2)
        if(verbose)
        {
            #message("Initial value of lambda in LP ",j," was set to default value (",LT$lambda,")")
        }
    }else{
    if(LT$lambda<0) stop("lambda should be positive");
        LT$lambda2=LT$lambda^2
    }
    
    # Checking lambda-type
    if(is.null(LT$type))
    {
        LT$type="gamma"
        if(verbose)
        {
            #message("By default, the prior density of lambda^2 in the LP ",j,"  was set to gamma")
        }
    }else{
        if(!LT$type%in%c("gamma","beta","FIXED")) stop("The prior for lambda^2 should be gamma, beta or a point of mass (i.e., fixed lambda)")
    }
    if(LT$type=="gamma")
    {
        if(is.null(LT$shape))
                {
             LT$shape=1.1
             if(verbose)
             {
                #message("shape parameter in LP ",j," was missing and was set to ",LT$shape)
             }
        }
        
        if(is.null(LT$rate))
                {
             LT$rate=(LT$shape-1)/LT$lambda2
             if(verbose)
             {
                #message("rate parameter in LP ",j," was missing and was set to ",LT$rate)
             }
        }   
    }
    
    if(LT$type=="beta")
    {
                if(is.null(LT$probIn))
        {
                LT$probIn=0.5
            if(verbose)
            {
                    #message("probIn in LP ",j," was missing and was set to ",LT$probIn)
            }
        }

        if(is.null(LT$counts))
        {
                LT$counts=2
            if(verbose)
            {
                    #message("Counts in LP ",j," was missing and was set to ",LT$counts)
            }
        } 

                LT$shape1=LT$probIn*LT$counts;
                LT$shape2=(1-LT$probIn)*LT$counts;

        if(is.null(LT$max))
        {
            LT$max=10*LT$lambda
            if(verbose)
            {
                #message("max parameter in LP ",j," was missing and was set to ",LT$max)
            }
        }
    }

    #Objects to storing information for MCMC iterations

    LT$b=rep(0,LT$p)
    LT$post_b=rep(0,LT$p)
    LT$post_b2=rep(0,LT$p)
    
    tmp=((var(y,na.rm=TRUE)*R2/nLT)/(LT$MSx))
    LT$tau2=rep(tmp,LT$p)
    LT$post_tau2=0  
    LT$post_lambda=0
    
    fname=paste(saveAt,LT$Name,"_lambda.dat",sep="");
    
    if(rmExistingFiles)
    { 
       unlink(fname) 
    }
    
    LT$NamefileOut=fname
    LT$fileOut=file(description=fname,open="w")

    LT$X=as.vector(LT$X)
    
    #*#
    if(is.null(LT$saveEffects)){LT$saveEffects=FALSE}
    if(LT$saveEffects){
        if(is.null(LT$storageMode)){LT$storageMode="double"}
        if(!LT$storageMode%in%c("single","double")) {
            stop("storageMode of LP ",j," can either be 'single' or 'double' (default)")
        }
        if(is.null(LT$thin)){ LT$thin=thin }
        fname=paste(saveAt,LT$Name,"_b.bin",sep="")
        if(rmExistingFiles){ unlink(fname) }
        LT$fileEffects=file(fname,open='wb')
        nRow=floor((nIter-burnIn)/LT$thin)
        writeBin(object=c(nRow,LT$p),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
    }#*#
    
    return(LT)
}


#Reproducing kernel Hilbert spaces
#This function simply sets hyperparameters and prepares inputs 
#for Reproducing Kernel Hilbert Spaces. The algorithm used here is 
#Fully described in de los Campos et al (2010).

setLT.RKHS=function(LT,y,n,j,weights,saveAt,R2,nLT,rmExistingFiles,verbose)  
{

    #Checking inputs
    if(is.null(LT$V))
    {
        if(is.null(LT$K)) stop("Kernel for linear term ",j, " was not provided, specify it with list(K=?,model='RKHS'), where ? is the kernel matrix")

        if(!is.matrix(LT$K)) stop("Kernel for linear term ",j, " should be a matrix, the kernel provided is of class ", class(LT$K))
        
        LT$K = as.matrix(LT$K)

        if(nrow(LT$K)!=ncol(LT$K)) stop("Kernel for linear term ",j, " is not a square matrix")

        #This code was rewritten to speed up computations
        #T = diag(weights)   
        #LT$K = T %*% LT$K %*% T 
        
        #Weight kernels
        #for(i in 1:nrow(LT$K))
        #{
        #   #j can not be used as subindex because its value is overwritten
        #   for(m in i:ncol(LT$K))
        #    {    
        #           LT$K[i,m]=LT$K[i,m]*weights[i]*weights[m];
        #            LT$K[m,i]=LT$K[i,m]
        #   }
        #}
        
        #Added January 10/2020
        #This is faster than the for loop
        
        LT$K=sweep(sweep(LT$K,1L,weights,"*"),2L,weights,"*")
    
        tmp =eigen(LT$K,symmetric=TRUE)
        LT$V =tmp$vectors
        LT$d =tmp$values
        rm(tmp)
    
    }else{
        if(any(weights!=1))
        { 
            warning("Eigen decomposition for LT",j," was provided and the model involves weights. Note: You should have weighted the kernel before computing eigen(K)") 
        }
    }
    
    #Defaul value for tolD
    #Only those eigenvectors whose eigenvalues> tolD are kept.
    if (is.null(LT$tolD)) 
    {
       LT$tolD = 1e-10
       if(verbose)
       {
            #message("Default value of minimum eigenvalue in LP ",j," was set to ",LT$tolD)
       }
    }
    
    #Removing elements whose eigenvalues < tolD
    tmp= LT$d > LT$tolD
    LT$levelsU = sum(tmp)
    LT$d = LT$d[tmp]
    LT$V = LT$V[, tmp]
    
    #Default degrees of freedom and scale parameter associated with the variance component for marker effect
    if (is.null(LT$df0)) 
    {
      LT$df0 = 5
      if(verbose)
      {
        #message("default value of df0 in LP ",j," was missing and was set to ",LT$df0)
      }
    }
   
    if(is.null(LT$R2))
    { 
           LT$R2=R2/nLT
    }

    if (is.null(LT$S0)) 
    {
          if(LT$df0<=0) stop("df0>0 in RKHS in order to set S0");

          LT$S0=((var(y,na.rm=TRUE)*LT$R2)/(mean(LT$d)))*(LT$df0+2)

      if(verbose)
      {
             #message("default value of S0 in LP ",j," was missing and was set to ",LT$S0)
      }
    }
    
    LT$u=rep(0,nrow(LT$V))
    
    LT$varU=LT$S0/(LT$df0+2)
       
    LT$uStar=rep(0, LT$levelsU)
    
    #Output files
    fname=paste(saveAt,LT$Name,"_varU.dat",sep="")
    LT$NamefileOut=fname; 

    if(rmExistingFiles)
    { 
       unlink(fname) 
    }

    #Objects for storing information for MCMC iterations
    LT$fileOut=file(description=fname,open="w")
    LT$post_varU=0
    LT$post_varU2=0
    LT$post_uStar = rep(0, LT$levelsU)
    LT$post_u = rep(0, nrow(LT$V))
    LT$post_u2 = rep(0,nrow(LT$V))
    
    #return object
    return(LT)
}

###Bayes B and C########################################################################################################################################                 

setLT.BayesBandC=function(LT,y,n,j,weights,saveAt,R2,nLT,rmExistingFiles, groups, nGroups,verbose,thin,nIter,burnIn)
{

  model=LT$model
 
  if(is.null(LT$X)) LT$X=set.X(LT)

  #Be sure that your X is a matrix
  LT$X=as.matrix(LT$X)  
  LT$p=ncol(LT$X)
  LT$colNames=colnames(LT$X)

  #Weight inputs if necessary
  LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights

  if(!is.null(groups))
  {

        x2=matrix(NA,nrow=nGroups,ncol=ncol(LT$X))
        for(g in 1:nGroups)
        {
                x2[g,]=apply(LT$X[groups==g,,drop=FALSE],2L,function(x) sum(x^2)) #the sum of the square of each of the columns for each group
        }
        LT$x2=x2;
    }else{
        LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
  }

  sumMeanXSq = sum((apply(LT$X,2L,mean))^2)
  LT$MSx=sum(LT$x2)/n-sumMeanXSq

  
  if(any(is.na(LT$X))){ stop("LP ",j," has NAs in X") }
  if(nrow(LT$X)!=n){ stop("Number of rows of LP ",j,"  not equal to the number of phenotypes") }
  
  if(is.null(LT$R2))
  {
    LT$R2=R2/nLT
    if(verbose)
    {
      #message("R2 in LP ",j," was missing and was set to ",LT$R2)
    }
  }

  #Default value for the degrees of freedom associated with the distribution assigned to the variance
  #of marker effects
  if(is.null(LT$df0))
  {
    LT$df0= 5
    if(verbose)
    {
        #message("DF in LP ",j," was missing and was set to ",LT$df0)
    }
  }


  #Default value for a marker being "in" the model
  if(is.null(LT$probIn))
  {
    LT$probIn=0.5
    if(verbose)
    {   
       #message("probIn in LP ",j," was missing and was set to ",LT$probIn)
    }
  } 


   #Default value for prior counts
  if(is.null(LT$counts))
  {
    LT$counts=10
    if(verbose)
    {
       #message("Counts in LP ",j," was missing and was set to ",LT$counts)
    }
  }

  LT$countsIn=LT$counts * LT$probIn
  LT$countsOut=LT$counts - LT$countsIn

  #Default value for the scale parameter associated with the distribution assigned to the variance of 
  #marker effects
  if(is.null(LT$S0))
  {
     if(LT$df0<=0) stop("df0>0 in ",model," in order to set S0");
     LT$S0=var(y, na.rm = TRUE)*LT$R2/(LT$MSx)*(LT$df0+2)/LT$probIn
     if(verbose)
     {
        #message("Scale parameter in LP ",j," was missing and was set to ",LT$S0)
     }
  }
 
  LT$b=rep(0, LT$p)
  LT$d=rbinom(n = LT$p, size = 1, prob = LT$probIn)

  if(model=="BayesB")
  {
        if(is.null(LT$shape0))
    {
            LT$shape0=1.1
    }
    if(is.null(LT$rate0)){
            LT$rate0=(LT$shape0-1)/LT$S0
    }
        LT$S=LT$S0
    LT$varB = LT$varB=rep(LT$S0/(LT$df0+2),LT$p)
        fname=paste(saveAt,LT$Name,"_parBayesB.dat",sep="")
  }else{
    LT$varB = LT$S0
    fname=paste(saveAt,LT$Name,"_parBayesC.dat",sep="")
        
  }

  LT$X=as.vector(LT$X)
  LT$x2=as.vector(LT$x2)

  if(rmExistingFiles)
  { 
      unlink(fname) 
  }

  LT$fileOut=file(description=fname,open="w")
  LT$NamefileOut=fname;
  
  if(model=="BayesB")
  {
    tmp=c('probIn','scale')
    write(tmp, ncolumns = LT$p, file = LT$fileOut, append = TRUE)
  }

  #Objects for storing MCMC information 
  LT$post_varB=0
  LT$post_varB2=0
  LT$post_d=0
  LT$post_probIn=0
  LT$post_probIn2=0
  LT$post_b=rep(0,LT$p)
  LT$post_b2=rep(0,LT$p)

  if(model=="BayesB")
  {
     LT$post_S=0
     LT$post_S2=0
  }
  
    #*#
    if(is.null(LT$saveEffects)){LT$saveEffects=FALSE}
    if(LT$saveEffects){
        if(is.null(LT$storageMode)){LT$storageMode="double"}
        if(!LT$storageMode%in%c("single","double")) {
            stop("storageMode of LP ",j," can either be 'single' or 'double' (default)")
        }
        if(is.null(LT$thin)){ LT$thin=thin }
        fname=paste(saveAt,LT$Name,"_b.bin",sep="")
        if(rmExistingFiles){ unlink(fname) }
        LT$fileEffects=file(fname,open='wb')
        nRow=floor((nIter-burnIn)/LT$thin)
        writeBin(object=c(nRow,LT$p),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
    }#*#
  
  #return object
  return(LT) 
}

#Bayes A, Mewissen et al. (2001).
#Prediction of Total Genetic Value Using Genome-Wide Dense Marker Maps
#Genetics 157: 1819-1829, Modified so that the Scale parameter is estimated from data (a gamma prior is assigned)

setLT.BayesA=function(LT,y,n,j,weights,saveAt,R2,nLT,rmExistingFiles,verbose,thin,nIter,burnIn)
{

  #Ckecking inputs 
  if(is.null(LT$X)) LT$X=set.X(LT)
 
  LT$X=as.matrix(LT$X)
  LT$p=ncol(LT$X)
  LT$colNames=colnames(LT$X)

  #Weight inputs if necessary
  LT$X=sweep(LT$X,1L,weights,FUN="*")  #weights
  LT$x2=apply(LT$X,2L,function(x) sum(x^2))  #the sum of the square of each of the columns
  sumMeanXSq = sum((apply(LT$X,2L,mean))^2)
  LT$MSx=sum(LT$x2)/n-sumMeanXSq

  #Default degrees of freedom for the prior assigned to the variance of markers
  if(is.null(LT$df0))
  {
     LT$df0 = 5 
     if(verbose)
     { 
        #message("DF in LP ",j," was missing and was set to ",LT$df0)
     }
  }
  if(is.null(LT$R2))
  {
    LT$R2=R2/nLT
    if(verbose)
    {
        #message("R2 in LP ",j," was missing and was set to ",LT$R2)
    }
  }

  #Defuault scale parameter for the prior assigned to the variance of markers
  if(is.null(LT$S0))
  {
     if(LT$df0<=0) stop("df0>0 in BayesA in order to set S0")
     LT$S0 = var(y, na.rm = TRUE)*LT$R2/(LT$MSx)*(LT$df0+2)
     if(verbose)
     {
        #message("Scale parameter in LP ",j," was missing and was set to ",LT$S0)
     }
  }

  # Improvement: Treat Scale as random, assign a gamma density 
  if(is.null(LT$shape0))
  {
     LT$shape0=1.1
  }

  if(is.null(LT$rate0))
  {    
     LT$rate0=(LT$shape0-1)/LT$S0
  }
  LT$S=LT$S0
  
  LT$b=rep(0,LT$p)
   
  LT$varB=rep(LT$S0/(LT$df0+2),LT$p)
  
  # Add one file when S0 is treated as random.
  fname=paste(saveAt,LT$Name,"_ScaleBayesA.dat",sep="") 
  if(rmExistingFiles)
  { 
    unlink(fname) 
  }

  LT$fileOut=file(description=fname,open="w")
  LT$NamefileOut=fname;
    
  LT$X=as.vector(LT$X)
  
  #Objects for storing information generated during MCMC iterations
  LT$post_varB=0
  LT$post_varB2=0
  
  LT$post_b=rep(0,LT$p)
  LT$post_b2=rep(0,LT$p)
  LT$post_S=0
  LT$post_S2=0

    #*#
    if(is.null(LT$saveEffects)){LT$saveEffects=FALSE}
    if(LT$saveEffects){
        if(is.null(LT$storageMode)){LT$storageMode="double"}
        if(!LT$storageMode%in%c("single","double")) {
            stop("storageMode of LP ",j," can either be 'single' or 'double' (default)")
        }
        if(is.null(LT$thin)){ LT$thin=thin }
        fname=paste(saveAt,LT$Name,"_b.bin",sep="")
        if(rmExistingFiles){ unlink(fname) }
        LT$fileEffects=file(fname,open='wb')
        nRow=floor((nIter-burnIn)/LT$thin)
        writeBin(object=c(nRow,LT$p),con=LT$fileEffects,size=ifelse(LT$storageMode=="single",4,8))
    }#*#
  
  #return object
  return(LT)
}

##################################################################################################
#Just the welcome function that will appear every time that your run the program
welcome=function()
{
#  message("\n");
#  message("#--------------------------------------------------------------------#");
#  message("#        _\\\\|//_                                                     #");
#  message("#       (` o-o ')      BGLR v1.1.0                                   #");
#  message("#------ooO-(_)-Ooo---------------------------------------------------#");
#  message("#                      Bayesian Generalized Linear Regression        #");
#  message("#                      Gustavo de los Campos, gdeloscampos@gmail.com #");
#  message("#    .oooO     Oooo.   Paulino Perez-Rodriguez, perpdgo@gmail.com    #");
#  message("#    (   )     (   )   April, 2022                                   #");
#  message("#_____\\ (_______) /_________________________________________________ #");
#  message("#      \\_)     (_/                                                   #");
#  message("#                                                                    #");
#  message("#------------------------------------------------------------------- #");
#  message("\n");
}
##################################################################################################

##################################################################################################
#The density of a scaled inverted chi-squered distribution
#df: degrees of freedom, S: Scale parameter
dScaledInvChisq=function (x, df, S)
{
    tmp = dchisq(S/x, df = df)/(x^2)
    return(tmp)
}


##################################################################################################
#The density function for lambda
#Density function for Regularization parameter in Bayesian LASSO
#Rate: rate parameter, shape: the value for the shape parameter
dLambda=function (rate, shape, lambda) 
{
    tmp = dgamma(x = I(lambda^2), rate = rate, shape = shape) * 2 * lambda
    return(tmp)
}

##################################################################################################
#Metropolis sampler for lambda in the Bayesian LASSO
metropLambda=function (tau2, lambda, shape1 = 1.2, shape2 = 1.2, max = 200, ncp = 0)
{
    lambda2 = lambda^2
    l2_new = rgamma(rate = sum(tau2)/2, shape = length(tau2),
        n = 1)
    l_new = sqrt(l2_new)
    logP_old = sum(dexp(x = tau2, log = TRUE, rate = (lambda2/2))) +
        dbeta(x = lambda/max, log = TRUE, shape1 = shape1, shape2 = shape2) -
        dgamma(shape = sum(tau2)/2, rate = length(tau2), x = (2/lambda2),
            log = TRUE)
    logP_new = sum(dexp(x = tau2, log = TRUE, rate = (l2_new/2))) +
        dbeta(x = l_new/max, log = TRUE, shape1 = shape1, shape2 = shape2) -
        dgamma(shape = sum(tau2)/2, rate = length(tau2), x = (2/l2_new),
            log = TRUE)
    accept = (logP_new - logP_old) > log(runif(1))
    if (accept) {
        lambda = l_new
    }
    return(lambda)
}

##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  if(interactive())
  {
    packageStartupMessage("# Gustavo de los Campos & Paulino Perez-Rodriguez")
    packageStartupMessage("# Support provided by the U.S., National Institutes of Health (NIH)")
    packageStartupMessage("# (Grant: R01GM101219, NIGMS)")
    packageStartupMessage("# and by the International Maize and Wheat Improvement Center (CIMMyT).")                        
    packageStartupMessage("# Type 'help(BGLR)' for summary information")
  }
  invisible()
}


##################################################################################################
#rtrun draws from a truncated univariate normal distribution using the inverse CDF algorithm
#Arguments:
#mu: mean
#sigma: standard deviation
#a: lower bound
#b: upper bound
#NOTES: 1) This routine was taken from bayesm package, December 18, 2012
#       2) The inputs are not checked, 
#It is assumed that are ok.
#rtrun=function (mu, sigma, a, b) 
#{
#    FA = pnorm(((a - mu)/sigma))
#    FB = pnorm(((b - mu)/sigma))
#    return(mu + sigma * qnorm(runif(length(mu)) * (FB - FA) + FA))
#}

#Using the rtruncnorm function in the truncnorm package
rtrun=function(mu,sigma,a,b)
{
        n=max(c(length(mu),length(sigma),length(a),length(b)))
    rtruncnorm(n,a,b,mu,sigma)
}


#Extract the values of z such that y[i]=j
#z,y vectors, j integer
#extract=function(z,y,j) subset(as.data.frame(z,y),subset=(y==j))
extract=function(z,y,j) z[y==j]

#This routine was adapted from rinvGauss function from S-Plus
# Random variates from inverse Gaussian distribution
# Reference:
#      Chhikara and Folks, The Inverse Gaussian Distribution,
#      Marcel Dekker, 1989, page 53.
# GKS  15 Jan 98
#n: Number of samples
#nu: nu parameter
#lambda: lambda parameter

rinvGauss=function(n, nu, lambda)
{
    if(any(nu<=0)) stop("nu must be positive")
    if(any(lambda<=0)) stop("lambda must be positive")
    if(length(n)>1) n = length(n)
    if(length(nu)>1 && length(nu)!=n) nu = rep(nu,length=n)
    if(length(lambda)>1 && length(lambda)!=n) lambda = rep(lambda,length=n)
        tmp = rnorm(n)
    y2 = tmp*tmp
    u = runif(n)
    r1 = nu/(2*lambda) * (2*lambda + nu*y2 - sqrt(4*lambda*nu*y2 + nu*nu*y2*y2))
    r2 = nu*nu/r1
    ifelse(u < nu/(nu+r1), r1, r2)
}


#log-likelihood for ordinal data
#y: response vector
#predicted response vector, yHat=X%*%beta
#threshold
loglik_ordinal=function(y,yHat,threshold)
{
    sum=0
        n=length(y)
    for(i in 1:n)
        {
           sum=sum + log(pnorm(threshold[y[i] + 1]-yHat[i])-pnorm(threshold[y[i]]-yHat[i]))
        }
        return(sum)
}

##################################################################################################

#Arguments:
#y: data vector, NAs allowed
#response_type: can be "gaussian", "ordinal",
#ETA: The linear predictor
#nIter: Number of MCMC iterations
#burnIn: burnIn
#thin: thin
#saveAt: string, where to save the information
#Se: Scale parameter for the prior for varE
#dfe: Degrees of freedom for the prior for varE
#weights: 
#R2
#Note: The function was designed to work with gaussian responses, some changes were made to deal binary and ordinal responses


#To add new method: 
#(a) create setLT, 
#(b) add it to the switch statement,
#(c) add code to update parameters in the Gibbs sampler, 
#(d) add code to save samples
#(e) add code to compute posterior means
#(f) Test:
#(f1) Test simple example without hyper-paramaeters, evaluate how  
#        default values were set
    #(f2)  Check posterior means and files
    #(f3)  Test simple example with some hyper-parameters give and 
    #         some set by default
#(f4) Run an example with a few missing values, compare to BLR 
#       example, check: (i) residual variance, (ii) plot of effects, (iii) plot 
#        of predictions in trn, (iv) plot of prediction in tst.


BGLR=function (y, response_type = "gaussian", a = NULL, b = NULL, 
    ETA = NULL, nIter = 1500, burnIn = 500, thin = 5, saveAt = "", 
    S0 = NULL, df0 = 5, R2 = 0.5, weights = NULL, 
    verbose = TRUE, rmExistingFiles = TRUE, groups=NULL) 
{
   
    if(verbose)
    {
    welcome()
    }

    IDs=names(y)
    if (!(response_type %in% c("gaussian", "ordinal")))  stop("Only gaussian and ordinal responses are allowed")

    if (saveAt == "") {
        saveAt = paste(getwd(), "/", sep = "")
    }

    y=as.vector(y)
    y0=y
    a = as.vector(a)
    b = as.vector(b)
    n = length(y)

    nGroups=1
    if(!is.null(groups))
    {
        groups<-as.character(groups)  #Groups as character and then as factor to avoid dummy levels
        groups<-as.factor(groups)
        #Number of records by group
        countGroups=table(groups)
        nGroups=length(countGroups)
        groupLabels=names(countGroups)
                groups=as.integer(groups)
                ggg=as.integer(groups-1);  #In C we begin to count in 0
        if(sum(countGroups)!=n) stop("length of groups and y differs, NA's not allowed in groups");
    }

    if(response_type=="ordinal")
    {

        y=factor(y,ordered=TRUE)
        lev=levels(y)
        nclass=length(lev)
        if(nclass==n) stop("The number of classes in y must be smaller than the number of observations");

        y=as.integer(y)
        z=y  

        fname = paste(saveAt, "thresholds.dat", sep = "")
        fileOutThresholds = file(description = fname, open = "w")
    }

    
    if (is.null(weights)) 
    {
        weights = rep(1, n)
    }

    if(!is.null(groups))
    {
      sumW2=tapply(weights^2,groups,"sum")
    }else{
      sumW2 = sum(weights^2)
    }

    nSums = 0

    whichNa = which(is.na(y))
    nNa = length(whichNa)

    Censored = FALSE
     
    if (response_type == "gaussian") 
    {
        if ((!is.null(a)) | (!is.null(b))) 
        {
            Censored = TRUE
            if ((length(a) != n) | (length(b) != n)) stop("y, a and b must have the same dimension")
            if (any(weights != 1)) stop("Weights are only implemented for Gausian uncensored responses")
        }
        mu = weighted.mean(x = y, w = weights, na.rm = TRUE)
    }
    post_mu = 0
    post_mu2 = 0

    fname = paste(saveAt, "mu.dat", sep = "")
    if (rmExistingFiles) 
    {
        unlink(fname)
    }
    else {
        #message("Note: samples will be appended to existing files")
    }

    fileOutMu = file(description = fname, open = "w")

    if (response_type == "ordinal") {
        if(verbose){message("Prior for residual is not necessary, if you provided it, it will be ignored")}
        if (any(weights != 1)) stop("Weights are not supported")
       
        countsZ=table(z)

        if (nclass <= 1) stop("Data vector y has only ", nclass, " differente values, it should have at least 2 different values")
        threshold=qnorm(p=c(0,cumsum(as.vector(countsZ)/n)))
          
        y = rtrun(mu =0, sigma = 1, a = threshold[z], b = threshold[ (z + 1)])

        mu=0
        #posterior for thresholds
        post_threshold = 0
        post_threshold2 = 0
        
    post_prob=matrix(nrow=n,ncol=nclass,0)
        post_prob2=post_prob
    }

    post_logLik = 0

    # yStar & yHat
    yStar = y * weights
    yHat = mu * weights
    
    if (nNa > 0) {
        yStar[whichNa] = yHat[whichNa]
    }

    post_yHat = rep(0, n)
    post_yHat2 = rep(0, n)

    # residual and residual variance
    e = (yStar - yHat)

    varE = var(e, na.rm = TRUE) * (1 - R2)

    if (is.null(S0)) {
        S0 = varE * (df0 + 2)
    }

    if(!is.null(groups))
    {
        varE=rep(varE/nGroups,nGroups)
        names(varE)=groupLabels
    }

    sdE = sqrt(varE)


    post_varE = 0
    post_varE2 = 0

    #File for storing sample for varE 

    fname = paste(saveAt, "varE.dat", sep = "")

    if (rmExistingFiles) {
        unlink(fname)
    }

    fileOutVarE = file(description = fname, open = "w")

    nLT = ifelse(is.null(ETA), 0, length(ETA))
    

    #Setting the linear terms
    if (nLT > 0) {
    
    if(is.null(names(ETA)))
        { 
             names(ETA)<-rep("",nLT)
        }

        for (i in 1:nLT) {  

        if(names(ETA)[i]=="")
        {
            ETA[[i]]$Name=paste("ETA_",i,sep="")
        }else{
               ETA[[i]]$Name=paste("ETA_",names(ETA)[i],sep="")
        }

            if (!(ETA[[i]]$model %in% c("FIXED", "BRR", "BL", "BayesA", "BayesB","BayesC", "RKHS","BRR_sets"))) 
            {
                stop("Error in ETA[[", i, "]]", " model ", ETA[[i]]$model, " not implemented (note: evaluation is case sensitive)")
                
            }

            if(!is.null(groups))
            {
        if(!(ETA[[i]]$model %in%  c("BRR","FIXED","BayesB","BayesC"))) stop("Error in ETA[[", i, "]]", " model ", ETA[[i]]$model, " not implemented for groups")
            }

            ETA[[i]] = switch(ETA[[i]]$model, 
                  FIXED = setLT.Fixed(LT = ETA[[i]],  n = n, j = i, weights = weights, y = y, nLT = nLT, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups), 
                              BRR = setLT.BRR(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn,lower_tri=ETA[[i]]$lower_tri),#*#
                              BL = setLT.BL(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn), 
                              RKHS = setLT.RKHS(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose), 
                              BayesC = setLT.BayesBandC(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn),
                              BayesA = setLT.BayesA(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn),
                              BayesB = setLT.BayesBandC(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,groups=groups,nGroups=nGroups,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn),
                              BRR_sets = setLT.BRR_sets(LT = ETA[[i]], n = n, j = i, weights = weights, y = y, nLT = nLT, R2 = R2, saveAt = saveAt, rmExistingFiles = rmExistingFiles,verbose=verbose,thin=thin,nIter=nIter,burnIn=burnIn)
                              )
        }
    }

    # Gibbs sampler

    time = proc.time()[3]

    for (i in 1:nIter) {
        # intercept
    if(!is.null(groups))
    {
        e = e + weights * mu 
                varEexpanded=varE[groups]
        #rhs = sum(tapply(e*weights,groups,"sum")/varE)
                rhs = as.numeric(crossprod(e/varEexpanded,weights));
                C = sum(sumW2/varE)
                sol = rhs/C
                mu = rnorm(n = 1, sd = sqrt(1/C)) + sol;
    }else{
            e = e + weights * mu
            rhs = sum(weights * e)/varE
            C = sumW2/varE
            sol = rhs/C
            mu = rnorm(n = 1, sd = sqrt(1/C)) + sol
    }
        if (response_type == "ordinal") {
            mu=0 
        }

        e = e - weights * mu
        
        #deltaSS and deltadf for updating varE
        deltaSS = 0
        deltadf = 0

        if (nLT > 0) {
            for (j in 1:nLT) {
                ## Fixed effects ####################################################################
                if (ETA[[j]]$model == "FIXED") {
                  #cat("varB=",ETA[[j]]$varB,"\n");
                  varBj = rep(ETA[[j]]$varB, ETA[[j]]$p)
                  if(!is.null(groups)){
                        ans = .Call("sample_beta_groups", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b,
                                             e, varBj, varE, 1e-9,ggg,nGroups)
          }else{
                    ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                         e, varBj, varE, 1e-9)
          }
                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]
                }#End of fixed effects

                ## Ridge Regression ##################################################################
                if (ETA[[j]]$model == "BRR") {
                  varBj = rep(ETA[[j]]$varB, ETA[[j]]$p)

                  if(!is.null(groups))
          {
                        ans = .Call("sample_beta_groups",n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                    e, varBj, varE, 1e-9,ggg,nGroups)
              }else{
            if(!(ETA[[j]]$lower_tri))
            {
                        ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                                 e, varBj, varE, 1e-9)
            }else{
                ans = .Call("sample_beta_lower_tri", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b,
                                                     e, ETA[[j]]$varB, varE, 1e-9)
            }
          }
                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]

                  DF = ETA[[j]]$df0 + ETA[[j]]$p
                  SS = sum(ETA[[j]]$b^2) + ETA[[j]]$S0
                  ETA[[j]]$varB = SS/rchisq(df = DF, n = 1)
                }# END BRR
                
                if(ETA[[j]]$model=="BRR_sets"){
                    ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b,
                                             e, ETA[[j]]$varB, varE, 1e-9)
                   ETA[[j]]$b = ans[[1]]
                    e = ans[[2]]
            SS=tapply(X=ETA[[j]]$b^2,INDEX=ETA[[j]]$sets,FUN=sum)+ETA[[j]]$S0
                   
                    ETA[[j]]$varSets=SS/rchisq(df=ETA[[j]]$DF1,n=ETA[[j]]$n_sets)
                    ETA[[j]]$varB=ETA[[j]]$varSets[ETA[[j]]$sets]
                }


                ## Bayesian LASSO ####################################################################
                if (ETA[[j]]$model == "BL") {
                  
                   varBj = ETA[[j]]$tau2 * varE
                   ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                e, varBj, varE, ETA[[j]]$minAbsBeta)

                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]

                  nu = sqrt(varE) * ETA[[j]]$lambda/abs(ETA[[j]]$b)
                  tmp = NULL
                  try(tmp <- rinvGauss(n = ETA[[j]]$p, nu = nu, lambda = ETA[[j]]$lambda2))
                  if (!is.null(tmp) && !any(tmp<0)) {
                    if (!any(is.na(sqrt(tmp)))) {
                      ETA[[j]]$tau2 = 1/tmp
                    }
                    else {
                      warning(paste("tau2 was not updated in iteration",i, "due to numeric problems with beta\n",sep=" "),immediate. = TRUE)
                    }
                  }
                  else {
                    warning(paste("tau2 was not updated  in iteration",i,"due to numeric problems with beta\n",sep=" "),immediate. = TRUE)
                  }

                  #Update lambda 
                  if (ETA[[j]]$type == "gamma") {
                    rate = sum(ETA[[j]]$tau2)/2 + ETA[[j]]$rate
                    shape = ETA[[j]]$p + ETA[[j]]$shape
                    ETA[[j]]$lambda2 = rgamma(rate = rate, shape = shape, n = 1)
                    if (!is.na(ETA[[j]]$lambda2)) {
                      ETA[[j]]$lambda = sqrt(ETA[[j]]$lambda2)
                    }
                    else {
                      warning(paste("lambda was not updated in iteration",i, "due to numeric problems with beta\n",sep=" "),immediate. = TRUE)
                    }
                  }

                  if (ETA[[j]]$type == "beta") {
                    ETA[[j]]$lambda = metropLambda(tau2 = ETA[[j]]$tau2, 
                                                   lambda = ETA[[j]]$lambda, shape1 = ETA[[j]]$shape1, shape2 = ETA[[j]]$shape2, 
                                                   max = ETA[[j]]$max)
                    ETA[[j]]$lambda2 = ETA[[j]]$lambda^2
                  }

                  deltaSS = deltaSS + sum((ETA[[j]]$b/sqrt(ETA[[j]]$tau2))^2)
                  deltadf = deltadf + ETA[[j]]$p
                }#END BL

                ## RKHS ####################################################################
                if (ETA[[j]]$model == "RKHS") {
                  #error
                  e = e + ETA[[j]]$u
                  rhs = crossprod(ETA[[j]]$V, e)/varE
                  varU = ETA[[j]]$varU * ETA[[j]]$d
                  C = as.numeric(1/varU + 1/varE)
                  SD = 1/sqrt(C)
                  sol = rhs/C
                  tmp = rnorm(n = ETA[[j]]$levelsU, mean = sol, sd = SD)
                  ETA[[j]]$uStar = tmp
                  ETA[[j]]$u = as.vector(ETA[[j]]$V %*% tmp)
          
                  #update error
                  e = e - ETA[[j]]$u
                   
                  #update the variance
                  tmp = ETA[[j]]$uStar/sqrt(ETA[[j]]$d)
                  SS = as.numeric(crossprod(tmp)) + ETA[[j]]$S0
                  DF = ETA[[j]]$levelsU + ETA[[j]]$df0
                  ETA[[j]]$varU = SS/rchisq(n = 1, df = DF)
                }#END RKHS

                ## BayesA ##############################################################################
                if (ETA[[j]]$model == "BayesA") {
                  varBj = ETA[[j]]$varB
                  ans = .Call("sample_beta", n, ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, 
                                             e, varBj, varE, 1e-9)
                  ETA[[j]]$b = ans[[1]]
                  e = ans[[2]]
                   
                  #Update variances
                  SS = ETA[[j]]$S + ETA[[j]]$b^2
                  DF = ETA[[j]]$df0 + 1
                  ETA[[j]]$varB = SS/rchisq(n = ETA[[j]]$p, df = DF)

                  tmpShape=ETA[[j]]$p*ETA[[j]]$df0/2+ETA[[j]]$shape0
                  tmpRate=sum(1/ETA[[j]]$varB)/2+ETA[[j]]$rate0
                  ETA[[j]]$S=rgamma(shape=tmpShape,rate=tmpRate,n=1)

                }#End BayesA

        #BayesB and BayesC
        if(ETA[[j]]$model %in% c("BayesB","BayesC"))
        {
            #Update marker effects
                        mrkIn=ETA[[j]]$d==1
                        pIn=sum(mrkIn)
                
                        if(ETA[[j]]$model=="BayesB")
                        {
                          if(!is.null(groups))
                          {
                             ans=.Call("sample_beta_BB_BCp_groups",n,ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, ETA[[j]]$d, e, ETA[[j]]$varB, varE, 1e-9, ETA[[j]]$probIn,ggg,nGroups);
                          }else{
                             ans=.Call("sample_beta_BB_BCp",n,ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, ETA[[j]]$d, e, ETA[[j]]$varB, varE, 1e-9, ETA[[j]]$probIn);
                          }
                        }else{
                          if(!is.null(groups))
                          {
                             ans=.Call("sample_beta_BB_BCp_groups",n,ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, ETA[[j]]$d, e, rep(ETA[[j]]$varB,ETA[[j]]$p), varE, 1e-9, ETA[[j]]$probIn,ggg,nGroups);
                          }else{   
                             ans=.Call("sample_beta_BB_BCp",n,ETA[[j]]$p, ETA[[j]]$X, ETA[[j]]$x2, ETA[[j]]$b, ETA[[j]]$d, e, rep(ETA[[j]]$varB,ETA[[j]]$p), varE, 1e-9, ETA[[j]]$probIn);
                          }
                        }

                        ETA[[j]]$d=ans[[1]]
                        e=ans[[2]]
                        ETA[[j]]$b=ans[[3]] 

            #Update the variance component associated with the markers
            if(ETA[[j]]$model=="BayesB")
            {
                SS = ETA[[j]]$b^2 + ETA[[j]]$S
                            DF = ETA[[j]]$df0+1
                            ETA[[j]]$varB = SS/rchisq(df=DF, n = ETA[[j]]$p)

                                # Update scale
                                tmpShape=ETA[[j]]$p*ETA[[j]]$df0/2+ETA[[j]]$shape0
                                tmpRate=sum(1/ETA[[j]]$varB)/2+ETA[[j]]$rate0
                                ETA[[j]]$S=rgamma(shape=tmpShape,rate=tmpRate,n=1)

            }else{
                SS = sum(ETA[[j]]$b^2) + ETA[[j]]$S0
                DF = ETA[[j]]$df0 + ETA[[j]]$p
                        ETA[[j]]$varB = SS/rchisq(df = DF, n = 1)
            }
                        mrkIn = sum(ETA[[j]]$d)
                        ETA[[j]]$probIn = rbeta(shape1 = (mrkIn + ETA[[j]]$countsIn + 1),
                                                shape2 = (ETA[[j]]$p - mrkIn + ETA[[j]]$countsOut + 1), n = 1)
        }
            }#Loop for
        }#nLT
        
        # yHat
        yHat = yStar - e
        
       #4#
         # residual variance # missing values
        if (response_type == "gaussian") {
        
            if(!is.null(groups))
        {   
            for(g in 1:nGroups)
            {
            SS=sum(e[groups==g]^2)+ S0 + deltaSS
                        DF=countGroups[g]+df0+deltadf
            varE[g]=SS/rchisq(n=1,df=DF)
        }
         }else{
                    SS = sum(e * e) + S0 + deltaSS
                    DF = n + df0 + deltadf
                    varE = SS/rchisq(n = 1, df = DF)
        }
            sdE = sqrt(varE)
        
          if (nNa > 0) {
            if (Censored) {
                if(!is.null(groups))
                {
                   #FIXME: Double check this, I was testing it and is ok
                   sdEexpanded=sdE[groups]
                   yStar[whichNa] = rtrun(mu = yHat[whichNa], a = a[whichNa], b = b[whichNa], sigma = sdEexpanded)

                }else{
                  yStar[whichNa] = rtrun(mu = yHat[whichNa], a = a[whichNa], b = b[whichNa], sigma = sdE)
                }
            }
            else{
                 if(!is.null(groups))
                 {
                    #FIXME: Double check this, I was testing it and is ok
                    sdEexpanded=sdE[groups]
                    yStar[whichNa] = yHat[whichNa] + rnorm(n = nNa, sd = sdEexpanded)
                 }else{
                    yStar[whichNa] = yHat[whichNa] + rnorm(n = nNa, sd = sdE)
                 }
            }
            e[whichNa] = yStar[whichNa] - yHat[whichNa]
          }
        }else{  #ordinal
            varE = 1
            sdE = 1
            
            #Update yStar, this is the latent variable
            if(nNa==0){
               yStar=rtrun(mu = yHat, sigma = 1, a = threshold[z], b = threshold[(z + 1)])
            }else{
               yStar[-whichNa]=rtrun(mu = yHat[-whichNa], sigma = 1, a = threshold[z[-whichNa]], b = threshold[(z[-whichNa] + 1)])
               yStar[whichNa]=yHat[whichNa] + rnorm(n = nNa, sd = sdE)           
            }

            #Update thresholds           
            if(nNa==0){ 
              for (m in 2:nclass) {
            
                lo = max(max(extract(yStar, z, m - 1)), threshold[m - 1])
                hi = min(min(extract(yStar, z, m)), threshold[m + 1])
                threshold[m] = runif(1, lo, hi)
              }
            }else{

              for (m in 2:nclass) {
                tmpY=yStar[-whichNa]
                tmpZ=z[-whichNa]
                lo = max(max(extract(tmpY, tmpZ, m - 1)), threshold[m - 1])
                hi = min(min(extract(tmpY, tmpZ, m)), threshold[m + 1])
                threshold[m] = runif(1, lo, hi)
              }
            }
            
            #Update error
            e = yStar - yHat
        }

        # Saving samples and computing running means
        if ((i%%thin == 0)) {
            if (nLT > 0) {
                for (j in 1:nLT) {

                  if (ETA[[j]]$model == "FIXED") {
                    write(ETA[[j]]$b,ncolumns=ETA[[j]]$p, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BRR") {
                    write(ETA[[j]]$varB, file = ETA[[j]]$fileOut, append = TRUE)
                  }
                  if (ETA[[j]]$model == "BRR_sets") {
                    write(ETA[[j]]$varSets, ncolumns=ETA[[j]]$n_sets,file = ETA[[j]]$fileOut, append = TRUE)
                  }
                  if (ETA[[j]]$model == "BL") {
                    write(ETA[[j]]$lambda, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "RKHS") {
                    write(ETA[[j]]$varU, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BayesC") {
                    tmp = c(ETA[[j]]$probIn, ETA[[j]]$varB)
                    write(tmp, ncolumns = 2, file = ETA[[j]]$fileOut, append = TRUE)
                  }

                  if (ETA[[j]]$model == "BayesA") {
                    tmp=ETA[[j]]$S
                    write(tmp, ncolumns = 1, file = ETA[[j]]$fileOut, append = TRUE)
                  }
                  
                  if(ETA[[j]]$model=="BayesB")
                  {
                        tmp=c(ETA[[j]]$probIn,ETA[[j]]$S)
                        write(tmp, ncolumns = 2, file = ETA[[j]]$fileOut, append = TRUE)
                  }
                }
            }

            #Output files
            write(x = mu, file = fileOutMu, append = TRUE)
            write(x = varE, ncolumns=nGroups,file = fileOutVarE, append = TRUE)

        if (response_type == "ordinal") {
          write(x=threshold[2:nclass],ncolumns=nclass-1,file=fileOutThresholds,append=TRUE)
        }


            if (i > burnIn) {
                nSums = nSums + 1
                k = (nSums - 1)/(nSums)
                if (nLT > 0) {
                  for (j in 1:nLT) {
                    if (ETA[[j]]$model == "FIXED") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                    }

                    if (ETA[[j]]$model == "BRR") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                          writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                      }#*#
                    }

                    if (ETA[[j]]$model == "BRR_sets") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      ETA[[j]]$post_varSets<-ETA[[j]]$post_varSets*k+ETA[[j]]$varSets/nSums
                      ETA[[j]]$post_varSets2<-ETA[[j]]$post_varSets2*k+(ETA[[j]]$varSets^2)/nSums
                      if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                          writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                      }#*#
                    }
                    
                    if (ETA[[j]]$model == "BL") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_tau2 = ETA[[j]]$post_tau2 * k + (ETA[[j]]$tau2)/nSums
                      ETA[[j]]$post_lambda = ETA[[j]]$post_lambda * k + (ETA[[j]]$lambda)/nSums
                      if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                          writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                      }#*#
                    }

                    if (ETA[[j]]$model == "RKHS") {
                      ETA[[j]]$post_varU = ETA[[j]]$post_varU * k + ETA[[j]]$varU/nSums
                      ETA[[j]]$post_varU2 = ETA[[j]]$post_varU2 * k + (ETA[[j]]$varU^2)/nSums
                      ETA[[j]]$post_uStar = ETA[[j]]$post_uStar * k + ETA[[j]]$uStar/nSums
                      ETA[[j]]$post_u = ETA[[j]]$post_u * k + ETA[[j]]$u/nSums
                      ETA[[j]]$post_u2 = ETA[[j]]$post_u2 * k + (ETA[[j]]$u^2)/nSums
                    }

                    if (ETA[[j]]$model == "BayesC") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      ETA[[j]]$post_d = ETA[[j]]$post_d * k + (ETA[[j]]$d)/nSums
                      ETA[[j]]$post_probIn = ETA[[j]]$post_probIn * k + (ETA[[j]]$probIn)/nSums
                      ETA[[j]]$post_probIn2 = ETA[[j]]$post_probIn2 * k + (ETA[[j]]$probIn^2)/nSums
                      if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                          writeBin(object=ETA[[j]]$b*ETA[[j]]$d,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                      }#*#
                    }

                    if (ETA[[j]]$model == "BayesA") {
                      ETA[[j]]$post_b = ETA[[j]]$post_b * k + ETA[[j]]$b/nSums
                      ETA[[j]]$post_b2 = ETA[[j]]$post_b2 * k + (ETA[[j]]$b^2)/nSums
                      ETA[[j]]$post_varB = ETA[[j]]$post_varB * k + (ETA[[j]]$varB)/nSums
                      ETA[[j]]$post_varB2 = ETA[[j]]$post_varB2 * k + (ETA[[j]]$varB^2)/nSums
                      ETA[[j]]$post_S = ETA[[j]]$post_S * k + (ETA[[j]]$S)/nSums
                      ETA[[j]]$post_S2 = ETA[[j]]$post_S2 * k + (ETA[[j]]$S^2)/nSums
                      if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                          writeBin(object=ETA[[j]]$b,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                      }#*#
                    }

                    if(ETA[[j]]$model=="BayesB")
                    {
                        ETA[[j]]$post_b=ETA[[j]]$post_b*k+ETA[[j]]$b/nSums
                        ETA[[j]]$post_b2=ETA[[j]]$post_b2*k+(ETA[[j]]$b^2)/nSums
                        ETA[[j]]$post_varB=ETA[[j]]$post_varB*k+(ETA[[j]]$varB)/nSums
                        ETA[[j]]$post_varB2=ETA[[j]]$post_varB2*k+(ETA[[j]]$varB^2)/nSums
                        ETA[[j]]$post_d = ETA[[j]]$post_d * k + (ETA[[j]]$d)/nSums
                        ETA[[j]]$post_probIn = ETA[[j]]$post_probIn * k + (ETA[[j]]$probIn)/nSums
                        ETA[[j]]$post_probIn2 = ETA[[j]]$post_probIn2 * k + (ETA[[j]]$probIn^2)/nSums
                        ETA[[j]]$post_S = ETA[[j]]$post_S * k + (ETA[[j]]$S)/nSums
                        ETA[[j]]$post_S2 = ETA[[j]]$post_S2 * k + (ETA[[j]]$S^2)/nSums
                        if(ETA[[j]]$saveEffects&&(i%%ETA[[j]]$thin)==0){
                            writeBin(object=ETA[[j]]$b*ETA[[j]]$d,con=ETA[[j]]$fileEffects,size=ifelse(ETA[[j]]$storageMode=="single",4,8))
                        }#*#
                    }
                  }
                }

                post_mu = post_mu * k + mu/nSums
                post_mu2 = post_mu2 * k + (mu^2)/nSums

                post_yHat = post_yHat * k + yHat/nSums
                post_yHat2 = post_yHat2 * k + (yHat^2)/nSums

                post_varE = post_varE * k + varE/nSums
                post_varE2 = post_varE2 * k + (varE^2)/nSums

                if (response_type == "ordinal") {
                  post_threshold = post_threshold * k + threshold/nSums
                  post_threshold2 = post_threshold2 * k + (threshold^2)/nSums

                  TMP=matrix(nrow=n,ncol=nclass,0)

                  TMP[,1]=pnorm(threshold[2]-yHat)

                  if(nclass>2){
                     for(m in 2:(nclass-1)){
                       TMP[,m]=pnorm(threshold[(m+1)]-yHat)-rowSums(as.matrix(TMP[,1:(m-1)]))
                     }
                  }
                  TMP[,nclass]=1-rowSums(TMP)

                  post_prob=post_prob*k+TMP/nSums
                  post_prob2=post_prob2*k+(TMP^2)/nSums
                 
                  if(nNa==0){
                    logLik=loglik_ordinal(z,yHat,threshold)
                  }else{
                    logLik=loglik_ordinal(z[-whichNa],yHat[-whichNa],threshold)
                  }
                }

                if(response_type == "gaussian") {
                  
                  tmpE = e/weights
                  if(!is.null(groups)){
                    tmpSD=rep(NA,n)
                    for(g in 1:nGroups)
                    {
                       index=(groups==g)
                       tmpSD[index]=sqrt(varE[g])/weights[index]
                    }
                  }else{
                    tmpSD = sqrt(varE)/weights
                  }

                  if (nNa > 0) {
                    tmpE = tmpE[-whichNa]
                    tmpSD = tmpSD[-whichNa]
                  }
                  
              logLik = sum(dnorm(tmpE, sd = tmpSD, log = TRUE))

                }#end gaussian

                post_logLik = post_logLik * k + logLik/nSums
            }
        }#end of saving samples and computing running means

        if (verbose) {
            #message("---------------------------------------")
            tmp = proc.time()[3]
            #message("Iter=",i," Time/Iter=",round(tmp-time,3))
            #message("VarE=",round(varE,3))
        #In the case of variance by groups
        #message("varE=",paste(round(varE,3),collapse=", "))
            time = tmp
        }
    }#end of Gibbs sampler

    #Closing files
    close(fileOutVarE)
    close(fileOutMu)
    if(response_type == "ordinal") close(fileOutThresholds)

    if (nLT > 0) {
        for (i in 1:nLT) {
            if (!is.null(ETA[[i]]$fileOut)) {
                flush(ETA[[i]]$fileOut)                
                close(ETA[[i]]$fileOut)
                ETA[[i]]$fileOut = NULL
            }
            if(!is.null(ETA[[i]]$fileEffects)){
                flush(ETA[[i]]$fileEffects)
                close(ETA[[i]]$fileEffects)
                if(!is.null(ETA[[i]]$compressEffects)&&ETA[[i]]$compressEffects==TRUE){
                    compressFile(paste0(saveAt,ETA[[i]]$Name,"_b.bin"))
                }
                ETA[[i]]$fileEffects = NULL
            }
            
        }
    }
    
    #return goodies

    out = list(y = y0, a=a,b=b,whichNa = whichNa, saveAt = saveAt, nIter = nIter, 
               burnIn = burnIn, thin = thin, 
               weights = weights, verbose = verbose, 
               response_type = response_type, df0 = df0, S0 = S0)

    out$yHat = post_yHat

    names(out$yHat)=IDs
    names(out$y)=IDs

    out$SD.yHat = sqrt(post_yHat2 - (post_yHat^2))
    out$mu = post_mu
    out$SD.mu = sqrt(post_mu2 - post_mu^2)
    out$varE = post_varE
    out$SD.varE = sqrt(post_varE2 - post_varE^2)
    
    #goodness of fit 
    out$fit = list()
    
    if(response_type=="gaussian")
    {
        tmpE = (yStar - post_yHat)/weights

        if(!is.null(groups))
        {
                    tmpSD=rep(NA,n)
                    for(g in 1:nGroups)
                    {
                       index=(groups==g)
                       tmpSD[index]=sqrt(varE[g])/weights[index]
                    }
         }else{
            tmpSD = sqrt(post_varE)/weights
     }
    
        if (nNa > 0) {
            tmpE = tmpE[-whichNa]
            tmpSD = tmpSD[-whichNa]
        }
        out$fit$logLikAtPostMean = sum(dnorm(tmpE, sd = tmpSD, log = TRUE))

    if (Censored) {
            cdfA = pnorm(q = a[whichNa], sd = sqrt(post_varE), mean = post_yHat[whichNa])
            cdfB = pnorm(q = b[whichNa], sd = sqrt(post_varE), mean = post_yHat[whichNa])
            out$fit$logLikAtPostMean = out$fit$logLikAtPostMean + sum(log(cdfB - cdfA))
        }
    }

    if(response_type=="ordinal")
    {
         out$probs=post_prob
         out$SD.probs=sqrt(post_prob2-post_prob^2)
         colnames(out$probs)=lev
         colnames(out$SD.probs)=lev
         out$threshold = post_threshold[-c(1, nclass + 1)]
         out$SD.threshold = sqrt(post_threshold2 - post_threshold^2)[-c(1, nclass + 1)] 
         
         #out$fit$logLikAtPostMean = loglik_ordinal(y,post_yHat,post_threshold)#*#
        tmp=0
         for(i in 1:nclass){
            tmp=tmp+sum(ifelse(y0==lev[i],log(out$probs[,i]),0))
         }
         
         out$fit$logLikAtPostMean=tmp
         out$levels=lev
         out$nlevels=nclass
    }

    out$fit$postMeanLogLik = post_logLik
    out$fit$pD = -2 * (post_logLik - out$fit$logLikAtPostMean)
    out$fit$DIC = out$fit$pD - 2 * post_logLik

    # Renaming/removing objects in ETA and appending names
    if (nLT > 0) {
        for (i in 1:nLT) {

            if (ETA[[i]]$model != "RKHS") {
                ETA[[i]]$b = ETA[[i]]$post_b
                ETA[[i]]$SD.b = sqrt(ETA[[i]]$post_b2 - ETA[[i]]$post_b^2)
                names(ETA[[i]]$b)=ETA[[i]]$colNames
                names(ETA[[i]]$SD.b)=ETA[[i]]$colNames
                tmp = which(names(ETA[[i]]) %in% c("post_b", "post_b2","X","x2"))
                ETA[[i]] = ETA[[i]][-tmp]
            }
            
            if(ETA[[i]]$model=="RKHS")
            {
               ETA[[i]]$SD.u=sqrt(ETA[[i]]$post_u2 - ETA[[i]]$post_u^2)
               ETA[[i]]$u=ETA[[i]]$post_u
               ETA[[i]]$uStar=ETA[[i]]$post_uStar
               ETA[[i]]$varU=ETA[[i]]$post_varU
               ETA[[i]]$SD.varU=sqrt(ETA[[i]]$post_varU2 - ETA[[i]]$post_varU^2)
               tmp=which(names(ETA[[i]])%in%c("post_varU","post_varU2","post_uStar","post_u","post_u2"))
               ETA[[i]]=ETA[[i]][-tmp]
            }

            if (ETA[[i]]$model %in% c("BRR","BRR_sets", "BayesA", "BayesC","BayesB")) {
                ETA[[i]]$varB = ETA[[i]]$post_varB
                ETA[[i]]$SD.varB = sqrt(ETA[[i]]$post_varB2 - (ETA[[i]]$post_varB^2))
                tmp = which(names(ETA[[i]]) %in% c("post_varB", "post_varB2"))
                ETA[[i]] = ETA[[i]][-tmp]
            }
            if(ETA[[i]]$model=="BRR_sets"){
                ETA[[i]]$varSets=ETA[[i]]$post_varSets
                ETA[[i]]$SD.varSets=sqrt(ETA[[i]]$post_varSets2-(ETA[[i]]$post_varSets^2))
                tmp<-which(names(ETA[[i]])%in%c("post_varSets","post_varSets2"))
                ETA[[i]]=ETA[[i]][-tmp]
            }

            if(ETA[[i]]$model %in% c("BayesB","BayesC"))
            {
                ETA[[i]]$d=ETA[[i]]$post_d
                ETA[[i]]$probIn=ETA[[i]]$post_probIn
                ETA[[i]]$SD.probIn=sqrt(ETA[[i]]$post_probIn2 - (ETA[[i]]$post_probIn^2))
                tmp = which(names(ETA[[i]]) %in% c("post_d", "post_probIn","post_probIn2"))
                ETA[[i]] = ETA[[i]][-tmp]
            }
            
            if(ETA[[i]]$model %in% c("BayesA","BayesB"))
            {
                ETA[[i]]$S=ETA[[i]]$post_S
                ETA[[i]]$SD.S=sqrt( ETA[[i]]$post_S2 - (ETA[[i]]$post_S^2))
                tmp=which(names(ETA[[i]])%in%c("post_S","post_S2"))
                ETA[[i]]=ETA[[i]][-tmp]
            }

            if(ETA[[i]]$model=="BL")
            {
                ETA[[i]]$tau2=ETA[[i]]$post_tau2
                ETA[[i]]$lambda=ETA[[i]]$post_lambda
                tmp = which(names(ETA[[i]]) %in% c("post_tau2", "post_lambda","lambda2"))
                ETA[[i]] = ETA[[i]][-tmp]
            }
        }
        out$ETA = ETA
    }
    class(out) = "BGLR"
    return(out)
}

#This function will be a wrapper for BGLR
#the idea is to maintain the compatibility with the function BLR in 
#the package BLR that was released in 2010, updated in 2011 and 2012

#NOTE: thin2 parameter is missing in BGLR, so it will be removed

BLR=function (y, XF = NULL, XR = NULL, XL = NULL, GF = list(ID = NULL, 
    A = NULL), prior = NULL, nIter = 1100, burnIn = 100, thin = 10, 
    thin2 = 1e+10, saveAt = "", minAbsBeta = 1e-09, weights = NULL) 
{

    ETA = NULL
    ETA = list()
    nLT = 0

    #message("This implementation is a simplified interface for the more general")
    #message("function BGLR, we keep it for backward compatibility with our package BLR")

    warning("thin2 parameter is not used any more and will be deleted in next releases\n",immediate. = TRUE);
   
    #message("Setting parameters for BGLR...")
    if (is.null(prior)) {
        #message("===============================================================")
        #message("No prior was provided, BGLR will be running with improper priors.")
        #message("===============================================================")
        prior = list(varE = list(S = NULL, df = 1), varBR = list(S = 0, 
            df = 0), varU = list(S = 0, df = 0), lambda = list(shape = 0, 
            rate = 0, type = "random", value = 50))
    }
    if (!is.null(XF)) {
        nLT = nLT + 1
        ETA[[nLT]] = list(X = XF, model = "FIXED")
    }
    if (!is.null(XR)) {
        nLT = nLT + 1
        ETA[[nLT]] = list(X = XR, model = "BRR", df0 = prior$varBR$df, 
            S0 = prior$varBR$S)
    }
    if (!is.null(XL)) {
        nLT = nLT + 1
        if (prior$lambda$type == "random") {
            if (is.null(prior$lambda$rate)) {
                #message("Setting prior for lambda^2 to beta")
                prior$lambda$type = "beta"
                prior$lambda$shape = NULL
                prior$lambda$rate = NULL
            }
            else {
                #message("Setting prior for lambda^2 to gamma")
                prior$lambda$type = "gamma"
                prior$lambda$max = NULL
                prior$lambda$shape1 = NULL
                prior$lambda$shape2 = NULL
            }
        }
        ETA[[nLT]] = list(X = XL, model = "BL", type = prior$lambda$type, 
                          rate = prior$lambda$rate, shape = prior$lambda$shape, 
                          max = prior$lambda$max, shape1 = prior$lambda$shape1, 
                          shape2 = prior$lambda$shape2, lambda = prior$lambda$value,
                          minAbsBeta=minAbsBeta)
    }

    #NOTE: In original BLR IDS are used to buid A matrix Z,
    #and then run the model y=Zu+e, u~MN(0,varU*A), and then using the Cholesky factorization
    #it was possible to fit the model. The algorithm used here is different (Orthogonal variables)
    #and it may be that the IDS are not longer necessary

    if (!is.null(GF[[1]])) {
        nLT = nLT + 1
        ETA[[nLT]] = list(K = GF$A, model = "RKHS", df0 = prior$varU$df, 
            S0 = prior$varU$S)
        warning("IDs are not used any more and will be deleted in next releases...\n",immediate. = TRUE) 
    }

    #message("Finish setting parameters for BGLR")
    #message("Fitting model using BGLR...")
    out = BGLR(y = y, ETA = ETA, df0 = prior$varE$df, S0 = prior$varE$S, 
               nIter = nIter, burnIn = burnIn, thin = thin, saveAt = saveAt, 
               weights = weights)

    #Backward compatibility with BLR
    if (nLT > 0) {
        for (j in 1:nLT) {
            if (ETA[[j]]$model == "FIXED") {
                out$bF = out$ETA[[j]]$b
                out$SD.bF = out$ETA[[j]]$SD.b
            }
            if (ETA[[j]]$model == "BL") {
                out$bL = out$ETA[[j]]$b
                out$SD.bL = out$ETA[[j]]$SD.b
                out$lambda = out$ETA[[j]]$lambda
            }
            if (ETA[[j]]$model == "BRR") {
                out$bR = out$ETA[[j]]$b
                out$SD.bR = out$ETA[[j]]$SD.b
                out$varBR = out$ETA[[j]]$varB
                out$SD.bR = out$ETA[[j]]$SD.varB
            }
            if (ETA[[j]]$model == "RKHS") {
                out$u = out$ETA[[j]]$u
                out$SD.u =out$ETA[[j]]$SD.u
                out$varU = out$ETA[[j]]$varU
            }
        }
    }
    out$ETA = NULL
    class(out) = "BLR"
    return(out)
}
