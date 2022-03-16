`G_E.Simulation` <-
function(GD,ha2,rg,re,NE,NQTN=NULL,NETN=NULL,effectunit=1,file.output=NULL){
    

    #############################
    #GD
    
    #############################
    if(is.null(NQTN))
    {print("Must give NQTN")}
    if(length(ha2)!=NE) 
    {
        #print("simulation paremeter error in h2")
        stop ("simulation paremeter error in h2")
    }
    if(is.null(file.output))
    {file.output=FALSE
     print("This project does not output phenotype files because file.output is null")
     print("The default is FALSE")
    }

    Y=matrix(0,nrow(GD),NE)
    X=GD[,-1]
    taxa=as.character(GD[,1])


    m=ncol(X)
    n=nrow(X)
    #duplicate taxa
    cov_g=rg
    mu_g=array(0,NE)
    if(!is.matrix(rg)) as.matrix(rg)
    matrix_cov_g=rg





    ####### simulate QTN effect and variance in each environment
    if(matrix_cov_g[1,2]==matrix_cov_g[2,1]&matrix_cov_g[1,2]!=1)
    {
     gv<-mvrnorm(NQTN,t(mu_g),matrix_cov_g)
    }else{
        gv=matrix(0,NQTN,NE)
        for(i in 1:NE)
        {
         gv[,i]=sort(rnorm(NQTN))
        }
        gv=gv[sample(1:NQTN,NQTN,replace=F),]
    }
    QTN.position=sample(1:m,NQTN,replace=F)
    SNPQ=as.matrix(X[,(QTN.position)])
    bv=SNPQ%*%gv


    bv_var=diag(var(bv))
    
    ####### calculate P A I E variance
    taxa_name=NULL
    cov_e=diag(1,NE)
    for(i in 1:NE)
    {   
        taxa_name=append(taxa_name,paste("trait_",i,sep=""))
        ve_i=bv_var[i]*(1/ha2[i]-1)
        for(j in 1:NE)
        {
            ve_j=bv_var[j]*(1/ha2[j]-1)
            cov_e[i,j]=cov_e[j,i]=sqrt(ve_i*ve_j)


        }
    }
    if(!is.null(NETN))
    {
        add_eff=rnorm(NETN,bv_var)
        add_posi=sample(c(1:m)[-c(QTN.position)],NQTN,replace=F)
        SNP_add=as.matrix(X[,(add_posi)])
        be=SNPQ%*%add_eff
        be=cbind(be,be)
        bv=bv+be
        QTN.position=c(QTN.position,add_posi)
    }

    residual=mvrnorm(n,t(mu_g),cov_e)
    Y=round((bv+residual),6)
    myY=cbind(as.data.frame(taxa),Y)
    colnames(myY)=c("taxa",taxa_name)
    
    #MeanY and GbyEY conversion
    Mean.Y=myY[,-3]
    Mean.Y[,2]=(myY[,2]+myY[,3])/2
    myY.1=cbind(paste(myY[,1],"-1",sep=""),myY[,2])
    myY.2=cbind(paste(myY[,1],"-2",sep=""),myY[,3])
    GbyE.Y=rbind(myY.1,myY.2)
    colnames(GbyE.Y)=c("taxa","trait")
    GbyE.Y=as.data.frame(as.matrix(data.table(GbyE.Y)))
    GbyE.Y[,2]=as.numeric(GbyE.Y[,2])
    
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
    #GbyE.QTN.Position conversion
    GbyE.QTN.Position=c(QTN.position,(QTN.position+ncol(GD[,-1])))
    
    return(list(Y=myY,GbyE.Y=GbyE.Y,Mean.Y=Mean.Y,u=bv,QTN.Position=QTN.position,GbyE.QTN.Position=GbyE.QTN.Position,residual=residual))
} #enf of phenotype simulation function
#=============================================================================================


