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