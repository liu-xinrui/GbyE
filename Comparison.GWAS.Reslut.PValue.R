Comparison.GWAS.Reslut.PValue <- function(GM,GWASresult){
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