PhenoComp<-
  function(expdata,label,gene,freq,method,freq1,outfile1,outfile2){
    control_exp<-expdata[,label==0];
    case_exp<-expdata[,label==1];
    cat("differential expressed analysis... \n");
    outlier_dir<-NULL;
    outlier_pvalue<-NULL;
    Lgene<-dim(gene)[1];
    for (k in 1:Lgene){
      cat(k)
      cat("\n")
      Nnorm=control_exp[k,]
      colN<-dim(control_exp)[2]
      colC<-dim(case_exp)[2]
      N_tmp=matrix(rep(Nnorm,Lgene),ncol=colN,byrow=T)-control_exp
      Nloc_up=which((rowSums(N_tmp>0)/colN)>freq)
      Nloc_down=which((rowSums(N_tmp<0)/colN)>freq)
      reverse=matrix(0,colC,4)
      reverse[,1]=rep(length(Nloc_up),colC)
      reverse[,2]=rep(length(Nloc_down),colC)
      Tcanc=case_exp[k,]
      if (length(Nloc_up)>0){
        N_tmp=matrix(rep(Tcanc,length(Nloc_up)),ncol=colC,byrow=T)-case_exp[Nloc_up,]
        case_p=colSums(N_tmp<0)
        reverse[,3]=case_p
      }
      if (length(Nloc_down)>0){
        N_tmpp=matrix(rep(Tcanc,length(Nloc_down)),ncol=colC,byrow=T)-case_exp[Nloc_down,]
        case_pp=colSums(N_tmpp>0)
        reverse[,4]=case_pp;
      }
      GenePair_sig=NULL
      GenePair=rep(0,colC)
      GenePair[which(reverse[,3]>reverse[,4])]<--1
      GenePair[which(reverse[,3]<reverse[,4])]<-1
      tmp=matrix(c(reverse[,1],reverse[,2],reverse[,1]-reverse[,3]+reverse[,4], reverse[,2]-reverse[,4]+reverse[,3]),ncol=4)
      GenePair_sig<-apply(tmp,1,function(x) fisher.test(matrix(x,ncol=2,byrow=T))$p.value)
      outlier_dir=rbind(outlier_dir,GenePair)
      outlier_pvalue=rbind(outlier_pvalue,GenePair_sig)
    }
    fdr<-apply(outlier_pvalue,2,function(x) p.adjust(x,method="fdr",length(x)))

	individual_level_DEGs<-matrix(0,nrow(outlier_dir),ncol(outlier_dir))
	for(i in 1:ncol(outlier_dir)){
	a<-which(fdr[,i]<0.05)
	individual_level_DEGs[a,i]<-outlier_dir[a,i]
	}
    #population-level up-regulation genes
      p_diff<-c()
      for(i in 1:ncol(individual_level_DEGs)){
         want_col_DEG<-individual_level_DEGs[,i]
         up_num<-sum(want_col_DEG==1)
         down_num<-sum(want_col_DEG==-1)
         all_num<-nrow(individual_level_DEGs)
         p_diff1<-(up_num+down_num)/all_num
         p_diff<-c(p_diff,p_diff1)
          }
	  p_up_result<-c()
      for(i in 1:ncol(individual_level_DEGs)){
         want_col_DEG<-individual_level_DEGs[,i]
         up_num<-sum(want_col_DEG==1)
         down_num<-sum(want_col_DEG==-1)
         p_up_ratio<-up_num/(up_num+down_num)
         p_up_result<-c(p_up_result,p_up_ratio)
          }
	  if(method==1){
        p0_up<-(median(p_diff))*(median(p_up_result))
		print(p0_up)}
	  if(method==2){
        p0_up<-(mean(p_diff))*(mean(p_up_result))
		print(p0_up)}
	  up_result<-matrix(0,dim(individual_level_DEGs)[1],4)
      z<-1
      for(i in 1:dim(individual_level_DEGs)[1]){
        k<-sum(individual_level_DEGs[i,]==1)
        s<-dim(individual_level_DEGs)[2]
        up_result[z,2]<-k/s
        up_result[z,3]<-1-pbinom(k-1,s,p0_up)
        z<-z+1
      }
      up_result[,1]<-gene
      up_result[,4]=p.adjust(up_result[,3],method="BH")
      colnames(up_result)<-c("geneid","up_ratio","p.value","FDR")
      up_DEG<-up_result[which(up_result[,4]<freq1),1]
      print(length(up_DEG))
    #population-level down-regulation genes
	  p_down_result<-c()
      for(i in 1:ncol(individual_level_DEGs)){
        want_col_DEG<-individual_level_DEGs[,i]
        up_num<-sum(want_col_DEG==1)
        down_num<-sum(want_col_DEG==-1)
        p_down_ratio<-down_num/(up_num+down_num)
        p_down_result<-c(p_down_result,p_down_ratio)
        }
	 if(method==1){
        p0_down<-(median(p_diff))*(median(p_down_result))
		print(p0_down)}
	  if(method==2){
        p0_down<-(mean(p_diff))*(mean(p_down_result))
		print(p0_down)}
      down_result<-matrix(0,dim(individual_level_DEGs)[1],4)
      z<-1
      for(i in 1:dim(individual_level_DEGs)[1]){
        k<-sum(individual_level_DEGs[i,]==-1)
        s<-dim(individual_level_DEGs)[2]
        down_result[z,2]<-k/s
        down_result[z,3]<-1-pbinom(k-1,s,p0_down)
        z<-z+1
      }
      down_result[,1]<-gene
      down_result[,4]=p.adjust(down_result[,3],method="BH")
      colnames(down_result)<-c("geneid","down_ratio","p.value","FDR")
      down_DEG<-down_result[which(down_result[,4]<freq1),1]
	  print(length(down_DEG))
      #Overlap of population-level up-regulated genes and population-level down regulated genes
      overlap_gene<-intersect(up_DEG,down_DEG)
	  overlap_gene_num<-length(overlap_gene)
      print(overlap_gene_num)
	  if(overlap_gene_num==0){
	  print(length(up_DEG)+length(down_DEG))
	  write.table(up_DEG,file=outfile1,sep="\t",row.names=FALSE,col.names=FALSE)
	  write.table(down_DEG,file=outfile2,sep="\t",row.names=FALSE,col.names=FALSE)
	  }
	  if(overlap_gene_num>0){
	  up_gene<-up_DEG[-match(overlap_gene,up_DEG)]
      down_gene<-down_DEG[-match(overlap_gene,down_DEG)]
      print(length(up_gene))
	  print(length(down_gene))
	  print(length(up_gene)+length(down_gene))

	  write.table(up_gene,file=outfile1,sep="\t",row.names=FALSE,col.names=FALSE)
	  write.table(down_gene,file=outfile2,sep="\t",row.names=FALSE,col.names=FALSE)
      }
  }



