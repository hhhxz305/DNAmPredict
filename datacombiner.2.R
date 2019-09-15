
#set working directory
setwd('d:/combodat')

#ram requirement very high
memory.limit(size=100000)
library(data.table)
setDTthreads(6)

#select cancer types to combine
cancers=c('BRCA','COAD','KIRC', 'CESC','LUAD')[c(4,5)]   #Currently set to CESC and LUAD as per Gini/Gap run

CESC12=NULL
types=NULL
if(sum(cancers=='BRCA')==1){
  if(is.null(CESC12)){
    holder=as.matrix(fread('brca.csv',stringsAsFactors = F))
    nam=holder[,1]
    nam=strsplit(nam,split='\\.')
    nam2=rep('',length(nam))
    for(i in 1:length(nam)){
      print(i)
      if(length(nam[[i]])>0){
        nam2[i]=nam[[i]][[length(nam[[i]])]]
      }
    }
    rownames(holder)=nam2
    rm(nam)
    rm(nam2)
    holder=t(holder[,-1])[,-1]
    #next need to filter race
    holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
    CESC12=holder
    types=c(types,rep('BRCA',dim(CESC12)[1]))
  }else{
    holder=as.matrix(fread('brca.csv',stringsAsFactors = F))
    nam=holder[,1]
    nam=strsplit(nam,split='\\.')
    nam2=rep('',length(nam))
    for(i in 1:length(nam)){
      print(i)
      if(length(nam[[i]])>0){
        nam2[i]=nam[[i]][[length(nam[[i]])]]
      }
    }
    rownames(holder)=nam2
    rm(nam)
    rm(nam2)
    holder=t(holder[,-1])[,-1] 
    holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
    cols=intersect(colnames(CESC12),colnames(holder))
    indx1=match(cols,colnames(CESC12))
    indx2=match(cols,colnames(holder))
    CESC12=rbind(CESC12[,indx1],holder[,indx2])
    #print(paste('Combining',cancers[i]))
    types=c(types,rep('BRCA',dim(holder)[1]))
  }
  rm(holder)
  gc()
}
if(sum(cancers=='CESC')==1){
  if(is.null(CESC12)){
    holder=as.matrix(fread('cesc.csv',stringsAsFactors = F))
    nam=holder[,1]
    nam=strsplit(nam,split='\\.')
    nam2=rep('',length(nam))
    for(i in 1:length(nam)){
      print(i)
      if(length(nam[[i]])>0){
        nam2[i]=nam[[i]][[length(nam[[i]])]]
      }
    }
    rownames(holder)=nam2
    rm(nam)
    rm(nam2)
    holder=t(holder[,-1])[,-1]
    holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
    CESC12=holder
    types=c(types,rep('CESC',dim(CESC12)[1]))
  }else{
    holder=as.matrix(fread('cesc.csv',stringsAsFactors = F))
    nam=holder[,1]
    nam=strsplit(nam,split='\\.')
    nam2=rep('',length(nam))
    for(i in 1:length(nam)){
      print(i)
      if(length(nam[[i]])>0){
        nam2[i]=nam[[i]][[length(nam[[i]])]]
      }
    }
    rownames(holder)=nam2
    rm(nam)
    rm(nam2)
    holder=t(holder[,-1])[,-1]
    holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
    cols=intersect(colnames(CESC12),colnames(holder))
    indx1=match(cols,colnames(CESC12))
    indx2=match(cols,colnames(holder))
    CESC12=rbind(CESC12[,indx1],holder[,indx2])
    #print(paste('Combining',cancers[i]))
    types=c(types,rep('CESC',dim(holder)[1]))
  }
  rm(holder)
  gc()
}
if(sum(cancers=='COAD')==1){
  if(is.null(CESC12)){
    holder=as.matrix(fread('coad.csv',stringsAsFactors = F))
    nam=holder[,1]
    nam=strsplit(nam,split='\\.')
    nam2=rep('',length(nam))
    for(i in 1:length(nam)){
      print(i)
      if(length(nam[[i]])>0){
        nam2[i]=nam[[i]][[length(nam[[i]])]]
      }
    }
    rownames(holder)=nam2
    rm(nam)
    rm(nam2)
    holder=t(holder[,-1])[,-1]
    holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
    CESC12=holder
    types=c(types,rep('COAD',dim(CESC12)[1]))
  }else{
    holder=as.matrix(fread('coad.csv',stringsAsFactors = F))
    nam=holder[,1]
    nam=strsplit(nam,split='\\.')
    nam2=rep('',length(nam))
    for(i in 1:length(nam)){
      print(i)
      if(length(nam[[i]])>0){
        nam2[i]=nam[[i]][[length(nam[[i]])]]
      }
    }
    rownames(holder)=nam2
    rm(nam)
    rm(nam2)
    holder=t(holder[,-1])[,-1] 
    holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
    cols=intersect(colnames(CESC12),colnames(holder))
    indx1=match(cols,colnames(CESC12))
    indx2=match(cols,colnames(holder))
    CESC12=rbind(CESC12[,indx1],holder[,indx2])
    #print(paste('Combining',cancers[i]))
    types=c(types,rep('COAD',dim(holder)[1]))
  }
  rm(holder)
  gc()
}
if(sum(cancers=='KIRC')==1){
  if(is.null(CESC12)){
    holder=as.matrix(fread('kirc.csv',stringsAsFactors = F))
    nam=holder[,1]
    nam=strsplit(nam,split='\\.')
    nam2=rep('',length(nam))
    for(i in 1:length(nam)){
      print(i)
      if(length(nam[[i]])>0){
        nam2[i]=nam[[i]][[length(nam[[i]])]]
      }
    }
    rownames(holder)=nam2
    rm(nam)
    rm(nam2)
    holder=t(holder[,-1])[,-1]
    holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
    CESC12=holder
    types=c(types,rep('KIRC',dim(CESC12)[1]))
  }else{
    holder=as.matrix(fread('kirc.csv',stringsAsFactors = F))
    nam=holder[,1]
    nam=strsplit(nam,split='\\.')
    nam2=rep('',length(nam))
    for(i in 1:length(nam)){
      print(i)
      if(length(nam[[i]])>0){
        nam2[i]=nam[[i]][[length(nam[[i]])]]
      }
    }
    rownames(holder)=nam2
    rm(nam)
    rm(nam2)
    holder=t(holder[,-1])[,-1] 
    holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
    cols=intersect(colnames(CESC12),colnames(holder))
    indx1=match(cols,colnames(CESC12))
    indx2=match(cols,colnames(holder))
    CESC12=rbind(CESC12[,indx1],holder[,indx2])
    #print(paste('Combining',cancers[i]))
    types=c(types,rep('KIRC',dim(holder)[1]))
  }
  rm(holder)
  gc()
}
if(sum(cancers=='LUAD')==1){
  if(is.null(CESC12)){
    holder=as.matrix(fread('luad.csv',stringsAsFactors = F))
    nam=holder[,1]
    nam=strsplit(nam,split='\\.')
    nam2=rep('',length(nam))
    for(i in 1:length(nam)){
      print(i)
      if(length(nam[[i]])>0){
        nam2[i]=nam[[i]][[length(nam[[i]])]]
      }
    }
    rownames(holder)=nam2
    rm(nam)
    rm(nam2)
    holder=t(holder[,-1])[,-1]
    holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
    CESC12=holder
    types=c(types,rep('LUAD',dim(CESC12)[1]))
  }else{
    holder=as.matrix(fread('luad.csv',stringsAsFactors = F))
    nam=holder[,1]
    nam=strsplit(nam,split='\\.')
    nam2=rep('',length(nam))
    for(i in 1:length(nam)){
      print(i)
      if(length(nam[[i]])>0){
        nam2[i]=nam[[i]][[length(nam[[i]])]]
      }
    }
    rownames(holder)=nam2
    rm(nam)
    rm(nam2)
    holder=t(holder[,-1])[,-1] 
    holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
    cols=intersect(colnames(CESC12),colnames(holder))
    indx1=match(cols,colnames(CESC12))
    indx2=match(cols,colnames(holder))
    CESC12=rbind(CESC12[,indx1],holder[,indx2])
    #print(paste('Combining',cancers[i]))
    types=c(types,rep('LUAD',dim(holder)[1]))
  }
  rm(holder)
  gc()
}


#types=types[(CESC12[,'race']=='white'|CESC12[,'race']=='black or african american')&is.na(CESC12[,'race'])]
#next need to filter race
#CESC12=CESC12[(CESC12[,'race']=='white'|CESC12[,'race']=='black or african american')&is.na(CESC12[,'race']),]

gc()



#cancer types
table(types)

dim(CESC12)

#separating genetic measures from clinical
CESC12_clin=CESC12[,-(1:505039)]
CESC12=CESC12[,(1:505039)]
storage.mode(CESC12)='numeric'

#set .7 as training and .3 as testing
set.seed(1314)
train.idx <- sample(1:dim(CESC12)[1],round(.7*dim(CESC12)[1]))
test.idx <- (1:dim(CESC12)[1])[-train.idx]

lapply(22619:505039,function(i){
a=log(CESC12[,i]/(1-CESC12[,i]),2)
a[a==Inf|a==-Inf]=NA
CESC12[,i]<<-a
NULL
})

#Pre-processing data and filtering data

###Based on Training dataset
all.var <- apply(CESC12[train.idx,1:505039], 2, var)
######CESC:482421 MET, 22618 SNP; LUAD:482421 MET, 22618 SNP;


met.var <- all.var[22619:505039]
snp.var <- all.var[1:22618]

hist(met.var,breaks=50)
hist(met.var,breaks=50,ylim=c(0,50))
abline(v=quantile(met.var, 0.9999, na.rm=T),col='red')

hist(snp.var,breaks=50)
hist(snp.var,breaks=50,ylim=c(0,50))
abline(v=quantile(snp.var, 0.995, na.rm=T),col='red')

met.data <- CESC12[,22619:505039][, met.var >= quantile(met.var, 0.9999, na.rm=T) & !is.na(met.var)]
snp.data <- CESC12[,1:22618][, snp.var >= quantile(snp.var, 0.995, na.rm=T) & !is.na(snp.var)]
race.data <-CESC12_clin[,'race']
type.data <-types
rm(CESC12)
gc()

#Saving final combined dataset
save.image(paste(paste(cancers,collapse = '_'),'_combined.rdata',sep=''))



