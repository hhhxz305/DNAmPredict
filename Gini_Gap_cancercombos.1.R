
#make sure methylation data is transposed: rows are individuals, colomns are methylations

if(dim(brca.met$merged.dat)[1]>dim(brca.met$merged.dat)[2]){
  brca.met$merged.dat=t(brca.met$merged.dat)
}
if(dim(coad.met$merged.dat)[1]>dim(coad.met$merged.dat)[2]){
  coad.met$merged.dat=t(coad.met$merged.dat)
}
if(dim(kirc.met$merged.dat)[1]>dim(kirc.met$merged.dat)[2]){
  kirc.met$merged.dat=t(kirc.met$merged.dat)
}
if(dim(cesc.met$merged.dat)[1]>dim(cesc.met$merged.dat)[2]){
  cesc.met$merged.dat=t(cesc.met$merged.dat)
}
if(dim(luad.met$merged.dat)[1]>dim(luad.met$merged.dat)[2]){
  luad.met$merged.dat=t(luad.met$merged.dat)
}

#may need to set working directory

#loading imported data
load('brca_met_raw.rdata')
load('kirc_met_raw.rdata')
load('coad_met_raw.rdata')
load('luad_met_raw.rdata')
load('cesc_met_raw.rdata')

#filter to only white and black

exam_dat=list()
length(exam_dat)=5


for(i in 1:5){
  exam_dat[[i]]=list()
}
names(exam_dat)=c('BRCA','COAD','KIRC', 'CESC','LUAD')

brcaids=rownames(brca.met$clinical)[which(brca.met$clinical[,grep('race',colnames(brca.met$clinical))]=='white'|brca.met$clinical[,grep('race',colnames(brca.met$clinical))]=='black or african american')]
brcaids=as.numeric(na.omit(match(brcaids,brca.met$merged.dat[,1])))
exam_dat$BRCA$met=as.matrix(brca.met$merged.dat[brcaids,-(1:3)])
#exam_dat$BRCA$snp=NULL
rm(brca.met)

coadids=rownames(coad.met$clinical)[which((coad.met$clinical[,grep('race',colnames(coad.met$clinical))]=='white'|coad.met$clinical[,grep('race',colnames(coad.met$clinical))]=='black or african american')&coad.met$clinical[,which(colnames(coad.met$clinical)=='gender')]=='female')]
coadids=as.numeric(na.omit(match(coadids,coad.met$merged.dat[,1])))
exam_dat$COAD$met=as.matrix(coad.met$merged.dat[coadids,-(1:3)])
#exam_dat$COAD$snp=NULL
rm(coad.met)

kircids=rownames(kirc.met$clinical)[which((kirc.met$clinical[,grep('race',colnames(kirc.met$clinical))]=='white'|kirc.met$clinical[,grep('race',colnames(kirc.met$clinical))]=='black or african american')&kirc.met$clinical[,which(colnames(kirc.met$clinical)=='gender')]=='female')]
kircids=as.numeric(na.omit(match(kircids,kirc.met$merged.dat[,1])))
exam_dat$KIRC$met=as.matrix(kirc.met$merged.dat[kircids,-(1:3)])
#exam_dat$KIRC$snp=NULL
rm(kirc.met)


cesc.met$clinical=cesc$clinical
cescids=rownames(cesc.met$clinical)[which((cesc.met$clinical[,grep('race',colnames(cesc.met$clinical))]=='white'|cesc.met$clinical[,grep('race',colnames(cesc.met$clinical))]=='black or african american')&cesc.met$clinical[,which(colnames(cesc.met$clinical)=='gender')]=='female')]
cescids=as.numeric(na.omit(match(cescids,rownames(cesc.met$merged.dat))))
exam_dat$CESC$met=as.matrix(cesc.met$merged.dat)[cescids,]
#exam_dat$CESC$snp=NULL
rm(cesc.met)


luad.met$clinical=luad$clinical
luadids=rownames(luad.met$clinical)[which((luad.met$clinical[,grep('race',colnames(luad.met$clinical))]=='white'|luad.met$clinical[,grep('race',colnames(luad.met$clinical))]=='black or african american')&luad.met$clinical[,which(colnames(luad.met$clinical)=='gender')]=='female')]
luadids=as.numeric(na.omit(match(luadids,rownames(luad.met$merged.dat))))
exam_dat$LUAD$met=as.matrix(luad.met$merged.dat)[luadids,]
#exam_dat$LUAD$snp=NULL
rm(luad.met)

#Create function for naive measure combining gini and Gap

dat=exam_dat
incl_cancer = 'CESC'  #Cancer type needs to keep
snp_var_thresh=0.995   #variance threshold for SNPs
met_var_thresh=0.9999  #variance threshold for methylation
seed=166
d.power=1 
samp=.7  #training proportion


library(cluster)

gini.cluster=function(x,clus,data){
  #library(reldist)
  gini=function (x, weights = rep(1, length = length(x))) 
  {
    ox <- order(x)
    x <- x[ox]
    weights <- weights[ox]/sum(weights)
    p <- cumsum(weights)
    nu <- cumsum(weights * x)
    n <- length(nu)
    nu <- nu/nu[n]
    sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
  }
  df=data.frame(data)
  x=names(df)[1]
  #gini.each=as.numeric(aggregate(as.formula(paste(x,'~clus')),data=data.frame(data),gini)[2])
  #sum.each=as.numeric(aggregate(as.formula(paste(x,'~clus')),data=data.frame(data),length)[2])
  gini.each=NULL
  sum.each=NULL
  clustypes=unique(clus)
  for(i in 1:length(unique(clus))){
    sum.each=c(sum.each,sum(clus==clustypes[i]))
    gini.each=c(gini.each,gini(data[clus==clustypes[i]]))
  }
  gini.cluster=sum(gini.each*sum.each)/length(data)
  return(gini.cluster)
}

clus_extract=function (x, method = "firstSEmax", SE.factor = 1, ...) 
{
  method <- match.arg(method, choices = eval(formals(maxSE)$method))
  stopifnot((K <- nrow(T <- x$Tab)) >= 1, SE.factor >= 0)
  nc <- maxSE(f = T[, "gap"], SE.f = T[, "SE.sim"], method = method, 
              SE.factor = SE.factor)
  nc
}

ent=function (truth, cluster, mode = c("entropy", "gini")){
  mode <- match.arg(mode, c("entropy", "gini"))
  k = length(unique(truth))
  tab = table(truth, cluster)
  clustersizes = colSums(tab)
  clustersizesnorm = clustersizes/sum(clustersizes)
  tabprop = tab
  lapply(1:(dim(tab)[2]), function(i) {
    tabprop[, i] <<- tabprop[, i]/clustersizes[i]
    NULL
  })
  if (mode == "entropy") {
    measure = 0
    maxmeasure = 0
    lapply(1:(dim(tabprop)[2]), function(i) {
      clustermeasure = 0
      maxclustermeasure = 0
      lapply(1:(dim(tabprop)[1]), function(j) {
        if (tabprop[j, i] != 0) {
          clustermeasure <<- clustermeasure + -tabprop[j,i] * log2(tabprop[j, i])
        }
        maxclustermeasure <<- maxclustermeasure + -1/k *log2(1/k)										   
        NULL
      })
      measure <<- measure + clustersizesnorm[i] * clustermeasure
      maxmeasure <<- maxmeasure + clustersizesnorm[i] *maxclustermeasure
    })
  }
  if (mode == "gini") {
    measure = 0
    maxmeasure = 0
    lapply(1:(dim(tabprop)[2]), function(i) {
      clustermeasure = 1
      maxclustermeasure = 1
      lapply(1:(dim(tabprop)[1]), function(j) {
        if (tabprop[j, i] != 0) {
          clustermeasure <<- clustermeasure + (-tabprop[j,i]^2)
        }
        maxclustermeasure <<- maxclustermeasure + (-1/k^2)
        NULL
      })
      measure <<- measure + clustersizesnorm[i] * clustermeasure
      maxmeasure <<- maxmeasure + clustersizesnorm[i] * maxclustermeasure
      NULL
    })
  }
  list(result = measure, measure = mode, normalized_measure = measure/maxmeasure)
}


cancertypes=names(dat)

#creating all combinations
combs=list()
lapply(1:length(cancertypes),function(i){
  combs<<-c(combs,combn(1:length(cancertypes), i, simplify = F))
  NULL
})
incl_cancer_indx=match(incl_cancer,cancertypes)
remvindx=NULL
if (!is.null(incl_cancer)){
  for(i in 1:length(combs)){
    if((sum(is.na(match(incl_cancer_indx,combs[[i]])))==length(incl_cancer_indx))){
      remvindx=c(remvindx,i)
    }
  }
  
}
for(i in length(remvindx):1){
  combs[[remvindx[i]]]=NULL
}
ginicollector=list()
gapcol=list()
entcol=list()
clustginicol=list()
overallgini=list()
tabcollect=list()
lapply(1:length(combs),function(i){
  
  
  
#combining datasets
  if(length(combs[[i]])==1){
    set.seed(1315)
    indx=sample(1:dim(dat[[combs[[i]]]]$met)[1],dim(dat[[combs[[i]]]]$met)[1]*samp)
    testdat_met=dat[[combs[[i]]]]$met[indx,]
    types=rep(1,length(indx))
    #types=rep(1,dim(dat[[combs[[i]]]]$met)[1])
    #testdat_snp=dat[[combs[[i]]]]$snp
  }else{
    set.seed(1315)
    indx=sample(1:dim(dat[[combs[[i]][1]]]$met)[1],dim(dat[[combs[[i]][1]]]$met)[1]*samp)
    testdat_met=dat[[combs[[i]][1]]]$met[indx,]
    types=rep(1,length(indx))
    #types=rep(1,dim(dat[[combs[[i]][1]]]$met)[1])
    #testdat_snp=dat[[combs[[i]][1]]]$snp
    for(j in 2:length(combs[[i]])){
      set.seed(1315)
      indx=sample(1:dim(dat[[combs[[i]][j]]]$met)[1],dim(dat[[combs[[i]][j]]]$met)[1]*samp)
      colnames(testdat_met)=NULL
      colnames(dat[[combs[[i]][j]]]$met)=NULL
      testdat_met=as.matrix(rbind(as.matrix(testdat_met),as.matrix(dat[[combs[[i]][j]]]$met)[indx,]))
      types=c(types,rep(j,length(indx)))
      #types=c(types,rep(j,dim(dat[[combs[[i]][j]]]$met)[1]))
      #testdat_snp=rbind(testdat_snp,dat[[combs[[i]][j]]]$snp)
    }
  }
  
  print('pre filter')
  print(dim(testdat_met))    
  print(class(testdat_met))
  print(class(testdat_met[1,1]))
  
  storage.mode(testdat_met) <- "numeric"
  #filtering to threshold
  met.var=apply(testdat_met,2,var)
  print(paste('length', length(met.var), 'NAs', sum(is.na(met.var))))
  
  #snp.var=apply(testdat_snp,2,var)
  testdat_met <- testdat_met[, met.var >= quantile(met.var, met_var_thresh, na.rm=T) & !is.na(met.var),drop=F]
  #testdat_snp <- testdat_snp[, snp.var >= quantile(snp.var, snp_var_thresh, na.rm=T) & !is.na(snp.var),drop=F]
  print('post filter')
  print(dim(testdat_met))    
  if(dim(testdat_met)[2]!=0){
    print('running clusters')
    pam1 <- function(x,k) list(cluster = pam(x,k, cluster.only=TRUE))
    set.seed(seed)
    gap2 <- clusGap(testdat_met,pam1,12)
    clusts=clus_extract(gap2)
    set.seed(seed)
    groups=kmeans(testdat_met,clusts)
    ginicol=list()
    ginicols=list()
    for(j in 1:dim(testdat_met)[2]){
      ginicols[[j]]=gini.cluster(colnames(testdat_met)[j],groups$cluster,testdat_met[,j])
    }
    clustginicollector=list()
    for(j in 1:clusts){
      clustginicollector[[j]]=ent(cluster=rep(j,sum(groups$cluster==j)),truth=types[groups$cluster==j],mode='gini')
    }
    ginicol$gini=ginicols
    ginicol$methyl_names=names(testdat_met)
    #names(ginicol)=c('Gini','Methyl_names')
   
 #calculating G/R
    #gapcol<<-c(gapcol,gap2$Tab[clusts,3])
    ginicollector[[i]]<<-ginicol #now collecting all gini for each feature
    gapcol[[i]]<<-gap2$Tab
    entcol[[i]]<<-ent(types,groups$cluster)
    overallgini[[i]]<<-ent(types,groups$cluster,mode='gini')
    tabcollect[[i]]<<-table(types,groups$cluster)
    clustginicol[[i]]<<-clustginicollector
  }else{
    #gapcol<<-c(gapcol,NA) #if p=0 then is unable to run
    gapcol[[i]]<<-NA
    entcol[[i]]<<-NA
    ginicollector[[i]]<<-NA
  }
  save(gapcol,ginicollector,combs,file='temp_gini_backup.RData')
})

combos=NULL
for(i in 1:length(combs)){
  combos=c(combos,paste(cancertypes[combs[[i]]],collapse = ','))
}

#Collecting resulting gini and Gaps
df=data.frame(matrix(NA,nrow=length(unlist(combos)),ncol=13))
df
for (i in 1:length(unlist(combos))) {
   gcol = NULL
   for (j in 1:length(clustginicol[[i]])) {
     gcol = c(gcol,clustginicol[[i]][[j]][[1]])
   }
   gcol=c(gcol,rep(NA,12-length(gcol)))
   
   df[i,] = c(combos[[i]],gcol)
}
names(df)=c('cancer_combo',paste('cluster',1:12))
df
str(df)
for (i in 2:dim(df)[2]){
  df[,i] = as.numeric(df[,i])
}
tabcollect[[5]]
tabcollect[[10]]
library(xtable)
library(clipr)
xtable(df)
xtable(tabcollect[[5]])
xtable(tabcollect[[10]])


prinfun=function(){
  for(i in 1:length(tabcollect)){
  print(xtable(tabcollect[[i]],caption = combos[i]))
  }
}
write_clip(prinfun())


#plotting resulting gini and Gaps
mgini=NULL
for(i in 1:length(gapcol)){
  mgini=c(mgini,gapcol[[i]][maxSE(gapcol[[i]][,3],gapcol[[i]][,4],method = "firstSEmax", SE.factor = 1),3])
}
par(mar=c(9.1, 4.1, 4.1, 2.1))
plot(1:dim(df)[1],apply(df[,2:dim(df)[2]],1,max,na.rm=T),type='l',xaxt='n',col=1,ylim=c(0,2),ylab='Measure',xlab='')
points(1:dim(df)[1],mgini,type='l',col=2)
points(1:dim(df)[1],apply(df[,2:dim(df)[2]],1,max,na.rm=T)+mgini,type='l',col=3)
axis(1,labels = unlist(combos),at=1:dim(df)[1],las=2,cex.axis=.6)
legend('topleft',col=1:3,legend=c('Gini', 'Gap','Gap+Gini'),lty=1)







