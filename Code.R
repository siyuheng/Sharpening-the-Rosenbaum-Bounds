#Input data (please adjust the file address to your own)
wls<-read.delim("/Users/heng6/Desktop/Biometrika Submission/Code and Data/Data.txt", sep=",")

#discard units with missing anger/hostility tendency
sub_data<-wls[which(wls$mu025rer>=0),] 
#discard missing outcomes
sub_data<-sub_data[which(sub_data$gx351re>0),]
#we focus on male respondents
sub_data<-sub_data[which(sub_data$sexrsp==1),]
#discard missing childhood maltreatment indicators
sub_data<-sub_data[which(sub_data$iv219rer>0 | sub_data$iv213rer>0 | sub_data$iv210rer>0 | sub_data$iv216rer>0 | sub_data$iv218rer>0 | sub_data$iv215rer>0),]
#age
sub_data<-sub_data[which(sub_data$ga003re>=0),]
#educational attainment
sub_data<-sub_data[which(sub_data$rb037re>0),]
#body mass index
sub_data<-sub_data[which(sub_data$mx011rec>=0),]
#smoking regularly or not 
sub_data<-sub_data[which(sub_data$ixt06rer==-2 | sub_data$ixt06rer>0),]
#drinking alcohol or not 
sub_data<-sub_data[which(sub_data$ru025re>0),]


ID<-sub_data$idpub

#define the treated unit and control
anger<-rep(-1,nrow(sub_data))
anger[sub_data$mu025rer>=3]=1
anger[sub_data$mu025rer==0]=0
table(anger)

#define the outcome
table(sub_data$gx351re)
heart<-rep(0, nrow(sub_data))
heart[sub_data$gx351re==1]=1
table(heart)

#covariates
sex<-sub_data$sexrsp
table(sex)
age<-sub_data$ga003re
table(age)
edu<-sub_data$rb037re
table(edu)
bodymass<-rep(0,nrow(sub_data))
bodymass[sub_data$mx011rec<25]=1  #normal
bodymass[sub_data$mx011rec>=25 & sub_data$mx011rec<30]=2  #overweight
bodymass[sub_data$mx011rec>=30 & sub_data$mx011rec<35]=3  #obese
bodymass[sub_data$mx011rec>=35]=4  #very obese
table(bodymass)
smoke<-as.numeric(sub_data$ixt06rer!=-2)
table(smoke)
alcohol<-sub_data$ru025re
table(alcohol)

#define the childhood maltreatment indicator
abuse<-as.numeric(sub_data$iv219rer>1 | sub_data$iv213rer>1 | sub_data$iv210rer>1 | sub_data$iv216rer>1 | sub_data$iv218rer>1 | sub_data$iv215rer>1)
table(abuse)

real_data<-data.frame(ID, anger, heart, sex, age, edu, bodymass, smoke, alcohol, abuse)
#We discard units with middle-level anger/hostility tendency
real_data<-real_data[which(real_data$anger!=-1),]
table(real_data$anger)

# Fit a propensity score using logistic regression with each covariate entering  
# linearly into the logistic link function 
# Put x=TRUE in order to have model object include design matrix 

propscore.model=glm(anger~age+edu+bodymass+smoke+alcohol+abuse,family=binomial,x=TRUE,y=TRUE,data=real_data); 
treated=propscore.model$y 


library(MASS) 

# Function for computing  

# rank based Mahalanobis distance.  Prevents an outlier from 
# inflating the variance for a variable, thereby decreasing its importance. 
# Also, the variances are not permitted to decrease as ties  
# become more common, so that, for example, it is not more important 
# to match on a rare binary variable than on a common binary variable 
# z is a vector, length(z)=n, with z=1 for treated, z=0 for control 
# X is a matrix with n rows containing variables in the distance 

smahal=function(z,X){ 
  X<-as.matrix(X) 
  n<-dim(X)[1] 
  rownames(X)<-1:n 
  k<-dim(X)[2] 
  m<-sum(z) 
  for (j in 1:k) X[,j]<-rank(X[,j]) 
  cv<-cov(X) 
  vuntied<-var(1:n) 
  rat<-sqrt(vuntied/diag(cv)) 
  cv<-diag(rat)%*%cv%*%diag(rat) 
  out<-matrix(NA,m,n-m) 
  Xc<-X[z==0,] 
  Xt<-X[z==1,] 
  rownames(out)<-rownames(X)[z==1] 
  colnames(out)<-rownames(X)[z==0] 
  library(MASS) 
  icov<-ginv(cv) 
  for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T) 
  out 
} 


# Function for adding a propensity score caliper to a distance matrix dmat 
# calipersd is the caliper in terms of standard deviation of the logit propensity scoe 

addcaliper=function(dmat,z,logitp,calipersd=.2,penalty=1000){ 
  sd.logitp=sd(logitp) 
  adif=abs(outer(logitp[z==1],logitp[z==0],"-")) 
  adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp)) 
  dmat=dmat+adif*penalty 
  dmat 
} 

# Matrix of covariates, excluding intercept 
Xmat=propscore.model$x[,-1] 
# Rank based Mahalanobis distance 
distmat=smahal(treated,Xmat) 
# Add caliper 
logit.propscore=predict(propscore.model) 
distmat2=addcaliper(distmat,treated,logit.propscore) 

### Create a subject index and name the rows and columns of distance matrix by ### this subject index 

subject.index=seq(1,length(treated),1) 
rownames(distmat2)=subject.index[treated==1] 
colnames(distmat2)=subject.index[treated==0] 

# Pair Matching 

library(optmatch) 

matchvec=pairmatch(distmat2) 

# Note: Can ignore warning message from matching 

# Create vectors of the subject indices of the treatment units ordered by 

# their matched set and corresponding control unit 

treated.subject.index=rep(0,sum(treated==1)) 
matched.control.subject.index=rep(0,length(treated.subject.index)) 
matchedset.index=substr(matchvec,start=3,stop=10) 
matchedset.index.numeric=as.numeric(matchedset.index) 
subjects.match.order=as.numeric(names(matchvec)) # The subject indices in  

# the order of matchvec 

for(i in 1:length(treated.subject.index)){ 
  
  matched.set.temp=which(matchedset.index.numeric==i) 
  matched.set.temp.indices=subjects.match.order[matched.set.temp] 
  
  if(treated[matched.set.temp.indices[1]]==1){ 
    treated.subject.index[i]=matched.set.temp.indices[1] 
    matched.control.subject.index[i]=matched.set.temp.indices[2] 
  } 
  
  if(treated[matched.set.temp.indices[2]]==1){ 
    treated.subject.index[i]=matched.set.temp.indices[2] 
    matched.control.subject.index[i]=matched.set.temp.indices[1] 
  } 
  
} 

### Check balance 

# Calculate standardized differences  

# Covariates used in propensity score model 

Xmat=propscore.model$x; 
treatedmat=Xmat[treated==1,]; 

# Standardized differences before matching 

controlmat.before=Xmat[treated==0,]; 
controlmean.before=apply(controlmat.before,2,mean); 
treatmean=apply(treatedmat,2,mean); 
treatvar=apply(treatedmat,2,var); 
controlvar=apply(controlmat.before,2,var); 
stand.diff.before=(treatmean-controlmean.before)/sqrt((treatvar+controlvar)/2); 

# Standardized differences after matching 

controlmat.after=Xmat[matched.control.subject.index,]; 
controlmean.after=apply(controlmat.after,2,mean); 

# Standardized differences after matching 

stand.diff.after=(treatmean-controlmean.after)/sqrt((treatvar+controlvar)/2)

library(ggplot2)   

abs.stand.diff.before=abs(stand.diff.before[-1]) 
abs.stand.diff.after=abs(stand.diff.after[-1]) 
covariates=names(stand.diff.before[-1]) 
plot.dataframe=data.frame(abs.stand.diff=c(abs.stand.diff.before,abs.stand.diff.after),covariates=rep(covariates,2),type=c(rep("Before",length(covariates)),rep("After",length(covariates)))) 

#The loveplot for checking the covariate balance
ggplot(plot.dataframe,aes(x=abs.stand.diff,y=covariates))+geom_point(size=5,aes(shape=factor(type)))+scale_shape_manual(values=c(4,1)) 

#Construct matched data
matched_data<-NULL

for (i in 1:length(treated.subject.index)){
  matched_data<-rbind(matched_data, real_data[treated.subject.index[i],])
  matched_data<-rbind(matched_data, real_data[matched.control.subject.index[i],])
}

#The childhood maltreatment indicator and smoking indicator were been exactly matched
sum(abs(real_data$abuse[treated.subject.index]-real_data$abuse[matched.control.subject.index])) 
sum(real_data$smoke[treated.subject.index]-real_data$smoke[matched.control.subject.index])

#The outcome vector of treated and controls
R_t=matched_data$heart[matched_data$anger==1]
R_c=matched_data$heart[matched_data$anger==0]
R=cbind(R_t, R_c)
sum(R_t)
sum(R_c)

#McNemar's test under random experiment
E_null=(1/2)*(sum(R_t)+sum(R_c))
V_null=E_null-(1/4)*sum((R_t+R_c)^2)
T_mcnemar=sum(R_t)
q_null=(T_mcnemar-E_null)/sqrt(V_null)
1-pnorm(q_null)

mcnemartest<-function(R, GAMMA){
  I=nrow(R)
  prob_bar<-rep(0, I)
  for (i in 1:I){
    if (R[i,1]==0 & R[i,2]==0){
      prob_bar[i]=0
    } 
    if (R[i,1]==1 & R[i,2]==1){
      prob_bar[i]=1
    }
    if (R[i,1]!=R[i,2]){
      prob_bar[i]=(GAMMA)/(1+GAMMA)
    }
  }
  E_null=sum(prob_bar)
  V_null=sum(prob_bar*(1-prob_bar))
  Test=sum(R[,1])
  pvalue=1-pnorm((Test-E_null)/sqrt(V_null))
  return(pvalue)
}


######Lambda>1##############

table(matched_data$abuse)
#The outcomes for treated and controls given childhood matreatment indicator=1
R_t_1=matched_data$heart[matched_data$anger==1 & matched_data$abuse==1]
R_c_1=matched_data$heart[matched_data$anger==0 & matched_data$abuse==1]
R_1<-cbind(R_t_1, R_c_1)
#The outcomes for treated and controls given childhood maltreatment indicator=0
R_t_0=matched_data$heart[matched_data$anger==1 & matched_data$abuse==0]
R_c_0=matched_data$heart[matched_data$anger==0 & matched_data$abuse==0]
R_0<-cbind(R_t_0, R_c_0)

#return the expectation and variance for McNemar's test in a sensitivity analysis
E_mcnemar<-function(R, GAMMA){
  I=nrow(R)
  prob_bar<-rep(0, I)
  for (i in 1:I){
    if (R[i,1]==0 & R[i,2]==0){
      prob_bar[i]=0
    } 
    if (R[i,1]==1 & R[i,2]==1){
      prob_bar[i]=1
    }
    if (R[i,1]!=R[i,2]){
      prob_bar[i]=(GAMMA)/(1+GAMMA)
    }
  }
  E_null=sum(prob_bar)
  V_null=sum(prob_bar*(1-prob_bar))
  E_g<-rep(0,2)
  E_g[1]=E_null
  E_g[2]=V_null
  return(E_g)
}

#Return the worst-case p-value from Corollary 2
mcnemar_extend<-function(R, R_1, R_0, Gamma, lambda){
  if (lambda==1){
    pv=mcnemartest(R, GAMMA = Gamma)
  }
  if (lambda>1){
    E_null=E_mcnemar(R_1, GAMMA = Gamma)[1]+E_mcnemar(R_0, GAMMA = (Gamma)^(1/lambda))[1]
    V_null=E_mcnemar(R_1, GAMMA = Gamma)[2]+E_mcnemar(R_0, GAMMA = (Gamma)^(1/lambda))[2]
    pv=1-pnorm((sum(R[,1])-E_null)/sqrt(V_null))
  }
  if (lambda<1){
    E_null=E_mcnemar(R_1, GAMMA = (Gamma)^lambda)[1]+E_mcnemar(R_0, GAMMA = Gamma)[1]
    V_null=E_mcnemar(R_1, GAMMA = (Gamma)^lambda)[2]+E_mcnemar(R_0, GAMMA = Gamma)[2]
    pv=1-pnorm((sum(R[,1])-E_null)/sqrt(V_null))
  }
  return(pv)
}


gamma=seq(from=1, to = 3, by = 0.01)
lambda=c(1/8, 1/4, 1/2, 1, 2, 4, 8)
p_value=matrix(0, nrow = length(gamma), ncol = length(lambda))
for (i in 1:nrow(p_value)){
  for (j in 1:ncol(p_value)){
    p_value[i,j]=mcnemar_extend(R, R_1, R_0, gamma[i], lambda[j])
  }
}

row.names(p_value)=gamma
colnames(p_value)=lambda 

#Report Table 1 in the article
p_value_report=p_value[c(32, 38, 43, 45, 53, 82, 112), ]

