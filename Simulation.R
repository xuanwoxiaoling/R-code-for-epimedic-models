# read in the data
setwd("/Users/qianwang/Desktop/Project/论文进度作业/R code for reference 2")
D = as.matrix((read.csv("DistMatrix.csv")))[18:20,18:20]
N = as.matrix.data.frame(read.csv("20pop.csv"))[,c(1,19:21)]
B = as.matrix.data.frame(read.csv("20births.csv"))[,c(1,19:21)]
y = data.frame(read.csv("20measles.csv"))[,c(1,19:21)]
year = y[,1]

tau = 4
M = 260 #times
J = 500 #particles
X = array(0,dim = c((M+1),J,3*5))
X.true = array(0,dim = c((M+1),3*5))
y.true = array(0,dim = c(M,3))

# Fixed parameters
sigma2 = 0.08
c = 0.4
start.year = 1944
enrty.age = tau
school.start.day = 251
tol = 10^(-6)

T = (tau+1) * 52 + 2 # start_year:1949

# holiday effect
p = 0.739
Beta.s = function(R0,v.ir,a) (1+2*(1-p)*a)*R0*v.ir
h = function(a) (1-2*p*a)/(1+2*(1-p)*a) #holiday coefficient
h_eff = function(t){
  d = (year[t]-start.year) %% 1 * 365.25
  if (d<7){1}
  else if (d>99 & d<116){1}
  else if (d>198 & d<253){1}
  else if (d>299 & d<309){1}
  else if (d>355){1}
  else {0}
}

# recruitment effect (birth)
r = function(k,t) { 
  d = function(t) (year[t]-start.year) %% 1 * 365.25
  if(d(t) < school.start.day && school.start.day < d(t+1)){
    c * as.numeric(B[round(year[t]-enrty.age-start.year+0.5),(k+1)])
  } else {(1-c) * as.numeric(B[round(year[t]-start.year-enrty.age+0.5),(k+1)] * 7/365.25)}
}

# measurement model 
dmeas <- function(y,rho,phi, H) {
  if(y>=0.5){
    pnorm(y+0.5,mean = rho*H,sd = sqrt(rho*(1-rho)*H + (phi*rho*H)^2+1)) - pnorm(y-0.5,mean = rho*H,sd = sqrt(rho*(1-rho)*H + (phi*rho*H)^2+1))
  }else {pnorm(0.5,mean = rho*H,sd = sqrt(rho*(1-rho)*H + (phi*rho*H)^2+1))}
}

rmeas <- function(rho,phi, H) {r=mean(round(rnorm(100,mean = rho*H,sd = sqrt(rho*(1-rho)*H + (phi*rho*H)^2+1))))
if(r<0){0
}else{r}}

# transmission model
sir.step <- function(t,R0,a,alpha,mu,v.ei,v.ir,sigma2,rho,phi,G,X){
  for (k in 1:3) {
    vkl = numeric(3)
    ss = numeric(3)
    for(l in 1:3){
      if(l==k){
        vkl[l] = 0
      }else{
        vkl[l] = round(G * (mean(D)/((mean(N[round(year[t]-start.year+0.5),]))^2)) *  
                         (as.numeric(N[round(year[t]-start.year+0.5),(k+1)]) * as.numeric(N[round(year[t]-start.year+0.5),(l+1)]))/D[k,l])
      }}
    tot_travel = sum(vkl)
    travel_prob = 1- exp(-tot_travel/as.numeric(N[round(year[t]-start.year+0.5),(k+1)]))
    
    for (l in 1:3) {
      if(l==k){
        ss[l] = 0
      }else{
        ss[l] = travel_prob * vkl[l]/tot_travel *
          ((X[m,j,3+(l-1)*5]/as.numeric(N[round(year[t]-start.year+0.5),(l+1)]))^alpha - 
             (X[m,j,3+(k-1)*5]/as.numeric(N[round(year[t]-start.year+0.5),(k+1)]))^alpha)
      }
    }
    
    S = X[m,j,1+(k-1)*5]
    E = X[m,j,2+(k-1)*5]
    I = X[m,j,3+(k-1)*5]
    N.ir = X[m,j,4+(k-1)*5]
    P = X[m,j,5+(k-1)*5]
    
    v.se = ifelse(h_eff(t)==1,h(a),1) * Beta.s(R0,v.ir,a)  * ((I /P)^alpha + sum(ss))
    
    gamma_SE = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_SD = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_EI = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_ED = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_IR = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_ID = rgamma(1,shape = 1/sigma2,scale = sigma2)
    
    S.Out = v.se * gamma_SE + mu * gamma_SD
    p.se = v.se * gamma_SE/  S.Out * (1-exp(-S.Out))
    p.sd = mu * gamma_SD/ S.Out * (1-exp(-S.Out))
    dN.se = rbinom(1, S, p.se)
    dN.sd = ifelse(S-dN.se==0,0,rbinom(1, S-dN.se, p.sd/(exp(-S.Out)+p.sd))) # because p.sd/(1-p.se) = p.sd/(exp(-S.Out)+p.sd)
    
    E.Out = v.ei * gamma_EI + mu * gamma_ED
    p.ei = v.ei * gamma_EI/  E.Out * (1-exp(-E.Out))
    p.ed = mu * gamma_ED/ E.Out * (1-exp(-E.Out))
    dN.ei = rbinom(1, E, p.ei)
    dN.ed = ifelse(E-dN.ei==0,0,rbinom(1, E-dN.ei, p.ed/(exp(-E.Out)+p.ed))) # because p.ed/(1-p.ei) = p.ed/(exp(-E.Out)+p.ed)
    
    I.Out = v.ir * gamma_IR + mu * gamma_ID
    p.ir = v.ir * gamma_IR/  I.Out * (1-exp(-I.Out))
    p.id = mu * gamma_ID/ I.Out * (1-exp(-I.Out))
    dN.ir = rbinom(1, I, p.ir)
    dN.id = ifelse(I-dN.ir==0,0,rbinom(1, I-dN.ir, p.id/(exp(-I.Out)+p.id))) # because p.id/(1-p.ir) = p.id/(exp(-I.Out)+p.id)
    
    S.new = S + round(r(k,t)) - dN.se - dN.sd
    E.new = E + dN.se - dN.ei - dN.ed
    I.new = I + dN.ei - dN.ir - dN.id
    N.ir.new = dN.ir
    P.new  = N[round(year[T]-start.year+0.5),(k+1)]
    
    X[m+1,j,1+(k-1)*5] = S.new
    X[m+1,j,2+(k-1)*5] = E.new
    X[m+1,j,3+(k-1)*5] = I.new
    X[m+1,j,4+(k-1)*5] = N.ir.new
    X[m+1,j,5+(k-1)*5] = P.new
  }
  X[m+1,j,]
}

# define parameters and initial parameters
R0 = 20
a = 0.163
G = 500   
alpha = 0.97 
v.ei = 1
v.ir = 1     
mu = 0.00032
rho = 0.5  
phi = 0.25 
s = 0.004
e = 0.00027
i = 0.00032

#generate the true trajectory of X and Y
#for (k in 1:3) {
#  X.true[1,1+(k-1)*5] = round(s*N[round(year[T]-start.year+0.5),(k+1)]) #S
#  X.true[1,2+(k-1)*5] = round(e*N[round(year[T]-start.year+0.5),(k+1)]) #E
#  X.true[1,3+(k-1)*5] = round(i*N[round(year[T]-start.year+0.5),(k+1)]) #I
#  X.true[1,4+(k-1)*5] = 0   #delta N.ir(weekly)
#  X.true[1,5+(k-1)*5] = N[round(year[T]-start.year+0.5),(k+1)] #P 
#}

#for (m in 1:M) {
#  X.true[m+1,] = sir.step(T + m -1,R0,a,alpha,mu,v.ei,v.ir,sigma2,rho,phi,G,X.true) # modify the sir.step function when we generate the true trajecory
#  y.true[m,] = sapply(1:3, function(k) round(rmeas(rho,phi,X.true[m+1,4+(k-1)*5])))
#}


set.seed(2018)
w = array(0,dim = c(M,J))
w.o = array(0,dim = c(M,J))
tolerance = 10^(-30)

# initial X
for (k in 1:3) {
  X[1,,1+(k-1)*5] = rep(round(s*N[round(year[T]-start.year+0.5),(k+1)]),J) #S
  X[1,,2+(k-1)*5] = rep(round(e*N[round(year[T]-start.year+0.5),(k+1)]),J) #E
  X[1,,3+(k-1)*5] = rep(round(i*N[round(year[T]-start.year+0.5),(k+1)]),J) #I
  X[1,,4+(k-1)*5] = rep(0,J)   #delta N.ir(weekly)
  X[1,,5+(k-1)*5] = rep(N[round(year[T]-start.year+0.5),(k+1)],J) #P 
}
# iteration
for(m in 1:M){ 
  for (j in 1:J) {
    t = T + m -1
    #Proposed X(S,E,I,N.ir,P) using proposed parameters
    X[m+1,j,] = sir.step(t,R0,a,alpha,mu,v.ei,v.ir,sigma2,rho,phi,G,X)
    w.o[m,j] = prod(sapply(1:3,function(k) dmeas(y.true[m,k],rho,phi, X[m+1,j,4+(k-1)*5])))}
  #calculate true weight
  for (j in 1:J) {
    w[m,j] = ifelse(w.o[m,j]==0,0,w.o[m,j]/sum(w.o[m,]))
    #Filtering failure: at some time point, the conditional likelihood of 
    if(w[m,j] > tolerance){
      w[m,j] = w[m,j]
    }else{
      w[m,j] = tolerance
    }}
  #Resample
  number = as.numeric(rmultinom(1,J,w[m,]))
  number[0] = 0
  for(ii in 1:J){
    if((number[ii]==0)==FALSE){
      for(iii in 1:number[ii]){
        n = sum(number[0:(ii-1)])
        for (k in 1:3) {
          X[m+1,iii+n,1+(k-1)*5]=X[m+1,ii,1+(k-1)*5]
          X[m+1,iii+n,2+(k-1)*5]=X[m+1,ii,2+(k-1)*5]
          X[m+1,iii+n,3+(k-1)*5]=X[m+1,ii,3+(k-1)*5]
          X[m+1,iii+n,4+(k-1)*5]=X[m+1,ii,4+(k-1)*5]
          X[m+1,iii+n,5+(k-1)*5]=X[m+1,ii,5+(k-1)*5]
        }
      }}}
}

log.likelihood = sum(log(rowMeans(w.o)))


par(mfrow=c(2,2))
plot(rowMeans(X[,,11]),type = 'l',main = 'London',xlab = 'Times(week)',ylab='Susceptible Population')
lines(sapply(1:M, function(m) median(X[m,,11])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,11],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,11],0.9)),col=3,lty=2)
lines(X.true[,11],type='l',col=4,lty=1)
legend("topright", bty="n",legend=c("Est. mean","Median",'10%,90% quantile',"True S"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

plot(rowMeans(X[,,12]),ylim=c(0,10000),type = 'l',main = 'London',xlab = 'Times(week)',ylab='Exposed Population')
lines(sapply(1:M, function(m) median(X[m,,12])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,12],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,12],0.9)),col=3,lty=2)
lines(X.true[,12],type='l',col=4,lty=1)
legend("topright", bty="n",legend=c("Est. mean","Median",'10%,90% quantile',"True E"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

plot(rowMeans(X[,,13]),ylim=c(0,9000),type = 'l',main = 'London',xlab = 'Times(week)',ylab='Infectious Population')
lines(sapply(1:M, function(m) median(X[m,,13])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,13],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,13],0.9)),col=3,lty=2)
lines(X.true[,8],type='l',col=4,lty=1)
legend("topright", bty="n",legend=c("Est. mean","Median",'10%,90% quantile',"True I"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

plot(rowMeans(X[,,14]),ylim=c(0,7000),type = 'l',main = 'London',xlab = 'Times(week)',ylab='Weekly Recovered')
lines(sapply(1:M, function(m) median(X[m,,14])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,14],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,14],0.9)),col=3,lty=2)
lines(X.true[,14],type='l',col=4,lty=1)
legend("topright", bty="n",legend=c("Est. mean","Median",'10%,90% quantile',"True N_ir"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)



par(mfrow=c(2,2))
plot(sapply(1:M,function(m) var(X[m,,11])),type='l',main='London',xlab = 'Times(week)',ylab='Variance of Susceptible Population')
plot(sapply(1:M,function(m) var(X[m,,12])),type='l',main='London',xlab = 'Times(week)',ylab='Variance of Exposed Population')
plot(sapply(1:M,function(m) var(X[m,,13])),type='l',main='London',xlab = 'Times(week)',ylab='Variance of Infectious Population')
plot(sapply(1:M,function(m) var(X[m,,14])),type='l',main='London',xlab = 'Times(week)',ylab='Variance of Weekly Recovered')

par(mfrow=c(1,1))
plot(y.true[,1],ylim=c(0,3500),type = 'l',main = 'Assigned Data',xlab = 'Times(week)',ylab='Weekly reported cases',col=1)
lines(y.true[,2],col=2)
lines(y.true[,3],col=3)
legend("topright", legend=c("Liverpool","Birmingham",'London'),col=c(1,2,3), lty=c(1,1,1), cex=0.7)

par(mfrow=c(1,1))
plot(y[,2],ylim=c(0,2000),type = 'l',main = 'Assigned Data',xlab = 'Times(week)',ylab='Weekly reported cases',col=1)
lines(y[,3],col=2)
lines(y[,4],col=3)
legend("topright", legend=c("Liverpool","Birmingham",'London'),col=c(1,2,3), lty=c(1,1,1), cex=0.7)




