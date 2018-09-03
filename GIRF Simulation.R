# read in the data
D = as.matrix((read.csv("DistMatrix.csv")))[18:20,18:20]
N = as.matrix.data.frame(read.csv("20pop.csv"))[,c(1,19:21)]
B = as.matrix.data.frame(read.csv("20births.csv"))[,c(1,19:21)]
y = data.frame(read.csv("20measles.csv"))[,c(1,19:21)]
year = y[,1]

tau = 4
M = 52 #times
J = 100 #particles
SSS = 15 #total interval points in a unit time
X = array(0,dim = c((M+1),(SSS+1),J,3*5))
X.true = array(0,dim = c((M+1),(SSS+1),3*5))
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
h_eff = function(t,sss){
  d = (year[t]-start.year) %% 1 * 365.25 + (sss-1)/SSS
  if (d<7){1}
  else if (d>99 & d<116){1}
  else if (d>198 & d<253){1}
  else if (d>299 & d<309){1}
  else if (d>355){1}
  else {0}
}

# recruitment effect (birth)
r = function(k,t,sss) { 
  d = function(t,sss) (year[t]-start.year) %% 1 * 365.25 + (sss-1)/SSS
  if(d(t,sss) < school.start.day && school.start.day < d(t,sss+1)){
    c * as.numeric(B[round(year[t]-start.year-enrty.age+0.5),(k+1)])
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

# assessment function :B = 2 in this case
u = function(m,sss,j) {
  if((sss-1)<SSS){
    inter1 = 1-(sss-1)/SSS
    Int1 = function(z,k) dnbinom(z,inter1/sigma2,1/(sigma2*v.ir*X[m,sss,j,3+(k-1)*5])) *  dmeas(y.true[m,k],rho,phi,z+sum(X[m,2:(sss),j,4+(k-1)*5]))
    b1.1=sum(sapply(0:10000, function(z) Int1(z,1)))
    b1.2=sum(sapply(0:10000, function(z) Int1(z,2)))
    b1.3=sum(sapply(0:10000, function(z) Int1(z,3)))
    
    inter2 = 2-(sss-1)/SSS
    Int2 = function(z,k) dnbinom(z,inter2/sigma2,1/(sigma2*v.ir*X[m,sss,j,3+(k-1)*5])) *  dmeas(y.true[(m+1),k],rho,phi,z+sum(X[m,2:(sss),j,4+(k-1)*5])) 
    b2.1=sum(sapply(0:10000, function(z) Int2(z,1)))
    b2.2=sum(sapply(0:10000, function(z) Int2(z,2)))
    b2.3=sum(sapply(0:10000, function(z) Int2(z,3)))
    
    (b1.1*b1.2*b1.3) * (b2.1*b2.2*b2.3)
  }else{
    b1=prod(sapply(1:3,function(k) dmeas(y.true[m,k],rho,phi,sum(X[m,2:(SSS+1),j,4+(k-1)*5])))) 
    
    inter2 = 2-(sss-1)/SSS
    Int2 = function(z,k) dnbinom(z,inter2/sigma2,1/(sigma2*v.ir*X[m,sss,j,3+(k-1)*5])) * dmeas(y.true[(m+1),k],rho,phi,z+sum(X[m,2:(sss),j,4+(k-1)*5])) 
    b2.1=sum(sapply(0:10000, function(z) Int2(z,1)))
    b2.2=sum(sapply(0:10000, function(z) Int2(z,2)))
    b2.3=sum(sapply(0:10000, function(z) Int2(z,3)))
    b1 * (b2.1*b2.2*b2.3)
  }
}

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
          ((X[m,sss,j,3+(l-1)*5]/as.numeric(N[round(year[t]-start.year+0.5),(l+1)]))^alpha - (X[m,sss,j,3+(k-1)*5]/as.numeric(N[round(year[t]-start.year+0.5),(k+1)]))^alpha)
      }
    }
    
    S = X[m,sss,j,1+(k-1)*5]
    E = X[m,sss,j,2+(k-1)*5]
    I = X[m,sss,j,3+(k-1)*5]
    N.ir = X[m,sss,j,4+(k-1)*5]
    P = X[m,sss,j,5+(k-1)*5]
    
    v.se = ifelse(h_eff(t,sss)==1,h(a),1) * Beta.s(R0,v.ir,a) * ((I /P)^alpha + sum(ss))
    
    gamma_SE = rgamma(1,shape = (1/SSS)/sigma2,scale = sigma2)
    gamma_SD = rgamma(1,shape = (1/SSS)/sigma2,scale = sigma2)
    gamma_EI = rgamma(1,shape = (1/SSS)/sigma2,scale = sigma2)
    gamma_ED = rgamma(1,shape = (1/SSS)/sigma2,scale = sigma2)
    gamma_IR = rgamma(1,shape = (1/SSS)/sigma2,scale = sigma2)
    gamma_ID = rgamma(1,shape = (1/SSS)/sigma2,scale = sigma2)
    
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
    
    S.new = S + round(r(k,t,sss)) - dN.se - dN.sd
    E.new = E + dN.se - dN.ei - dN.ed
    I.new = I + dN.ei - dN.ir - dN.id
    N.ir.new = dN.ir
    P.new  = N[round(year[T]-start.year+0.5),(k+1)]
    X[m,sss+1,j,1+(k-1)*5] = S.new
    X[m,sss+1,j,2+(k-1)*5] = E.new
    X[m,sss+1,j,3+(k-1)*5] = I.new
    X[m,sss+1,j,4+(k-1)*5] = N.ir.new
    X[m,sss+1,j,5+(k-1)*5] = P.new
  }
X[m,sss+1,j,]
}

# define parameters and initial parameters
R0 = 20 
a = 0.163
G = 500   
alpha = 0.97 
v.ei = 1
v.ir = 1     
mu = 0.016/52
rho = 0.5  
phi = 0.25 
s = 0.04
e = 0.00027
i = 0.00032

#generate the true trajectory of X and Y
#for (k in 1:3) {
#  X.true[1,1,1+(k-1)*5] = round(s*N[round(year[T]-start.year+0.5),(k+1)]) #S
#  X.true[1,1,2+(k-1)*5] = round(e*N[round(year[T]-start.year+0.5),(k+1)]) #E
#  X.true[1,1,3+(k-1)*5] = round(i*N[round(year[T]-start.year+0.5),(k+1)]) #I
#  X.true[1,1,4+(k-1)*5] = 0   #delta N.ir(weekly)
#  X.true[1,1,5+(k-1)*5] = N[round(year[T]-start.year+0.5),(k+1)] #P 
#}

#for (m in 1:M) {
#  for (sss in 1:SSS) {
#  X.true[m,sss+1,] = sir.step(T + m -1,R0,a,alpha,mu,v.ei,v.ir,sigma2,rho,phi,G,X.true) # modify the sir.step function when we generate the true trajecory
#  }
#  X.true[m+1,1,] = X.true[m,SSS+1,]
#  y.true[m,] = sapply(1:3, function(k) round(rmeas(rho,phi,X.true[m+1,1,4+(k-1)*5])))
#}

set.seed(2018)
w = array(0,dim = c(M,SSS,J))
w.o = array(0,dim = c(M,SSS,J))
tolerance = 10^(-30)
L = array(0,dim = c(M))
u.old = numeric(J)
u.new = numeric(J)

# initial X
for (k in 1:3) {
  X[1,1,,1+(k-1)*5] = rep(round(s*N[round(year[T]-start.year+0.5),(k+1)]),J) #S
  X[1,1,,2+(k-1)*5] = rep(round(e*N[round(year[T]-start.year+0.5),(k+1)]),J) #E
  X[1,1,,3+(k-1)*5] = rep(round(i*N[round(year[T]-start.year+0.5),(k+1)]),J) #I
  X[1,1,,4+(k-1)*5] = rep(0,J)   #delta N.ir(weekly)
  X[1,1,,5+(k-1)*5] = rep(N[round(year[T]-start.year+0.5),(k+1)],J) #P 
}


#Begin: iteratied filtering method 
for(m in 1:52){
  if(m==1){
    u.old = rep(1,J)
  }else{
    div = sapply(1:J,function(j) prod(sapply(1:3, function(k) dmeas(y.true[m-1,k],rho,phi,sum(X[m-1,2:(SSS+1),j,4+(k-1)*5])))))
  }
  for (sss in 1:SSS) {
    #Propsed particles
    for (j in 1:J) {
      t = T + m -1
      #Proposed X(S,E,I,N.ir,P) using proposed parameters
      X[m,sss+1,j,] = sir.step(t,R0,a,alpha,mu,v.ei,v.ir,sigma2,rho,phi,G,X)
      u.new[j] = u(m,sss+1,j)
      error.time = 0
      while (u.new[j]==0) {
        X[m,sss+1,j,] = sir.step(t,R0,a,alpha,mu,v.ei,v.ir,sigma2,rho,phi,G,X)
        u.new[j] = u(m,sss+1,j)
        error.time = error.time + 1
      }
      w.o[m,sss,j] = ifelse(m==1,u.new[j]/u.old[j],u.new[j]/(u.old[j]/div[j]))}
      #calculate true weight
      for (j in 1:J) {
        w[m,sss,j] = ifelse(w.o[m,sss,j]==0,0,w.o[m,sss,j]/sum(w.o[m,sss,]))
        #Filtering failure: at some time point, the conditional likelihood of 
        if(w[m,sss,j] > tolerance){
          w[m,sss,j] = w[m,sss,j]
        }else{
          w[m,sss,j] = tolerance
        }
      }
        #Resample
        number = as.numeric(rmultinom(1,J,w[m,sss,]))
        number[0] = 0
        for(ii in 1:J){
          if((number[ii]==0)==FALSE){
            for(iii in 1:number[ii]){
              n = sum(number[0:(ii-1)])
              for (k in 1:3) {
                X[m,(sss+1),iii+n,1+(k-1)*5]=X[m,(sss+1),ii,1+(k-1)*5]
                X[m,(sss+1),iii+n,2+(k-1)*5]=X[m,(sss+1),ii,2+(k-1)*5]
                X[m,(sss+1),iii+n,3+(k-1)*5]=X[m,(sss+1),ii,3+(k-1)*5]
                X[m,(sss+1),iii+n,4+(k-1)*5]=X[m,(sss+1),ii,4+(k-1)*5]
                X[m,(sss+1),iii+n,5+(k-1)*5]=X[m,(sss+1),ii,5+(k-1)*5]
              }
              u.old[iii+n] = u.new[ii]}}}}
  
  #loglikelihood for time m
  L[m] = sum(log(rowMeans(w.o[m,,]))) 
  
  for (k in 1:3) {
    X[m+1,1,,1+(k-1)*5] = X[m,SSS+1,,1+(k-1)*5]
    X[m+1,1,,2+(k-1)*5] = X[m,SSS+1,,2+(k-1)*5]
    X[m+1,1,,3+(k-1)*5] = X[m,SSS+1,,3+(k-1)*5]
    X[m+1,1,,4+(k-1)*5] = X[m,SSS+1,,4+(k-1)*5]
    X[m+1,1,,5+(k-1)*5] = X[m,SSS+1,,5+(k-1)*5]
  }
  
}


par(mfrow=c(2,2))
plot(rowMeans(X[1:53,1,,11]),ylim=c(50000,250000),type = 'l',main = 'London',xlab = 'Times(week)',ylab='Susceptible Population')
lines(sapply(1:M, function(m) median(X[m,1,,11])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,1,,11],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,1,,11],0.9)),col=3,lty=2)
lines(X.true[1:53,11],type='l',col=4,lty=1)
legend("topright", bty="n",legend=c("Est. mean","Median",'10%,90% quantile',"True S"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

plot(rowMeans(X[1:53,1,,12]),ylim=c(0,50000),type = 'l',main = 'London',xlab = 'Times(week)',ylab='Exposed Population')
lines(sapply(1:M, function(m) median(X[m,1,,12])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,1,,12],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,1,,12],0.9)),col=3,lty=2)
lines(X.true[1:53,12],type='l',col=4,lty=1)
legend("topright", bty="n",legend=c("Est. mean","Median",'10%,90% quantile',"True E"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

plot(rowMeans(X[1:53,1,,13]),ylim=c(0,25000),type = 'l',main = 'London',xlab = 'Times(week)',ylab='Infectious Population')
lines(sapply(1:M, function(m) median(X[m,1,,13])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,1,,13],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,1,,13],0.9)),col=3,lty=2)
lines(X.true[1:53,8],type='l',col=4,lty=1)
legend("topright", bty="n",legend=c("Est. mean","Median",'10%,90% quantile',"True I"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

Nir = sapply(1:(M+1),function(m) sapply(1:J,function(j) sum(X[m,2:(SSS+1),j,14]))) 
plot(rowMeans(t(Nir)),ylim=c(0,7000),type = 'l',main = 'London',xlab = 'Times(week)',ylab='Weekly Recovered')
lines(sapply(1:M, function(m) median(t(Nir)[m,])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(t(Nir)[m,],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(t(Nir)[m,],0.9)),col=3,lty=2)
lines(X.true[1:53,14],type='l',col=4,lty=1)
legend("topright", bty="n",legend=c("Est. mean","Median",'10%,90% quantile',"True N_ir"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)


log.likelihood = sum(L)
log.likelihood
