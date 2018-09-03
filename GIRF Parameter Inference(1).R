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
SSS = 15 #total interval points in a unit time
X = array(0,dim = c((M+1),(SSS+1),J,3*5))

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


set.seed(2018)

#theta.1 = array(0,dim = c((MM+1),(M+1),(SSS+1),J))
#theta.2 = array(0,dim = c((MM+1),(M+1),(SSS+1),J))
#theta.3 = array(0,dim = c((MM+1),(M+1),(SSS+1),J))
#theta.4 = array(0,dim = c((MM+1),(M+1),(SSS+1),J))
#theta.5 = array(0,dim = c((MM+1),(M+1),(SSS+1),J))
#theta.6 = array(0,dim = c((MM+1),(M+1),(SSS+1),J))
#theta.7 = array(0,dim = c((MM+1),(M+1),(SSS+1),J))
#theta.8 = array(0,dim = c((MM+1),(M+1),(SSS+1),J))
#theta.9 = array(0,dim = c((MM+1),(M+1),(SSS+1),J))
theta.10 = array(0,dim = c((MM+1),(M+1),(SSS+1),J))
theta.11 = array(0,dim = c((MM+1),(M+1),(SSS+1),J))
theta.12 = array(0,dim = c((MM+1),(M+1),(SSS+1),J))

#theta.1[1,1,1,] = runif(J,0,1)
#theta.2[1,1,1,] = runif(J,0,1)
#theta.3[1,1,1,] = runif(J,0,1)
#theta.4[1,1,1,] = runif(J,0,1)
#theta.5[1,1,1,] = runif(J,0,1)
#theta.6[1,1,1,] = runif(J,0,1)
#theta.7[1,1,1,] = runif(J,0,1)
#theta.8[1,1,1,] = runif(J,0,1)
#theta.9[1,1,1,] = runif(J,0,1)
theta.10[1,1,1,] = rnorm(J,-3,1)
theta.11[1,1,1,] = sapply(1:J,function(j) min(rnorm(1,-8,1),-theta.10[1,1,j]))
theta.12[1,1,1,] = sapply(1:J,function(j) min(rnorm(1,-8,1),
                                            (log(1-exp(theta.10[1,1,j]+theta.11[1,1,j])) - log(exp(theta.10[1,1,j])+exp(theta.11[1,1,j])+2*exp(theta.10[1,1,j]+theta.11[1,1,j])) )
))

logit = function(x) 1/(1+exp(-x))
# define parameters and initial parameters
# paras = numeric(15)
R0 = 20   #exp(theta.1)
a = 0.163 #logit(theta.2)
G = 500   #exp(theta.3)
alpha = 0.97 #theta.4
v.ei = 1     #1/(exp(theta.5) * (1-logit(theta.6)))
v.ir = 1     #1/(exp(theta.5) * logit(theta.6))
mu = 0.00032 #theta.7
rho = 0.5    #logit(theta.8)
phi = 0.25  #exp(theta.9)

w = array(0,dim = c(M,SSS,J))
w.o = array(0,dim = c(M,SSS,J))
tolerance = 10^(-30)
u.old = numeric(J)
u.new = numeric(J)

for (mm in 1:MM) {
  s = sapply(1:J,function(j) logit(theta.10[mm,1,1,j]))
  e = sapply(1:J,function(j) logit(theta.11[mm,1,1,j]))
  i = sapply(1:J,function(j) logit(theta.12[mm,1,1,j]))
  #Begin: iteratied filtering method 
  for(m in 1:M){
    if(m==1){
      u.old = rep(1,J)
    }else{
      div = sapply(1:J,function(j) prod(sapply(1:3, function(k) dmeas(y.true[m-1,k],rho,phi,sum(X[m-1,2:(SSS+1),j,4+(k-1)*5])))))
    }
    for (sss in 1:SSS) {
      #Propsed particles
      for (j in 1:J) {
        t = T + m -1
        
        if(m==1 & sss==1){
          # initial X
          for (k in 1:3) {
            X[1,1,j,1+(k-1)*5] = round(s[j]*N[round(year[T]-start.year+0.5),(k+1)]) #S
            X[1,1,j,2+(k-1)*5] = round(e[j]*N[round(year[T]-start.year+0.5),(k+1)]) #E
            X[1,1,j,3+(k-1)*5] = round(i[j]*N[round(year[T]-start.year+0.5),(k+1)]) #I
            X[1,1,j,4+(k-1)*5] = 0   #delta N.ir(weekly)
            X[1,1,j,5+(k-1)*5] = N[round(year[T]-start.year+0.5),(k+1)] #P 
          }
        }
        
        #perturbation density for parameters is Norm with sd based on the geometric cooling scheme
        theta.10[mm,m,sss+1,] = sapply(1:J,function(j) theta.10[mm,m,sss,j] = theta.10[mm,m,sss,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(mm-1)*M)/(50*M)) * ifelse(mm > 20,0.002,0.007))))
        theta.11[mm,m,sss+1,] = sapply(1:J,function(j) theta.11[mm,m,sss,j] = min(theta.11[mm,m,sss,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(mm-1)*M)/(50*M)) * ifelse(mm>20,0.0003,0.005))),-theta.10[mm,m,sss,j]))
        theta.12[mm,m,sss+1,] = sapply(1:J,function(j) theta.12[mm,m,sss,j] = min(theta.12[mm,m,sss,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(mm-1)*M)/(50*M)) * ifelse(mm>20,0.0003,0.005))),
                                                                            (log(1-exp(theta.10[mm,m+1,sss,j]+theta.11[mm,m,sss,j])) - 
                                                                               log(exp(theta.10[mm,m+1,sss,j])+exp(theta.11[mm,m,sss,j]) + 2*exp(theta.10[mm,m,sss,j]+theta.11[mm,m,sss,j])))))
        
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
            n = sum(number[0:(i-1)])
            for (k in 1:3) {
              X[m,(sss+1),iii+n,1+(k-1)*5]=X[m,(sss+1),ii,1+(k-1)*5]
              X[m,(sss+1),iii+n,2+(k-1)*5]=X[m,(sss+1),ii,2+(k-1)*5]
              X[m,(sss+1),iii+n,3+(k-1)*5]=X[m,(sss+1),ii,3+(k-1)*5]
              X[m,(sss+1),iii+n,4+(k-1)*5]=X[m,(sss+1),ii,4+(k-1)*5]
              X[m,(sss+1),iii+n,5+(k-1)*5]=X[m,(sss+1),ii,5+(k-1)*5]
            }
            u.old[iii+n] = u.new[ii]
            theta.10[mm,m,sss+1,iii+n]=theta.10[mm,m,sss,ii]
            theta.11[mm,m,sss+1,iii+n]=theta.11[mm,m,sss,ii]
            theta.12[mm,m,sss+1,iii+n]=theta.12[mm,m,sss,ii]}}}}
    
    for (k in 1:3) {
      X[m+1,1,,1+(k-1)*5] = X[m,SSS+1,,1+(k-1)*5]
      X[m+1,1,,2+(k-1)*5] = X[m,SSS+1,,2+(k-1)*5]
      X[m+1,1,,3+(k-1)*5] = X[m,SSS+1,,3+(k-1)*5]
      X[m+1,1,,4+(k-1)*5] = X[m,SSS+1,,4+(k-1)*5]
      X[m+1,1,,5+(k-1)*5] = X[m,SSS+1,,5+(k-1)*5]
    }
    theta.10[mm,m+1,1,]=theta.10[mm,m,SSS+1,]
    theta.11[mm,m+1,1,]=theta.11[mm,m,SSS+1,]
    theta.12[mm,m+1,1,]=theta.12[mm,m,SSS+1,]
  }
  #Initial parameters for next iteration
  theta.10[mm+1,1,1,] = sapply(1:J,function(j) theta.10[mm+1,1,1,j] = theta.10[mm,M+1,1,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(mm-1)*M)/(50*M)) * ifelse(mm > 20,0.001,0.01))))
  theta.11[mm+1,1,1,] = sapply(1:J,function(j) theta.11[mm+1,1,1,j] = min(theta.11[mm,M+1,1,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(mm-1)*M)/(50*M)) * ifelse(mm > 20,0.0001,0.01))),-theta.10[mm,M+1,1,j]))
  theta.12[mm+1,1,1,] = sapply(1:J,function(j) theta.12[mm+1,1,1,j] = min(theta.12[mm,M+1,1,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(mm-1)*M)/(50*M)) * ifelse(mm > 20,0.0001,0.01))),
                                                                      (log(1-exp(theta.10[mm,M+1,1,j]+theta.11[mm,M+1,1,j])) - 
                                                                         log(exp(theta.10[mm,M+1,1,j])+exp(theta.11[mm,M+1,1,j]) + 2*exp(theta.10[mm,M+1,1,j]+theta.11[mm,M+1,1,j])))))
  
}

par(mfrow=c(1,3))
plot(rowMeans(logit(theta.10[,1,1,])),ylim=c(min(logit(theta.10)),0.2),type='l',col=1,main = 'Initial proportion of Susceptible',xlab = 'Iteraion times',ylab = 'I(S)')
lines(sapply(1:MM,function(mm) median(logit(theta.10[mm,1,1,]),0.1)),col=2,lty=1)
lines(sapply(1:MM,function(mm) quantile(logit(theta.10[mm,1,1,]),0.1)),col=3,lty=2)
lines(sapply(1:MM,function(mm) quantile(logit(theta.10[mm,1,1,]),0.9)),col=3,lty=2)
lines(sapply(1:MM,function(mm) 0.04),col=4,lty=2)
legend("topright", bty="n",legend=c("Est. I(S):mean","Est. median",'Est. 10%,90% quantile',"True Value"),col=c(1,2,3,4), lty=c(1,1,2,2), cex=0.7)

plot(rowMeans(logit(theta.11[,1,1,])),ylim=c(min(logit(theta.11)),0.0015),type='l',col=1,main = 'Initial proportion of Exposed',xlab = 'Iteraion times',ylab = 'I(E)')
lines(sapply(1:MM,function(mm) median(logit(theta.11[mm,1,1,]),0.1)),col=2,lty=1)
lines(sapply(1:MM,function(mm) quantile(logit(theta.11[mm,1,1,]),0.1)),col=3,lty=2)
lines(sapply(1:MM,function(mm) quantile(logit(theta.11[mm,1,1,]),0.9)),col=3,lty=2)
lines(sapply(1:MM,function(mm) 0.00027),col=4,lty=2)
legend("topright", bty ="n",legend=c("Est. I(S):mean","Est. median",'Est. 10%,90% quantile',"True Value"),col=c(1,2,3,4), lty=c(1,1,2,2), cex=0.7)

plot(rowMeans(logit(theta.12[,1,1,])),ylim=c(min(logit(theta.12)),0.0015),type='l',col=1,main = 'Initial proportion of Infectious',xlab = 'Iteraion times',ylab = 'I(I)')
lines(sapply(1:MM,function(mm) median(logit(theta.12[mm,1,1,]),0.1)),col=2,lty=1)
lines(sapply(1:MM,function(mm) quantile(logit(theta.12[mm,1,1,]),0.1)),col=3,lty=2)
lines(sapply(1:MM,function(mm) quantile(logit(theta.12[mm,1,1,]),0.9)),col=3,lty=2)
lines(sapply(1:MM,function(mm) 0.00032),col=4,lty=2)
legend("topright", bty ="n",legend=c("Est. I(S):mean","Est. median",'Est. 10%,90% quantile',"True Value"),col=c(1,2,3,4), lty=c(1,1,2,2), cex=0.7)
