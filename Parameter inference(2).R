# read in the data
setwd("/Users/qianwang/Desktop/Project/论文进度作业/R code for reference 2")
D = as.matrix((read.csv("DistMatrix.csv")))[18:20,18:20]
N = as.matrix.data.frame(read.csv("20pop.csv"))[,c(1,19:21)]
B = as.matrix.data.frame(read.csv("20births.csv"))[,c(1,19:21)]
y = data.frame(read.csv("20measles.csv"))[,c(1,19:21)]
year = y[,1]

tau = 4
MM = 4 #total iteration
M = 260 #times
J = 500 #particles
X = array(0,dim = c((MM+1),(M+1),J,3*5))

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
sir.step <- function(mm,t,R0,a,alpha,mu,v.ei,v.ir,sigma2,rho,phi,G,X){
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
          ((X[mm,m,j,3+(l-1)*5]/as.numeric(N[round(year[t]-start.year+0.5),(l+1)]))^alpha - 
           (X[mm,m,j,3+(k-1)*5]/as.numeric(N[round(year[t]-start.year+0.5),(k+1)]))^alpha)
      }
    }
    
    S = X[mm,m,j,1+(k-1)*5]
    E = X[mm,m,j,2+(k-1)*5]
    I = X[mm,m,j,3+(k-1)*5]
    N.ir = X[mm,m,j,4+(k-1)*5]
    P = X[mm,m,j,5+(k-1)*5]
    
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
    
    X[mm,m+1,j,1+(k-1)*5] = S.new
    X[mm,m+1,j,2+(k-1)*5] = E.new
    X[mm,m+1,j,3+(k-1)*5] = I.new
    X[mm,m+1,j,4+(k-1)*5] = N.ir.new
    X[mm,m+1,j,5+(k-1)*5] = P.new
  }
  X[mm,m+1,j,]
}


set.seed(2018)

theta.1 = array(0,dim = c((MM+1),(M+1),J))
theta.2 = array(0,dim = c((MM+1),(M+1),J))
theta.3 = array(0,dim = c((MM+1),(M+1),J))
theta.4 = array(0,dim = c((MM+1),(M+1),J))
theta.5 = array(0,dim = c((MM+1),(M+1),J))
theta.6 = array(0,dim = c((MM+1),(M+1),J))
theta.7 = array(0,dim = c((MM+1),(M+1),J))
theta.8 = array(0,dim = c((MM+1),(M+1),J))
theta.9 = array(0,dim = c((MM+1),(M+1),J))

theta.1[1,1,] = rnorm(J,3,0.5)
theta.2[1,1,] = rnorm(J,-1,1)
theta.3[1,1,] = rnorm(J,6,0.5)
theta.4[1,1,] = rnorm(J,0,1)
theta.5[1,1,] = rnorm(J,0.7,0.5)
theta.6[1,1,] = rnorm(J,0,0.5)
theta.7[1,1,] = rnorm(J,-6.8,0.5)
theta.8[1,1,] = rnorm(J,0,1)
theta.9[1,1,] = rnorm(J,-1,1)
s = 0.04
e = 0.00027
i = 0.00032

# initial X
for (k in 1:3) {
  X[mm,1,,1+(k-1)*5] = rep(round(s*N[round(year[T]-start.year+0.5),(k+1)]),J) #S
  X[mm,1,,2+(k-1)*5] = rep(round(e*N[round(year[T]-start.year+0.5),(k+1)]),J) #E
  X[mm,1,,3+(k-1)*5] = rep(round(i*N[round(year[T]-start.year+0.5),(k+1)]),J) #I
  X[mm,1,,4+(k-1)*5] = rep(0,J)   #delta N.ir(weekly)
  X[mm,1,,5+(k-1)*5] = rep(N[round(year[T]-start.year+0.5),(k+1)],J) #P 
}

logit = function(x) 1/(1+exp(-x))
w = array(0,dim = c(MM,M,J))
w.o = array(0,dim = c(MM,M,J))
tolerance = 10^(-30)

for (mm in 1:MM) {
  #Begin: iteratied filtering method 
  for(m in 1:M){ #when estimate s.e.i only use the data from initial couples of weeks we choose firt year 
    #perturbation density for parameters is Norm with sd based on the geometric cooling scheme
    theta.1[mm,m+1,] = sapply(1:J,function(j) theta.1[mm,m+1,j] = theta.1[mm,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.002)))
    theta.2[mm,m+1,] = sapply(1:J,function(j) theta.2[mm,m+1,j] = theta.2[mm,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.008)))
    theta.3[mm,m+1,] = sapply(1:J,function(j) theta.3[mm,m+1,j] = theta.3[mm,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.001)))
    theta.4[mm,m+1,] = sapply(1:J,function(j) theta.4[mm,m+1,j] = theta.4[mm,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.005)))
    theta.5[mm,m+1,] = sapply(1:J,function(j) theta.5[mm,m+1,j] = theta.5[mm,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.0003)))
    theta.6[mm,m+1,] = sapply(1:J,function(j) theta.6[mm,m+1,j] = theta.6[mm,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.0003)))
    theta.7[mm,m+1,] = sapply(1:J,function(j) theta.7[mm,m+1,j] = theta.7[mm,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.001)))
    theta.8[mm,m+1,] = sapply(1:J,function(j) theta.8[mm,m+1,j] = theta.8[mm,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.001)))
    theta.9[mm,m+1,] = sapply(1:J,function(j) theta.9[mm,m+1,j] = theta.9[mm,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.0005)))
    #Propsed particles
    for (j in 1:J) {
      t = T + m -1
      #Proposed X(S,E,I,N.ir,P) using proposed parameters
      # define parameters and initial parameters
      # paras = numeric(15)
      R0 = exp(theta.1[mm,m+1,j]) #20 
      a = -1/(2*(1-p))+logit(theta.2[mm,m+1,j])*1/(2*p)-(-1/(2*(1-p))) #0.163
      G = exp(theta.3[mm,m+1,j])#500   
      alpha = exp(theta.4[mm,m+1,j]) #0.97 
      v.ei = 1/(exp(theta.5[mm,m+1,j]) * (1-logit(theta.6[mm,m+1,j]))) #1
      v.ir = 1/(exp(theta.5[mm,m+1,j]) * logit(theta.6[mm,m+1,j])) #1     
      mu = logit(theta.7[mm,m+1,j])  #0.016/52
      rho = logit(theta.8[mm,m+1,j]) #0.5  
      phi = exp(theta.9[mm,m+1,j]) #0.25 
      X[mm,m+1,j,] = sir.step(mm,t,R0,a,alpha,mu,v.ei,v.ir,sigma2,rho,phi,G,X)
      w.o[mm,m,j] = prod(sapply(1:3,function(k) dmeas(y.true[m,k],rho,phi, X[mm,m+1,j,4+(k-1)*5])))}
    #calculate true weight
    for (j in 1:J) {
      w[mm,m,j] = ifelse(w.o[mm,m,j]==0,0,w.o[mm,m,j]/sum(w.o[mm,m,]))
      #Filtering failure: at some time point, the conditional likelihood of 
      if(w[mm,m,j] > tolerance){
        w[mm,m,j] = w[mm,m,j]
      }else{
        w[mm,m,j] = tolerance
      }}
    #Resample
    number = as.numeric(rmultinom(1,J,w[mm,m,]))
    number[0] = 0
    for(ii in 1:J){
      if((number[ii]==0)==FALSE){
        for(iii in 1:number[ii]){
          n = sum(number[0:(ii-1)])
          for (k in 1:3) {
            X[mm,m+1,iii+n,1+(k-1)*5]=X[mm,m+1,ii,1+(k-1)*5]
            X[mm,m+1,iii+n,2+(k-1)*5]=X[mm,m+1,ii,2+(k-1)*5]
            X[mm,m+1,iii+n,3+(k-1)*5]=X[mm,m+1,ii,3+(k-1)*5]
            X[mm,m+1,iii+n,4+(k-1)*5]=X[mm,m+1,ii,4+(k-1)*5]
            X[mm,m+1,iii+n,5+(k-1)*5]=X[mm,m+1,ii,5+(k-1)*5]
          }
          theta.1[mm,m+1,iii+n]=theta.1[mm,m+1,ii]
          theta.2[mm,m+1,iii+n]=theta.2[mm,m+1,ii]
          theta.3[mm,m+1,iii+n]=theta.3[mm,m+1,ii]
          theta.4[mm,m+1,iii+n]=theta.4[mm,m+1,ii]
          theta.5[mm,m+1,iii+n]=theta.5[mm,m+1,ii]
          theta.6[mm,m+1,iii+n]=theta.6[mm,m+1,ii]
          theta.7[mm,m+1,iii+n]=theta.7[mm,m+1,ii]
          theta.8[mm,m+1,iii+n]=theta.8[mm,m+1,ii]
          theta.9[mm,m+1,iii+n]=theta.9[mm,m+1,ii]
        }}}
  }
  #Initial parameters for next iteration
  theta.1[mm+1,1,] = sapply(1:J,function(j) rnorm(1,mean(theta.1[mm,M+1,j]),sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.02)))
  theta.2[mm+1,1,] = sapply(1:J,function(j) rnorm(1,mean(theta.2[mm,M+1,j]),sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.008)))
  theta.3[mm+1,1,] = sapply(1:J,function(j) rnorm(1,mean(theta.3[mm,M+1,j]),sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.001)))
  theta.4[mm+1,1,] = sapply(1:J,function(j) rnorm(1,mean(theta.4[mm,M+1,j]),sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.005)))
  theta.5[mm+1,1,] = sapply(1:J,function(j) rnorm(1,mean(theta.5[mm,M+1,j]),sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.0003)))
  theta.6[mm+1,1,] = sapply(1:J,function(j) rnorm(1,mean(theta.6[mm,M+1,j]),sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.0003)))
  theta.7[mm+1,1,] = sapply(1:J,function(j) rnorm(1,mean(theta.7[mm,M+1,j]),sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.001)))
  theta.8[mm+1,1,] = sapply(1:J,function(j) rnorm(1,mean(theta.8[mm,M+1,j]),sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.001)))
  theta.9[mm+1,1,] = sapply(1:J,function(j) rnorm(1,mean(theta.9[mm,M+1,j]),sd = sqrt(0.92^(((m)+(mm-1)*M)/(50*M)) * 0.001)))
  
  #Initial X with new parameters, prepare for next iteration 
  for(k in 1:3){
    X[mm+1,1,,1+(k-1)*5] = rep(round(s*N[round(year[T]-start.year+0.5),(k+1)]),J) #S
    X[mm+1,1,,2+(k-1)*5] = rep(round(e*N[round(year[T]-start.year+0.5),(k+1)]),J) #E
    X[mm+1,1,,3+(k-1)*5] = rep(round(i*N[round(year[T]-start.year+0.5),(k+1)]),J) #I
    X[mm+1,1,,4+(k-1)*5] = rep(0,J)   #delta N.ir(weekly)
    X[mm+1,1,,5+(k-1)*5] = rep(N[round(year[T]-start.year+0.5),(k+1)],J) #P
  }
}

plot(rowMeans(exp(theta.1[mm,,])),ylim=c(0,100),type='l',col=1,main = 'R0',xlab = 'Iteraion times',ylab = 'Iterated R0')
lines(sapply(1:M,function(m) median(exp(theta.1[mm,m,]))),col=2,lty=1)
lines(sapply(1:M,function(m) quantile(exp(theta.1[mm,m,]),0.1)),col=3,lty=2)
lines(sapply(1:M,function(m) quantile(exp(theta.1[mm,m,]),0.9)),col=3,lty=2)
legend("topright", bty="n",legend=c("Est. R0:mean","Est. median",'Est. 10%,90% quantile'),col=c(1,2,3), lty=c(1,1,2), cex=0.7)

plot(rowMeans(-1/(2*(1-p))+logit(theta.2[mm,,])*1/(2*p)-(-1/(2*(1-p)))),ylim=c(0.05,0.75),type='l',col=1,main = 'a',xlab = 'Iteraion times',ylab = 'Iterated a')
lines(sapply(1:M,function(m) median(-1/(2*(1-p))+logit(theta.2[mm,m,])*1/(2*p)-(-1/(2*(1-p))))),col=2,lty=1)
lines(sapply(1:M,function(m) quantile(-1/(2*(1-p))+logit(theta.2[mm,m,])*1/(2*p)-(-1/(2*(1-p))),0.1)),col=3,lty=2)
lines(sapply(1:M,function(m) quantile(-1/(2*(1-p))+logit(theta.2[mm,m,])*1/(2*p)-(-1/(2*(1-p))),0.9)),col=3,lty=2)
legend("topright", bty ="n",legend=c("Est. a:mean","Est. median",'Est. 10%,90% quantile'),col=c(1,2,3), lty=c(1,1,2), cex=0.7)

plot(rowMeans(exp(theta.3[mm,,])),ylim=c(200,1200),type='l',col=1,main = 'G',xlab = 'Iteraion times',ylab = 'Iterated G')
lines(sapply(1:M,function(m) median(exp(theta.3[mm,m,]),0.1)),col=2,lty=1)
lines(sapply(1:M,function(m) quantile(exp(theta.3[mm,m,]),0.1)),col=3,lty=2)
lines(sapply(1:M,function(m) quantile(exp(theta.3[mm,m,]),0.9)),col=3,lty=2)
legend("topright", bty ="n",legend=c("Est. G:mean","Est. median",'Est. 10%,90% quantile'),col=c(1,2,3), lty=c(1,1,2), cex=0.7)

plot(rowMeans(exp(theta.4[mm,,])),ylim=c(0,7),type='l',col=1,main = 'alpha',xlab = 'Iteraion times',ylab = 'Iterated alpha')
lines(sapply(1:M,function(m) median(exp(theta.4[mm,m,]),0.1)),col=2,lty=1)
lines(sapply(1:M,function(m) quantile(exp(theta.4[mm,m,]),0.1)),col=3,lty=2)
lines(sapply(1:M,function(m) quantile(exp(theta.4[mm,m,]),0.9)),col=3,lty=2)
legend("topright", bty ="n",legend=c("Est. alpha:mean","Est. median",'Est. 10%,90% quantile'),col=c(1,2,3), lty=c(1,1,2), cex=0.7)

plot(rowMeans(1/(exp(theta.5[mm,,]) * (1-logit(theta.6[mm,,])))),ylim=c(0.2,4),type='l',col=1,main = 'v.ei',xlab = 'Iteraion times',ylab = 'Iterated v.ei')
lines(sapply(1:M,function(m) median(1/(exp(theta.5[mm,m,]) * (1-logit(theta.6[1,m,]))))),col=2,lty=1)
lines(sapply(1:M,function(m) quantile(1/(exp(theta.5[mm,m,]) * (1-logit(theta.6[1,m,]))),0.1)),col=3,lty=2)
lines(sapply(1:M,function(m) quantile(1/(exp(theta.5[mm,m,]) * (1-logit(theta.6[1,m,]))),0.9)),col=3,lty=2)
legend("topright", bty ="n",legend=c("Est. v.ei:mean","Est. median",'Est. 10%,90% quantile'),col=c(1,2,3), lty=c(1,1,2), cex=0.7)

plot(rowMeans(1/(exp(theta.5[mm,,]) * (logit(theta.6[mm,,])))),ylim=c(0.2,4.2),type='l',col=1,main = 'v.ir',xlab = 'Iteraion times',ylab = 'Iterated v.ir')
lines(sapply(1:M,function(m) median(1/(exp(theta.5[mm,m,]) * (logit(theta.6[1,m,]))))),col=2,lty=1)
lines(sapply(1:M,function(m) quantile(1/(exp(theta.5[mm,m,]) * (logit(theta.6[1,m,]))),0.1)),col=3,lty=2)
lines(sapply(1:M,function(m) quantile(1/(exp(theta.5[mm,m,]) * (logit(theta.6[1,m,]))),0.9)),col=3,lty=2)
legend("topright", bty ="n",legend=c("Est. v.ir:mean","Est. median",'Est. 10%,90% quantile'),col=c(1,2,3), lty=c(1,1,2), cex=0.7)

plot(rowMeans(logit(theta.7[mm,,])),ylim=c(0,0.01),type='l',col=1,main = 'mu',xlab = 'Iteraion times',ylab = 'Iterated mu')
lines(sapply(1:M,function(m) median(logit(theta.7[mm,m,]))),col=2,lty=1)
lines(sapply(1:M,function(m) quantile(logit(theta.7[mm,m,]),0.1)),col=3,lty=2)
lines(sapply(1:M,function(m) quantile(logit(theta.7[mm,m,]),0.9)),col=3,lty=2)
legend("topright", bty ="n",legend=c("Est. mu:mean","Est. median",'Est. 10%,90% quantile'),col=c(1,2,3), lty=c(1,1,2), cex=0.7)

plot(rowMeans(logit(theta.8[mm,,])),ylim=c(0,1),type='l',col=1,main = 'rho',xlab = 'Iteraion times',ylab = 'Iterated rho')
lines(sapply(1:M,function(m) median(logit(theta.8[mm,m,]))),col=2,lty=1)
lines(sapply(1:M,function(m) quantile(logit(theta.8[mm,m,]),0.1)),col=3,lty=2)
lines(sapply(1:M,function(m) quantile(logit(theta.8[mm,m,]),0.9)),col=3,lty=2)
legend("topright", bty ="n",legend=c("Est. rho:mean","Est. median",'Est. 10%,90% quantile'),col=c(1,2,3), lty=c(1,1,2), cex=0.7)

plot(rowMeans(exp(theta.9[mm,,])),ylim=c(0,2),type='l',col=1,main = 'psi',xlab = 'Iteraion times',ylab = 'Iterated psi')
lines(sapply(1:M,function(m) median(exp(theta.9[mm,m,]))),col=2,lty=1)
lines(sapply(1:M,function(m) quantile(exp(theta.9[mm,m,]),0.1)),col=3,lty=2)
lines(sapply(1:M,function(m) quantile(exp(theta.9[mm,m,]),0.9)),col=3,lty=2)
legend("topright", bty="n",legend=c("Est. psi:mean","Est. median",'Est. 10%,90% quantile'),col=c(1,2,3), lty=c(1,1,2), cex=0.7)
