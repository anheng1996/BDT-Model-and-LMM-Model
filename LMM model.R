t1=proc.time()

setwd("~/R/TSM")
library(readxl)

dat1<-read_xlsx("Data Input for R.xlsx",sheet=1)
dat2<-read_xlsx("Data Input for R.xlsx",sheet=2)

####Calibrate caplet vols to get a,b,c,d, and phai
list_black=dat1$`Black Vols`[2:16]^2*dat1$`Maturity, T`[2:16]

error_cap<-function(x) {
  if (x[1]>x[2]/x[3] || x[1]+x[4]<0 || x[2]<0 || x[4]<0) {
    return(999999999)
  } else {
    list_fun7<-rep(0,15)
    for (i in 1:15) {
      T=dat1$`Maturity, T`[i+1]
      sigma_squared<-function(t) {
        sigma_i=(x[1]+x[2]*(T-t))*exp(-x[3]*(T-t))+x[4]
        return(sigma_i^2)
      }
      list_fun7[i]=integrate(sigma_squared,0,T)$value
    }
    return(sum((list_black-list_fun7)^2))
  }
}

guess1<-c(-0.02,0.3,2,0.14)
output_1=optim(guess1,error_cap)

a=output_1$par[1]
b=output_1$par[2]
c=output_1$par[3]
d=output_1$par[4]

list_fun7<-rep(0,15)
for (i in 1:15) {
  T=dat1$`Maturity, T`[i+1]
  sigma_squared<-function(t) {
    sigma_i=(a+b*(T-t))*exp(-c*(T-t))+d
    return(sigma_i^2)
  }
  list_fun7[i]=integrate(sigma_squared,0,T)$value
}

phai=sqrt(list_black/list_fun7)

####Calibrate swaption to get correlation

cor_gerenation<-function(alpha1,yita1) {
  cor_matrix<-matrix(nrow = 15,ncol = 15)
  for (i in 1:15){
    for (j in 1:i){
      #cor_matrix[i,j]=exp(-abs(i-j)/14*(-log(alpha1)+yita1*(15-i-j+1)/(13)))
      cor_matrix[i,j]=exp(-abs(i-j)/14*(-log(alpha1)+yita1*(i^2+j^2+i*j-3*15*i-3*15*j+3*i+3*j+2*15^2-15-4)/(13*12)))
    }
  }

  for (i in 1:14) {
    for (j in (i+1):15){
      cor_matrix[i,j]=cor_matrix[j,i]
    }

  }
  return(cor_matrix)
}

list1<-dat2$`Black Vol`^2*dat1$`Maturity, T`[(dat2$T0/0.25+1)]#Black vols squared * T0

error_sum<-function(x) {
  if (x[1]<=0 || x[2]<0 || x[2]>-log(x[1])) {
    return(9999999)
  }else {
    list2<-rep(0,11)
    cor<-cor_gerenation(x[1],x[2])
    
    for (m in 1:11) {
      T0=dat2$T0[m]
      tenor=dat2$tenor[m]
      swap_rate=dat2$`Swap Rate at t=0`[m]
      
      for (i in 1:(tenor/0.25)) {
        for (j in 1:(tenor/0.25)) {
          fi=dat1$`Forward Rate(%)`[((T0/0.25)+i)]
          fj=dat1$`Forward Rate(%)`[((T0/0.25)+j)]
          
          denominator=sum(dat1$`t x P(0,T)`[(T0/0.25+1):(T0/0.25+tenor/0.25)])
          wi=dat1$`t x P(0,T)`[((T0/0.25)+i)]/denominator
          wj=dat1$`t x P(0,T)`[((T0/0.25)+j)]/denominator
          
          phai_i=phai[((T0/0.25)+i-1)]
          phai_j=phai[((T0/0.25)+j-1)]
          
          Ti=dat1$`Maturity, T`[((T0/0.25)+i)]
          Tj=dat1$`Maturity, T`[((T0/0.25)+j)]
          
          rho_ij=cor[(T0/0.25+i-1),(T0/0.25+j-1)]
          
          rho_sigma_integral<-function(t) {
            sigmai=phai_i*((b*(Ti-t)+a)*exp(-c*(Ti-t))+d)
            sigmaj=phai_j*((b*(Tj-t)+a)*exp(-c*(Tj-t))+d)
            return(rho_ij*sigmai*sigmaj)
          }
          
          list2[m]=list2[m]+((wi*wj*fi*fj)/swap_rate^2)*integrate(rho_sigma_integral,0,dat1$`Maturity, T`[(T0/0.25+1)])$value
        }
      }
    }
    
    return(sum((list2-list1)^2))
  }
}

guess<-c(0.2,1)
output_1=optim(guess,error_sum)

alpha11=output_1$par[1]
yita11=output_1$par[2]

alpha11
yita11
cor_matrix=cor_gerenation(alpha11,yita11)

####Simulate libors
libor_sim<-function() {
  l_sim<-array(0,dim = c(15,16))
  l_sim[,1]<-dat1$`Forward Rate(%)`[2:16]
  for (k in 1:14) {
    C_matrix<-array(0,dim = c((16-k),(16-k)))
    
    for (i in 1:(16-k)) {
      for (j in 1:(16-k)) {
        phai_i=phai[(i+k-1)]
        phai_j=phai[(j+k-1)]
        
        Ti=dat1$`Maturity, T`[(i+k)]
        Tj=dat1$`Maturity, T`[(j+k)]
        
        rho_ij=cor_matrix[(i+k-1),(j+k-1)]
        
        rho_sigma_integral<-function(t) {
          sigmai=phai_i*((b*(Ti-t)+a)*exp(-c*(Ti-t))+d)
          sigmaj=phai_j*((b*(Tj-t)+a)*exp(-c*(Tj-t))+d)
          return(rho_ij*sigmai*sigmaj)
        }
        C_matrix[i,j]=integrate(rho_sigma_integral,dat1$`Maturity, T`[k],dat1$`Maturity, T`[k+1])$value
      }
    }
    
    sigma_matrix=t(chol(C_matrix))
    random_num<-rnorm((16-k),0,1)
    l_sim[15,(k+1)]=l_sim[15,k]*exp(-0.5*C_matrix[(16-k),(16-k)]+sum(random_num*sigma_matrix[16-k,]))
    
    for (i in 14:k) {
      drift=0
      for (j in (i+1):15) {
        tar_j=dat1$`Year fraction, t`[(j+2)]
        drift=drift-0.5*((l_sim[j,k]*tar_j)/(1+l_sim[j,k]*tar_j)+
                           (l_sim[j,(k+1)]*tar_j)/(1+l_sim[j,(k+1)]*tar_j))*C_matrix[(j-k+1),(i-k+1)]
      }
      l_sim[i,(k+1)]=l_sim[i,k]*exp(drift-0.5*C_matrix[(1-k+i),(1-k+i)]+sum(random_num*sigma_matrix[(1-k+i),]))
    }
  }
  
  sigma_integral<-function(t) {
    Ti=dat1$`Maturity, T`[16]
    sigmai=phai[15]*((b*(Ti-t)+a)*exp(-c*(Ti-t))+d)
    return(sigmai^2)
  }
  C=integrate(sigma_integral,dat1$`Maturity, T`[15],dat1$`Maturity, T`[16])$value
  l_sim[15,16]<-l_sim[15,15]*exp(-0.5*C+sqrt(C)*rnorm(1,0,1))
  return(l_sim)
}

num_sim=10000
libor_matrix<-array(0,dim = c(15,16,num_sim))
for (i in 1:num_sim){
  libor_matrix[,,i]=libor_sim()
}

####Value caplets
caplet_value<-function(x) {
  caplet_value_sum=0
  cap_strike=dat1$`Cap Strike`[(x+1)]
  
  if(x != 15) {
    for (i in 1:num_sim) {
      pam=max(0,libor_matrix[x,(x+1),i]-cap_strike)*dat1$`Year fraction, t`[(x+2)]        #payment at the caplet maturity
      vam=pam*prod(1+libor_matrix[(x+1):15,(x+1),i]*dat1$`Year fraction, t`[(x+3):17])    #value at the maturity
      va0=vam*dat1$`P(0,T)`[17]                            #value at t=0
      caplet_value_sum=caplet_value_sum+va0
    }
  } else {
    for (i in 1:num_sim) {
      caplet_value_sum=caplet_value_sum+max(0,libor_matrix[x,(x+1),i]-cap_strike)*dat1$`Year fraction, t`[17]
    }
  }

  return(caplet_value_sum/num_sim)
}

capletvalues<-array(dim = c(15,1))
for (i in 1:15) {
  capletvalues[i,1]=caplet_value(i)
}

####Value European Swaptions

euro_swaption<-function(T0,tenor,K) {
  euswap_valuesum=0
  
  for (i in 1:num_sim) {
    fix_leg=dat1$`Year fraction, t`[(T0/0.25+2):(T0/0.25+1+tenor/0.25)]*K
    df_fix<-rep(0,tenor/0.25)  #discount factors for the fixed legs
    df_fix[1]=1/(1+libor_matrix[(T0/0.25),(T0/0.25+1),i]*dat1$`Year fraction, t`[(T0/0.25+2)])
    
    for (j in 2:(tenor/0.25)) {
      df_fix[j]=df_fix[(j-1)]/(1+libor_matrix[((T0/0.25)+j-1),(T0/0.25+1),i]*dat1$`Year fraction, t`[(T0/0.25+1+j)])
    }
    fix_leg_sum=sum(fix_leg*df_fix)
    
    float_leg=1-df_fix[(tenor/0.25)]
    
    value_e=max(0,float_leg-fix_leg_sum)  #swaption value at exercise date
    value_m=value_e*prod(1+libor_matrix[(T0/0.25):15,(T0/0.25+1),i]*dat1$`Year fraction, t`[(T0/0.25+2):17])  #swaption value at maturity
    value_0=value_m*dat1$`P(0,T)`[17] #swaption value at t=0
    euswap_valuesum=euswap_valuesum+value_0
  }
  return(euswap_valuesum/num_sim)
}

euswaption<-array(dim = c(11,1))
for (i in 1:11) {
  euswaption[i,1]=euro_swaption(T0=dat2$T0[i],tenor=dat2$tenor[i],K=dat2$`Swap Rate at t=0`[i])
}


####Value the Product
dat3<-read_xlsx("Data Input for R.xlsx",sheet=3)

m=14 #number of exercise date
cfe=1000+dat3$`accrued interest if redempted`  #cash flow if redemption
coupon=dat3$`accrued interest if redempted`
coupon[c(1,3,5,7,9,11,13)]=0

scale1<-function(x) {
  xmin=min(x)
  xmax=max(x)
  a=2/(xmax-xmin)
  b=1-a*xmax
  return(a*x+b)
}

product<-function( ) {
  value<-matrix(nrow = num_sim,ncol = 14)
  value[,14]=cfe[14]
  for (i in 1:13) {
    Y=as.matrix(value[,(15-i)]/(1+libor_matrix[(14-i),(15-i),]*dat1$`Year fraction, t`[(16-i)])+coupon[(15-i)])
    
    p=array(dim = c((i+2),num_sim))
    
    p[1,]=1/(1+libor_matrix[(14-i),(15-i),]*dat1$`Year fraction, t`[(16-i)])
    for (j in 2:(i+2)) {
      p[j,]=p[(j-1),]/(1+libor_matrix[(13-i+j),(15-i),]*dat1$`Year fraction, t`[(15-i+j)])
    }
    
    x1=libor_matrix[(14-i),(15-i),]
    x2<-rep(0,num_sim)
    for (j in 1:num_sim) {
      x2[j]=(1-p[(i+2),j])/sum(p[1:(i+2),j]*dat1$`Year fraction, t`[(16-i):17])
    }
    X=matrix(nrow = num_sim,ncol = 6)
    X[,1]=1
    X[,2]=scale1(x1)
    X[,3]=scale1(x2)
    X[,4]=X[,2]^2
    X[,5]=X[,3]^2
    X[,6]=X[,2]*X[,3]
    
    theta1=solve(t(X)%*%X+0.00001*diag(6))%*%t(X)%*%Y
    ECV=X%*%theta1   #expected continuation value
    
    for(j in 1:num_sim){
      value[j,(14-i)]=ifelse(ECV[j]<Y[j,1],cfe[(14-i)],Y[j,1])
    }
  }
  
  value_m=rep(0,num_sim)   #value at maturity
  for (i in 1:num_sim) {
    value_m[i]=value[i,1]*prod(1+libor_matrix[1:15,2,i]*dat1$`Year fraction, t`[3:17])
  }
  value_0=(sum(value_m)*dat1$`P(0,T)`[17])/num_sim
  return(value_0)
}

pro_value=product()

####Value Bermudan Swaption
Bermudan<-function() {
  euro_swaption1<-function(T0,tenor,K) {
    euswap_value<-rep(0,num_sim)
    
    for (i in 1:num_sim) {
      fix_leg=dat1$`Year fraction, t`[(T0/0.25+2):(T0/0.25+1+tenor/0.25)]*K
      df_fix<-rep(0,tenor/0.25)  #discount factors for the fixed legs
      df_fix[1]=1/(1+libor_matrix[(T0/0.25),(T0/0.25+1),i]*dat1$`Year fraction, t`[(T0/0.25+2)])
      
      for (j in 2:(tenor/0.25)) {
        df_fix[j]=df_fix[(j-1)]/(1+libor_matrix[((T0/0.25)+j-1),(T0/0.25+1),i]*dat1$`Year fraction, t`[(T0/0.25+1+j)])
      }
      fix_leg_sum=sum(fix_leg*df_fix)
      
      float_leg=1-df_fix[(tenor/0.25)]
      
      value_e=max(0,float_leg-fix_leg_sum)  #swaption value at exercise date
      value_m=value_e*prod(1+libor_matrix[(T0/0.25):15,(T0/0.25+1),i]*dat1$`Year fraction, t`[(T0/0.25+2):17])  #swaption value at maturity
      value_0=value_m*dat1$`P(0,T)`[17] #swaption value at t=0
      euswap_value[i]=value_0
    }
    return(euswap_value)
  }
  exer2=euro_swaption1(2,2,0.015)
  exer1=euro_swaption1(1,2,0.015)
  
  Y<-matrix(nrow = num_sim,ncol = 1)
  for (i in 1:num_sim) {
    Y[i,1]=exer2[i]/((1+dat1$`Year fraction, t`[6]*libor_matrix[4,5,i])*(1+dat1$`Year fraction, t`[7]*libor_matrix[5,6,i])*
                       (1+dat1$`Year fraction, t`[8]*libor_matrix[6,7,i])*(1+dat1$`Year fraction, t`[9]*libor_matrix[7,8,i]))
  }
  
  X<-matrix(nrow = num_sim,ncol = 6)
  X[,1]=1
  X[,2]=scale1(libor_matrix[4,5,])
  
  p<-array(dim = c(12,num_sim))
  p[1,]=1/(1+libor_matrix[4,5,]*dat1$`Year fraction, t`[6])
  for (i in 2:12) {
    p[i,]=p[(i-1),]/(1+libor_matrix[(3+i),5,]*dat1$`Year fraction, t`[(5+i)])
  }
  
  for (i in 1:num_sim) {
    X[i,3]=(1-p[12,i])/sum(p[1:12,i]*dat1$`Year fraction, t`[6:17])
  }
  
  X[,3]=scale1(X[,3])
  X[,4]=X[,2]^2
  X[,4]=scale1(X[,4])
  X[,5]=X[,3]^2
  X[,5]=scale1(X[,5])
  X[,6]=X[,2]*X[,3]
  X[,6]=scale1(X[,6])
  
  theta2=solve(t(X)%*%X+0.00001*diag(6))%*%t(X)%*%Y
  ECV=X%*%theta2   #expected continuation value
  
  valuea_4<-rep(0,num_sim)  #value at t=4
  for (i in 1:num_sim) {
    valuea_4[i]=ifelse(exer1[i]<ECV[i,1],Y[i,1],exer1[i])*prod(1+dat1$`Year fraction, t`[6:17]*libor_matrix[4:15,5,i])
  }
  return(mean(valuea_4)*dat1$`P(0,T)`[17])
}

ber_swaption=Bermudan( )

t2=proc.time()
t=t2-t1
print(paste0('Runtime£º',t[3][[1]],'s'))