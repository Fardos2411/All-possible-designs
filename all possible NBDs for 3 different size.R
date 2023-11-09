################################################################################
# CBND_3diffsize: Circular balance neighbor design for block of three different 
# sizes (K1,k2 and k3)

# Algorithm from paper:
# Akbar Firdos,Mahmood Ul Hassan,Farrukh Jamal,Hurria Ali,Khadija Noreen and Rashid Ahmed Algorithms to Construct Minimal Circular Strongly. 
# Coded by Fardos Akbar et al., 2023
################################################################################
################################################################################
# Selection of i group of size K1 from adjusted A. The set of remaining 
# (Unselected) elements are saved in object B2. 
################################################################################
grouping1<-function(A,k,v,i){
  bs<-c()
  z=0;f=1
  A1=A
  while(f<=i){
    for(y in 1:5000){
      comp<-sample(1:length(A1),k)
      com<-A1[comp]
      cs<-sum(com)
      if(cs%%v==0){
        bs<-rbind(bs,com)
        A1<-A1[-comp]
        z<-z+1
        f=f+1
      }
      if(z==i) break
    }
    if(z<i) {bs<-c();z=0;f=1;A1=A}  
  }
  list(B1=bs,B2=A1)
}
################################################################################
# Selection of i group of size K1 from adjusted A and selection of required 
# number of groups of size K2 from B2. The set of remaining (Unselected) 
# elements are saved in B3.
################################################################################
grouping2<-function(A,k,v,i,sk2){
  bs1<-c()
  j=i+sk2
  z=0;f=1
  A1=A
  while(f<=j){
    s<-grouping1(A1,k[1],v,i)
    A2<-s$B2
    z=i;f=f+i
    for(y in 1:2000){
      comp<-sample(1:length(A2),k[2])
      com<-A2[comp]
      cs<-sum(com)
      if(cs%%v==0){
        bs1<-rbind(bs1,com)
        A2<-A2[-comp]
        z<-z+1
        f=f+1
      }
      if(z==j) break
    }
    if(z<j) {bs1<-c();z=0;f=1;A1=A}
  }
  list(B1=s$B1,B2=bs1,B3=A2)
}
################################################################################
# Selection of i group of size K1 from adjusted A, selection of required number 
# of groups of size K2 from B2 and division of required number of groups of size
# K3 from B3.
################################################################################
grouping3<-function(A,k,v,i,sk2,sk3){
  bs1<-c()
  j=i+sk2+sk3
  z=0;f=1
  A1=A
  while(f<=j){
    s<-grouping2(A1,k,v,i,sk2)
    A3<-s$B3
    z=i+sk2;f=f+i+sk2
    for(y in 1:1000){
      comp<-sample(1:length(A3),k[3])
      com<-A3[comp]
      cs<-sum(com)
      if(cs%%v==0){
        bs1<-rbind(bs1,com)
        A3<-A3[-com]
        z<-z+1
        f=f+1
      }
      if(z==j) break
    }
    if(z<j) {bs1<-c();z=0;f=1;A1=A}  
  }
  gs1<-t(apply(s$B1,1,sort))
  gs1<-cbind(gs1,rowSums(gs1),rowSums(gs1)/v)
  rownames(gs1)<-paste("G",1:i, sep="")
  colnames(gs1)<-c(paste(1:k[1], sep=""),"sum" ,"sum/v")
  gs2<-t(apply(s$B2,1,sort))
  gs2<-cbind(gs2,rowSums(gs2),rowSums(gs2)/v)
  rownames(gs2)<-paste("G",(i+1):(i+sk2), sep="")
  colnames(gs2)<-c(paste(1:k[2], sep=""),"sum" ,"sum/v")
  gs3<-t(apply(bs1,1,sort))
  gs3<-cbind(gs3,rowSums(gs3),rowSums(gs3)/v)
  rownames(gs3)<-paste("G",(i+sk2+1):(i+sk2+sk3), sep="")
  colnames(gs3)<-c(paste(1:k[3], sep=""),"sum" ,"sum/v")
  fs1<-t(apply(s$B1,1,sort))
  fs1<-delmin(fs1)
  rownames(fs1)<-paste("S",1:i, sep="")
  colnames(fs1)<-rep("",(k[1])-1)
  fs2<-t(apply(s$B2,1,sort))
  fs2<-delmin(fs2)
  rownames(fs2)<-paste("S",(i+1):(i+sk2), sep="")
  colnames(fs2)<-rep("",(k[2])-1)
  fs3<-t(apply(bs1,1,sort))
  fs3<-delmin(fs3)
  rownames(fs3)<-paste("S",(i+sk2+1):(i+sk2+sk3), sep="")
  colnames(fs3)<-rep("",(k[3]-1))
  list(B1=list(fs1,fs2,fs3),B4=list(gs1,gs2,gs3),B5=A3)
}
#######################################################################
# Obtain set(s) of shifts by deleting smallest value of each group
#######################################################################
delmin<-function(z){
  fs<-c()
  n<-nrow(z)
  c<-ncol(z)-1
  for(i in 1:n){
    z1<-z[i,]
    z2<-z1[z1!=min(z1)]
    fs<-rbind(fs,z2)
  }
  return(fs)
}
################################################################################
# Selection of adjusted A and the set(s) of shifts to obtain Circular  
# balance neighbor design for three different block size.
################################################################################
# D=2: Circular Strongly Balanced Neighbor Designs
# D=1: Circular Balanced Neighbor Designs
#   K: Vector of three different block sizes 
#   i: Number of sets of shifts for K1
# Sk2: Number of sets of shifts for K2
# Sk3: Number of sets of shifts for K3
CBND_3diffsize<-function(k,i,D,sk2,sk3){
  if(length(k)>3 | length(k)<3){stop("length(k)=3")}
  if(any(k<=2)!=0) stop("k=Block size: Each block size must be greater than 3")
  if(i<=0) stop("i=must be a positive integer")
  if(k[1]<k[2] | k[2]<k[3] |  k[1]<k[3]  ) stop("k1>K2>K3")
  
  setClass( "stat_test", representation("list"))
  
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 52), collapse = "")
    cat(row, "\n")
    
    if(D==2){
      cat("Following are required sets of shifts to obtain the 
          minimal CSBND for", "v=" ,object$R[1], ",","k1=",object$R[2],
          ",","k2=",object$R[3],"and","k3=",object$R[4],"\n")
      row <- paste(rep("=", 52), collapse = "")}
    if(D==1){
      cat("Following are required sets of shifts to obtain the 
          minimal CBND for", "v=" ,object$R[1], ",","k1=",object$R[2],
          ",","k2=",object$R[3],"and","k3=",object$R[4],"\n")
      row <- paste(rep("=", 52), collapse = "")}
    
    cat(row, "\n")
    print(object$S[[1]])
    cat("\n")
    print(object$S[[2]])
    cat("\n")
    print(object$S[[3]])
  })
  
  
  if(D==2 & sk2==1 & sk3==1){  
    v=2*i*k[1]+2*k[2]+2*k[3]-1 ; m=(v-1)/2
    if(m%%8==0){
      j=m/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(j-1),(j+1):m,(v-j))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==1){
      j=(m-1)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==2){
      j=(m-2)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==3){
      j=(m-3)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(m-j-1),(m-j+1):m,(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    
    if(m%%8==4){
      j=(m-4)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:j,(j+2):(m-1),(m+1),(v-(j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==5){
      j=(m-5)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==6){
      j=(m-6)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==7){
      j=(m-7)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    } 
    
    
  }  
  
  if(D==1 & sk2==1 & sk3==1){
    v= 2*i*k[1]+2*k[2]+2*k[3]+1; m=(v-1)/2
    
    if(m%%8==0){
      j=m/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(j-1),(j+1):m,(v-j))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==1){
      j=(m-1)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==2){
      j=(m-2)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==3){
      j=(m-3)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(m-j-1),(m-j+1):m,(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    
    if(m%%8==4){
      j=(m-4)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:j,(j+2):(m-1),(m+1),(v-(j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==5){
      j=(m-5)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==6){
      j=(m-6)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==7){
      j=(m-7)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }   
    
    
  }   
  
  
  if(D==2 & sk2==2 & sk3==2){  
    v=2*i*k[1]+4*k[2]+4*k[3]-1 ; m=(v-1)/2
    
    if(m%%8==0){
      j=m/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(j-1),(j+1):m,(v-j))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==1){
      j=(m-1)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==2){
      j=(m-2)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==3){
      j=(m-3)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(m-j-1),(m-j+1):m,(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    
    if(m%%8==4){
      j=(m-4)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:j,(j+2):(m-1),(m+1),(v-(j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==5){
      j=(m-5)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==6){
      j=(m-6)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==7){
      j=(m-7)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    } 
  }  
  
  
  if(D==1 & sk2==2 & sk3==2){
    v= 2*i*k[1]+4*k[2]+4*k[3]+1  ;m=(v-1)/2
    
    if(m%%8==0){
      j=m/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(j-1),(j+1):m,(v-j))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==1){
      j=(m-1)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==2){
      j=(m-2)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==3){
      j=(m-3)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(m-j-1),(m-j+1):m,(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    
    if(m%%8==4){
      j=(m-4)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:j,(j+2):(m-1),(m+1),(v-(j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==5){
      j=(m-5)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==6){
      j=(m-6)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==7){
      j=(m-7)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }   
    
  } 
  
  
  if(D==2 & sk2==1 & sk3==2){  
    v=2*i*k[1]+2*k[2]+4*k[3]-1 ; m=(v-1)/2
    
    if(m%%8==0){
      j=m/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(j-1),(j+1):m,(v-j))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==1){
      j=(m-1)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==2){
      j=(m-2)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==3){
      j=(m-3)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(m-j-1),(m-j+1):m,(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    
    if(m%%8==4){
      j=(m-4)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:j,(j+2):(m-1),(m+1),(v-(j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==5){
      j=(m-5)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==6){
      j=(m-6)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==7){
      j=(m-7)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    } 
  }  
  
  
  if(D==1 & sk2==1 & sk3==2){
    v= 2*i*k[1]+2*k[2]+4*k[3]+1  ;m=(v-1)/2
    
    
    if(m%%8==0){
      j=m/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(j-1),(j+1):m,(v-j))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==1){
      j=(m-1)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==2){
      j=(m-2)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==3){
      j=(m-3)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(m-j-1),(m-j+1):m,(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    
    if(m%%8==4){
      j=(m-4)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:j,(j+2):(m-1),(m+1),(v-(j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==5){
      j=(m-5)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==6){
      j=(m-6)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==7){
      j=(m-7)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    } 
  } 
  
  if(D==2 & sk2==2 & sk3==1){  
    v=2*i*k[1]+4*k[2]+2*k[3]-1; m=(v-1)/2
    if(m%%8==0){
      j=m/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(j-1),(j+1):m,(v-j))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==1){
      j=(m-1)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==2){
      j=(m-2)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==3){
      j=(m-3)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(m-j-1),(m-j+1):m,(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    
    if(m%%8==4){
      j=(m-4)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:j,(j+2):(m-1),(m+1),(v-(j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==5){
      j=(m-5)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==6){
      j=(m-6)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==7){
      j=(m-7)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    } 
  }  
  
  
  if(D==1 & sk2==2 & sk3==1){
    v= 2*i*k[1]+4*k[2]+2*k[3]+1  ;m=(v-1)/2
    
    
    if(m%%8==0){
      j=m/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(j-1),(j+1):m,(v-j))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==1){
      j=(m-1)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==2){
      j=(m-2)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==3){
      j=(m-3)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(m-j-1),(m-j+1):m,(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    
    if(m%%8==4){
      j=(m-4)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:j,(j+2):(m-1),(m+1),(v-(j+1)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==5){
      j=(m-5)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==6){
      j=(m-6)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    
    if(m%%8==7){
      j=(m-7)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    } 
    
  } 
  new("stat_test", x)  
}
###################################################################
# Generation of design using sets of cyclical shifts
###################################################################
# H is an output object from CBND_3diffsize
# The output is called using the design_CBND to generate design
design_CBND<-function(H){
  setClass( "CBND_design", representation("list"))
  setMethod("show", "CBND_design", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following is minimal CBND for", "v=" ,object$R[1], "and","k=",object$R[2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    for(i in 1:length(ss)){
      W<-ss[[i]]
      nr<-dim(W)[1]
      for(j in 1:nr){
        print(object$Design[[i]][[j]])
        cat("\n\n")
      }}
  })  
  v<-H$R[1]
  k<-H$R[2]
  ss<-H$S  
  treat<-(1:v)-1
  fn<-(1:v)
  G<-list()
  for(j in 1:length(ss)){ 
    W<-ss[[j]]
    nr<-dim(W)[1]
    nc<-dim(W)[2]
    D<-list()
    for(i in 1:nr){
      dd<-c()
      d1<-matrix(treat,(nc+1),v,byrow = T)
      ss1<-cumsum(c(0,W[i,]))
      dd2<-d1+ss1
      dd<-rbind(dd,dd2)
      rr<-dd[which(dd>=v)]%%v
      dd[which(dd>=v)]<-rr
      colnames(dd)<-paste("B",fn, sep="")
      rownames(dd)<-rep("",(nc+1))
      fn<-fn+v
      D[[i]]<-dd
    }
    G[[j]]<-D
  }
  x<-list(Design=G,R=H$R)
  new("CBND_design",x)
}
#################################################################################
# Selection of adjusted A and the set(s) of shifts to obtain Circular Generalized 
# neighbor design for three different block size.
##################################################################################
# D=3: Minimal CPNDs in which v/2 unordered pairs do not appear
# D=4: Minimal CPNDs in which 3v/2 unordered pairs do not appear 
#   K: Vector of three different block sizes
#   i: Number of sets of shifts for K1
# Sk2: Number of sets of shifts for K2
# Sk3: Number of sets of shifts for K3
CPND_3diffsize<-function(k,i,D,sk2,sk3){
  if(length(k)>3 | length(k)<3){stop("length(k)=3")}
  if(any(k<=2)!=0) stop("k=Block size: Each block size must be greater than 3")
  if(i<=0) stop("i=must be a positive integer")
  if(k[1]<k[2] | k[2]<k[3] |  k[1]<k[3]  ) stop("k1>K2>K3")
  setClass( "stat_test", representation("list"))
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 52), collapse = "")
    cat(row, "\n")
    cat("Following are required sets of shifts to obtain the 
        minimal CPND for", "v=" ,object[[5]][1], ",","k1=",object[[5]][2],
        ",","k2=",object[[5]][3],"and","k3=",object[[5]][4],"\n")
    row <- paste(rep("=", 52), collapse = "")
    cat(row, "\n")
    print(object[[1]])
    cat("\n")
    print(object[[2]])
    cat("\n")
    print(object[[3]])
  })
  if(D==3 & sk2==1 & sk3==1){  
    v=2*i*k[1]+2*k[2]+2*k[3]+2 ; m=(v-2)/2
    if(m%%4==0){
      A<-1:m
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1 | m%%4==2){return("The minimal CPNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+2k3+2 and k1,k2 k3 block sizes")}
  } 
  if(D==4 & sk2==1 & sk3==1){
    v= 2*i*k[1]+2*k[2]+2*k[3]+4 
    m<-(v-2)/2
    if(m%%4==0){
      A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(2:(m-2),m,(2*m+1))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
  }
  if(D==3 & sk2==2 & sk3==2){  
    v=2*i*k[1]+4*k[2]+4*k[3]+2 ; m=(v-2)/2
    if(m%%4==0){
      A<-1:m
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1 | m%%4==2){return("The minimal CPNDs in which v/2 unordered pair cannot be constructed for V=2ik1+4k2+4k3+2 and k1,k2 k3 block sizes")}
  }  
  if(D==4 & sk2==2 & sk3==2){
    v= 2*i*k[1]+4*k[2]+4*k[3]+4  ;m<-(v-2)/2
    if(m%%4==0){
      A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(2:(m-2),m,(2*m+1))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
  } 
  if(D==3 & sk2==1 & sk3==2){  
    v=2*i*k[1]+2*k[2]+4*k[3]+2 ; m=(v-2)/2
    if(m%%4==0){
      A<-1:m
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1 | m%%4==2){return("The minimal CPNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+4k3+2 and k1,k2 k3 block sizes")}
  }  
  if(D==4 & sk2==1 & sk3==2){
    v= 2*i*k[1]+2*k[2]+4*k[3]+4  ;m<-(v-2)/2
    if(m%%4==0){
      A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(2:(m-2),m,(2*m+1))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
  } 
  if(D==3 & sk2==2 & sk3==1){  
    v=2*i*k[1]+4*k[2]+2*k[3]+2 ; m=(v-2)/2
    if(m%%4==0){
      A<-1:m
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1 | m%%4==2){return("The minimal CPNDs in which v/2 unordered pair cannot be constructed for V=2ik1+4k2+2k3+2 and k1,k2 k3 block sizes")}
  }  
  if(D==4 & sk2==2 & sk3==1){
    v= 2*i*k[1]+4*k[2]+2*k[3]+4  ;m<-(v-2)/2
    if(m%%4==0){
      A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(2:(m-2),m,(2*m+1))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S1=A1$B1,S2=A1$B2,S3=A1$B3,G=A1$B4,R=A2,A=A)
    }
  } 
  new("stat_test", x)  
}
#################################################################################
# Selection of adjusted A and the set(s) of shifts to obtain Circular Weakly 
# balance neighbour design for three different block size.
##################################################################################
# D=5: Minimal CWBNDs in which v/2 of the ordered pairs appear twice as neighbors
# D=6: Minimal CWBNDs in which 3v/2 of the unordered pairs appear twice as neighbors 
#   K: Vector of three different block sizes 
#   i: Number of sets of shifts for K1
# Sk2: Number of sets of shifts for K2
# Sk3: Number of sets of shifts for K3
CWBND_3diffsize<-function(k,i,D,sk2,sk3){
  if(length(k)>3 | length(k)<3){stop("length(k)=3")}
  if(any(k<=2)!=0) stop("k=Block size: Each block size must be greater than 3")
  if(i<=0) stop("i=must be a positive integer")
  if(k[1]<k[2] | k[2]<k[3] |  k[1]<k[3]  ) stop("k1>K2>K3")
  setClass( "stat_test", representation("list"))
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 52), collapse = "")
    cat(row, "\n")
    cat("Following are required sets of shifts to obtain the 
        minimal CWBND for", "v=" ,object$R[1], ",","k1=",object$R[2],
        ",","k2=",object$R[3],"and","k3=",object$R[4],"\n")
    row <- paste(rep("=", 52), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
    cat("\n")
    print(object$S[[2]])
    cat("\n")
    print(object$S[[3]])
  })
  if(D==5 & sk2==1 & sk3==1){  
    v=2*i*k[1]+2*k[2]+2*k[3] ; m<-(v-2)/2
    if(m%%4==2){
      A=1:(m+1)
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(1:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==0 | m%%4==1){return("The minimal CWBNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+2k2+2k3 and and k1,k2 k3 block sizes")}
  }  
  if(D==6 & sk2==1 & sk3==1){
    v= 2*i*k[1]+2*k[2]+2*k[3]-2; m<-(v-2)/2
    if(m%%4==0){
      A=c(1:(m-1),(m-1),(m+1),(m+2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(1:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(1:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(1:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
  }   
  if(D==5 & sk2==2 & sk3==2){  
    v=2*i*k[1]+4*k[2]+4*k[3] ; m=(v-2)/2
    if(m%%4==2){
      A=1:(m+1)
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(1:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==0 | m%%4==1){return("The minimal CWBNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+4k2+4k3 and and k1,k2 k3 block sizes")}
  }  
  if(D==6 & sk2==2 & sk3==2){
    v= 2*i*k[1]+4*k[2]+4*k[3]-2  ;m<-(v-2)/2
    if(m%%4==0){
      A=c(1:(m-1),(m-1),(m+1),(m+2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(1:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(1:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(1:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
  } 
  if(D==5 & sk2==1 & sk3==2){  
    v=2*i*k[1]+2*k[2]+4*k[3] ; m=(v-2)/2
    if(m%%4==2){
      A=1:(m+1)
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(1:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==0 | m%%4==1){return("The minimal CWBNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+2k2+4k3 and and k1,k2 k3 block sizes")}
  }  
  if(D==6 & sk2==1 & sk3==2){
    v= 2*i*k[1]+2*k[2]+4*k[3]-2  ;m<-(v-2)/2
    if(m%%4==0){
      A=c(1:(m-1),(m-1),(m+1),(m+2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(1:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(1:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(1:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
  } 
  if(D==5 & sk2==2 & sk3==1){  
    v=2*i*k[1]+4*k[2]+2*k[3] ; m=(v-2)/2
    if(m%%4==2){
      A=1:(m+1)
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(1:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==0 | m%%4==1){return("The minimal CWBNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+4k2+2k3 and and k1,k2 k3 block sizes")}
  }  
  if(D==6 & sk2==2 & sk3==1){
    v= 2*i*k[1]+4*k[2]+2*k[3]-2  ;m<-(v-2)/2
    if(m%%4==0){
      A=c(1:(m-1),(m-1),(m+1),(m+2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(1:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(1:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(1:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
  } 
  new("stat_test", x)  
}
###################################################################
# Generation of design using sets of cyclical shifts
###################################################################
# H is an output object from CWBND_equalsize
# The output is called using the design_CWBND to generate design
design_CWBND<-function(H){
  setClass( "CWBND_design", representation("list"))
  setMethod("show", "CWBND_design", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following is minimal CWBND for", "v=" ,object$R[1], "and","k=",object$R[2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    for(i in 1:length(ss)){
      W<-ss[[i]]
      nr<-dim(W)[1]
      for(j in 1:nr){
        print(object$Design[[i]][[j]])
        cat("\n\n")
      }}
  })  
  v<-H$R[1]
  k<-H$R[2]
  ss<-H$S  
  treat<-(1:v)-1
  fn<-(1:v)
  G<-list()
  for(j in 1:length(ss)){ 
    W<-ss[[j]]
    nr<-dim(W)[1]
    nc<-dim(W)[2]
    D<-list()
    for(i in 1:nr){
      dd<-c()
      d1<-matrix(treat,(nc+1),v,byrow = T)
      ss1<-cumsum(c(0,W[i,]))
      dd2<-d1+ss1
      dd<-rbind(dd,dd2)
      rr<-dd[which(dd>=v)]%%v
      dd[which(dd>=v)]<-rr
      colnames(dd)<-paste("B",fn, sep="")
      rownames(dd)<-rep("",(nc+1))
      fn<-fn+v
      D[[i]]<-dd
    }
    G[[j]]<-D
  }
  x<-list(Design=G,R=H$R)
  new("CWBND_design", x)
}
#################################################################################
# Selection of adjusted A and the set(s) of shifts to obtain Circular Strongly 
# Partially Balance Neighbor Design for three different block size.

#################################################################################
# D=7: Minimal CSPBNDs in which v/2 un-ordered pairs do not appear
# D=8: Minimal CSPBNDs in which 3v/2 un-ordered pairs do not appear 
# K: Vector of three different block sizes
# i: Number of sets of shifts for K1
# Sk2: Number of sets of shifts for K2
# Sk3: Number of sets of shifts for K3
CSPBND_3diffsize<-function(k,i,D,sk2,sk3){
  if(length(k)>3 | length(k)<3){stop("length(k)=3")}
  if(any(k<=2)!=0) stop("k=Block size: Each block size must be greater than 3")
  if(i<=0) stop("i=must be a positive integer")
  if(k[1]<k[2] | k[2]<k[3] | k[1]<k[3] | k[3]<3 ) stop("k1>K2>K3>2")
  setClass( "stat_test", representation("list"))
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 52), collapse = "")
    cat(row, "\n")
    cat("Following are required sets of shifts to obtain the 
        minimal CSPBNDs for", "v=" ,object$R[1], ",","k1=",object$R[2],
        ",","k2=",object$R[3],"and","k3=",object$R[4],"\n")
    row <- paste(rep("=", 52), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
    cat("\n")
    print(object$S[[2]])
    cat("\n")
    print(object$S[[3]])
  })
  if(D==7 & sk2==1 & sk3==1){ 
    v=2*i*k[1]+2*k[2]+2*k[3] ; m=(v-2)/2
    if(m%%4==0){
      A<-0:m
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(0:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1 | m%%4==2){return("The minimal CSPBNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+2k3+2 and k1,k2 k3 block sizes")}
  } 
  if(D==8 & sk2==1 & sk3==1){
    v= 2*i*k[1]+2*k[2]+2*k[3]+2
    m<-(v-2)/2
    if(m%%4==0){
      A=c(0:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(0:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(0,2:(m-2),m,(2*m+1))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(0:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
  } 
  if(D==7 & sk2==2 & sk3==2){ 
    v=2*i*k[1]+4*k[2]+4*k[3] ; m=(v-2)/2
    if(m%%4==0){
      A<-0:m
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(0:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1 | m%%4==2){return("The minimal CSPBNDs in which v/2 unordered pair cannot be constructed for V=2ik1+4k2+4k3+2 and k1,k2 k3 block sizes")}
  }
  if(D==8 & sk2==2 & sk3==2){
    v= 2*i*k[1]+4*k[2]+4*k[3]+2 ;m<-(v-2)/2
    if(m%%4==0){
      A=c(0:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(0:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(0,2:(m-2),m,(2*m+1))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(0:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
  } 
  if(D==7 & sk2==1 & sk3==2){ 
    v=2*i*k[1]+2*k[2]+4*k[3] ; m=(v-2)/2
    if(m%%4==0){
      A<-0:m
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(0:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1 | m%%4==2){return("The minimal CSPBNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+4k3+2 and k1,k2 k3 block sizes")}
  } 
  if(D==8 & sk2==1 & sk3==2){
    v= 2*i*k[1]+2*k[2]+4*k[3]+2 ;m<-(v-2)/2
    if(m%%4==0){
      A=c(0:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(0:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(0,2:(m-2),m,(2*m+1))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(0:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
  } 
  if(D==7 & sk2==2 & sk3==1){ 
    v=2*i*k[1]+4*k[2]+2*k[3] ; m=(v-2)/2
    if(m%%4==0){
      A<-0:m
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(0:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1 | m%%4==2){return("The minimal CSPBNDs in which v/2 unordered pair cannot be constructed for V=2ik1+4k2+2k3+2 and k1,k2 k3 block sizes")}
  } 
  if(D==8 & sk2==2 & sk3==1){
    v= 2*i*k[1]+4*k[2]+2*k[3]+2 ;m<-(v-2)/2
    if(m%%4==0){
      A=c(0:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(0:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(0,2:(m-2),m,(2*m+1))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(0:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
  } 
  new("stat_test", x) 
}
##################################################################
# Generation of design using sets of cyclical shifts
###################################################################
# H is an output object from CSPBND_3diffsize
# The output is called using the design_CSPBND to generate design
design_CSPBND<-function(H){
  setClass( "CSPBND_design", representation("list"))
  setMethod("show", "CSPBND_design", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following is minimal CSPBNDs for", "v=" ,object$R[1], "and","k=",object$R[2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    for(i in 1:length(ss)){
      W<-ss[[i]]
      nr<-dim(W)[1]
      for(j in 1:nr){
        print(object$Design[[i]][[j]])
        cat("\n\n")
      }}
  }) 
  v<-H$R[1]
  k<-H$R[2]
  ss<-H$S 
  treat<-(1:v)-1
  fn<-(1:v)
  G<-list()
  for(j in 1:length(ss)){ 
    W<-ss[[j]]
    nr<-dim(W)[1]
    nc<-dim(W)[2]
    D<-list()
    for(i in 1:nr){
      dd<-c()
      d1<-matrix(treat,(nc+1),v,byrow = T)
      ss1<-cumsum(c(0,W[i,]))
      dd2<-d1+ss1
      dd<-rbind(dd,dd2)
      rr<-dd[which(dd>=v)]%%v
      dd[which(dd>=v)]<-rr
      colnames(dd)<-paste("B",fn, sep="")
      rownames(dd)<-rep("",(nc+1))
      fn<-fn+v
      D[[i]]<-dd
    }
    G[[j]]<-D
  }
  x<-list(Design=G,R=H$R)
  new("CSPBND_design", x)
}
#################################################################################
# Selection of adjusted A and the set(s) of shifts to obtain Circular Strongly 
# Balance Generalized Neighbor Design for three different block size.
##################################################################################
# D=9: Minimal CSBGNDs in which v/2 of the ordered pairs appear twice as neighbors
# D=10: Minimal CSBGNDs in which 3v/2 of the unordered pairs appear twice as neighbors 
#   K: Vector of three different block sizes 
#   i: Number of sets of shifts for K1
# Sk2: Number of sets of shifts for K2
# Sk3: Number of sets of shifts for K3
CSBGND_3diffsize<-function(k,i,D,sk2,sk3){
  if(length(k)>3 | length(k)<3){stop("length(k)=3")}
  if(any(k<=2)!=0) stop("k=Block size: Each block size must be greater than 3")
  if(i<=0) stop("i=must be a positive integer")
  if(k[1]<k[2] | k[2]<k[3] |  k[1]<k[3]  ) stop("k1>K2>K3")
  setClass( "stat_test", representation("list"))
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 52), collapse = "")
    cat(row, "\n")
    cat("Following are required sets of shifts to obtain the 
        minimal CSGBND for", "v=" ,object$R[1], ",","k1=",object$R[2],
        ",","k2=",object$R[3],"and","k3=",object$R[4],"\n")
    row <- paste(rep("=", 52), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
    cat("\n")
    print(object$S[[2]])
    cat("\n")
    print(object$S[[3]])
  })
  if(D==9 & sk2==1 & sk3==1){  
    v=2*i*k[1]+2*k[2]+2*k[3]-2 ; m<-(v-2)/2
    if(m%%4==2){
      A=0:(m+1)
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+2k2+2k3 and and k1,k2 k3 block sizes")}
  }
  if(D==10 & sk2==1 & sk3==1){
    v= 2*i*k[1]+2*k[2]+2*k[3]-4; m<-(v-2)/2
    if(m%%4==0){
      A=c(0:(m-1),(m-1),(m+1),(m+2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }  
  }  
  if(D==9 & sk2==2 & sk3==2){  
    v=2*i*k[1]+4*k[2]+4*k[3]-2 ; m=(v-2)/2
    if(m%%4==2){
      A=0:(m+1)
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+4k2+4k3 and and k1,k2 k3 block sizes")}
  }
  if(D==10 & sk2==2 & sk3==2){
    v= 2*i*k[1]+4*k[2]+4*k[3]-4  ;m<-(v-2)/2
    if(m%%4==0){
      A=c(0:(m-1),(m-1),(m+1),(m+2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
  }
  if(D==9 & sk2==1 & sk3==2){  
    v=2*i*k[1]+2*k[2]+4*k[3]-2 ; m=(v-2)/2
    if(m%%4==2){
      A=0:(m+1)
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+2k2+4k3 and and k1,k2 k3 block sizes")}
  } 
  if(D==10 & sk2==1 & sk3==2){
    v= 2*i*k[1]+2*k[2]+4*k[3]-4  ;m<-(v-2)/2
    if(m%%4==0){
      A=c(0:(m-1),(m-1),(m+1),(m+2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
  }
  if(D==9 & sk2==2 & sk3==1){  
    v=2*i*k[1]+4*k[2]+2*k[3]-2 ; m=(v-2)/2
    if(m%%4==2){
      A=0:(m+1)
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+4k2+2k3 and and k1,k2 k3 block sizes")}
  }
  if(D==10 & sk2==2 & sk3==1){
    v= 2*i*k[1]+4*k[2]+2*k[3]-4  ;m<-(v-2)/2
    if(m%%4==0){
      A=c(0:(m-1),(m-1),(m+1),(m+2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
      A1<-grouping3(A,k,v,i,sk2,sk3)
      A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
      x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
    }
  } 
  new("stat_test", x)  
}
###################################################################
# Generation of design using sets of cyclical shifts
###################################################################
# H is an output object from CSBGND_equalsize
# The output is called using the design_CSBGND to generate design
design_CSBGND<-function(H){
  setClass( "CSBGND_design", representation("list"))
  setMethod("show", "CSBGND_design", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following is minimal CSBGND for", "v=" ,object$R[1], "and","k=",object$R[2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    for(i in 1:length(ss)){
      W<-ss[[i]]
      nr<-dim(W)[1]
      for(j in 1:nr){
        print(object$Design[[i]][[j]])
        cat("\n\n")
      }}
  })  
  v<-H$R[1]
  k<-H$R[2]
  ss<-H$S  
  treat<-(1:v)-1
  fn<-(1:v)
  G<-list()
  for(j in 1:length(ss)){ 
    W<-ss[[j]]
    nr<-dim(W)[1]
    nc<-dim(W)[2]
    D<-list()
    for(i in 1:nr){
      dd<-c()
      d1<-matrix(treat,(nc+1),v,byrow = T)
      ss1<-cumsum(c(0,W[i,]))
      dd2<-d1+ss1
      dd<-rbind(dd,dd2)
      rr<-dd[which(dd>=v)]%%v
      dd[which(dd>=v)]<-rr
      colnames(dd)<-paste("B",fn, sep="")
      rownames(dd)<-rep("",(nc+1))
      fn<-fn+v
      D[[i]]<-dd
    }
    G[[j]]<-D
  }
  x<-list(Design=G,R=H$R)
  new("CSBGND_design", x)
}




Create_Designs<-function(v,k1,k2,k3,i){  
  if((i=(v+1-(2*k2)-(2*k3))/(2*k1)) & (v+1-(2*k2)-(2*k3))%%(2*k1)==0){
    
    
    cat("MCSBNDs is possible  one set for k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3="
        ,k3,"i=",i,"put k1,k2,K3 and i value in example # 2")
    (H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=2,sk2=1,sk3=1))
    
  }else
  
  if(v%%2!=0){if((i=(v+1-(4*k2)-(2*k3))/(2*k1)) & (v+1-(4*k2)-(2*k3))%%(2*k1)==0){
    
    
    cat("MSBNDs is possible with two set of k2 and one set for k3","v=",v,
        "k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1,k2,K3 and i value in example # 2")
    (H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=2,sk2=2,sk3=1))
  }else
    if((i=(v+1-(2*k2)-(4*k3))/(2*k1)) & (v+1-(2*k2)-(4*k3))%%(2*k1)==0){
      
      
      cat("MSBNDs is possible with one set of k2 and two set for k3","v=",v,
          "k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1,k2,K3 and i value in example # 2")
      (H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=1,sk2=1,sk3=2))
    }else
      if((i=(v+1-(4*k2)-(4*k3))/(2*k1)) & (v+1-(4*k2)-(4*k3))%%(2*k1)==0){
        
        
        cat("MSBNDs is possible with two set of k2 and k3","v=",v,"k1=",k1,
            "k2=",k2,"k3=",k3,"i=",i,"put k1,k2,K3 and i value in example # 2")
        (H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=2,sk2=2,sk3=2))
      }else
    if((i=(v-1-(2*k2)-(2*k3))/(2*k1)) & (v-1-(2*k2)-(2*k3))%%(2*k1)==0){
      
     
      cat("MCBNDs is possible equale for both k2 and k3 ","v=",v,"k1=",k1,"k2=",k2,
          "k3=",k3, "i=", i,"put k1,k2,K3 and i value in example # 1")
     (H<-CBND_3diffsize(k=c(k1, k2, k3), i=i,D=1, sk2=1, sk3=1))
    }else
      if((i=(v-1-(4*k2)-(2*k3))/(2*k1)) & (v-1-(4*k2)-(2*k3))%%(2*k1)==0){
        
        
        cat("MCBNDs is possible with two sets for  k2 and one set for k3 ","v=",v,"k1="
            ,k1,"k2=",k2,"k3=",k3, "i=", i,"put k1,k2,K3 and i value in example # 1")
        (H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=1,sk2=2,sk3=1))
      }else
        if((i=(v-1-(2*k2)-(4*k3))/(2*k1)) & (v-1-(2*k2)-(4*k3))%%(2*k1)==0){
         
        
          cat("MCBNDs is possible with one set for  k2 and  two sets for k3 ","v=",v,"k1=",k1,
              "k2=",k2,"k3=",k3, "i=", i,"put k1,k2,K3 and i value in example # 1")
          (H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=1,sk2=1,sk3=2))
        }else
        
          if((i=(v-1-(4*k2)-(2*k3))/(2*k1)) & ((v-1-(4*k2)-(4*k3))%%(2*k1)==0) & i%%1==0){
            
           
            cat("MCBNDs is possible with two sets for k2 and k3 ","v=",v,"k1=",k1,
                "k2=",k2,"k3=",k3,"k3=",k3,"i=",i,"put k1,k2,K3 and i value in example # 1")
            (H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=1,sk2=2,sk3=2))
            
          }else
          
              
                    
                    cat("MBNDs and MSBNDs are not possible","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i)
    
    
    
  }else
    if(v%%2==0){
      
      if( ( i=(v-(2*k2)-(2*k3)+2)/(2*k1)) & (i %% 1==0) & (v-(2*k2)-(2*k3)+2)%%(2*k1)==0 & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
        
        cat("MCWBNDs-II and MCSGBNDS-I are possible with one set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 6 and example 9 respectivly ")
        (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=6,sk2=1,sk3=1))
        (H<-CSBGND_3diffsize(k=c(k1,k2,k3),i=i,D=9,sk2=1,sk3=1))
      }else
        if( (i=(v-(4*k2)-(2*k3)+2)/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(2*k3)+2)%%(2*k1)==0 & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
          
          cat("MCWBNDs-II and MCSGBNDS-I are possible with two set of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 6 and example 9 respectivly ")
          (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=6,sk2=2,sk3=1))
          (H<-CSBGND_3diffsize(k=c(k1,k2,k3),i=i,D=9,sk2=2,sk3=1))
        }else
          if( ( i=(v-(2*k2)-(4*k3)+2)/(2*k1)) & (i %% 1==0) & (v-(2*k2)-(4*k3)+2)%%(2*k1)==0 & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
            
            cat("MCWBNDs-II and MCSGBNDS-I are possible with one set of k2 and two set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 6 and example 9 respectivly ")
            (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=6,sk2=1,sk3=2))
            (H<-CSBGND_3diffsize(k=c(k1,k2,k3),i=i,D=9,sk2=1,sk3=2))
          }else
            if( (i=(v-(4*k2)-(4*k3)+2)/(2*k1))  & (i %% 1==0) & (v-(4*k2)-(4*k3)+2)%%(2*k1)==0 & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
              
              cat("MCWBNDs-II and MCSGBNDS-I are possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 6 and example 9 respectivly ")
              
              (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=6,sk2=2,sk3=2))
              (H<-CSBGND_3diffsize(k=c(k1,k2,k3),i=i,D=9,sk2=2,sk3=2))
              
              
              
            }else
      
      
      
      
      if((i=(v-(2*k2)-(2*k3)-2)/(2*k1)) & (i %% 1==0) & (v-(2*k2)-(2*k3)-2)%%(2*k1)==0 & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
        
      
        cat("MCPBNDs-I and MCSPBNDS-II are possible with one set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 3 and example 8 respectivly")
        (H<-CPND_3diffsize(k=c(k1,k2,k3),i=i,D=3,sk2=1,sk3=1))
        (H<-CSPBND_3diffsize(k=c(k1,k2,k3),i=i,D=8,sk2=1,sk3=1))
      }else
        if((i=(v-(4*k2)-(2*k3)-2)/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(2*k3)-2)%%(2*k1)==0 & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
          
          
          cat("MCPBNDs-I and MCSPBNDS-II are possible with two sets of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 3 and example 8 respectivly")
          (H<-CPND_3diffsize(k=c(k1,k2,k3),i=i,D=3,sk2=2,sk3=1))
          (H<-CSPBND_3diffsize(k=c(k1,k2,k3),i=i,D=8,sk2=2,sk3=1))
        }else
          if((i=(v-(2*k2)-(4*k3)-2)/(2*k1))& (i %% 1==0) & (v-(2*k2)-(4*k3)-2)%%(2*k1)==0 & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
            
           
            cat("MCPBNDs-I and MCSPBNDS-II are possible with one set of k2 and two sets of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 3 and example 8 respectivly")
            (H<-CPND_3diffsize(k=c(k1,k2,k3),i=i,D=3,sk2=1,sk3=2))
            (H<-CSPBND_3diffsize(k=c(k1,k2,k3),i=i,D=8,sk2=1,sk3=2))
          }else
            if((i=(v-(4*k2)-(4*k3)-2)/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(4*k3)-2)%%(2*k1)==0 & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
              
                           cat("MCPBNDs-I and MCSPBNDS-II are possible with two sets of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 3 and example 8 respectivly")
              (H<-CPND_3diffsize(k=c(k1,k2,k3),i=i,D=3,sk2=2,sk3=2))
              (H<-CSPBND_3diffsize(k=c(k1,k2,k3),i=i,D=8,sk2=2,sk3=2))
              
            }else
              if( (i=(v-(2*k2)-(2*k3)-2)/(2*k1)) & (i %% 1==0) & (v-(2*k2)-(2*k3)-2)%%(2*k1)==0){
               
               
                cat("MCSPBNDS-II is possible with one set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 8 ")
                (H<-CSPBND_3diffsize(k=c(k1,k2,k3),i=i,D=8,sk2=1,sk3=1))
              }else
                if((i=(v-(4*k2)-(2*k3)-2)/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(2*k3)-2)%%(2*k1)==0){
                  
                                 cat("MCSPBNDS-II is possible with two sets of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 8 ")
                  (H<-CSPBND_3diffsize(k=c(k1,k2,k3),i=i,D=8,sk2=2,sk3=1))
                }else
                  if((i=(v-(2*k2)-(4*k3)-2)/(2*k1))& (i %% 1==0) & (v-(2*k2)-(4*k3)-2)%%(2*k1)==0){
                    
                                      cat("MCSPBNDS-II is possible with one set of k2 and two sets of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 8 ")
                    (H<-CSPBND_3diffsize(k=c(k1,k2,k3),i=i,D=8,sk2=1,sk3=2))
                  }else
                    if((i=(v-(4*k2)-(4*k3)-2)/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(4*k3)-2)%%(2*k1)==0){
                      
                     
                      cat("MCSPBNDS-II is possible with two sets of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 8 ")
                      (H<-CSPBND_3diffsize(k=c(k1,k2,k3),i=i,D=8,sk2=2,sk3=2))
                      
                      
                      
                      
                      
                      
                    }else
                      
                      if((i=(v-(2*k2)-(2*k3)-4)/(2*k1))& (i %% 1==0) & (v-(2*k2)-(2*k3)-4)%%(2*k1)==0){
                        
                        
                        cat("MCPBNDs-II is possible with one sets of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 4 ")
                        (H<-CPND_3diffsize(k=c(k1,k2,k3),i=i,D=4,sk2=1,sk3=1))
                      }else
                        if( (i=(v-(4*k2)-(2*k3)-4)/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(2*k3)-4)%%(2*k1)==0){
                          
                          
                          cat("MCPBNDs-II is possible with two sets of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 4 ")
                          (H<-CPND_3diffsize(k=c(k1,k2,k3),i=i,D=4,sk2=2,sk3=1))
                        }else
                          if( (i=(v-(2*k2)-(4*k3)-4)/(2*k1)) & (i %% 1==0) & (v-(2*k2)-(4*k3)-4)%%(2*k1)==0){
                            
                            
                            cat("MCPBNDs-II is possible with one set of k2 and two sets of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 4 ")
                            (H<-CPND_3diffsize(k=c(k1,k2,k3),i=i,D=4,sk2=1,sk3=2))
                          }else
                            if( (i=(v-(4*k2)-(4*k3)-4)/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(4*k3)-4)%%(2*k1)==0){
                              
                              cat("MCPBNDs-II is possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 4 ")
                              (H<-CPND_3diffsize(k=c(k1,k2,k3),i=i,D=4,sk2=2,sk3=2))
                            }else
                              
                              
                              
                              
                              
                              if( (i=(v-(2*k2)-(2*k3))/(2*k1)) & (i %% 1==0) & (v-(2*k2)-(2*k3))%%(2*k1)==0 & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                                
                                
                                cat("MCWBNDs-I and MCSPBNDS-I are possible with one set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 5 and example 7 respectivly ")
                                (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=5,sk2=1,sk3=1))
                                (H<-CSPBND_3diffsize(k=c(k1,k2,k3),i=i,D=7,sk2=1,sk3=1))
                              }else
                                if( (i=(v-(4*k2)-(2*k3))/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(2*k3))%%(2*k1)==0 & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                                  
                                    cat("MCWBNDs-I and MCSPBNDS-I are possible with two sets of k2 and one set 
                                      of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 5 and example 7 respectivly ")
                                  (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=5,sk2=2,sk3=1))
                                  (H<-CSPBND_3diffsize(k=c(k1,k2,k3),i=i,D=7,sk2=2,sk3=1))
                                }else
                                  if( (i=(v-(2*k2)-(4*k3))/(2*k1)) & (i %% 1==0) &(v-(2*k2)-(4*k3))%%(2*k1)==0 & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                                    
                                  
                                    cat("MCWBNDs-I and MCSPBNDS-I are possible with one set of k2 and two sets of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 5 and example 7 respectivly ")
                                    (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=5,sk2=1,sk3=2))
                                    (H<-CSPBND_3diffsize(k=c(k1,k2,k3),i=i,D=7,sk2=1,sk3=2))
                                  }else
                                    if( (i=(v-(4*k2)-(4*k3))/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(4*k3))%%(2*k1)==0 & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                                      
                                      cat("MCWBNDs-I and MCSPBNDS-I are possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 5 and example 7 respectivly ")
                                      (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=5,sk2=2,sk3=2))
                                      (H<-CSPBND_3diffsize(k=c(k1,k2,k3),i=i,D=7,sk2=2,sk3=2))
                                    }else
                                      if( (i=(v-(2*k2)-(2*k3))/(2*k1)) & (i %% 1==0) & (v-(2*k2)-(2*k3))%%(2*k1)==0 & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                                        
                                       
                                        cat("MCWBNDs-I is possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example # 5")
                                        (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=5,sk2=2,sk3=2))
                                      }else
                                        if( (i=(v-(4*k2)-(2*k3))/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(2*k3))%%(2*k1)==0 & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                                          
                                          cat("MCWBNDs-I is possible with two sets of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example # 5")
                                          (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=5,sk2=2,sk3=1))
                                        }else
                                          if( (i=(v-(2*k2)-(4*k3))/(2*k1)) & (i %% 1==0) & (v-(2*k2)-(4*k3))%%(2*k1)==0 & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                                            
                                            
                                            cat("MCWBNDs-I is possible with one set of k2 and two sets of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example # 5")
                                            (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=5,sk2=1,sk3=2))
                                            }else
                                            if( (i=(v-(2*k2)-(4*k3))/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(4*k3))%%(2*k1)==0 & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                                              
                                             
                                              cat("MCWBNDs-I is  possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example # 5")
                                           
                                              (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=5,sk2=2,sk3=2))
                                               }else
                                              
                                             
                                                      
                                                      if( (i=(v-(2*k2)-(2*k3)+2)/(2*k1)) & (i %% 1==0) & (v-(2*k2)-(2*k3)+2)%%(2*k1)==0){
                                                        
                                                       
                                                        cat("MCWBNDs-II is possible with one set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 6 ")
                                                        (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=6,sk2=1,sk3=1))
                                                        }else
                                                        if( (i=(v-(4*k2)-(2*k3)+2)/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(2*k3)+2)%%(2*k1)==0){
                                                          
                                                         
                                                          cat("MCWBNDs-II is possible with two set of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 6 ")
                                                          (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=6,sk2=2,sk3=1))
                                                          }else
                                                          if( (i=(v-(2*k2)-(4*k3)+2)/(2*k1)) & (i %% 1==0) & (v-(2*k2)-(4*k3)+2)%%(2*k1)==0 ){
                                                            
                                                            
                                                            cat("MCWBNDs-II is possible with one set of k2 and two set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 6 ")
                                                            (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=6,sk2=1,sk3=2))  
                                                             }else
                                                            if( (i=(v-(4*k2)-(4*k3)+2)/(2*k1)) & (i %% 1==0) &(v-(4*k2)-(4*k3)+2)%%(2*k1)==0 & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                                                              
                                                              
                                                              cat("MCWBNDs-II is possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example 6 ")
                                                              (H<-CWBND_3diffsize(k=c(k1,k2,k3),i=i,D=6,sk2=2,sk3=2))
                                                              
                                                              
                                                              
                                                              
                                                            }else
                                                              
                                                              
                                                              
                                                              if( (i=(v-(2*k2)-(2*k3)+4)/(2*k1)) & (i %% 1==0) & (v-(2*k2)-(2*k3)+4)%%(2*k1)==0 ){
                                                                
                                                               
                                                                cat("MCGBNDs-II is possible with one set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example # 10")
                                                                (H<-CSBGND_3diffsize(k=c(k1,k2,k3),i=i,D=10,sk2=1,sk3=1))
                                                                 }else
                                                                if( (i=(v-(4*k2)-(2*k3)+4)/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(2*k3)+4)%%(2*k1)==0){
                                                                  
                                                                  
                                                                  cat("MCGBNDs-II is possible with two set of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example # 10")
                                                                  (H<-CSBGND_3diffsize(k=c(k1,k2,k3),i=i,D=10,sk2=2,sk3=1)) 
                                                                  }else
                                                                  if( (i=(v-(2*k2)-(4*k3)+4)/(2*k1)) & (i %% 1==0) & (v-(2*k2)-(4*k3)+4)%%(2*k1)==0 ){
                                                                    
                                                                  
                                                                    cat("MCGBNDs-II is possible with one set of k2 and two sets of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example # 10")
                                                                    (H<-CSBGND_3diffsize(k=c(k1,k2,k3),i=i,D=10,sk2=1,sk3=2))
                                                                    }else
                                                                    if( (i=(v-(4*k2)-(4*k3)+4)/(2*k1)) & (i %% 1==0) & (v-(4*k2)-(4*k3)+4)%%(2*k1)==0){
                                                                      
                                                                     
                                                                      cat("MCGBNDs-II is possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i,"put k1, k2,K3 and i value in example # 10")
                                                                      (H<-CSBGND_3diffsize(k=c(k1,k2,k3),i=i,D=10,sk2=2,sk3=2))
                                                                      }else
                                                                      cat("MCPBNDs-I,MCPBNDs-II,MCSPBNDs-II, MCWBNDs-I,MCWBNDs-II,MCSBGNDs-I and MCSBGNDs-II are not possible")
    }
  
}


###############################################################################
# Examples: Using CBND_3diffsize function to obtain the set(s) of shifts
# for construction of circular weakly balance neighbour design for block of 
# three different sizes (k1, K2 and k3)
###############################################################################
# D=1: Circular Balanced Neighbor Designs
# D=2: Circular Strongly Balanced Neighbor Designs
#MCBNDs and MCSBNDs, say here D1 & D2 respectively, for v odd

(H<-Create_Designs(v=99,k1=9,k2=8,k3=6))

# MCBNDs#1
(H<-CBND_3diffsize(k=c(10,9,5),i=3,D=2,sk2=1,sk3=1))
#MCSBNDs#2
(H<-CBND_3diffsize(k=c(7,5,4),i=1,D=2,sk2=1,sk3=1))
# Example#3
(H<-CPND_3diffsize(k=c(10,7,3),i=3,D=3,sk2=1,sk3=1))
#Example#4
(H<-CPND_3diffsize(k=c(10,7,3),i=3,D=4,sk2=1,sk3=1))
# Example#5
(H<-CWBND_3diffsize(k=c(12,7,4),i=1,D=5,sk2=1,sk3=1))
# Example#6
(H<-CWBND_3diffsize(k=c(6,5,3),i=4,D=6,sk2=2,sk3=2))
# Example#7
(H<-CSPBND_3diffsize(k=c(12,6,4),i=3,D=7,sk2=2,sk3=1))
#Example#8
(H<-CSPBND_3diffsize(k=c(6,5,3),i=5,D=8,sk2=2,sk3=2))
# Example#9
(H<-CSBGND_3diffsize(k=c(10,8,5),i=3,D=9,sk2=1,sk3=2))
#Example#10
(H<-CSBGND_3diffsize(k=c(8,5,4),i=6,D=10,sk2=1,sk3=1))
(design_CSBGND(H))




