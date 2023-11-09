
################################################################################
# CQRNDs2_equalsize: Circular Quasi Rees neighbor designs for block of equal size(K)
################################################################################
# Algorithm from paper:
# Akbar Firdos1, Mahmood Ul Hassan2, Farrukh Jamal1, Hurria Ali1, Khadija Noreen1 and Rashid Ahmed1. Efficient Balanced and Strongly 
# Balanced Neighbor Designs through  CQRNDs2 -designs Generated with R 
# Coded by Ali et al., 01-03-2022 to 05-07-2022
# Version 2.0  (2022-07-05)
################################################################
# Division of adjusted A in i-1 groups of size k and one
# of size K-2 to get the set(s) of shifts
################################################################
grouping1<-function(A,k,v,i){
  bs<-c()
  z=0;f=1
  A1=A
  while(f<=((i+1)-1)){
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
  bs1<-t(apply(bs,1,sort))
  bs1<-cbind(bs1,rowSums(bs),rowSums(bs)/v)
  rownames(bs1)<-paste("G",1:i, sep="")
  colnames(bs1)<-c(paste(1:k, sep=""),"sum" ,"sum/v")
  bs2<-t(apply(bs,1,sort))
  bs2<-delmin(bs2)
  list(B1=list(bs2),B2=list(bs1),B3=A1)
}
#######################################################################
# Obtaing set(s) of shifts by deleting smallest value of each group
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
  rownames(fs)<-paste("S",1:n, sep="")
  colnames(fs)<-rep("",c)
  return(fs)
}
################################################################################
# Selection of adjusted A and the set(s) of shifs to obtain Circular  
# balance neighbor design for block of equal size. 
################################################################################
# D=1: Circular Balanced Neighbor Designs
# D=2: Circular Strongly Balanced Neighbor Designs
#   K: Block sizes
#   i: Number of set of shifts for K
CBND_equalsize<-function(k,i,D=1){
  if(k<=2) stop("k= Block size: Block size must be greater than 2")
  if(i<=0) stop("i= Must be a positive integer")
  setClass( "stat_test", representation("list"))
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    if(D==2){
      cat("Following are required sets of shifts to obtain the 
minimal CSBND for", "v=" ,object[[3]][1], "and","k=",object[[3]][2], "\n")}
    if(D==1){
      cat("Following are required sets of shifts to obtain the 
minimal CBND for", "v=" ,object[[3]][1], "and","k=",object[[3]][2], "\n")}
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
  })
  if(D==2){ 
    v=2*i*k-1; m=(v-1)/2
    if(m%%8==0){
      j=m/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(j-1),(j+1):m,(v-j))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==1){
      j=(m-1)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==2){
      j=(m-2)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==3){
      j=(m-3)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(m-j-1),(m-j+1):m,(v-(m-j)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==4){
      j=(m-4)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:j,(j+2):(m-1),(m+1),(v-(j+1)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==5){
      j=(m-5)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==6){
      j=(m-6)/8
      if(j<0) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==7){
      j=(m-7)/8
      if(j<1) {return("Conditions are not satisfied for CSBND")}
      A=c(0:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
  }
  if(D==1){
    v=2*i*k+1; m=(v-1)/2
    if(m%%8==0){
      j=m/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(j-1),(j+1):m,(v-j))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==1){
      j=(m-1)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==2){
      j=(m-2)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==3){
      j=(m-3)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(m-j-1),(m-j+1):m,(v-(m-j)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==4){
      j=(m-4)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:j,(j+2):(m-1),(m+1),(v-(j+1)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==5){
      j=(m-5)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==6){
      j=(m-6)/8
      if(j<0) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%8==7){
      j=(m-7)/8
      if(j<1) {return("Conditions are not satisfied for CBNDs")}
      A=c(1:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
  }
  new("stat_test", x)
}
##################################################################
# Generation of design using sets of cyclical shifts
###################################################################
# H is an output object from CWBND_equalsize
# The output is called using the design_CWBND to generate design
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
  new("CBND_design", x)
}
# D=3: minimal CGNDs in which v/2 unordered pairs do not appear
# D=4: minimal CGNDs in which 3v/2 unordered pairs do not appear 
#   K: Block sizes
#   i: Number of sets of shifts for K
CGND_equalsize<-function(k,i,D=3){
  if(k<=3) stop("k= Block size: Block size must be greater than 3")
  if(i<=0) stop("i= Must be a positive integer")
  setClass( "stat_test", representation("list"))
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following are required sets of shifts to obtain the 
minimal CPND for", "v=" ,object[[3]][1], "and","k=",object[[3]][2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    print(object[[1]])
  })
  if(D==3){  
    v=2*i*k+2; m=(v-2)/2
    if(m%%4==0){
      A=1:m
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==1 |  m%%4==2){return("The minimal CPNDs in which v/2 unordered pairs cannot be constructed for v=2ik+2 and k=block size")}
  }
  if(D==4){
    v=2*i*k+4; m=(v-2)/2
    if(m%%4==0){
      A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(2:(m-2),m,(2*m+1))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
  }
  new("stat_test", x)
}
# D=5: Minimal CWBNDs in which v/2 of the ordered pairs appear twice as neighbors
# D=6: Minimal CWBNDs in which 3v/2 of the unordered pairs appear twice as neighbors
#   K: Block sizes
#   i: Number of set of shifts for K
CWBND_equalsize<-function(k,i,D=5){
  if(k<3) stop("k= Block size: Block size must be greater than 2")
  if(i<=0) stop("i= Must be a positive integer")
  setClass( "stat_test", representation("list"))
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following are required sets of shifts to obtain the 
minimal CWBND for", "v=" ,object[[3]][1], "and","k=",object[[3]][2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
  })
  if(D==5){ 
    v=2*i*k; m=(v-2)/2
    if(m%%4==2){
      A=1:(m+1)
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(1:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==1 |  m%%4==0){return("The minimal CWBNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik and k=block size")}
  }
  if(D==6){
    v=2*i*k-2; m=(v-2)/2
    if(m%%4==0){
      A=c(1:(m-1),(m-1),(m+1),(m+2))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(1:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(1:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(1:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
  }
  new("stat_test", x)
}
##################################################################
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
CSPBND_equalsize<-function(k,i,D=1){
  if(k<=3) stop("k= Block size: Block size must be greater than 3")
  if(i<=0) stop("i= Must be a positive integer")
  setClass( "stat_test", representation("list"))
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following are required sets of shifts to obtain the 
        minimal CSPBNDs for", "v=" ,object[[3]][1], "and","k=",object[[3]][2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
  })
  if(D==7){
    v=2*i*k; m=(v-2)/2
    if(m%%4==0){
      A=0:m
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(0:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==1 | m%%4==2){return("The minimal CSPBNDs in which v/2 unordered pairs cannot be constructed for v=2ik+2 and k=block size")}
  }
  if(D==8){
    v=2*i*k+2; m=(v-2)/2
    if(m%%4==0){
      A=c(0:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(0:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(0,2:(m-2),m,(2*m+1))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(0:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
  }
  new("stat_test", x)
}
##################################################################
# Generation of design using sets of cyclical shifts
###################################################################
# H is an output object from CSPBND_equalsize
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

# D=9: Minimal CSBGNDs in which v/2 of the ordered pairs appear twice as neighbors
# D=10: Minimal CSBGNDs in which 3v/2 of the unordered pairs appear twice as neighbors
#   K: Block sizes
#   i: Number of set of shifts for K
CSBGND_equalsize<-function(k,i,D=9){
  if(k<=3) stop("k= Block size: Block size must be greater than 3")
  if(i<=0) stop("i= Must be a positive integer")
  setClass( "stat_test", representation("list"))
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following are required sets of shifts to obtain the 
minimal CSBGNDs for", "v=" ,object[[3]][1], "and","k=",object[[3]][2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
  })
  if(D==9){
    v=(2*i*k)-2; m=(v-2)/2
    if(m%%4==2){
      A=0:(m+1)
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==3){
      A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==1 |  m%%4==0){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik and k=block size")}
  }
  if(D==10){
    v=2*i*k-4; m=(v-2)/2
    if(m%%4==0){
      A=c(0:(m-1),(m-1),(m+1),(m+2))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==1){
      A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==2){
      A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(m%%4==3){
      A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
      A1<-grouping1(A,k,v,i)
      A2<-c(v,k);names(A2)<-c("V","K")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
  }
  new("stat_test", x)
}
##################################################################
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

Create_Designs<-function(v,k,i){
  if(v%%2!=0){
    if((v-1)%%(2*k)==0){
      i=(v-1)/(2*k)
      cat("MCBNDs is possible for","v=",v,"k=",k,"i=",i,"put k and i value in example # 1",rep("", 30))
      (H<-CBND_equalsize(k=k,i=i,D=1))
      }else
        if((v+1)%%(2*k)==0){
          i=(v+1)/(2*k)
          cat("MCSBNDs is possible for","v=",v,"k=",k,"i=",i,"put k and i value in example # 2",rep("", 30))
          (H<-CBND_equalsize(k=k,i=i,D=2))
        }else
          cat("MCBNDs and MCSBNDs are not possible and change value ","v=",v,"and", "k=",k, " run again",rep("", 30))
    
  }else
    if(v%%2==0){
      
      if((v-2)%%(2*k)==0 & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
        i=(v-2)/(2*k)
        cat("MCPBNDs-I and MCSPBNDS-II are possible for","v=",v,"k=",k,"i=",i,"put k and i value in example 3 and example 8 respectivly  ",rep("",15))
       (H<-CSPBND_equalsize(k=k,i=i,D=8))
        (H<-CGND_equalsize(k=k,i=i,D=3))
      }else
        
        if((v-2)%%(2*k)==0){
          i=(v-2)/(2*k)
          cat(" MCSPBNDs-II is possible for","v=",v,"k=",k,"i=",i,"put k and i value in example # 8",rep("",15))
          (H<-CSPBND_equalsize(k=k,i=i,D=8))
        }else
          
          if((v-4)%%(2*k)==0){
            i=(v-4)/(2*k)
            cat("MCPBNDs-II is possible for","v=",v,"k=",k,"i=",i,"put k and i value in example # 4",rep("",15))
            (H<-CGND_equalsize(k=k,i=i,D=4))
          }else
            if(v%%(2*k)==0 & ((v-2)/2)%%4==3){
              i=v/(2*k)
              cat("MCWBNDs-I and MCSPBNDs-I is possible for","v=",v,"k=",k,"i=",i,"put k and i value in example 5 and example 7 respectivly ",rep("",15))
              (H<-CSPBND_equalsize(k=k,i=i,D=7))
              (H<-CWBND_equalsize(k=k,i=i,D=5))
              }else
              if(v%%(2*k)==0& ((v-2)/2)%%4==2){
                i=v/(2*k)
                cat("MCWBNDs-I  is possible for","v=",v,"k=",k,"i=",i,"put k and i value in example # 5",rep("",15))
                (H<-CWBND_equalsize(k=k,i=i,D=5))
              }else
                if(v%%(2*k)==0& ((v-2)/2)%%4==0 ){
                  i=v/(2*k)
                  cat("MCSPBNDs-I is possible for","v=",v,"k=",k,"i=",i,"put k and i value in example # 7",rep("",15))
                  (H<-CSPBND_equalsize(k=k,i=i,D=7))
                  }else
                  if((v+2)%%(2*k)==0 & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                    i=(v+2)/(2*k)
                    cat("MCWBNDs-II and MCSBGNDs-I are possible for","v=",v,"k=",k,"i=",i,
                        "put k and i value in example # 6 and example # 9 respictively ",rep("",15))
                    (H<-CSBGND_equalsize(k=k,i=i,D=9))
                    (H<-CWBND_equalsize(k=k,i=i,D=6))
                  
                  }else
                    if((v+2)%%(2*k)==0){
                      i=(v+2)/(2*k)
                      cat("MCWBNDs-II  is possible for","v=",v,"k=",k,"i=",i,"put k and i value in example # 6",rep("",15))
                      (H<-CWBND_equalsize(k=k,i=i,D=6))
                    
                    }else
                      
                      if((v+4)%%(2*k)==0){
                        i=(v+4)/(2*k)
                        cat("MCSGNDs-II is possible for","v",v,"k=",k,"i=",i,
                            "put k and i value in example # 10",rep("",15))
                        (H<-CSBGND_equalsize(k=k,i=i,D=10))
                        
                      }else
                        cat("MCSPBNDs-II,MCPBNDs,MCPBNDs-II,MCWBNDs-I,MCSPBNDs-I,MCSBGNDs-I and MCSBGNDs-II  are not possible for given  value","v=",v,"and", "k=",k,rep("",15))
    }
  
}

cat("Run the program to create designs and put the v and k value in the Create_Designs Function")

################################################################################
# Examples: Using CBND_equal size function to obtain the set(s) of shifts
# for construction of circular balance neighbor design for equal block  
# sizes (k)
################################################################################

#Run this program, check possible design and get i value
(H<-Create_Designs(v=71, k=4))
(H<-CGND_equalsize(k=6, i=2, D=4))
(H<-CWBND_equalsize(k=5, i=1, D=6))
(H<-CSPBND_equalsize(k=4, i=3, D=7))
(H<-CSPBND_equalsize(k=8, i=1, D=8))
(H<-CSBGND_equalsize(k=8,i=1,D=9))
#MCBNDs and MCSBNDs, say here D1 & D2 respectively, for v odd
# example#1
#put k and i value in example#1
(H<-CBND_equalsize(k=4,i=3,D=1))
(D<-design_CBND(H))
# example#2
#put k and i value in example#2
(H<-CBND_equalsize(k=4,i=3,D=2))
(D<-design_CBND(H))
#MCPBNDs-I and MCPBNDs-II, say D3 & D4 for v even
#put k and i value in example#3
# example#3
(H<-CGND_equalsize(k=4,i=2,D=3))
#put k and i value in example#4
# example #4
(H<-CGND_equalsize(k=6,i=2,D=4))

#	MCWBNDs-I and MCWBNDs-II, say D5 & D6 for v even
#put k and i value in example#5
# example#5
(H<-CWBND_equalsize(k=4,i=2,D=5))
(D<-design_CWBND(H))
#put k and i value in example#6
# example#6
(H<-CWBND_equalsize(k=6, i=2,D=6))
(D<-design_CWBND(H))
#	MCSPBNDs-I and MCSPBNDs-II, say D7 & D8 for v even
#put k and i value in example#7
# example#7
(H<-CSPBND_equalsize(k=9,i=4,D=7))
(D<-design_CSPBND(H))
#put k and i value in example#8
# example#8
(H<-CSPBND_equalsize(k=7,i=5,D=8))
(D<-design_CSPBND(H))
#	MCSBGNDs-I and MCSBGNDs-II, say D9 & D10 for v even 
#put k and i value in example#9
# example#9
(H<-CSBGND_equalsize(k=4,i=2,D=9))
(D<-design_CSBGND(H))
#put k and i value in example#10
# example #10
(H<-CSBGND_equalsize(k=5,i=3,D=10))
design_CSBGND(H)



