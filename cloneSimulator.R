#R code

#The error rate of E.coli polymerase system is 10^-9~10^-11
#DNA replication fidelity in Escherichia coli: a multi-DNA polymerase affair.

MutRate <- 10^-4
Gen <- 14

seed <- as.double(1)
RANDU <- function() {
  seed <<- ((2^16 + 3) * seed) %% (2^31)
  seed/(2^31)
}

seq <- matrix(nrow=Gen,ncol=88*2^Gen)
for(k in 1:88) {
  seq[1,k] <- 0
}
flag <- 1
for(i in 2:Gen) {
  for(j in 1:2^(i-1)) {
    for(k in 1:88) {
      if (flag==1) {
        seq[i,k+(j-1)*88] <- seq[i-1,k+floor((j-1)/2)*88]
      }
      else if ((RANDU() < MutRate) && (flag==0)) {
        seq[i,k+(j-1)*88] <- 1
      }
      else seq[i,k+(j-1)*88] <- seq[i-1,k+floor((j-1)/2)*88]
    }
    if (flag==1) { flag <- 0 }
    else if (flag==0) { flag <- 1 }
  }
}
error <- vector()
for(k in 1:88) {
  error[k] <- 0
}
for(j in 1:2^(Gen-1)) {
  for(k in 1:88) {
    error[k] <- error[k] + seq[Gen,k+(j-1)*88]
  }
}
error <- error/2^(Gen-1)
plot(error,log='y')
