
# setwd("C:/Users/jisaa/Documents/400")



# 1A
TrafficRaw <- read.delim("./webtraffic.txt",sep = "\t")
row1<-as.vector(colSums(TrafficRaw)[1:9])
row2<-as.vector(colSums(TrafficRaw)[10:18])
row3<-as.vector(colSums(TrafficRaw)[19:27])
row4<-as.vector(colSums(TrafficRaw)[28:36])
row5<-as.vector(colSums(TrafficRaw)[37:45])
row6<-as.vector(colSums(TrafficRaw)[46:54])
row7<-as.vector(colSums(TrafficRaw)[55:63])
row8<-as.vector(colSums(TrafficRaw)[64:72])
row9<-as.vector(colSums(TrafficRaw)[73:81])
Traffic<-rbind(row1,row2,row3,row4,row5,row6,row7,row8,row9)
Traffic[9,1]=1000
print(Traffic)

#1B
# This Markov chain is irreduciable because all states can communicate with one another.
# Without artificial line of code ("Traffic[9,1]=1000") the chain would not be irreduciable, for 9 wouldn't be able to communicate with 1
# And because all states can be communicated two ways, the chain is recurrent, not transient.
# Also the chain is aperiodic because from state 2-8 can return to itself, and the greatest divisior of state 1 - 9 to itself is one.
# And since the chain is a periodic and recurrent, the chain is ergodic.

# 1C
row1<-as.vector(colSums(TrafficRaw)[1:9])/sum(colSums(TrafficRaw)[1:9])
row2<-as.vector(colSums(TrafficRaw)[10:18])/sum(colSums(TrafficRaw)[10:18])
row3<-as.vector(colSums(TrafficRaw)[19:27])/sum(colSums(TrafficRaw)[19:27])
row4<-as.vector(colSums(TrafficRaw)[28:36])/sum(colSums(TrafficRaw)[28:36])
row5<-as.vector(colSums(TrafficRaw)[37:45])/sum(colSums(TrafficRaw)[37:45])
row6<-as.vector(colSums(TrafficRaw)[46:54])/sum(colSums(TrafficRaw)[46:54])
row7<-as.vector(colSums(TrafficRaw)[55:63])/sum(colSums(TrafficRaw)[55:63])
row8<-as.vector(colSums(TrafficRaw)[64:72])/sum(colSums(TrafficRaw)[64:72])
row9<-as.vector(colSums(TrafficRaw)[73:81])
p<-rbind(row1,row2,row3,row4,row5,row6,row7,row8,row9)
p[9,1]<-1000/1000
print(p)

# 1D
a <- c(1,rep(0,8))
five_click_prob <- a %*% p %*% p %*% p %*% p %*% p
page5prob<-five_click_prob[5]
print(page5prob)

# 1E
Q<-t(p)-diag(9)
Q[9,]<-rep(1,9)
rhs=c(rep(0,8),1)
Pi=solve(Q,rhs)
print(Pi)

# 1F
B=p[1:8,1:8]
Q=diag(8)-B
rhs=rep(c(0.1,2,3,5,5,3,3,2))
m=solve(Q,rhs)
print(m[1])


#2a
samplefunc<-function(lambda)
   {
   var_px=1/(lambda^2)
   tolerance=10^(-3)
   n=var_px/tolerance^2*100
   return (n)
   }
print(samplefunc(1))

#2b
lambdalist<-c(1,2,4)
for (lambda in lambdalist)
   {
   n<- samplefunc(lambda)
   x<- runif(n,0, 1)
   y<- -log(x)/lambda
   z= sin(y)/lambda
   integral<-sum(z)/n
   print(integral)
   }


#3a
#Metropolis-Hastings is a better algorithm to use for drawing from gamma
#distribution because its shape is not symmetric and is not too large for the alogorithm

#3b
qfunction<-function(x){
   return (rexp(1,rate = x))
}
pfunction<- function(x,shape,scale){
   return (x^(shape-1)*exp(-1*x/scale))
}

n<-15000
x<-c(1)
for (i in seq(1, n))
   {
      x_p<-tail(x,n=1)
      x_c=qfunction(x_p)
      a<-(pfunction(x_c,shape=2,scale=2)*(qfunction(x_p)))/
             (pfunction(x_p,shape=2,scale=2)*qfunction(x_c))
      if (runif(1)<min(1,a))
         {
         x[i]<-x_c  
         }
      else
         {
         x[i]<-x_p
         }
   }
output<-x[seq(1, length(x), 100)][50:150]

print(output)

#3c
hist(output,20,prob=TRUE)
curve(dgamma(x, shape = 2, scale = 2), add=TRUE,col="red")

observation = seq(1,100)
cor(output,observation)
#According to the histogram, the output seems to be aligned with random samples
#from gamma distribution. While there are more samples in the right end of the distribution
#then the gamma distribution suggests, I believe the shape of the sample will follow the distribution
#more closely when the number of sample increase.
#The low correlation between observation also suggests randomness in the sample as well.