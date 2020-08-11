##########################################################
####  Model for schattering pattern of Pseudanabaena #####
####  Based on Sugiら(Nat. Comm., 2019)              #####
##########################################################
# 2020.2.25 Hisamoto 

#########################################################
# Parameter setting
#########################################################
v0 <- 0.1           # Mean velocity 0.1
rr <- 0.2           # Range of repulsive force 0.2
l <- 1.0            # Length of filament 1.0

tmax <- 2000        # Maximum steps 2000
Areasize <- 2       # Initial area range 2
N<- 500             # Particle number 500

ka <- 0.1           # Coefficient of attractive force 0.1
kr <- 10            # Coefficient of repulsive force 10

w0 <- 0.0005        # Mean rotation rate 0.0005
sigw0 <- 0.2        # Standard deviation of rotation rate 0.2
tau <- 500          # Correlation time of rotation rate 500
#########################################################
# Initialization
#########################################################
X <- numeric(0)
Y <- numeric(0)
theata <- numeric(0)
Omega <- numeric(0)
Xt <- numeric(0)
Yt <- numeric(0)
theatat <- numeric(0)
Omegat <- numeric(0)


# particles randomly distributed in initial
X <- runif(N,min = 0, max = Areasize )
Y <- runif(N,min = 0, max = Areasize )


# Initial direction
theata <- rnorm(N,mean = 0, sd = pi)
Omega <- rnorm(N,mean = w0, sd = sigw0)

Log <- list(cbind(X,Y))
Xt <- X
Yt <- Y
theatat <- theata
Omegat <- Omega


#########################################################
# Calculation
#########################################################
for(t in 2:tmax){
  
  # Distance between particles
  Distance <- dist(cbind(Xt,Yt))
  D <- as.matrix(Distance)  
  
  for (i in 1:N){
    
    # Distance between particle i and others
    d <- D[,i]
    
    # Relationship between particles
    Interaction_distance <- cbind(theatat,d,Xt,Yt) 
    
    # Range which attraction force works
    At_force <-subset(Interaction_distance, d<l & d>rr )
    
    # Range which repulsion force works
    Rp_force <-subset(Interaction_distance, d<rr & d>0 )
    
    # Calculation of rotating ratio (equation (3))
    Omega[i] <- Omegat[i] - (Omegat[i] - w0) / tau + sqrt(2 / tau) * sigw0 * rnorm(1,mean = 0, sd = 1)
    
    if(length(At_force[,1]) > 0){ 
      
      # Calculation of moving direction (equation (2))
      theata[i] <- theatat[i] + Omega[i] +  sum(sin(2 * (At_force[,1] - theatat[i]))) / length(At_force[,1])
     
      # Calculation of attraction force (equation (5))
      Fax <- sum(ka * ( At_force[,3] - Xt[i]) /  At_force[,2]^2)/ length(At_force[,1])
      Fay <- sum(ka * ( At_force[,4] - Yt[i]) /  At_force[,2]^2)/ length(At_force[,1])
  
      
    # In case no particle is in the range of d < l
    }else{
      theata[i] <- theatat[i] +  Omega[i]
      Fax <- 0
      Fay <- 0 
      
    }
    
    if(length(Rp_force[,1]) > 0){ 
      # Calculation of repulsive force (equation (4))
      Frx <- sum(kr*(Rp_force[,2]-rr)*(Rp_force[,3] - Xt[i]) / Rp_force[,2])
      Fry <- sum(kr*(Rp_force[,2]-rr)*(Rp_force[,4] - Yt[i]) / Rp_force[,2])
      
    }else{
      Frx <- 0
      Fry <- 0 
      
    }
    
    # Calculation of the position (equation (1))
    X[i] <- Xt[i] + v0 * cos(theata[i])  + v0 *Fax + v0 *Frx
    Y[i] <- Yt[i] + v0 * sin(theata[i])  + v0 *Fay + v0 *Fry
    
    
  }
  

  Log[[t]] <- cbind(X,Y)
  
  Xt <- numeric(0)
  Yt <- numeric(0)
  theatat <- numeric(0)
  Omegat <- numeric(0)
  
  Xt <- X
  Yt <- Y 
  theatat <- theata
  Omegat <- Omega
}


#########################################################
####  Output                                 
####  Warning: tmax JPEG pictures will saved on disk.          
#########################################################

 for(t in 1:tmax){
   png(paste("Commet_",t,".png",sep=""), width = 1000, height = 1000)               # 描画デバイスを開く
   plot(Log[[t]][,1],Log[[t]][,2],xlim=c(-50,50),ylim=c(-50,50),pch=16,cex=0.4)      # グラフを描く
   dev.off()
 }

# Writing out the csv file of position x and y
tempx <- matrix(0,nrow=N,ncol=tmax)
for(p in 1:N){
  tempx[p,] <- sapply(Log,"[[",p,1)
}
write.table(tempx,"X.csv",sep=",")

tempy <- matrix(0,nrow=N,ncol=tmax)
for(p in 1:N){
  tempy[p,] <- sapply(Log,"[[",p,2)
}
write.table(tempy,"Y.csv",sep=",")