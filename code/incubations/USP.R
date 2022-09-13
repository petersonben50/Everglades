USP<-function(y,x,C=1000){

T.A<-array(0,dim=c((C+1),1))
T.B<-array(0,dim=c((C+1),1))
T.AB.a<-array(0,dim=c((C+1),1))
T.AB.b<-array(0,dim=c((C+1),1))
T.AB<-array(0,dim=c((C+1),1))

NI<-array(0,dim=c(C,1))

##########PRE-PROCESSING#########

lA<-unique(x[,1])
lB<-unique(x[,2])

A<-length(lA)
B<-length(lB)

n<-length(y)/(A*B)

Y<-array(0,dim=c(A,B,n))



#### Set up the blocks of results ####
for(i in 1:A){
for(j in 1:B){

Y[i,j,] <- y[x[,1] == lA[i] & x[,2] == lB[j]]
}
}

##########COMPUTING THE OBSERVED STATISTICS########################

#Factor A

for(i in 1:(A-1)){
for(s in (i+1):A){
T.A[1] <- T.A[1] + (sum(Y[i,,])-sum(Y[s,,]))^2
}
}

#Factor B

for(j in 1:(B-1)){
for(h in (j+1):B){
T.B[1]<-T.B[1]+(sum(Y[,j,])-sum(Y[,h,]))^2
}
}

#Interaction (a)

for(i in 1:(A-1)){
for(s in (i+1):A){

for(j in 1:(B-1)){
for(h in (j+1):B){

T.AB.a[1]<-T.AB.a[1] + (sum(Y[i,j, ])-sum(Y[s,j, ])-sum(Y[i,h, ])+sum(Y[s,h, ]))^2

}
}
}
}

#Interaction (b)

T.AB.b[1]<-T.AB.a[1]



########COMPUTING THE PROBABILITY OF ni EXCHANGES

Pr.A<-array()
Pr.B<-array()

if((n %% 2) > 0){
for(ni in 0:((n-1)/2)){

Pr.A[ni+1]=choose(n,ni)^(B*A*(A-1))
Pr.B[ni+1]=choose(n,ni)^(A*B*(B-1))

}
}

if((n %% 2) == 0){

for(ni in 0:(n/2-1)){
Pr.A[ni+1]<-choose(n,ni)^(B*A*(A-1))
Pr.B[ni+1]<-choose(n,ni)^(A*B*(B-1))
}
Pr.A<-c(Pr.A,((choose(n,(n/2))^(2*B))/2)^(A*(A-1)/2))
Pr.B<-c(Pr.B,((choose(n,(n/2))^(2*A))/2)^(B*(B-1)/2))
}

####
#Pr.A<-c(1,729)
#Pr.B<-c(1,32)
####

Pr.A<-Pr.A/sum(Pr.A)
Pr.B<-Pr.B/sum(Pr.B)

for(ni in 2:length(Pr.A)){
Pr.A[ni]<-Pr.A[ni-1]+Pr.A[ni]
Pr.B[ni]<-Pr.B[ni-1]+Pr.B[ni]
}

######OBTAINING THE SYNCHRONIZED PERMUTATION DISTRIBUTION#######




for(cc in 2:(C+1)){
#print(cc)

Y.perm<-Y


u<-runif(1)

ni=0
while(u>Pr.A[ni+1]){
ni=ni+1
}


NI[cc]<-ni
###Column permutations


for(i in 1:(A-1)){
for(s in (i+1):A){




if(ni>0){
for(j in 1:B){


pool1<-sample(Y[i,j,])
pool2<-sample(Y[s,j,])


Y.perm[i,j,]<-c(pool2[c(1:ni)],pool1[-c(1:ni)])
Y.perm[s,j,]<-c(pool1[c(1:ni)],pool2[-c(1:ni)])
}
}

if(ni==0){Y.perm<-Y}

T.A[cc]<-T.A[cc]+(sum(Y.perm[i,,])-sum(Y.perm[s,,]))^2


for(j in 1:(B-1)){
for(h in (j+1):B){

T.AB.a[cc]<-T.AB.a[cc] + (sum(Y.perm[i,j, ])-sum(Y.perm[s,j, ])-sum(Y.perm[i,h, ])+sum(Y.perm[s,h, ]))^2

}
}
}
}


###Row permutations

u<-runif(1)
ni=0
while(u>Pr.B[ni+1]){
ni=ni+1
}



Y.perm<-Y

for(j in 1:(B-1)){
for(h in (j+1):B){





if(ni>0){
for(i in 1:A){


pool1<-sample(Y[i,j,])
pool2<-sample(Y[i,h,])



Y.perm[i,j,]<-c(pool2[c(1:ni)],pool1[-c(1:ni)])
Y.perm[i,h,]<-c(pool1[c(1:ni)],pool2[-c(1:ni)])
}
}

if(ni==0){Y.perm<-Y}

T.B[cc]<-T.B[cc]+(sum(Y.perm[,j,])-sum(Y.perm[,h,]))^2


for(i in 1:(A-1)){
for(s in (i+1):A){

T.AB.b[cc]<-T.AB.b[cc] + (sum(Y.perm[i,j, ])-sum(Y.perm[s,j, ])-sum(Y.perm[i,h, ])+sum(Y.perm[s,h, ]))^2

}
}
}
}


}#end cc

T.A<-round(T.A,digits=8)
T.B<-round(T.B,digits=8)
T.AB.a<-round(T.AB.a,digits=8)
T.AB.b<-round(T.AB.b,digits=8)

T.AB<-apply(cbind(T.AB.a,T.AB.b),1,sum)

#pa<-round(sum(T.A[-1]>=T.A[1])/C,digits=log(C,10))
#pb<-round(sum(T.B[-1]>=T.B[1])/C,digits=log(C,10))

# pa<-round(sum(T.A>=T.A[1])/(C+1),
#           digits=log(C,10))
# pb<-round(sum(T.B>=T.B[1])/(C+1),
#           digits=log(C,10))
pab<-round(sum(T.AB>=T.AB[1])/(C+1),digits=log(C,10))

pab.a<-round(sum(T.AB.a>=T.AB.a[1])/(C+1),digits=log(C,10))
pab.b<-round(sum(T.AB.b>=T.AB.b[1])/(C+1),digits=log(C,10))



pa<-round(sum(T.A>=T.A[1])/(C+1),
          digits = 8)
pb<-round(sum(T.B>=T.B[1])/(C+1),
          digits = 8)
pab <- round(sum(T.AB >= T.AB[1])/(C+1),
             digits = 8)



######RESULTS

min.sig<-c(Pr.A[1],Pr.B[1])

return(list(pa=pa,pb=pb,pab=pab,pab.a=pab.a,pab.b=pab.b,TA=T.A[1],TB=T.B[1],TAB.a=T.AB.a[1],TAB.b=T.AB.b[1],type="Unconstrained",C=C,min.sig=min.sig,exact=FALSE))

#return(list(T.A=T.A,T.B=T.B,TAB.a=T.AB.a,TAB.b=T.AB.b))
}
