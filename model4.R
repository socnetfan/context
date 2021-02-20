
library(RSiena)

rm(list=ls())

# load network and behavior data
load("sch058.Rdata")
z <- nrow(friend.data.w1)

#############################################################################

# generate siena data sets
friendship <- sienaNet(
     array( c( friend.data.w1, friend.data.w2, friend.data.w3),
     dim = c( z, z, 3 ) ) )
smkebeh <- sienaNet( cbind(attr[,8],attr[,30],attr[,39]), type = "behavior" )
# coCovar: Constant actor covariates
# varCovar: Time-varying actor covariates
# coDyadCovar: Constant dyadic covariates
# varDyadCovar: Time-varying dyadic covariates
homesmke <- coCovar(attr$homesmke, centered=FALSE)
grade <- varCovar( cbind(attr$grade1, attr$grad2), centered=FALSE)
female <- coCovar( attr$female, centered=FALSE)
paredu <- coCovar( attr$paredu, centered=FALSE )
psupport <- coCovar( attr$psupport, centered=FALSE)
pmonitor <- coCovar (attr$pmonitor, centered=FALSE)
depression <- coCovar (attr$depression, centered=FALSE)
limitnomi <- varCovar( cbind(attr$limitnomi12,attr$limitnomi23), centered=FALSE)
rm(attr)

# load club data
load('sch058club.Rdata')

# sport teams
v <- cbind(attr[,21],attr[,25:35])
x <- matrix(rep(0,z*z),ncol=z)
for (o in 1:z)
{
  for (p in 1:z)
  {
    for (q in 1:ncol(v)) 
    {
      if (v[o,q]==v[p,q] & v[o,q]==1 & o!=p) {x[o,p] <- x[o,p]+1}
    }
  }
}
# number of sport teams each pairs of individuals share
coclubs1 <- coDyadCovar(x)
# number of sport teams
nclubs11 <- coCovar (rowSums(v), centered=FALSE)
# number of sport teams squared
nclubs12 <- coCovar (rowSums(v)^2, centered=FALSE)
# each pairs of individuals share any sport team
x[which(x>=1)] <- 1
coclubs12 <- coDyadCovar(x)
rm(v,x,o,p,q)

# drama club
v <- attr[,15]
x <- matrix(rep(0,z*z),ncol=z)
for (o in 1:z)
{
  for (p in 1:z)
  {
    if (v[o]==v[p] & v[o]==1 & o!=p) {x[o,p] <- 1}
  }
}
# each pairs of individuals share drama club membership
codrama <- coDyadCovar(x)
# each individual in a drama club
drama <- coCovar (v, centered=FALSE)
rm(v,x,o,p)

# band
v <- attr[,20]
x <- matrix(rep(0,z*z),ncol=z)
for (o in 1:z)
{
  for (p in 1:z)
  {
    if (v[o]==v[p] & v[o]==1 & o!=p) {x[o,p] <- 1}
  }
}
# each pairs of individuals share band membership
coband <- coDyadCovar(x)
# each individual in a band
band <- coCovar (v, centered=FALSE)
rm(v,x,o,p)

# chorus
v <- attr[,22]
x <- matrix(rep(0,z*z),ncol=z)
for (o in 1:z)
{
  for (p in 1:z)
  {
    if (v[o]==v[p] & v[o]==1 & o!=p) {x[o,p] <- 1}
  }
}
# each pairs of individuals share chorus membership
cochorus <- coDyadCovar(x)
# each individual in a chorus
chorus <- coCovar (v, centered=FALSE)
rm(v,x,o,p)

# orchestra
v <- attr[,23]
x <- matrix(rep(0,z*z),ncol=z)
for (o in 1:z)
{
  for (p in 1:z)
  {
    if (v[o]==v[p] & v[o]==1 & o!=p) {x[o,p] <- 1}
  }
}
# each pairs of individuals share orchestra  membership
coorch <- coDyadCovar(x)
# each individual in an orchestra
orch <- coCovar (v, centered=FALSE)
rm(v,x,o,p)

# academic clubs
v <- cbind(attr[,8:14],attr[,17:19],attr[,37:40])
x <- matrix(rep(0,z*z),ncol=z)
for (o in 1:z)
{
  for (p in 1:z)
  {
    for (q in 1:ncol(v)) 
    {
      if (v[o,q]==v[p,q] & v[o,q]==1 & o!=p) {x[o,p] <- x[o,p]+1}
    }
  }
}
# number of academic clubs each pairs of individuals share
coclubs3 <- coDyadCovar(x)
# number of academic clubs
nclubs31 <- coCovar (rowSums(v), centered=FALSE)
# number of academic clubs squared
nclubs32 <- coCovar (rowSums(v)^2, centered=FALSE)
# each pairs of individuals share any academic club
x[which(x>=1)] <- 1
coclubs32 <- coDyadCovar(x)
rm(v,x,o,p,q)

mymodel <- sienaModelCreate(projname = 'model4')
mydata <- sienaDataCreate(friendship,limitnomi,coclubs1,coclubs12,codrama,drama,coband,band,cochorus,chorus,coorch,orch,coclubs3,coclubs32,nclubs11,nclubs12,nclubs31,nclubs32,female,psupport,depression,paredu,homesmke,grade,pmonitor,smkebeh)

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,transTrip,transRecTrip,cycle3,inPop,inInAss)
myeff <- includeEffects(myeff, X, interaction1 = "coclubs1")
myeff <- includeEffects(myeff, X, interaction1 = "coclubs3")
myeff <- includeEffects(myeff, egoX, interaction1 = "psupport")
myeff <- includeEffects(myeff, egoX, interaction1 = "pmonitor")
myeff <- includeEffects(myeff, egoX, interaction1 = "homesmke")
myeff <- includeEffects(myeff, simX, interaction1 = "female")
myeff <- includeEffects(myeff, simX, interaction1 = "grade")
myeff <- includeEffects(myeff, simX, interaction1 = "paredu")
myeff <- includeEffects(myeff, egoX, interaction1 = "limitnomi")
myeff <- includeEffects(myeff, egoX, altX, simX, interaction1 = "nclubs11")
myeff <- includeEffects(myeff, egoX, altX, simX, interaction1 = "drama")
myeff <- includeEffects(myeff, egoX, altX, simX, interaction1 = "band")
myeff <- includeEffects(myeff, egoX, altX, simX, interaction1 = "chorus")
myeff <- includeEffects(myeff, egoX, altX, simX, interaction1 = "orch")
myeff <- includeEffects(myeff, egoX, altX, simX, interaction1 = "nclubs31")
myeff <- includeEffects(myeff, egoX, altX, simX, interaction1 = "smkebeh")
myeff <- includeInteraction( myeff, simX, egoX, interaction1 = c("smkebeh", "chorus"))
myeff <- includeEffects(myeff,name="smkebeh", avNegAbsDiffX, interaction1 = "friendship")
myeff <- includeInteraction(myeff, avNegAbsDiffX, effFrom, name="smkebeh", interaction1=c("friendship","orch"))
myeff <- includeEffects(myeff,name="smkebeh", avNegAbsDiff2m, interaction1 = "friendship", interaction2="coclubs12")
myeff <- includeEffects(myeff,name="smkebeh", avNegAbsDiff2m, interaction1 = "friendship", interaction2="codrama")
myeff <- includeEffects(myeff,name="smkebeh", avNegAbsDiff2m, interaction1 = "friendship", interaction2="coband")
myeff <- includeEffects(myeff,name="smkebeh", avNegAbsDiff2m, interaction1 = "friendship", interaction2="cochorus")
myeff <- includeEffects(myeff,name="smkebeh", avNegAbsDiff2mm, interaction1 = "friendship", interaction2="coorch")
myeff <- includeEffects(myeff,name="smkebeh", avNegAbsDiff2m, interaction1 = "friendship", interaction2="coclubs32")
myeff <- includeEffects(myeff,name="smkebeh",indeg, interaction1 = "friendship")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "nclubs11")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "nclubs12")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "drama")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "band")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "chorus")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "orch")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "nclubs31")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "nclubs32")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "psupport")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "homesmke")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "depression")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "pmonitor")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "female")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "grade")
myeff

ans <- siena07(mymodel, data=mydata, effects=myeff, returnDeps=TRUE, batch=T)
save(ans,file='model4.Rdata')
while (summary(ans)$tconv.max[1,1] >=0.25 | max(abs(ans$tconv))>=0.1) {
  ans <- siena07(mymodel, data=mydata, effects=myeff, prevAns=ans, returnDeps=TRUE, batch=T)
  save(ans,file='model4.Rdata')
}
rm(list=ls())


