
library(RSiena)

rm(list=ls())

# load network and behavior data
load("sch058.Rdata")
z <- nrow(friend.data.w1)

# generate siena data sets
friendship <- sienaNet(
     array( c( friend.data.w1, friend.data.w2, friend.data.w3),
     dim = c( z, z, 3 ) ) )
smkebeh <- sienaNet( cbind(attr$smke1,attr$smke2,attr$smke3), type = "behavior" )
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
v <- cbind(club[,8:15],club[,17:23],club[,25:35],club[,37:40])
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
# number of clubs each pairs of individuals share
coclubs <- coDyadCovar(x)
# number of clubs
nclubs <- coCovar (rowSums(v), centered=FALSE)
# number of clubs squared
nclubs2 <- coCovar (rowSums(v)^2, centered=FALSE)
# each pairs of individuals share any club
x[which(x>=1)] <- 1
coclubs2 <- coDyadCovar(x)
rm(v,x,o,p,q)

mymodel <- sienaAlgorithmCreate(projname = 'model1', behModelType = NULL)
mydata <- sienaDataCreate(friendship,limitnomi,coclubs,coclubs2,nclubs,nclubs2,female,psupport,depression,paredu,homesmke,grade,pmonitor,smkebeh)

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,transTrip,transRecTrip,cycle3,inPop,inInAss)
myeff <- includeEffects(myeff, X, interaction1 = "coclubs")
myeff <- includeEffects(myeff, egoX, interaction1 = "psupport")
myeff <- includeEffects(myeff, egoX, interaction1 = "pmonitor")
myeff <- includeEffects(myeff, egoX, interaction1 = "homesmke")
myeff <- includeEffects(myeff, simX, interaction1 = "female")
myeff <- includeEffects(myeff, simX, interaction1 = "grade")
myeff <- includeEffects(myeff, simX, interaction1 = "paredu")
myeff <- includeEffects(myeff, egoX, interaction1 = "limitnomi")
myeff <- includeEffects(myeff, egoX, altX, simX, interaction1 = "nclubs")
myeff <- includeEffects(myeff, egoX, altX, simX, interaction1 = "smkebeh")
myeff <- includeEffects(myeff,name="smkebeh", avNegAbsDiffX, interaction1 = "friendship")
myeff <- includeEffects(myeff,name="smkebeh", avNegAbsDiff2m, interaction1 = "friendship", interaction2="coclubs2")
myeff <- includeEffects(myeff,name="smkebeh",indeg, interaction1 = "friendship")
myeff <- includeEffects(myeff,name="smkebeh",outdeg, interaction1 = "friendship")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "nclubs")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "nclubs2")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "psupport")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "homesmke")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "depression")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "pmonitor")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "female")
myeff <- includeEffects(myeff,name="smkebeh",effFrom, interaction1 = "grade")
myeff

ans <- siena07(mymodel, data=mydata, effects=myeff, returnDeps=TRUE, batch=T)
save(ans,file='model1.Rdata')
while (summary(ans)$tconv.max[1,1] >=0.25 | max(abs(ans$tconv))>=0.1) {
  ans <- siena07(mymodel, data=mydata, effects=myeff, prevAns=ans, returnDeps=TRUE, batch=T)
  save(ans,file='model1.Rdata')
}
rm(list=ls())


