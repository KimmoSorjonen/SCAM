
############################# BIENVENUE #############################

############# I SOLEMNLY SWEAR THAT I AM UP TO NO GOOD ##############

#################################
## Loading packages

library(lavaan)

#################################
## Data

#### 1. Xiong et al. (2026) ####

## X = Psychological resilience (PR)
## M1 = Perceived social support (PSS)
## M2 = Career decision-making self-efficacy (CSE)
## Y = Employability (EB)

sn1 <- "Xiong et al. (2026)"

n1 <- 1709 ## sample size

ind1 <- c(0.048,0.390,0.140,0.578,0.236,0.814) ## effects

low1 <- c(0.018,0.339,0.106,0.532,0.198,0.784) ## lower CI

upp1 <- c(0.079,0.439,0.175,0.622,0.273,0.844) ## upper CI

se1 <- (upp1-low1)/3.92 ## standard error

rm1 <- matrix(c( ## correlations
  
  1.00, 0.72, 0.79, 0.80,
  0.72, 1.00, 0.71, 0.69,
  0.79, 0.71, 1.00, 0.88,
  0.80, 0.69, 0.88, 1.00),nrow=4)

colnames(rm1) <- rownames(rm1) <- c("X","M1","M2","Y") ## order


#### 2. Wang et al. (2025) ####

## X = Physical activity (PA)
## M1 = Self-efficacy (SE)
## M2 = Emotion regulation (ER)
## Y = Mental health (MH)

sn2 <- "Wang et al. (2025)"

n2 <- 1395 ## sampe size

ind2 <- c(-0.0111,-0.0156,-0.0026,-0.0294,-0.0534,-0.0828) ## effects

se2 <- c(0.0041,0.0078,0.0012,0.0088,0.0254,0.0267) ## SE

rm2 <- matrix(c( ## correlations
  
   1.0000,  0.1069,  0.0616, -0.0828,
   0.1069,  1.0000,  0.0890, -0.1362,
   0.0616,  0.0890,  1.0000, -0.3089
  -0.0828, -0.1362, -0.3089,  1.0000),nrow=4)

colnames(rm2) <- rownames(rm2) <- c("X","M1","M2","Y") ## order


#### 3. Wu et al. (2025) ####

## X = Internet addiction (IA)
## M1 = Self-control (SC)
## M2 = Anxiety (AN)
## Y = Academic procrastination (AP)

sn3 <- "Wu et al. (2026)"

n3 <- 2421 ## sample size

ind3 <- c(0.241,0.018,0.010,0.269,0.224,0.493) ## effects

low3 <- c(0.213,0.009,0.005,0.240,0.186,0.458) ## lower CI

upp3 <- c(0.270,0.028,0.015,0.299,0.262,0.527) ## upper CI

se3 <- (upp3-low3)/3.92 ## SE

rm3 <- matrix(c( ## correlations
  
   1.00,  0.48, -0.58,  0.33,
   0.48,  1.00, -0.38,  0.38,
  -0.58, -0.38,  1.00, -0.55,
   0.33,  0.38, -0.55,  1.00),nrow=4)

colnames(rm3) <- rownames(rm3) <- c("X","Y","M1","M2") ## order


#### 4. Yue et al. (2026) ####

## X = Physical activity (PA)
## M1 = Exercise self-efficacy (ESE)
## M2 = Psychological resilience (PR)
## Y = Negative emotions (NE)

sn4 <- "Yue et al. (2026)"

n4 <- 877 ## sample size

ind4 <- c(-0.061,-0.022,-0.008,-0.092,-0.192,-0.284) ## effects

se4 <- c(0.012,0.007,0.003,0.014,0.031,0.029) ## SE

rm4 <- matrix(c( ## correlations
  
   1.000,  0.415, 0.333, -0.382,
   0.415,  1.000, 0.356, -0.409,
   0.333,  0.356, 1.000, -0.319,
  -0.382, -0.409, -0.319, 1.000),nrow=4)

colnames(rm4) <- rownames(rm4) <- c("X","M1","M2","Y") ## order

#### lists of studies, sample size, effects, SE, and correlations

studyv <- c(sn1,sn2,sn3,sn4)
nv <- c(n1,n2,n3,n4)
indlist <- list(ind1,ind2,ind3,ind4)
selist <- list(se1,se2,se3,se4)
rmlist <- list(rm1,rm2,rm3,rm4)

#################################
## Mediated effects

## single shared cause model

ccmod <- "

CC =~ X+M1+M2+Y

## (Error) variances

CC ~~ 1*CC

X ~~ X
M1 ~~ M1
M2 ~~ M2
Y ~~ Y

"

## Mediation model, corresponding to model 
## in scrutinized studies 

medmod <- "

## Regression effects

M1 ~ a1*X
M2 ~ a2*X+mm*M1
Y ~ c*X+b1*M1+b2*M2

## (Error) variances

X ~~ X
M1 ~~ M1
M2 ~~ M2
Y ~~ Y

## Mediated, direct, 
## and total effects

ind1 := a1*b1
ind2 := a2*b2
ind3 := a1*mm*b2
totind := ind1+ind2+ind3
dir := c
tot := totind+dir

"

#################################
## Table results

restab <- vector() ## empty vector, filled below

for(i in 1:length(nv)){ ## for each study
  
  ######
  
  n <- nv[i] ## picking sample size from list
  rm <- rmlist[i] ## picking correlations
  ind.o <- indlist[[i]] ## picking effects
  se.o <- selist[[i]] ## picking SE
  study <- studyv[i] ## picking study name
  
  ######
  
  fit.cc <- lavaan(ccmod, sample.cov=rm, ## fitting cc-model to data 
                   sample.nobs=n, meanstructure=F)
  
  rm.pred <- lavInspect(fit.cc, what="cor.ov") ## implied correlations
  
  fit.med <- lavaan(medmod, sample.cov=rm.pred, ## fitting mediation model 
                    sample.nobs=n, meanstructure=F)
  
  pares <- parameterestimates(fit.med) ## parameter estimates
  
  ## Predicted indirect effects
  ind.p <- pares[pares$op==":=","est"] ## vector of effects
  se.p <- pares[pares$op==":=","se"] ## vector of SE
  
  #### Fisher transformations
  
  b.d.fish <- atanh(ind.o)-atanh(ind.p) ## difference
  se.d.fish <- sqrt(atanh(se.o)^2 + atanh(se.p)^2) ## SE of difference
  
  z.d <- b.d.fish/se.d.fish ## Z-value of difference
  p.d <- 2 * pnorm(-abs(z.d)) ## p-value of difference
  
  low.d.fish <- b.d.fish-1.96*se.d.fish ## lower CI
  upp.d.fish <- b.d.fish+1.96*se.d.fish ## upper CI
  
  #### Reverse transformations
  
  b.d <- tanh(b.d.fish)
  low.d <- tanh(low.d.fish)
  upp.d <- tanh(upp.d.fish)
  
  #### Results
  
  res <- data.frame(study = rep(study,length(b.d)),
                    indeff = 1:6,
                    ind.o=ind.o,
                    se.o=se.o,
                    ind.p=ind.p,
                    se.p=se.p,
                    b.d=b.d,
                    low.d=low.d,
                    upp.d=upp.d,
                    p.d=p.d)
  
  restab <- rbind(restab,res) ## adding results to table
}

restabnew <- cbind(restab[,1],round(restab[,2:10],3)) ## rounding

restab$low.o <- restab$ind.o-1.96*restab$se.o ## lower CI, observed
restab$upp.o <- restab$ind.o+1.96*restab$se.o ## upper CI, observed
restab$low.p <- restab$ind.p-1.96*restab$se.p ## lower CI, expected
restab$upp.p <- restab$ind.p+1.96*restab$se.p ## upper CI, expected

pick <- c("ind.o","low.o","upp.o", ## to be used
          "ind.p","low.p","upp.p","p.d")

res1 <- restab[1:6,pick] ## study 1
res2 <- restab[7:12,pick] ## study 2
res3 <- restab[13:18,pick] ## study 3
res4 <- restab[19:24,pick] ## study 4

reslist <- list(res1,res2,res3,res4) ## all in a list

#################################
## Figure 4

indeff <- rev(c("Ind.M1","Ind.M2","Ind.M1.M2", ## y-labels
                "Tot. Ind","Direct","Total"))

sf <- 0.8 ## sizing factor

if(dev.cur()==2) dev.off() ## removing earlier plots
par(mar=c(0,1,1.2,0), oma=c(4,4,0.5,1), mfrow=c(1,4)) ## margins and layout

for(j in 1:4){ ## for each study
  
  res <- reslist[[j]] ## picking results
  ind.o <- res[,1:3] ## picking observed effect
  ind.p <- res[,4:6] ## picking expected effect
  pvv <- res$p.d ## picking p-values
  
  mi <- min(c(res$low.o,res$low.p)) ## lowest value
  ma <- max(c(res$upp.o,res$upp.p)) ## highest value
  
  plot(c(mi,ma),c(0.5,6.5), type="n",yaxt="n", xlab="",
       ylab="", cex.axis=sf) ## empty plot
  
  mtext(studyv[j], side = 3, line=0.2, outer=F, cex=sf) ## study name on top
  if(j==2) mtext("Standardized Effect", side = 1, line=2.5, outer=T, cex=sf)
  axis(2,1:6, labels=F,las=2, cex.axis=sf, cex.lab=0.09) ## Y-axis
  if(j==1) axis(2,1:6, labels=indeff,las=2, cex.axis=sf, cex.lab=0.09) ## Y-axis
  abline(h=6.5, lwd=0.5, col="dark gray") ## gray line
  
  for(i in 1:6){ ## for each effect
    
    yvo <- 7-i+0.2 ## y-value for reported effect
    yvp <- 7-i-0.2 ## y-value for predicted effect
    
    pc <- ifelse(pvv[i]<0.05,"black","white") ## color of marker
    
    ## Confidence intervals
    arrows(ind.o[i,2],yvo,ind.o[i,3],yvo,angle=90,length=0.03,code=3,lwd=2)
    arrows(ind.p[i,2],yvp,ind.p[i,3],yvp,angle=90,length=0.03,code=3,lwd=2)
    
    ## Effects
    points(ind.o[i,1],yvo,pch=21,col="black",bg=pc, cex=1.5)
    points(ind.p[i,1],yvp,pch=21,col="black",bg="gray", cex=1.5)
    
    ## Gray line
    abline(h=7-i-0.5, lwd=0.5, col="dark gray")
  }
}

difftab <- restabnew[,c("b.d","low.d","upp.d","p.d")] ## table

cbind(difftab[1:6,],difftab[7:12,], ## reorganized table
      difftab[13:18,],difftab[19:24,])


########################## MISCHIEF MANAGED #########################

############################# AU REVOIR #############################

