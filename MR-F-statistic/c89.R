
##################start
### need data
dat_es2o <- read.csv("dat_eL2o-89.csv", header=TRUE, sep=",",stringsAsFactors=F, colClasses = NA)
exp_es2o <- read.csv("exp_eL2o-89.csv", header=TRUE, sep=",",stringsAsFactors=F, colClasses = NA)



###MR-PRESSO
press01 <-run_mr_presso(dat_es2o, NbDistribution = 5000, SignifThreshold = 0.05)
press01 <- as.data.frame(press01)

######################################################################################
#8. Regression dilution I2 - wellbeing
######################################################################################
I2<-c()
# I-squared function
Isq = function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}
#Read in expsoure data

wellbeing <- exp_es2o
#Rename required columns
wellbeing$BetaXG<-wellbeing$beta.exposure
wellbeing$seBetaXG<-wellbeing$se.exposure
BetaXG   = wellbeing$BetaXG
seBetaXG = wellbeing$seBetaXG 

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)

#unweightedIsq  I^2
unIsq = Isq(BXG,seBetaXG) #unweighted

#Save results
output<-cbind("89", mF, unIsq)
I2<-rbind(I2, output)
colnames(I2) <- c("Exposure", "mF", "Isq Unweighted")

c89 <- cbind(I2,press01)
function_savecsv('c89')

#######################################################################

