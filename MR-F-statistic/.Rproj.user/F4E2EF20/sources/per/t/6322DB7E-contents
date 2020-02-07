library("TwoSampleMR")
es2o <- extract_instruments(outcomes=c(1055,1001))#uili

exp_es2o <- clump_data(es2o)

out_es2o <- extract_outcome_data(
  snps = exp_es2o$SNP,
  outcomes = c(297), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3,access_token=NULL
) # 允许回文序列？

dat_es2o <- harmonise_data(
  exposure_dat = exp_es2o, 
  outcome_dat = out_es2o, action = 2
) # dat2, action = 2尝试对其回文序列？

res_es2o <-mr(dat_es2o, method_list=c("mr_ivw","mr_egger_regression", "mr_weighted_median", "mr_weighted_mode") )

##################start
### need data
dat_es2o <- read.csv("dat_eL2o-1001.csv", header=TRUE, sep=",",stringsAsFactors=F, colClasses = NA)
exp_es2o <- read.csv("exp_eL2o-1001.csv", header=TRUE, sep=",",stringsAsFactors=F, colClasses = NA)



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
urate <- exp_es2o
wellbeing <- urate
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
output<-cbind("1001", mF, unIsq)
I2<-rbind(I2, output)
colnames(I2) <- c("Exposure", "mF", "Isq Unweighted")

c1001 <- cbind(I2,press01)
function_savecsv('c1001')

#######################################################################
