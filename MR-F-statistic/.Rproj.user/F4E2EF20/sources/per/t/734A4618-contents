######################################################
#
######################################################
# save csv data. x <- 'xxx'
function_savecsv <- function(x) {
  mpg <- x
  #mpg <- as.character(x)
  # 去掉引号 eval(parse(text=mpg))
  
  csvname <- paste(mpg, ".csv", sep = "")
  write.csv(eval(parse(text=mpg)), file = csvname)
  
}

function_savecsv('gwas_catalog')



# read
mpg <- read.csv("mpg.csv", header=TRUE, sep=",",stringsAsFactors=F, colClasses = NA)
######################################################
