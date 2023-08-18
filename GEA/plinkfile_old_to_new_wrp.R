# Put ONLY the chromosome number in the first column of the .map file

library(tidyverse)

setwd("/0_data")

# File name to rename
filename <- 'plinkfile.map'

# Check if file exists.
if (file.exists(filename)) {
  # Rename file name
  file.rename(filename,'plinkfile_old.map')
}else{
  print('File Not found :')
}

old_map <- read.delim("plinkfile_old.map", sep = "", header = F)

old_map$V1 <- substring(old_map$V2, 3, 3)

write.table(old_map, 'plinkfile.map', sep = " ", row.names = F, col.names = F, quote = F)
