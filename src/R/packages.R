list.of.packages <-c("faraway","R2OpenBUGS","modeest")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.packages
if(length(new.packages)) install.packages(new.packages)
