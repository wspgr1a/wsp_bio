#### ANNA SPERNOL

source("http://bioconductor.org/biocLite.R")

#biocLite('affy')
#biocLite("ggplot2")
#biocLite('hgu95av2cdf')
#biocLite('nortest')
library('affy')
library('Biobase')
library("ggplot2")
library('hgu95av2cdf')
library('nortest')

setwd('D:/Magisterka/WSP')

load("D:/Magisterka/WSP/dane.RData")

dane1=data[,1:5]
dane2=data[,50:55]

ekspresja=selekcja(dane1, dane2, 0.3)
