# dodatkowe biblioteki
biocLite("MASS")
library(MASS)
biocLite("sortinghat")
library(sortinghat)

#------------------------------------------------------------------------------------
# klasyfikacja z ocena bledu klasyfikatora bootstrap 632
mojeda=exprs(expset)
mojeda=t(mojeda)
mojeda=as.data.frame(mojeda)
mojeda=mojeda[,1:100] # sprawdzanie dzia³ania na mniejszej ilosci danych, dla wszystkich danych zakomentowaæ
pomocnicza=expset@phenoData@data$CLASS

moje_x <- mojeda
moje_y <- pomocnicza


lda_wrapper <- function(object, newdata) { predict(object, newdata)$class }
set.seed(42)
apparent <- errorest_apparent(x = moje_x, y = moje_y, train = MASS:::lda,
                              classify = lda_wrapper)
set.seed(42)
loo_boot <- errorest_loo_boot(x = moje_x, y = moje_y, train = MASS:::lda,
                              classify = lda_wrapper)


errorest_632(x = moje_x, y = moje_y, train = MASS:::lda, classify = lda_wrapper,
             apparent = apparent, loo_boot = loo_boot, num_bootstraps = 50)

#-------------------------------------------------------------------------------------


























