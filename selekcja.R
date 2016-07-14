##### ANNA SPERNOL


####### SELEKCJA GENOW ########

#sprawdzic czy wariancje sa jednorodne i czy rozklad jest normalny. Jesli tak o zastosowac test t-Studenta, jesli nie test Wilcoxona. 
#Zastosowac jeden rodzaj testu dla wszystkich danych, zeby miec jak porownywac wyniki
#Wybrac geny do dalszej analizy, takie ktorych p-wartosc jest nizsza od p-wartocci progowej

# dane1- dane ekspresji grupy 1; dane2- dane ekspresji grupy 2, p_prog - wartosc progowa dla p-wartosci
selekcja <- function(dane1, dane2, p_prog) { 
  
###### WYKRESY ######
qqnorm(dane1, main = "Normal Q-Q Plot grupa 1")
hist(dane1,breaks=100, main='Dane wejsciowe grupa 1')
qqnorm(dane2, main = "Normal Q-Q Plot grupa 2")
hist(dane2,breaks=100, main='Dane wejsciowe grupa 2')
  
  
###### NORMALNOSC ROZKLADOW (Test Test W Shapiro-Wilka) ######
#czy rozklad normalny:  H0(p-value>0.5) - rozklad normalny;  Ha-rozklad rozni sie od normalnego  
# pierwsza grupa danych
  p_valueSh1=matrix(data=NA,nrow=nrow(dane1),ncol=1)
  for (i in 1:nrow(dane1)) {
    test1=shapiro.test(dane1[i,]) 
# p wartosc
    p_valueSh1[i,]=test1$p.value  
  }
# procent danych z rokladem normalnym
  norm_proc1=(length(which(p_valueSh1[]>=0.05))*100)/nrow(p_valueSh1)
  
# druga grupa danych
  p_valueSh2=matrix(data=NA,nrow=nrow(dane2),ncol=1)
  for (i in 1:nrow(dane2)) {
    test2=shapiro.test(dane2[i,])
# p-wartosci 
    p_valueSh2[i,]=test2$p.value 
  }
# procent danych z rokladem normalnym
  norm_proc2=(length(which(p_valueSh2[]>=0.05))*100)/nrow(p_valueSh2)
  

###### JEDNORODNOSC WARIANCJI (Test F) ######
# czy wariancje w grupach sa jednorodne H0- sa jednorodne; Ha- nie sa jednorodne
  macierz=matrix(data=NA,nrow=nrow(dane1),ncol=1) 
  macierz_p_value=matrix(data=NA,nrow=nrow(dane1),ncol=1) 
  for (i in 1:nrow(dane1)) {
    test3=var.test(dane1[i,],dane2[i,])
# p-wartosc
    macierz_p_value[i]=test3$p.value 
    if ((test3$statistic>=test3$conf.int[1])&&(test3$statistic<=test3$conf.int[2]))macierz[i,1]=1
    else macierz[i,1]=0
  }
# przepisanie nazw wierszy  
  rownames(macierz_p_value)=rownames(dane1)
  rownames(macierz)=rownames(dane1)
  
# procent danych z jednorodna wariancja
  var_proc=(length(which(macierz[]==1))*100)/nrow(macierz)
 
###### TESTY ######
## t-Student
# H0: srednie w grupach sa takie same; Ha: srednie w grupach nie sa takie same

  if ((norm_proc1>=0.5)&& (norm_proc2>=50)){ 
    if (var_proc>50) {
      var.equal=T
      print("Wykonano test t-Studenta")
    } else {
      var.equal=F
      print("Wykonano test t-Studenta dla danych z niejednorodna wariancja") }
    
    pval_test=matrix(data=NA,nrow=nrow(dane1),ncol=1)
    for (i in 1:nrow(dane1)) {
      test4=t.test(dane1[i,],dane2[i,],var.equal=var.equal)
# p-warto?ci
      pval_test[i]=test4$p.value 
    }
    
## Mann-Whitney-Wilcoxon 
# H0: mediany w grupach sa takie same, Ha: mediany sa rozne 

  } else { 
    print('Wykonano test Manna-Whitneya-Wilcoxona')
    pval_test=matrix(data=NA,nrow=nrow(dane1),ncol=1)
    for (i in 1:nrow(dane1)) {
      test5=wilcox.test(dane1[i,],dane2[i,])
# p-wartosci      
      pval_test[i]=test5$p.value 
    }
  }
# przepisanie nazw wierszy  
  rownames(pval_test)=rownames(dane1)
  
# roznica wartosci srednich
mean1=rowMeans(dane1)
mean2=rowMeans(dane2)
rozn=mean2-mean1
zmiany=cbind(rozn)
row.names(zmiany)=rownames(dane1)
  
zmiany=cbind(zmiany,pval_test)

###### GENY ROZNICUJACE ######

#geny z p-wartosci ponizej progu
geny=which(pval_test<p_prog) 
zmiany2=zmiany[geny[],]

# geny o zmianie ekspresji mniejszej/wiekszej niz 0.5
geny2=which(abs(zmiany2[,1])>=0.5)

zmiany3=zmiany2[geny2[],] 
# indeksy wybranych genow
geny3=geny[geny2]

###### EKPRESJA WYBRANYCH GENOW ######
exp_1=dane1[geny3[],] 
exp_2=dane2[geny3[],] 


###### VOLCANO PLOT ######

lista_genow=data.frame(rozn,pval_test)
lista_genow$threshold = as.factor((abs(lista_genow$rozn)>(0.5)) & (lista_genow$pval_test < p_prog))

volcano <- ggplot(data=lista_genow, aes(x=rozn, y=-log10(pval_test),colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-2, 2)) + ylim(c(0, 5)) +
  xlab("zmiany ekspresji") + ylab("-log10 p-wartosc")+
  ggtitle("Volcano plot koncowych danych") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"))
print(volcano)

###### WYBRANE GENY ######
# wybrane geny
wynik=cbind(zmiany3[,1], zmiany3[,2])

# utworzenie listy genow w pliku
write(row.names(wynik),file="Wybrane geny.txt")

result<-list("expresja1"=exp_1,"expresja2"=exp_2)

return (result)

}
    