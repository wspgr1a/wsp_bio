####### KLASTROWANIE - HIERARCHICZNE, DIANA, PAM - AGATA ¯ERTKA #############
#Do klastrowania wybiera siê zazwyczaj grupê genów z najwiêksz¹ wariancj¹ - top100 .
#Skraca to czas obliczeñ i redukuje szum w danych. 
#Ja przeprowadzê klasteryzacjê dla ca³ej macierzy z Expressionset, je¿eli bêdzie osoba, która
#robi selekcjê i wyselekcjonuje te 100 genów, mo¿na je tutaj w³o¿yæ, wystarczy podmieniæ zmienn¹ 
#mat, reszta bêdzie dzia³aæ.

mat <- exprs(expset)

##Klastrowanie - organizowanie danych które s¹ „blisko” w pewne grupy czyli klastry.

### KLASTROWANIE HIERARCHICZNE

## BOTTOM-UP (hclust) - podejœcie aglomeracyjne
#Zasada dzia³ania opiera siê na tym, ¿e znajdujemy 2 najbli¿ej po³o¿one punkty, ³¹czymy je ze sob¹
#w "super punkt" nastêpnie szukamy kolejnych najbli¿szych punktów (traktuj¹c ju¿
#po³¹czone jako jeden "super punkt")

d.cor=as.dist(1-cor(mat[,],method='spearman')) # korzystam ze Spearmana, bo mog¹ byæ wartoœci odstaj¹ce
h=hclust(d.cor, method = "average", members=NULL) #parametry: definicja odleg³oœci punktów, definicja ³¹czenia
#wyniki reprezentujemy w formie dendogramu
plot(as.dendrogram(h),
     main="Dendrogram otrzymany metod¹ 'hclust'")

## TOP-DOWN (diana) - DIvisive ANAlysis Clustering
#Podejœcie "top down" - wszystkie obserwacje startuj¹ z jednego klastra, a podzia³y s¹ wykonywane
#rekurencyjnie jako jeden ruch w dó³ hierarchi. Klastry s¹ dzielone tak d³ugo, a¿ ka¿dy klaster 
#zawiera tylko 1 obserwacjê. Klastry dzielone s¹ na zasadzie szukania najwiêkszej odmiennoœci w 
#pojedynczym duzym klastrze. W kolejnych krokach algorytm przypisuje obserwacje które s¹ "bli¿ej"
#od³amu ni¿ "starych" wartoœci. W rezultacie 1 du¿y klaster dzieli siê na 2 mniejsze.

diana=diana(t(d.cor))
plot(diana,
     main="Dendrogram otrzymany metod¹ 'diana'")
#zakomentowany plot, w którym wyœwietli nam siê kilka rysunków, mo¿na je wybieraæ wpisuj¹c cyfry, 
#aby wyjœæ, wciskamy 0, (dodany parametr ask=T)
#plot(diana,
#     main="Dendrogram otrzymany metod¹ 'diana'",ask=T)


### PAM - Partitioning around medoids

#Opis dzia³ania algorytmu: 
#Wybierz losowo k punktów jako medoidy. Medoidy - punkty centralne klastrów.
#Po³¹cz ka¿dy punkt z najbli¿szym punktem centralnym u¿ywaj¹c miary podobieñstwa
#Dla ka¿dego punktu centralnego zamieñ go z ka¿dym innym punktem nie bêd¹cym medoidem i 
#wyznacz koszt takiej konfiguracji
#Wybierz konfiguracjê z najni¿szym kosztem
#Powtarzaj kroki od po³¹czenia z najbli¿szym punktem centralnym do wybierania konfiguracji,
#dopóki nie ma ¿adnych zmian w medoidach


## WYBÓR LICZBY KLASTRÓW METOD¥ SILHOUETTE - podejœcie heurystyczne
#Maksymalna œrednia szerokoœæ silhouette dla wszystkich obserwacji
#mo¿e byæ u¿yta do wyboru liczby klastrów
#Dla ka¿dej obserwacji i, definiuje siê szerokoœæ silhouette s(i)
#za pomoc¹ algorytmu, Obserwacje z s(i) bliskim 1 mog¹ byæ rozwa¿ane jako dobrze
#klastrowane, a s(i)<0 jako z³e.
#
#Wybór liczby klastrów metod¹ Silhouette
asw=numeric(10)
for (k in 2:10) # testowanie liczby klastrów od 2 do 10
  asw[[k]]=pam(t(d.cor), k)$silinfo$avg.width # œrednia szerokoœæ silhouette
k.best=which.max(asw)+1 # liczba klastrów z najwiêksz¹ œredni¹ szerokoœcia silhouette
cat("Nalepsza liczba klastrów dla danych:", k.best, "\n")

## Klastrowanie PAM
pam=pam(t(d.cor), k.best)
plot(pam,
     main='Wyniki klastrowania metod¹ PAM.') # pam$clustering zawiera nazwy macierzy z numerem odpowiadaj¹cym klastrowi, do którego dany pacjent zosta³ przyporz¹dkowany

#zakomentowany plot, w którym wyœwietli nam siê kilka rysunków, mo¿na je wybieraæ wpisuj¹c cyfry, 
#aby wyjœæ, wciskamy 0
#plot(pam,ask=T,
#     main='Wyniki klastrowania metod¹ PAM.') # pam$clustering zawiera nazwy macierzy z numerem odpowiadaj¹cym klastrowi, do którego dany pacjent zosta³ przyporz¹dkowany
