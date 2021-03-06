####### KLASTROWANIE - HIERARCHICZNE, DIANA, PAM - AGATA �ERTKA #############
#Do klastrowania wybiera si� zazwyczaj grup� gen�w z najwi�ksz� wariancj� - top100 .
#Skraca to czas oblicze� i redukuje szum w danych. 
#Ja przeprowadz� klasteryzacj� dla ca�ej macierzy z Expressionset, je�eli b�dzie osoba, kt�ra
#robi selekcj� i wyselekcjonuje te 100 gen�w, mo�na je tutaj w�o�y�, wystarczy podmieni� zmienn� 
#mat, reszta b�dzie dzia�a�.

mat <- exprs(expset)

##Klastrowanie - organizowanie danych kt�re s� �blisko� w pewne grupy czyli klastry.

### KLASTROWANIE HIERARCHICZNE

## BOTTOM-UP (hclust) - podej�cie aglomeracyjne
#Zasada dzia�ania opiera si� na tym, �e znajdujemy 2 najbli�ej po�o�one punkty, ��czymy je ze sob�
#w "super punkt" nast�pnie szukamy kolejnych najbli�szych punkt�w (traktuj�c ju�
#po��czone jako jeden "super punkt")

d.cor=as.dist(1-cor(mat[,],method='spearman')) # korzystam ze Spearmana, bo mog� by� warto�ci odstaj�ce
h=hclust(d.cor, method = "average", members=NULL) #parametry: definicja odleg�o�ci punkt�w, definicja ��czenia
#wyniki reprezentujemy w formie dendogramu
plot(as.dendrogram(h),
     main="Dendrogram otrzymany metod� 'hclust'")

## TOP-DOWN (diana) - DIvisive ANAlysis Clustering
#Podej�cie "top down" - wszystkie obserwacje startuj� z jednego klastra, a podzia�y s� wykonywane
#rekurencyjnie jako jeden ruch w d� hierarchi. Klastry s� dzielone tak d�ugo, a� ka�dy klaster 
#zawiera tylko 1 obserwacj�. Klastry dzielone s� na zasadzie szukania najwi�kszej odmienno�ci w 
#pojedynczym duzym klastrze. W kolejnych krokach algorytm przypisuje obserwacje kt�re s� "bli�ej"
#od�amu ni� "starych" warto�ci. W rezultacie 1 du�y klaster dzieli si� na 2 mniejsze.

diana=diana(t(d.cor))
plot(diana,
     main="Dendrogram otrzymany metod� 'diana'")
#zakomentowany plot, w kt�rym wy�wietli nam si� kilka rysunk�w, mo�na je wybiera� wpisuj�c cyfry, 
#aby wyj��, wciskamy 0, (dodany parametr ask=T)
#plot(diana,
#     main="Dendrogram otrzymany metod� 'diana'",ask=T)


### PAM - Partitioning around medoids

#Opis dzia�ania algorytmu: 
#Wybierz losowo k punkt�w jako medoidy. Medoidy - punkty centralne klastr�w.
#Po��cz ka�dy punkt z najbli�szym punktem centralnym u�ywaj�c miary podobie�stwa
#Dla ka�dego punktu centralnego zamie� go z ka�dym innym punktem nie b�d�cym medoidem i 
#wyznacz koszt takiej konfiguracji
#Wybierz konfiguracj� z najni�szym kosztem
#Powtarzaj kroki od po��czenia z najbli�szym punktem centralnym do wybierania konfiguracji,
#dop�ki nie ma �adnych zmian w medoidach


## WYB�R LICZBY KLASTR�W METOD� SILHOUETTE - podej�cie heurystyczne
#Maksymalna �rednia szeroko�� silhouette dla wszystkich obserwacji
#mo�e by� u�yta do wyboru liczby klastr�w
#Dla ka�dej obserwacji i, definiuje si� szeroko�� silhouette s(i)
#za pomoc� algorytmu, Obserwacje z s(i) bliskim 1 mog� by� rozwa�ane jako dobrze
#klastrowane, a s(i)<0 jako z�e.
#
#Wyb�r liczby klastr�w metod� Silhouette
asw=numeric(10)
for (k in 2:10) # testowanie liczby klastr�w od 2 do 10
  asw[[k]]=pam(t(d.cor), k)$silinfo$avg.width # �rednia szeroko�� silhouette
k.best=which.max(asw)+1 # liczba klastr�w z najwi�ksz� �redni� szeroko�cia silhouette
cat("Nalepsza liczba klastr�w dla danych:", k.best, "\n")

## Klastrowanie PAM
pam=pam(t(d.cor), k.best)
plot(pam,
     main='Wyniki klastrowania metod� PAM.') # pam$clustering zawiera nazwy macierzy z numerem odpowiadaj�cym klastrowi, do kt�rego dany pacjent zosta� przyporz�dkowany

#zakomentowany plot, w kt�rym wy�wietli nam si� kilka rysunk�w, mo�na je wybiera� wpisuj�c cyfry, 
#aby wyj��, wciskamy 0
#plot(pam,ask=T,
#     main='Wyniki klastrowania metod� PAM.') # pam$clustering zawiera nazwy macierzy z numerem odpowiadaj�cym klastrowi, do kt�rego dany pacjent zosta� przyporz�dkowany
