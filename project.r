# Treść: Niech X_1,...X_n będzie próbą prostą z rozkładu potęgowego o gęstości f(x)=a/(x^{a+1})I(x>1) (I to jest indykator), 
# gdzie a>0 to nieznany parametr. Skonstruować test o rozmiarze Alpha =0.05 najmocniejszy z testów nieobciążonych hipotezy 
# H0: a=2, przeciwko alternatywie a ≠ 2. Narysować wykres mocy empirycznej dla n =40 w zależności od a na przedziale [1,3].

####################################################################################################
# Wyznaczenie c_a i c_b metoda bisekcji z układu rownan
# 1: F_20(c_b) - F_20(c_a) = 1 - alpha
# 2: F_22(c_b) - F_22(c_a) = 1 - alpha

alpha <- 0.05
fun <-function(c_a){
  pchisq(qchisq(pchisq(c_a, df = 20) + 1 - alpha, df = 20), df = 22) - pchisq(c_a, df = 22) - 1 + alpha
}

curve(fun, xlim=c(9,11), col='blue', lwd=1.5, lty=2)
abline(h=0)
abline(v=9)
# miejsce zerowe w okolicy 10

# gotowiec
#library(NLRoot)
#BFfzero(fun, 9, 11)

# funkcja bisekcji
bisection <- function(f, a, b, n = 1000, tol = 1e-15) { # argument tol - tolerancja (dokladnosc przyblizenia)
  
  # Znaki f(a) i f(b) musza sie roznic
  if(is.finite(f(a)) && is.finite(f(b)) && !(f(a)*f(b) < 0)) {
    stop('Znaki f(a) i f(b) musza sie roznic')
  }
  
  for (i in 1:n) {
    c <- (a + b) / 2 # Obliczamy srodek przedzialu
    
    # Jeśli wartosc funkcji w punkcie srodkowym wynosi 0 lub środek jest poniżej tolerancji zwracamy pierwiastek
    if ((f(c) == 0) || ((b - a) / 2) < tol) {
      return(c)
    }
    ifelse(sign(f(c)) == sign(f(a)), a <- c, b <- c) # przypisujemy srodek jako lewy lub prawy koniec przedzialu
    # sign() to funkcja signum
  }
  
  # Jesli osiagnieto maksymalna ilosc iteracji
  print('Przekroczono limit iteracji')
}

c_a <- bisection(fun, 9, 11)
c_a
c_b <- qchisq(pchisq(c_a, df = 20) + 1 - alpha, df = 20)
c_b
pchisq(c_a, df = 20) - pchisq(c_b, df = 20) + 1
pchisq(c_a, df = 22) - pchisq(c_b, df = 22) + 1

####################################################################################################

library(VGAM)
# ?ppareto

m <- 100000
k <- 400
n <- 40 # ilość obserwacji
alpha <- 0.05

a <- c() # wartosci z przedzialu 1,3
T <- c() # statystyki testowe
M <- c() # wektor wyników testu
Moc <- c() # wektor mocy testu

for (j in 1:k) { # krzywa empirycznej mocy
  a[j] <- j*0.005+1 # liczymy pkty siatki
  # dla takiego pktu siatki a[j] liczymy wartosc testu i sprawdzamy ile razy test odrzucil hipoteze
  for (i in 1:m) {
    X <- rpareto(n, shape = a[j], scale = 1) # generujemy próbkę z rozkładu potęgowego (Pareta)
    T[i] <- sum(log(X)) # wartosc statystyki testowej
    
    # liczymy czy test odrzucil hipoteze czy nie
    M[i]<-1
    if(T[i] > c_a & T[i] < c_b) {
      M[i] = 0
    }
  }
  
  # liczymy moc testu dla ustalonego a
  Moc[j]=mean(M)
}

plot(a, Moc) # pkt siatki vs moc empiryczna
