# Projekt 1 – Transformacje
## Opis programu:
Program pozwala użytkownikowi na przeliczanie współrzędnych umieszczonych w pliku na inny układ współrzędnych. Pozwala on także na wybranie elipsoidy do której przypisane będą obliczane współrzędne.  

## Wymagania do obsługi  programu
Projekt tworzony był za pomocą programu Python w wersji 3.11.5.
Zaimportowane biblioteki to‘numpy’ oraz ‘argparse’. 

## Opcje dostępne w programie:
### Dostępne elipsoidy:
* GRS80
* WGRS84
* KRASOWSKIEGO

### Dostępne transformacje: 
* XYZ - BLH
* BLH - XYZ
* XYZ - NEU
* BL- GK2000
* BL - GK1992

### Funkcje zdefiniowane w programie:
- `plik(sciezka, trans)` - wczytuje plik z danymi oraz przygotowuje je do dalszego przetwarzania;
- `hirvonen(X, Y, Z)` - konwertuje współrzędne X,Y,Z na współrzędne fi, lambda, H za pomocą metody Hirvonena.
- `flh2xyz(F, L, H)` - konwertuje współrzędne fi, lambda, h na X,Y,Z;
- `pl1992(f, l)` - konwertuje współrzędne fi, lambda, h na układ 1992;
- `pl2000(f, l)` - Konwertuje współrzędne fi, lambda, h na układ 2000;
- `Rneu(phi, lam)` - Oblicza macierz obrotu NEU;
- `xyz2neup(X, Y, Z, X0, Y0, Z0)` - konwertuje współrzędne X,Y,Z na współrzędne NEU względem punktu odniesienia;
- `licz(plik_input, trans)` - wykonuje odpowiednie przekształcenie na danych z pliku;

## Opis działania programu:
Program po wystartowaniu wymaga od użytkownika podania danych na podstawie kolejnych flag:
1.  Na jakiej elipsoidzie wykonywane będą obliczenia? (np.: grs80)
2.  Wklej ścieżkę do pliku txt z danymi: (np.: "C:\Users\User\Downloads\wsp_inp.txt")
3. Jaką transformację chcesz wykonać?: (np.: XYZ2BLH)
Transformacje oraz elipsoidy obsługiwane przez program podane zostały podane powyżej.
Po poprawnym wprowadzeniu wszystkich niezbędnych danych w folderze z programem pojawi się plik z przeliczonymi danymi – program wyświetla komunikat o końcu obliczeń i nazwie napisanego pliku.

### Przykładowy plik wejściowy i wyjściowy
* wejściowy: wsp.txt
* wyjściowy: resultsXYZ2BLH_GRS80.txt

### Przykładowe wywołanie pliku za pomocą wiersza poleceń:
* użytkownik musi znaleźć się w folderze w którym zapisany jest program
np.: C:\Users\User\Desktop\all\informatyka sem.4\git projekt 1>
* następnie w konsoli wpisuje polecenie wybierając odpowiednie opcje - py `nazwa projektu` -el `nazwa elipsoidy` -wsp `plik ze współrzędnymi` -t `nazwa transformacji` (nazwy zarówno transformacji jak i elipsoid możliwych do użycia w obliczeniach podane zostały wyżej)
np.: py projekt1.py -el grs80 -wsp wsp.txt -t xyz2blh

## Uwagi dotyczące pliku wejściowego
Plik ze współrzędnymi, który użytkownik wgrywa, nie może posiadać nagłówka, ponieważ program nie da rady go odczytać. (ze względu na brak flagi z zapytaniem o ilość linijek nagłówka)

