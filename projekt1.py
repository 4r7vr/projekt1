from math import sin, cos, tan, sqrt, atan, atan2, radians, pi
from numpy import array, savetxt, column_stack, vstack, dot, reshape
from argparse import ArgumentParser


class Transformacje:
    def __init__(self, model):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            e2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        """
        if model == "WGS84":
            self.a = 6378137.0
            self.b = 6356752.31424518
        elif model == "GRS80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "KRASOWSKI":
            self.a = 6378245.000
            self.b = 6356863.019
        else:
            raise NotImplementedError(f"{model} model not implemented")

        self.e = sqrt(self.a**2 - self.b**2)/self.a
        self.e2 = self.e**2

    def plik(self, sciezka, trans):
        """
        INPUT: przyjmuje plik zewnętrzny oraz parametr trans (string)
        OUTPUT: zwraca współrzędne z pliku w formie list dostosowywując je do 
            późniejszej transformacji (np. [X], [Y], [Z])

        """
        with open(sciezka, 'r') as file:
            lines = file.readlines()

        if len(lines) == 1:  # Sprawdzenie, czy plik zawiera tylko jedną linię
            # Usunięcie znaków białych z początku i końca linii
            line = lines[0].strip()
            parts = line.split()
            if trans == 'XYZ2NEUP':
                lists = {"1": [float(parts[0])], "2": [float(parts[1])], "3": [float(parts[2])],
                         "X0": [], "Y0": [], "Z0": []}
                lists["X0"].append(float(parts[0]))
                lists["Y0"].append(float(parts[1]))
                lists["Z0"].append(float(parts[2]))
            else:
                lists = {"1": [float(parts[0])], "2": [
                    float(parts[1])], "3": [float(parts[2])]}
        else:
            if trans == 'XYZ2BLH' or trans == 'BLH2XYZ':
                lists = {"1": [], "2": [], "3": []}
                for line in lines:
                    parts = line.split()
                    lists["1"].append(float(parts[0]))
                    lists["2"].append(float(parts[1]))
                    lists["3"].append(float(parts[2]))

            elif trans == 'BL2PL2000' or trans == 'BL2PL1992':
                lists = {"1": [], "2": []}
                for line in lines:
                    parts = line.split()
                    lists["1"].append(float(parts[0]))
                    lists["2"].append(float(parts[1]))

            elif trans == 'XYZ2NEUP':
                lists = {"1": [], "2": [], "3": [],
                         "X0": [], "Y0": [], "Z0": []}
                parts = lines[0].split()
                lists["X0"].append(float(parts[0]))
                lists["Y0"].append(float(parts[1]))
                lists["Z0"].append(float(parts[2]))

                for line in lines[1:]:
                    parts = line.split()
                    lists["1"].append(float(parts[0]))
                    lists["2"].append(float(parts[1]))
                    lists["3"].append(float(parts[2]))

        return lists

    def xyz2flh(self, X, Y, Z):
        """
        INPUT: przyjmuje współrzędne X, Y, Z - w formie list
        OUTPUT: zwraca współrzędne fi, lambda, h - w formie list

        """
        fi_list = []
        l_list = []
        h_list = []

        for x, y, z in zip(X, Y, Z):
            l = atan2(y, x)
            p = sqrt(x**2 + y**2)
            f = atan(z / (p * (1 - self.e2)))

            while True:
                N = self.a / sqrt(1 - self.e2 * sin(f)**2)
                h = p / cos(f) - N
                fs = f
                f = atan(z / (p * (1 - (self.e2 * (N / (N + h))))))
                if abs(fs - f) < (0.000001/206265):
                    break

            fi_list.append((f)*(180/pi))
            l_list.append((l)*(180/pi))
            h_list.append(h)

        return fi_list, l_list, h_list

    def flh2xyz(self, F, L, H):
        """
        INPUT: przyjmuje współrzędne fi, lambda, h - w formie list
        OUTPUT: zwraca współrzędne X, Y, Z - w formie list

        """
        x_list = []
        y_list = []
        z_list = []

        for f, l, h in zip(F, L, H):
            f = f * (pi/180)
            l = l * (pi/180)
            N = self.a / sqrt(1 - self.e2 * sin(f)**2)
            X = (N + h) * cos(f) * cos(l)
            Y = (N + h) * cos(f) * sin(l)
            Z = ((N * (1 - self.e2)) + h) * sin(f)
            x_list.append(X)
            y_list.append(Y)
            z_list.append(Z)
        return(x_list, y_list, z_list)

    def pl1992(self, f, l):
        """
        INPUT: przyjmuje współrzędne fi oraz lambda - w formie list
        OUTPUT: zwraca współrzędne X, Y w układzie 1992 - w formie list

        """
        results_x92 = []
        results_y92 = []

        for f, l in zip(f, l):
            f = f*(pi/180)
            l = l*(pi/180)
            l0 = radians(19)
            m = 0.9993
            N = self.a / sqrt(1 - self.e2 * sin(f)**2)

            b2 = self.a**2 * (1 - self.e2)
            ep2 = (self.a**2 - b2) / b2

            delta_l = l - l0
            t = tan(f)
            ni2 = ep2 * (cos(f)**2)
            A0 = 1 - (self.e2 / 4) - (3 * self.e2**2 / 64) - \
                (5 * self.e2**3 / 256)
            A2 = (3 / 8) * (self.e2 + (self.e2**2 / 4) + (15 * self.e2**3 / 128))
            A4 = (15 / 256) * (self.e2**2 + ((3 * self.e2**3) / 4))
            A6 = (35 * self.e2**3) / 3072

            sigma = self.a * (A0 * f - A2 * sin(2 * f) +
                              A4 * sin(4 * f) - A6 * sin(6 * f))

            xgk = sigma + (((delta_l**2 / 2) * N * sin(f) * cos(f)) * (1 + ((delta_l**2 / 12) * (cos(f)**2) * (5 - t**2 + 9 *
                           ni2 + 4 * ni2**2)) + ((delta_l**4 / 360) * (cos(f)**4) * (61 - 58 * t**2 + t**4 + 270 * ni2 - 330 * ni2 * t**2))))

            ygk = (delta_l * N * cos(f)) * (1 + ((delta_l**2 / 6) * (cos(f)**2) * (1 - t**2 + ni2)) +
                                            (((delta_l**4 / 120) * (cos(f)**4)) * (5 - (18 * t**2) + t**4 + (14 * ni2) - (58 * ni2 * t**2))))

            x92 = xgk * m - 5300000
            y92 = ygk * m + 500000
            results_x92.append(x92)
            results_y92.append(y92)

        return results_x92, results_y92

    def pl2000(self, f, l):
        """
        INPUT: przyjmuje współrzędne fi oraz lambda - w formie list
        OUTPUT: zwraca współrzędne X, Y w układzie 2000 - w formie list

        """
        results_x2000 = []
        results_y2000 = []

        for f, l in zip(f, l):
            f = f*(pi/180)
            l = l*(pi/180)
            m = 0.999923
            l0 = 0
            strefa = 0

            if l > radians(13.5) and l < radians(16.5):
                strefa = 5
                l0 = radians(15)
            elif l > radians(16.5) and l < radians(19.5):
                strefa = 6
                l0 = radians(18)
            elif l > radians(19.5) and l < radians(22.5):
                strefa = 7
                l0 = radians(21)
            elif l > radians(22.5) and l < radians(25.5):
                strefa = 8
                l0 = radians(24)
            else:
                print("Punkt poza strefami odwzorowawczymi układu PL-2000")
                continue

            b2 = self.a**2 * (1 - self.e2)
            ep2 = (self.a**2 - b2) / b2

            delta_l = l - l0
            t = tan(f)
            ni2 = ep2 * (cos(f)**2)
            N = self.a / sqrt(1 - self.e2 * sin(f)**2)

            A0 = 1 - (self.e2/4) - (3 * self.e2**2 / 64) - \
                (5 * self.e2**3 / 256)
            A2 = (3/8) * (self.e2 + (self.e2**2 / 4) + (15 * self.e2**3 / 128))
            A4 = (15/256) * (self.e2**2 + ((3 * self.e2**3) / 4))
            A6 = (35 * self.e2**3) / 3072

            sigma = self.a * (A0 * f - A2 * sin(2 * f) + A4 *
                              sin(4 * f) - A6 * sin(6 * f))

            xgk = sigma + (((delta_l**2 / 2) * N * sin(f) * cos(f)) * (1 + ((delta_l**2 / 12) * (cos(f)**2) * (5 - t**2 + 9 *
                           ni2 + 4 * ni2**2)) + ((delta_l**4 / 360) * (cos(f)**4) * (61 - 58 * t**2 + t**4 + 270 * ni2 - 330 * ni2 * t**2))))

            ygk = (delta_l * N * cos(f)) * (1 + ((delta_l**2 / 6) * (cos(f)**2) * (1 - t**2 + ni2)) +
                                            (((delta_l**4 / 120) * (cos(f)**4)) * (5 - (18 * t**2) + t**4 + (14 * ni2) - (58 * ni2 * t**2))))

            x2000 = xgk * m
            y2000 = ygk * m + (strefa * 1000000) + 500000
            results_x2000.append(x2000)
            results_y2000.append(y2000)

        return results_x2000, results_y2000

    def Rneu(self, phi, lam):
        """
        INPUT: przyjmuje współrzędne fi oraz lambda - w formie list
        OUTPUT: zwraca macierz obrotu NEU - w formie array

        """

        R = array([[-sin(phi)*cos(lam), -sin(lam), cos(phi)*cos(lam)],
                   [-sin(phi)*sin(lam), cos(lam), cos(phi)*sin(lam)],
                   [cos(phi), 0, sin(phi)]])
        return(R)

    def xyz2neup(self, X, Y, Z, X0, Y0, Z0):
        """
        INPUT: przyjmuje współrzędne X, Y, Z - w formie list, oraz X0, Y0,
        Z0 - również w formie list    
        OUTPUT: zwraca wyniki transformacji NEU w formie array

        """
        neu_results = []

        for X0, Y0, Z0 in zip(X0, Y0, Z0):
            lam = atan2(Y0, X0)
            p = sqrt(X0**2 + Y0**2)
            phi = atan(Z0 / (p * (1 - self.e2)))

            while True:
                N = self.a / sqrt(1 - self.e2 * sin(phi)**2)
                h = p / cos(phi) - N
                phi2 = phi
                phi = atan(Z0 / (p * (1 - (self.e2 * (N / (N + h))))))
                if abs(phi2 - phi) < (0.000001/206265):
                    break

        R = array([[-sin(phi) * cos(lam), -sin(lam), cos(phi) * cos(lam)],
                   [-sin(phi) * sin(lam),  cos(lam), cos(phi) * sin(lam)],
                   [cos(phi), 0, sin(phi)]])

        X_list = []
        Y_list = []
        Z_list = []
        for X, Y, Z in zip(X, Y, Z):
            X_list.append(X)
            Y_list.append(Y)
            Z_list.append(Z)

        X_list = array(X_list)
        Y_list = array(Y_list)
        Z_list = array(Z_list)

        xyz = column_stack([reshape(X_list, (len(X_list), 1)), reshape(
            Y_list, (len(Y_list), 1)), reshape(Z_list, (len(Z_list), 1))])

        xyz0 = array([X0, Y0, Z0]).T

        xyzt = xyz-xyz0

        neu = R.T @ xyzt.T
        neu_results.append(neu.T)

        return neu_results

    def licz(self, plik_input, trans):
        """
        INPUT: przyjmuje plik z zewnątrz oraz parametr trans (string)
        OUTPUT: zwraca wyniki transformacji wprowadzonej przez użytkownika w 
        formie listy, a poszczególne funkcje warunkowe zpisują je do pliku

        """
        dane = self.plik(plik_input, trans)  # Wczytanie danych z pliku

        if trans == 'XYZ2BLH':
            X = dane["1"]
            Y = dane["2"]
            Z = dane["3"]
            wyniki = self.xyz2flh(X, Y, Z)
            savetxt(f"results{trans}_{args.el}.txt",
                    column_stack(wyniki), fmt='%.6f', delimiter=' ')
            return wyniki

        if trans == 'BLH2XYZ':
            X = dane["1"]
            Y = dane["2"]
            Z = dane["3"]
            wyniki = self.flh2xyz(X, Y, Z)
            savetxt(f"results{trans}_{args.el}.txt",
                    column_stack(wyniki), fmt='%.4f', delimiter=' ')
            return wyniki

        if trans == 'XYZ2NEUP':
            X = dane["1"]
            Y = dane["2"]
            Z = dane["3"]
            X0 = dane["X0"]
            Y0 = dane["Y0"]
            Z0 = dane["Z0"]
            wyniki = self.xyz2neup(X, Y, Z, X0, Y0, Z0)
            savetxt(f"results{trans}_{args.el}.txt",
                    vstack(wyniki),fmt='%.4f', delimiter=' ')
            return wyniki

        if trans == 'BL2PL2000':
            X = dane["1"]
            Y = dane["2"]
            wyniki = self.pl2000(X, Y)
            savetxt(f"results{trans}_{args.el}.txt",
                    column_stack(wyniki),fmt='%.4f', delimiter=' ')
            return wyniki

        if trans == 'BL2PL1992':
            X = dane["1"]
            Y = dane["2"]
            wyniki = self.pl1992(X, Y)
            savetxt(f"results{trans}_{args.el}.txt",
                    column_stack(wyniki),fmt='%.4f', delimiter=' ')
            return wyniki


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument(
        '-el', type=str, help='Przyjmuje nazwe elipsoidy, na której ma zostac wykonana transformacja. Dostepne elipsoidy: WGS84, GRS80, KRASOWSKI')
    parser.add_argument(
        '-wsp', type=str, help='Przyjmuje sciezke do pliku ze wspolrzednymi; jezeli plik jest w tym samym folderze co program wystarczy wpisac nazwe pliku z rozszerzeniem')
    parser.add_argument(
        '-t', type=str, help='Przyjmuje nazwe transformacji, ktora ma zostac wykonana. Dostepne transformacje: XYZ2BLH, BLH2XYZ, XYZ2NEUP, BL2PL2000, BL2PL1992')
    args = parser.parse_args()

    try:
        if args.el is None:
            args.el = input(
                'Na jakiej elpisoidzie wykonywane będą obliczenia?: ').upper()
        if args.wsp is None:
            args.wsp = input('Wklej ścieżkę do pliku txt z danymi: ')
        if args.t is None:
            args.t = input('Jaką transformację chcesz wykonać?: ').upper()

        obiekt = Transformacje(args.el.upper())
        dane = obiekt.licz(args.wsp, args.t.upper())

    except FileNotFoundError:
        print('Podany plik nie istnieje.')
    except KeyError:
        print('Nieprawidlowa elipsoida lub transformacja.')
    except IndexError:
        print('Nieprawidlowy format danych w pliku.')
    except ValueError:
        print('Nieprawidlowy format danych w pliku.')
    finally:
        nazwa_pliku = f"results{args.t.upper()}_{args.el.upper()}.txt"
        print('Transformacja zakończona pomyślnie - dane zapisane do pliku:', 
              nazwa_pliku)