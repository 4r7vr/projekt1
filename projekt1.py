# - XYZ (geocentryczne) -> BLH (elipsoidalne fi, lambda, h)
# - BLH -> XYZ
# - XYZ -> NEUp
# - BL(GRS80, WGS84, ew. Krasowski) -> 2000
# - BL(GRS80, WGS84, ew. Krasowski) -> 1992

from math import sin, cos, tan, sqrt, atan, atan2, radians
from numpy import array, savetxt, column_stack
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
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
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
        with open(sciezka, 'r') as file:
            lines = file.readlines()

            if trans == 'XYZ2BLH' or trans == 'XYZ2NEUP' or trans == 'BLH2XYZ':
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
            return lists

    def hirvonen(self, X, Y, Z):
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

            fi_list.append(f)
            l_list.append(l)
            h_list.append(h)

        return fi_list, l_list, h_list

    def flh2xyz(f, l, h, self):
        N = self.a / sqrt(1 - self.e2 * sin(f)**2)
        X = (N * cos(f)+h) * cos(l)
        Y = (N * cos(f)+h) * sin(l)
        Z = ((N * (1 - self.e2)) + h) * sin(f)
        return(X, Y, Z)

    def pl1992(self, f, l):
        results_x92 = []
        results_y92 = []
        for f, l in zip(f, l):
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
        results_x2000 = []
        results_y2000 = []
        for f, l in zip(f, l):
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

        b2 = self.a**2 * (1 - self.e2)
        ep2 = (self.a**2 - b2) / b2

        delta_l = l - l0
        t = tan(f)
        ni2 = ep2 * (cos(f)**2)
        N = self.a / sqrt(1 - self.e2 * sin(f)**2)

        A0 = 1 - (self.e2/4) - (3 * self.e2**2 / 64) - (5 * self.e2**3 / 256)
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

    def Rneu(fa, la):
        R = array([[-sin(fa)*cos(la), -sin(la), cos(fa)*cos(la)],
                   [-sin(fa)*sin(la), cos(la), cos(fa)*sin(la)],
                   [cos(fa), 0, sin(fa)]])
        return(R)

    def xyz2neup(self, X, Y, Z, X0, Y0, Z0):
        neu = []
        p = sqrt(X0**2 + Y0**2)
        fi = atan(Z0 / (p*(1 - self.e2)))
        while True:
            N = self.a/sqrt(1-self.e2 * sin(fi)**2)
            h = (p / cos(fi)) - N
            fi_poprzednia = fi
            fi = atan((Z0 / p)/(1-((N * self.e2)/(N + h))))
            if abs(fi_poprzednia - fi) < (0.000001/206265):
                break
        N = self.a/sqrt(1-self.e2 * sin(fi)**2)
        h = p/cos(fi) - N
        lam = atan(Y0 / X0)

        R_neu = self.Rneu(fi, lam)
        X_sr = [X - X0, Y - Y0, Z - Z0]
        X_neu = R_neu.T@X_sr
        neu.append(X_neu)

        return(neu)

    def licz(self, plik_input, trans):
        dane = self.plik(plik_input, trans)  # Wczytanie danych z pliku
        X = dane["1"]
        Y = dane["2"]
        # Z = dane["3"]

        if trans == 'XYZ2BLH':
            Z = dane["3"]
            wyniki = self.hirvonen(X, Y, Z)
            savetxt(f"results{trans}_{args.el}.txt",
                    column_stack(wyniki), delimiter=' ')
            return wyniki
        if trans == 'BLH2XYZ':
            Z = dane["3"]
            wyniki = self.flh2xyz(X, Y, Z)
            savetxt(f"results{trans}_{args.el}.txt",
                    column_stack(wyniki), delimiter=' ')
            return wyniki
        # if trans == 'XYZ2NEUp':
        #     wyniki = self.hirvonen(X, Y, Z)
        if trans == 'BL2PL2000':
            wyniki = self.pl2000(X, Y)
            savetxt(f"results{trans}_{args.el}.txt",
                    column_stack(wyniki), delimiter=' ')
            return wyniki
        if trans == 'BL2PL1992':
            wyniki = self.pl1992(X, Y)
            savetxt(f"results{trans}_{args.el}.txt",
                    column_stack(wyniki), delimiter=' ')
            return wyniki


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-wsp', type=str)
    parser.add_argument('-el', type=str)
    parser.add_argument('-t', type=str)
    args = parser.parse_args()

    i = 0
    try:
        while i == 0:

            if args.el is None:
                args.el = input(
                    'Na jakiej elpisoidzie wykonywane będą obliczenia?: ').upper()
            if args.wsp is None:
                args.wsp = input('Wklej ścieżkę do pliku txt z danymi: ')
            if args.t is None:
                args.t = input('Jaką transformację chcesz wykonać?: ').upper()
            obiekt = Transformacje(args.el.upper())
            dane = obiekt.licz(args.wsp, args.t.upper())
            print(dane)
            i = i + 1

    finally:
        print('FIN :)')
