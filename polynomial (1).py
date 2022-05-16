import math


class Polynomial:
    def __init__(self, coefficients):
        self.coefficients = coefficients

    def value(self, z):
        z_power, result = 1, 0
        for i in range(len(self.coefficients)):
            result += self.coefficients[i] * z_power
            z_power *= z
        return result

    def derivative(self):
        return Polynomial([self.coefficients[i] * i for i in range(1, len(self.coefficients))])

    @staticmethod
    def psi_mu(n, mu):  # mu can be complex
        ei = math.cos(-math.pi / n) + math.sin(-math.pi / n) * 1j
        k1 = 1j * ei / (2 * n * math.sin(math.pi / n))
        k2 = (mu ** (1 / n)) / (2 * n ** 2 * (2 * math.sin(math.pi / n)) ** 2) - (n - 1) * k1 / (2 * n)
        p1 = Polynomial([1 - k1] + [0] * (n - 1) + [k1])
        p2 = Polynomial([1] + [0] * (n - 1) + [-1])
        return p1 + Polynomial([k2]) * p2 * p2

    @staticmethod
    def __li(xi, i):
        p = Polynomial.develop_roots(xi[:i] + xi[i + 1:])
        return p * Polynomial([1 / p.value(xi[i])])

    @staticmethod
    def __fi(xi, vi, di, d2i, i):
        li = Polynomial.__li(xi, i)
        dli = li.derivative().value(xi[i])
        if d2i is None:
            da = di[i] / vi[i] - 2 * dli
            di = Polynomial([1 - da * xi[i], da])
            return Polynomial([vi[i]]) * li * li * di
        else:
            d2li = li.derivative().derivative().value(xi[i])
            c1 = di[i]/vi[i] - 3 * dli
            c2 = d2i[i]/vi[i] - 3 * d2li - 6 * dli * (c1 + dli)
            delta = 4 * c1 ** 2 - 8 * c2
            da = (2 * c1 + delta ** .5)/4
            db = c1 - da
            ai, bi = Polynomial([1 - da * xi[i], da]), Polynomial([1 - db * xi[i], db])
            return Polynomial([vi[i]]) * li * li * li * ai * bi

    @staticmethod
    def interpolator_polynomial(xi, vi, di, d2i=None):
        sum_fi = Polynomial([0])
        for i in range(len(xi)):
            sum_fi += Polynomial.__fi(xi, vi, di, d2i, i)
        return sum_fi

    @staticmethod
    def __coproduct(numbers, order):
        if order == 0:
            return 1
        result = 0
        for i in range(len(numbers) - order + 1):
            result += numbers[i] * Polynomial.__coproduct(numbers[i + 1:], order - 1)
        return result

    @staticmethod
    def develop_roots(roots, alpha=1+0j):
        coefficients = []
        for k in range(len(roots) + 1):
            coefficients = [Polynomial.__coproduct(roots, k) * alpha * (-1) ** k] + coefficients
        return Polynomial(coefficients)

    def __add__(self, other):
        n1 = len(self.coefficients)
        n2 = len(other.coefficients)
        if n1 > n2:
            coefficients = [self.coefficients[k] + other.coefficients[k] for k in range(n2)] + self.coefficients[n2:]
        else:
            coefficients = [self.coefficients[k] + other.coefficients[k] for k in range(n1)] + other.coefficients[n1:]
        while coefficients[-1] == 0 and len(coefficients) > 1:
            del coefficients[-1]
        return Polynomial(coefficients)

    def __mul__(self, other):
        n1 = len(self.coefficients)
        n2 = len(other.coefficients)
        coefficients = [0] * (n1 + n2 - 1)
        for i in range(n1):
            for j in range(n2):
                coefficients[i+j] += self.coefficients[i] * other.coefficients[j]
        return Polynomial(coefficients)

    def __str__(self):
        result = ''
        for k in range(len(self.coefficients)):
            result += str(self.coefficients[k]) + ' Z^' + str(k) + ' ; '
        return result
