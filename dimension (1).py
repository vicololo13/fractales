from polynomial import Polynomial
import random
import math

poly = Polynomial([-1, 0, 0, 1])
poly_der = poly.derivative()
max_iter, epsilon = 100, .1


def convergence_index(z0):
    un = z0
    for i in range(max_iter):
        poly_der_un = poly_der.value(un)
        if poly_der_un == 0:
            return -1
        un, un_1 = un - poly.value(un) / poly_der_un, un
        if abs(un_1 - un) < epsilon:
            if un.real > 0:
                return 0
            elif un.imag > 0:
                return 1
            else:
                return 2
    return -1


n, tn = int(input('n = ')), 0
n_points = 50
e = 2 / n

for x in range(n):
    for y in range(n):
        z = 1 + 1j - 2 * (x / n + 1j * y / n)
        index = convergence_index(z)
        for p in range(n_points):
            delta = e * (random.random() + random.random() * 1j)
            if convergence_index(z + delta) != index:
                tn += 1
                break

print('estimated T_n :', tn)
