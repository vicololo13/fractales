from PIL import Image
from polynomial import Polynomial
import math


class NewtonRenderer:
    x_light, y_light, z_light = -math.sin(1) / math.sqrt(2), -math.cos(1) / math.sqrt(2), -1 / math.sqrt(2)
    max_iter, epsilon, alpha, beta, delta = 100, 1e-10, .1, 10, 1e-5
    color_list = [(255, 150, 100), (150, 255, 100), (100, 150, 255),
                  (255, 255, 150), (255, 150, 255), (150, 255, 255), (0, 0, 0)]

    def __init__(self, poly, center=0j, radius=1+0j, size=250, mode='full'):  # radius can be complex
        self.poly, self.center, self.radius, self.size, self.mode = poly, center, radius, size, mode
        self.poly_der = poly.derivative()
        self.roots = []
        if len(poly.coefficients) > len(NewtonRenderer.color_list):
            self.mode = 'no_color'

    def __root_index(self, z, epsilon):
        for i in range(len(self.roots)):
            if abs(z - self.roots[i]) < epsilon:
                return i
        self.roots.append(z)
        return len(self.roots) - 1

    def __iterate(self, z):
        poly_der_z = self.poly_der.value(z)
        if poly_der_z == 0:
            return 'undefined'
        return z - self.poly.value(z) / poly_der_z

    def __convergence_index(self, z):
        un = z
        for i in range(self.max_iter):
            un, un_1 = self.__iterate(un), un
            if un == 'undefined':
                return -1
            elif abs(un - un_1) < self.epsilon:
                return self.__root_index(un, 2 * self.epsilon)
        return -1

    def __v_alpha(self, z):
        root_index = self.__convergence_index(z)
        if root_index == -1:
            return 0
        limit, un, sum_distances = self.roots[root_index], z, 0
        while abs(un - limit) > self.epsilon:
            sum_distances, un = sum_distances + abs(un - limit), self.__iterate(un)
        return self.alpha / (self.alpha + sum_distances)

    def __orientation_coefficient(self, z0):
        if self.mode == 'color_only':
            return 1
        z1, z2 = z0 + self.delta, z0 + self.delta * 1j
        v_z0, v_z1, v_z2 = self.__v_alpha(z0), self.__v_alpha(z1), self.__v_alpha(z2)
        n_x, n_y = (v_z1 - v_z0) / self.delta, (v_z2 - v_z0) / self.delta
        n_norm = (self.beta ** -2 + n_x ** 2 + n_y ** 2) ** .5
        return (1 - (self.x_light * n_x + self.y_light * n_y + self.z_light / self.beta) / n_norm) / 2

    @staticmethod
    def __interpolate(x, a):
        return x ** a / (x ** a + (1 - x) ** a)

    def __point_color(self, z):
        cr, cg, cb = 255, 255, 255
        if self.mode != 'no_color':
            cr, cg, cb = self.color_list[self.__convergence_index(z)]
        oz = self.__orientation_coefficient(z)
        oz_r, oz_b = NewtonRenderer.__interpolate(oz, 2), NewtonRenderer.__interpolate(oz, .5)
        return int(cr * oz_r), int(cg * oz), int(cb * oz_b)

    def render(self):
        image = Image.new('RGB', (self.size,) * 2)
        for x in range(self.size):
            for y in range(self.size):
                z = self.center + self.radius * (2 * (x / self.size + (y / self.size) * 1j) - 1 - 1j)
                image.putpixel((x, y), self.__point_color(z))
        image.save('rendered_fractal.png', 'PNG')


# quick demonstration with psi_2^1

p = Polynomial.psi_mu(2, 1)

renderer = NewtonRenderer(p, radius=1+1j)
renderer.render()
