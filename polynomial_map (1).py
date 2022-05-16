from PIL import Image
from polynomial import Polynomial

center, radius, max_iter, epsilon, size = 4.525+0j, .1, 250, 1e-10, 250
image = Image.new('RGB', (size,) * 2)

for x in range(size):
    for y in range(size):
        third_root = center + radius * (2 * (x / size + (y / size) * 1j) - 1 - 1j)
        image.putpixel((x, y), (0,) * 3)
        poly = Polynomial.develop_roots([1j, -1j, third_root])
        poly_der = poly.derivative()
        un = third_root / 3

        for i in range(max_iter):
            poly_der_un = poly_der.value(un)
            if poly_der_un == 0:
                image.putpixel((x, y), (255,) * 3)
                break
            un, un_1 = un - poly.value(un) / poly_der_un, un
            if abs(un - un_1) < epsilon:
                image.putpixel((x, y), (255,) * 3)
                break

image.save('polynomial_map.png', 'PNG')
