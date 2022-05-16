class Rational:
    __displayed_digits = 10
    __max_size = 1000

    def __init__(self, num, den):
        d = Rational.__gcd(num, den)
        self.numerator = num // d
        self.denominator = den // d
        self.__trunc()

    def __trunc(self):  # maximum relative error is under 10 ** (1 - Rational.__max_size)
        if Rational.__max_size > 0:
            t = min(len(str(self.numerator)), len(str(self.denominator))) - Rational.__max_size
            if t > 0:
                self.numerator //= 10 ** t
                self.denominator //= 10 ** t

    @staticmethod
    def __gcd(a, b):
        a, b = abs(a), abs(b)
        while a != 0:
            if b < a:
                a, b = b, a
            else:
                a, b = b % a, a
        return b

    @staticmethod
    def set_parameters(displayed_digits=10, max_size=1000):  # 0, 0 means rational form display and exact mode
        Rational.__displayed_digits, Rational.__max_size = displayed_digits, max_size

    def __add__(self, other):
        num = self.numerator * other.denominator + self.denominator * other.numerator
        den = self.denominator * other.denominator
        return Rational(num, den)

    def __neg__(self):
        return Rational(-self.numerator, self.denominator)

    def __sub__(self, other):
        return -other + self

    def __mul__(self, other):
        num, den = self.numerator * other.numerator, self.denominator * other.denominator
        d = Rational.__gcd(num, den)
        return Rational(num // d, den // d)

    def __truediv__(self, other):
        signum = 2 * int(other.numerator > 0) - 1
        return Rational(signum * self.numerator * other.denominator, self.denominator * abs(other.numerator))

    def __pow__(self, power, modulo=None):
        if power < 0:
            return Rational(1, 1) / (self ** (-power))
        else:
            return Rational(self.numerator ** power, self.denominator ** power)

    def __abs__(self):
        return Rational(abs(self.numerator), self.denominator)

    def __eq__(self, other):
        return self.numerator * other.denominator == self.denominator * other.numerator

    def __ne__(self, other):
        return self.numerator * other.denominator != self.denominator * other.numerator

    def __le__(self, other):
        return self.numerator * other.denominator <= other.numerator * self.denominator

    def __lt__(self, other):
        return self.numerator * other.denominator < other.numerator * self.denominator

    def __ge__(self, other):
        return self.numerator * other.denominator >= other.numerator * self.denominator

    def __gt__(self, other):
        return self.numerator * other.denominator > other.numerator * other.denominator

    def __str__(self):
        if Rational.__displayed_digits == 0:
            return str(self.numerator) + '/' + str(self.denominator)
        else:
            return self.__decimal_form()

    def __decimal_form(self):
        result = ''
        a, b = self.numerator, self.denominator
        if a < 0:
            a = -a
            result += '-'
        exponent = len(str(a)) - len(str(b))
        if exponent > 0:
            b *= 10 ** exponent
        else:
            a *= 10 ** -exponent
        if a // b == 0 and a != 0:
            exponent -= 1
            a *= 10
        for i in range(self.__displayed_digits):
            if i == 1:
                result += '.'
            result += str((a // b) % 10)
            a *= 10
        return result + 'e' + str(exponent)


# quick test: calculate square root of 2 using Newton method

Rational.set_parameters(displayed_digits=100)

r2 = Rational(1, 1)

for i in range(10):
    r2 = (r2 ** 2 + Rational(2, 1))/(Rational(2, 1) * r2)

print(r2)
