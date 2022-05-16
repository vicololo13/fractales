class PRNGPrototype:
    def __init__(self, seed=1/2):
        self.seed = seed
        self.re_init = 0
        self.un = 1/4 + self.seed*1/2

    def __bit(self):
        self.un = (self.un - 1/self.un)/2
        if self.un == 0:
            self.re_init += 1
            self.un = self.seed + self.re_init
        return int(self.un > 0)

    def next_number(self, bits=32):
        n = 0
        for i in range(bits):
            n = 2*n + self.__bit()
        return n / (2 ** bits - 1)


# quick test: average and variance
# (it should be 1/2 and 1/12, respectively)

gen = PRNGPrototype()
sample_size = 10 ** 4
sum_numbers = 0
sum_squares = 0

for k in range(sample_size):
    d = gen.next_number()
    sum_numbers += d
    sum_squares += d ** 2

average_numbers = sum_numbers / sample_size
average_squares = sum_squares / sample_size

print('average :', average_numbers)
print('variance :', average_squares - average_numbers ** 2)
