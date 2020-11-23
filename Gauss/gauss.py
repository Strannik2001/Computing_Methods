import sys
import random
from math import sin


class Test:
    def __init__(self, mode, n, a, b=1):
        self.mode = mode
        self.n = n
        self.need_optimization = b
        self.matrix = []
        self.b = []
        self.q = 0
        self.determinant = 1
        self.invert_matrix = []
        for i in range(n):
            self.invert_matrix.append([])
            self.invert_matrix[i] = [0] * self.n
            self.invert_matrix[i][i] = 1
        self.prepare_matrix(a)
        self.matrix_copy = []
        for i in range(n):
            self.matrix_copy.append([])
            for j in range(n):
                self.matrix_copy[i].append(self.matrix[i][j])

    def prepare_matrix(self, a):
        if self.mode == 'static':
            with open(a, "r") as f:
                for i in range(self.n):
                    self.matrix.append(list(map(float, f.readline().split())))
                self.b = list(map(float, f.readline().split()))
        else:
            for i in range(self.n):
                self.matrix.append([])
                self.matrix[i] = [0] * self.n
            self.b = [0] * self.n
            self.q = 1.001 - 2 * a / 1000
            x = random.random()
            for i in range(self.n):
                self.b[i] = abs(x - self.n / 10) * (i + 1) * sin(x)
            for i in range(1, self.n + 1):
                for j in range(1, self.n + 1):
                    if i == j:
                        self.matrix[i - 1][j - 1] = pow(self.q - 1, i + j)
                    else:
                        self.matrix[i - 1][j - 1] = pow(self.q, i + j) + 0.1 * (j - i)

    def out_matrix(self, a='matrix'):
        if a == 'invert matrix':
            print("Invert matrix")
            m = self.invert_matrix
        else:
            print("Matrix")
            m = self.matrix
        out_matrix(m)

    def gaus(self):
        sign = 1
        for i in range(self.n):
            if self.need_optimization:
                maxi = i
                for j in range(i, self.n):
                    if abs(self.matrix[maxi][i]) < abs(self.matrix[j][i]):
                        maxi = j
                if maxi != i:
                    sign *= -1
                self.b[i], self.b[maxi] = self.b[maxi], self.b[i]
                self.matrix[i], self.matrix[maxi] = self.matrix[maxi], self.matrix[i]
                self.invert_matrix[i], self.invert_matrix[maxi] = self.invert_matrix[maxi], self.invert_matrix[i]
            cur = self.matrix[i][i]
            self.determinant *= cur
            self.matrix[i] = [z / cur for z in self.matrix[i]]
            self.b[i] /= cur
            self.invert_matrix[i] = [z / cur for z in self.invert_matrix[i]]
            for j in range(i + 1, self.n):
                cur = self.matrix[j][i]
                self.b[j] -= cur * self.b[i]
                self.matrix[j] = [self.matrix[j][z] - cur * self.matrix[i][z] for z in range(self.n)]
                self.invert_matrix[j] = [self.invert_matrix[j][z] - cur * self.invert_matrix[i][z] for z in range(self.n)]
        self.determinant *= sign
        print("Determinant = {}".format(self.determinant))
        for i in range(self.n)[::-1]:
            for j in range(0, i):
                self.invert_matrix[j] = [self.invert_matrix[j][z] - self.matrix[j][i] * self.invert_matrix[i][z] for z in range(self.n)]
                self.b[j] -= self.b[i] * self.matrix[j][i]
                self.matrix[j] = [self.matrix[j][z] - self.matrix[j][i] * self.matrix[i][z] for z in range(self.n)]
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        print("Answer")
        for j in range(self.n - 1):
            print("%015.10f" % self.b[j], end='  |^|  ')
        print("%015.10f" % self.b[self.n - 1])
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")

def out_matrix(a):
    n = len(a)
    print("##########################################")
    for i in range(n):
        for j in range(n - 1):
            print("%015.10f" % a[i][j], end='  |^|  ')
        print("%015.10f" % a[i][n - 1])
    print("##########################################")


def mul_matrix(a, b):
    n = len(a)
    res = []
    for i in range(n):
        res.append([])
        for j in range(n):
            res[i].append(0)
            for k in range(n):
                res[i][j] += a[i][k] * b[k][j]
    return res


def main():
    if sys.argv[1] == 'static':
        matrix = Test('static', int(sys.argv[2]), sys.argv[3])
    elif sys.argv[1] == 'dynamic':
        matrix = Test('dynamic', int(sys.argv[2]), float(sys.argv[3]))
    else:
        print("##Usage##")
        print("First argument - 'static' or 'dynamic'")
        print("Second argument - rang of matrix")
        print("Third argument - path_to_file if static, starting number if dynamic")
        exit(0)
    matrix.gaus()
    out_matrix(matrix.invert_matrix)


if __name__ == "__main__":
    main()
