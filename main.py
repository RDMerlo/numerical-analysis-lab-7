import numpy as np


def f(x, t):
    return alpha * (x ** 2 - 2 * t)


def print_matrix(a):
    a = a.transpose()
    for i in range(len(a) - 1, -1, -1):
        print("%1.2f" % (i * 0.01), end=' ')
        for j in range(len(a[i])):
            print("%2.4f" % (np.abs(a[i][j])), end=' ')
        print()
    print("T/X ", end=' ')
    for i in range(len(a[0])):
        print("%2.2f  " % (0.1 * i), end=' ')
        print("--------")
    print()


def tma(a, b, c, d):
    n = len(d)
    x = np.zeros(n)
    p = [-c[0] / b[0]]
    q = [d[0] / b[0]]
    for i in range(1, n):
        p.append(-c[i] / (b[i] + a[i] * p[i - 1]))
        q.append((d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]))
    x[-1] = q[-1]
    for i in reversed(range(n - 1)):
        x[i] = p[i] * x[i + 1] + q[i]
    return x


N = 4
alpha = 0.5 * N
n = 11
m = 6

x = np.empty(n)
t = np.empty(m)
u = np.zeros([n, m])

h = 0.1
tau = 0.01
s = tau / h ** 2

a = s
b = s
c = (1 - 2 * s)

for i in range(n):
    x[i] = i * h
for j in range(m):
    t[j] = j * tau
print('Явная')

for j in range(m):
    u[0, j] = 0
    u[n - 1, j] = alpha * t[j]

for i in range(1, n):
    u[i, 0] = 0

for j in range(0, m - 1):
    for i in range(1, n - 1):
        u[i, j + 1] = a * u[i - 1, j] + c * u[i, j] + b * u[i + 1, j] + tau * f(x[i], t[j])

print_matrix(u)

# print('Кранка - Николсона')
# print('teta = ', end = '')
# teta = float(input())
teta = 0.5
U = np.zeros((m, n))
for j in range(n):
    U[0, j] = 0
for i in range(1, m):
    aa = np.zeros(n)
    bb = np.zeros(n)
    cc = np.zeros(n)
    dd = np.zeros(n)
    dd[0] = 0
    dd[-1] = alpha * t[i]
    bb[0] = 1
    bb[-1] = 1
    for j in range(1, n - 1):
        aa[j] = teta * (2 * a - b * h) * tau * 0.5 / (h ** 2)
        bb[j] = c * tau * teta - teta * 2 * a * tau / (h ** 2) - 1
        cc[j] = teta * (2 * a + h * b) * tau * 0.5 / (h ** 2)
        dd[j] = -U[i - 1, j] - teta * tau * f(x[j], t[i]) - (1 - teta) * (
                    (a * tau / h ** 2) * (U[i - 1, j - 1] - 2 * U[i - 1, j] + U[i - 1, j + 1]) + (b * tau * 0.5 / h) * (
                        - U[i - 1, j - 1] + U[i - 1, j + 1]) + c * tau * U[i - 1, j] + tau * f(x[j], t[i - 1]))
    xx = tma(aa, bb, cc, dd)
    for j in range(n):
        U[i, j] = xx[j]
print_matrix(U.transpose())

