import numpy as np
import matplotlib.pyplot as plt
import functions as f


def Runge_Kutta(mode, x_l, x_r, y, n, f):
    if mode == 0:
        h = (x_r - x_l) / n
        graph = [y]
        last = y
        for i in range(0, n):
            cur = [0, 0]
            cur[0] = x_l + i * h
            cur[1] = last
            last = cur[1] + (h / 2) * (f(cur[0], cur[1]) + f(cur[0] + h, cur[1] + h * f(cur[0], cur[1])))
            graph.append(last)
        return graph
    elif mode == 1:
        h = (x_r - x_l) / n
        graph = [y]
        last = y
        for i in range(0, n):
            cur = [0, 0]
            coef = [0, 0, 0]
            cur[0] = x_l + i * h
            cur[1] = last
            coef[0] = f(cur[0] + (h / 2), cur[1] + (h / 2) * f(cur[0], cur[1]))
            coef[1] = f(cur[0] + (h / 2), cur[1] + (h / 2) * coef[0])
            coef[2] = f(cur[0] + h, cur[1] + h * coef[1])
            last = cur[1] + (h / 6) * (f(cur[0], cur[1]) + 2 * (coef[0] + coef[1]) + coef[2])
            graph.append(last)
        return graph
    else:
        c = len(y)
        h = (x_r - x_l) / n
        graph = [[], []]
        graph[0].append(y[0])
        graph[1].append(y[1])
        last = y
        for i in range(0, n):
            cur = [0, 0]
            coef = []
            for j in range(4):
                coef.append(np.empty(c))
            cur[0] = x_l + i * h
            cur[1] = last
            coef[0] = np.fromiter((f[j](cur[0], cur[1]) for j in range(c)), float)
            print(coef[0])
            coef[1] = np.fromiter((f[j](cur[0] + .5 * h, cur[1] + .5 * h * coef[0]) for j in range(c)), float)
            coef[2] = np.fromiter((f[j](cur[0] + .5 * h, cur[1] + .5 * h * coef[1]) for j in range(c)), float)
            coef[3] = np.fromiter((f[j](cur[0] + h, cur[1] + h * coef[2]) for j in range(c)), float)
            last = cur[1] + (h / 6) * (coef[0] + 2 * (coef[1] + coef[2]) + coef[3])
            graph[0].append(last[0])
            graph[1].append(last[1])
        return graph


def boundary_problem(n, l, r, y, funcs):
    p, q, func = funcs
    sigma, gamma, delta = y
    h = (r - l) / n
    graph = np.empty(n + 1)
    A = lambda x: 1 / (h * h) - .5 * p(l + x * h) / h
    B = lambda x: 1 / (h * h) + .5 * p(l + x * h) / h
    C = lambda x: -2 / (h * h) + q(l + x * h)
    a = np.empty(n + 1)
    b = np.empty(n + 1)
    a[1] = -gamma[0] / (sigma[0] * h - gamma[0])
    b[1] = delta[0] / (sigma[0] - gamma[0] / h)
    for i in range(1, n):
        a[i + 1] = -B(i) / (A(i) * a[i] + C(i))
        b[i + 1] = (func(l + h * i) - A(i) * b[i]) / (A(i) * a[i] + C(i))
    graph[n - 1], graph[n] = np.linalg.solve(np.array([[1, -a[n]], [gamma[1], -(sigma[1] * h + gamma[1])]]), np.array([b[n], -delta[1] * h]))
    for i in range(n - 2, -1, -1):
        graph[i] = a[i + 1] * graph[i + 1] + b[i + 1]
    return graph


def main():
    mode, test = map(int, input().split())
    npoints = 50
    preciseans = 1000
    plt.figure(figsize=(20, 10))
    if mode == 0:
        left = [0, 0]
        right = [1, 1]
        funcs = [f.f0, f.f1]
        funcs_solutions = [f.f0_solution, f.f1_solution]
        initializers = [10, 2]
        x = np.linspace(left[test], right[test], preciseans)
        plt.plot(x, funcs_solutions[test](x), color='blue', label='Solution')
        runge_res1 = Runge_Kutta(0, left[test], right[test], initializers[test], npoints - 1, funcs[test])
        plt.plot(np.linspace(left[test], right[test], npoints), runge_res1, 'o', color='red', label="SecondOrder")
        runge_res2 = Runge_Kutta(1, left[test], right[test], initializers[test], npoints - 1, funcs[test])
        plt.plot(np.linspace(left[test], right[test], npoints), runge_res2, 'o', color='green', label="FourthOrder")
        plt.grid()
        plt.legend()
        plt.show()
    elif mode == 1:
        left = [0, 0]
        right = [1, 1]
        funcs = [[f.system1_y1, f.system1_y2], [f.system4_y1, f.system4_y2]]
        funcs_solutions = [None, [f.system4_y1_solution, f.system4_y2_solution]]
        initializers = [[.5, 1], [1, 1]]
        if funcs_solutions[test] is not None:
            x = np.linspace(left[test], right[test], preciseans)
            plt.plot(x, funcs_solutions[test][0](x), color='blue', label='SolutionY1')
            plt.plot(x, funcs_solutions[test][1](x), color='green', label='SolutionY2')
        runge_res1 = Runge_Kutta(3, left[test], right[test], initializers[test], npoints, funcs[test])
        plt.plot(np.linspace(left[test], right[test], npoints + 1), runge_res1[0], 'o', label='Y1')
        plt.plot(np.linspace(left[test], right[test], npoints + 1), runge_res1[1], 'o', label='Y2')
        plt.grid()
        plt.legend()
        plt.show()
    else:
        left = [0, 0, 0]
        right = [1, 1, .5 * np.pi]
        funcs = [[f.p1, f.q1, f.func1], [f.p2, f.q2, f.func2], [f.p3, f.q3, f.func3]]
        funcs_solutions = [None, f.Bundle2_solution, f.Bundle3_solution]
        initializers = [[[0, .5], [1, 1], [1.3, 2]], [[1, -1], [0, 1], [-1, 2]], [[1, 1], [0, 0], [0, 0]]]
        if funcs_solutions[test] is not None:
            x = np.linspace(left[test], right[test], preciseans)
            plt.plot(x, funcs_solutions[test](x), color='blue', label="Solution")
        res = boundary_problem(npoints, left[test], right[test], initializers[test], funcs[test])
        plt.plot(np.linspace(left[test], right[test], npoints + 1), res, 'o', label='y')
        plt.grid()
        plt.legend()
        plt.show()
    #res = boundary_problem(npoints, left, right, , f.p1, f.q1, f.func1)
    #
    #res = boundary_problem(npoints, left, right, [1, .5], [0, -1], [2, 1], f.p2, f.q2, f.func2)

if __name__ == "__main__":
    main()

