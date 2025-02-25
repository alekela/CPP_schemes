import numpy as np
import matplotlib.pyplot as plt


def Boundary(P, u, uleft, pright):
    u[0] = uleft
    P[0] = P[1]
    P[-1] = pright
    u[-1] = u[-2]


def SIMPLE():
    xstart = 0
    xend = 1
    N = 10
    dx = (xend - xstart) / N

    rho0 = 1
    nu = 0.01
    u0 = 0
    P0 = 0
    uleft = 1
    pright = 0
    prec = 1e-5
    tau = 0.001

    grid = [xstart + i * dx for i in range(N + 1)]
    grid_centers = [(grid[i + 1] + grid[i]) / 2 for i in range(N)]

    P = [0 for i in range(N)];
    u = [0 for i in range(N)];
    rho = [0 for i in range(N)];
    I = [0 for i in range(N)]
    for j in range(N):
        P[j] = P0
        u[j] = u0
        rho[j] = rho0
    Boundary(P, u, uleft, pright)
    for t in range(1):
        iter = 0
        flag = True
        while (flag):
            tmp_u = [0 for _ in range(N)]
            coeffs_dp = [[0 for _ in range(N - 2)] for _ in range(N - 2)]
            res_dp = [0 for _ in range(N - 2)]
            for i in range(1, N - 1):
                tmp_u[i] = u[i] + tau * (-u[i] * (u[i] - u[i - 1]) / dx + nu * (u[i + 1] - 2 * u[i] + u[i - 1]) / dx / dx)
            tmp_u[0] = uleft

            for j in range(1, N - 1):
                if j == 1:
                    coeffs_dp[j - 1][j - 1] = -1
                    coeffs_dp[j - 1][j] = 1
                    res_dp[j - 1] = dx / tau / 2 * rho0 * (tmp_u[j + 1] - tmp_u[j - 1])
                elif j == N - 2:
                    coeffs_dp[j - 1][j - 1] = 1
                    res_dp[j - 1] = pright
                else:
                    coeffs_dp[j - 1][j - 1] = -2
                    coeffs_dp[j - 1][j] = 1
                    coeffs_dp[j - 1][j - 2] = 1
                    res_dp[j - 1] = dx / tau / 2 * rho0 * (tmp_u[j + 1] - tmp_u[j - 1])

            tmp_dP = np.linalg.solve(coeffs_dp, res_dp)
            tmp_dP = [tmp_dP[0]] + list(tmp_dP) + [pright]

            new_P = [0 for _ in range(N)];
            new_u = [0 for _ in range(N)]
            for j in range(1, N - 1):
                new_P[j] = P[j] + tmp_dP[j]
                new_u[j] = tmp_u[j] - tau / rho0 / 2 / dx * (tmp_dP[j + 1] - tmp_dP[j - 1])

            for j in range(1, N - 1):
                P[j] = new_P[j]
                u[j] = new_u[j]
            Boundary(P, u, uleft, pright)

            iter += 1
            R = max([abs(u[j + 1] - u[j - 1]) / 2 / dx for j in range(1, N - 1)])
            if R < prec:
                flag = False
        print("Сошлось за", iter, "итераций")
    print(*u)
    print(*P)
    plt.plot(grid_centers, P)
    #plt.plot(grid_centers, u)
    plt.show()


if __name__ == "__main__":
    SIMPLE()
