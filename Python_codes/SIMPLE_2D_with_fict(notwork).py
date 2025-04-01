import numpy as np
import matplotlib.pyplot as plt


def Boundary(Nx, Ny, P, ux, uy, uxleft, pright):
    for i in range(Ny):
        ux[i][0] = uxleft
        ux[i][1] = uxleft
        uy[i][0] = 0
        uy[i][1] = 0
        P[i][-1] = pright
        P[i][-2] = pright
        ux[i][-1] = uxleft # 0
        ux[i][-2] = uxleft # 0
        uy[i][-1] = 0
        uy[i][-2] = 0
    for j in range(Nx):
        ux[0][j] = uxleft # 0
        ux[1][j] = uxleft # 0
        uy[0][j] = 0
        uy[1][j] = 0
        ux[-1][j] = uxleft # 0
        ux[-2][j] = uxleft # 0
        uy[-1][j] = 0
        uy[-2][j] = 0


def SIMPLE():
    xstart = 0
    xend = 5
    ystart = 0
    yend = 1
    Nx = 25
    Ny = 25
    dx = (xend - xstart) / Nx
    dy = (yend - ystart) / Ny

    rho0 = 1
    rhoshar = 1
    nu = 0.01
    u0 = 0
    P0 = 0
    uxleft = 1
    uyleft = 0
    pright = 0
    prec = 1e-5
    tau = 0.001

    x = np.array([xstart + i * dx for i in range(Nx + 1)])
    y = np.array([ystart + i * dy for i in range(Ny + 1)])
    X, Y = np.meshgrid(x, y)
    grid = [[(X[i][j], Y[i][j]) for i in range(Nx + 1)] for j in range(Ny + 1)]
    grid_centers = [
        [((grid[j + 1][i + 1][0] + grid[j][i][0]) / 2, (grid[j + 1][i + 1][1] + grid[j][i][1]) / 2) for j in range(Ny)]
        for i in range(Nx)]

    P = [[P0 for _ in range(Nx)] for _ in range(Ny)]
    ux = [[u0 for _ in range(Nx)] for _ in range(Ny)]
    uy = [[u0 for _ in range(Nx)] for _ in range(Ny)]
    Boundary(Nx, Ny, P, ux, uy, uxleft, pright)

    for t in range(1):
        iter = 0
        flag = True
        while (flag and iter < 100):
            dux = [[0 for _ in range(Nx)] for _ in range(Ny)];
            duy = [[0 for _ in range(Nx)] for _ in range(Ny)];
            for i in range(1, Ny - 1):
                for j in range(1, Nx - 1):
                    dux[i][j] = ux[i][j] + tau * (
                                -ux[i][j] * (ux[i][j] - ux[i][j - 1]) / dx - uy[i][j] * (ux[i][j] - ux[i - 1][j]) / dy
                                + nu * (ux[i][j + 1] - 2 * ux[i][j] + ux[i][j - 1]) / dx / dx
                                + nu * (ux[i + 1][j] - 2 * ux[i][j] + ux[i - 1][j]) / dy / dy)
                    duy[i][j] = uy[i][j] + tau * (
                                -ux[i][j] * (uy[i][j] - uy[i][j - 1]) / dx - uy[i][j] * (uy[i][j] - uy[i - 1][j]) / dy
                                + nu * (uy[i][j + 1] - 2 * uy[i][j] + uy[i][j - 1]) / dx / dx
                                + nu * (uy[i + 1][j] - 2 * uy[i][j] + uy[i - 1][j]) / dy / dy)

            Boundary(Nx, Ny, P, dux, duy, uxleft, pright)

            coeffs_dp = [[0 for _ in range((Nx - 2) * (Ny - 2))] for _ in range((Nx - 2) * (Ny - 2))]
            res_dp = [0 for _ in range((Nx - 2) * (Ny - 2))]

            for i in range(Ny - 2):
                for j in range(Nx - 2):
                    ic = i + 1
                    jc = j + 1
                    # уравнение круга: (x - 1.5) ** 2 + (y - 0.5) ** 2 <= 0.5**2
                    if (grid_centers[ic][jc][0] - 1.5) ** 2 + (grid_centers[ic][jc][1] - 0.5) ** 2 <= 0.25 ** 2:
                        rho = rhoshar
                    else:
                        rho = rho0
                    if i == 0 and j == 0:
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j + 1] = 1 / dx / dx
                        coeffs_dp[i * (Nx - 2) + j][(i + 1) * (Nx - 2) + j] = 1 / dy / dy
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j] = -(1 / dx / dx + 1 / dy / dy)
                        res_dp[i * (Nx - 2) + j] = rho / tau * ((ux[ic][jc + 1] - ux[ic][jc - 1]) / 2 / dx)
                    elif i == Ny - 3 and j == Ny - 3:
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j - 1] = 1 / dx / dx
                        coeffs_dp[i * (Nx - 2) + j][(i - 1) * (Nx - 2) + j] = 1 / dy / dy
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j] = -(1 / dx / dx + 1 / dy / dy)
                        res_dp[i * (Ny - 2) + j] = rho / tau * ((ux[ic][jc + 1] - ux[ic][jc - 1]) / 2 / dx)
                    elif i == Ny - 3 and j == 0:
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j + 1] = 1 / dx / dx
                        coeffs_dp[i * (Nx - 2) + j][(i - 1) * (Nx - 2) + j] = 1 / dy / dy
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j] = -(1 / dx / dx + 1 / dy / dy)
                        res_dp[i * (Nx - 2) + j] = rho / tau * ((ux[ic][jc + 1] - ux[ic][jc - 1]) / 2 / dx)
                    elif i == 0 and j == Ny - 3:
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j - 1] = 1 / dx / dx
                        coeffs_dp[i * (Nx - 2) + j][(i + 1) * (Nx - 2) + j] = 1 / dy / dy
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j] = -(1 / dx / dx + 1 / dy / dy)
                        res_dp[i * (Nx - 2) + j] = 0
                    elif i == 0:
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j + 1] = 1 / dx / dx
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j - 1] = 1 / dx / dx
                        coeffs_dp[i * (Nx - 2) + j][(i + 1) * (Nx - 2) + j] = 1 / dy / dy
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j] = -(2 / dx / dx + 1 / dy / dy)
                        res_dp[i * (Nx - 2) + j] = rho / tau * ((ux[ic][jc + 1] - ux[ic][jc - 1]) / 2 / dx)
                    elif j == 0:
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j + 1] = 1 / dx / dx
                        coeffs_dp[i * (Nx - 2) + j][(i + 1) * (Nx - 2) + j] = 1 / dy / dy
                        coeffs_dp[i * (Nx - 2) + j][(i - 1) * (Nx - 2) + j] = 1 / dy / dy
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j] = -(1 / dx / dx + 2 / dy / dy)
                        res_dp[i * (Nx - 2) + j] = rho / tau * ((ux[ic][jc + 1] - ux[ic][jc - 1]) / 2 / dx + (
                                    uy[ic + 1][jc] - uy[ic - 1][jc]) / 2 / dy)

                    elif i == Ny - 3:
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j + 1] = 1 / dx / dx
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j - 1] = 1 / dx / dx
                        coeffs_dp[i * (Nx - 2) + j][(i - 1) * (Nx - 2) + j] = 1 / dy / dy
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j] = -(2 / dx / dx + 1 / dy / dy)
                        res_dp[i * (Nx - 2) + j] = rho / tau * ((ux[ic][jc + 1] - ux[ic][jc - 1]) / 2 / dx)
                    elif j == Nx - 3:
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j - 1] = 1 / dx / dx
                        coeffs_dp[i * (Nx - 2) + j][(i + 1) * (Nx - 2) + j] = 1 / dy / dy
                        coeffs_dp[i * (Nx - 2) + j][(i - 1) * (Nx - 2) + j] = 1 / dy / dy
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j] = -(1 / dx / dx + 2 / dy / dy)
                        res_dp[i * (Nx - 2) + j] = rho / tau * ((uy[ic + 1][jc] - uy[ic - 1][jc]) / 2 / dy)

                    elif (i != 0 and i != Ny - 3) and (j != 0 and j != Nx - 3):
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j + 1] = 1 / dx / dx
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j - 1] = 1 / dx / dx
                        coeffs_dp[i * (Nx - 2) + j][(i + 1) * (Nx - 2) + j] = 1 / dy / dy
                        coeffs_dp[i * (Nx - 2) + j][(i - 1) * (Nx - 2) + j] = 1 / dy / dy
                        coeffs_dp[i * (Nx - 2) + j][i * (Nx - 2) + j] = -2 * (1 / dx / dx + 1 / dy / dy)
                        res_dp[i * (Nx - 2) + j] = rho / tau * ((ux[ic][jc + 1] - ux[ic][jc - 1]) / 2 / dx + (uy[ic + 1][jc] - uy[ic - 1][jc]) / 2 / dy)

            tmp_dP = np.linalg.solve(coeffs_dp, res_dp)
            tmp_dP_matrix = [[0 for _ in range(Nx)] for _ in range(Ny)]
            for i in range(Ny - 2):
                for j in range(Nx - 2):
                    ic = i + 1
                    jc = j + 1
                    tmp_dP_matrix[ic][jc] = tmp_dP[i * (Ny - 2) + j]

            # for i in range(Ny):
            #    tmp_dP_matrix[i][0] = tmp_dP_matrix[i][1]
            #    tmp_dP_matrix[i][-1] = 0
            # здесь еще гр условия сверху и снизу, но давления нулевые, так что да

            new_P = [[0 for _ in range(Nx)] for _ in range(Ny)];
            new_ux = [[0 for _ in range(Nx)] for _ in range(Ny)];
            new_uy = [[0 for _ in range(Nx)] for _ in range(Ny)]
            for i in range(1, Ny - 1):
                for j in range(1, Nx - 1):
                    new_P[i][j] = P[i][j] + tmp_dP_matrix[i][j]
                    # уравнение круга: (x - 1.5) ** 2 + (y - 0.5) ** 2 <= 0.5**2
                    if (grid_centers[i][j][0] - 1.5) ** 2 + (grid_centers[i][j][1] - 0.5) ** 2 <= 0.25 ** 2:
                        rho = rhoshar
                    else:
                        rho = rho0

                    new_ux[i][j] = dux[i][j] - tau / rho / 2 / dx * (tmp_dP_matrix[i][j + 1] - tmp_dP_matrix[i][j - 1])
                    new_uy[i][j] = duy[i][j] - tau / rho / 2 / dy * (tmp_dP_matrix[i + 1][j] - tmp_dP_matrix[i - 1][j])

            Boundary(Nx, Ny, new_P, new_ux, new_uy, uxleft, pright)

            for i in range(1, Ny - 1):
                for j in range(1, Nx - 1):
                    P[i][j] = new_P[i][j]
                    ux[i][j] = new_ux[i][j]
                    uy[i][j] = new_uy[i][j]
            iter += 1
            R = [[abs(ux[i][j + 1] - ux[i][j - 1] / 2 / dx + uy[i + 1][j] - ux[i - 1][j] / 2 / dy) for i in
                  range(1, Ny - 1)] for j in range(1, Nx - 1)]
            R = max(list(map(lambda x: max(x), R)))
            if R < prec:
                flag = False
        print("Сошлось за", iter, "итераций")
    for i in range(Ny):
        for j in range(Nx):
            if (grid_centers[i][j][0] - 1.5) ** 2 + (grid_centers[i][j][1] - 0.5) ** 2 <= 0.25 ** 2:
                P[i][j] = -np.inf
                ux[i][j] = -np.inf
                uy[i][j] = -np.inf
    print("P:")
    for i in P[::-1]:
        print(*i)
    print()
    print("ux:")
    for i in ux[::-1]:
        print(*i)
    print("uy:")
    for i in uy[::-1]:
        print(*i)
    print()
    fig, ax = plt.subplots(1, 1)

    c = ax.pcolormesh(X, Y, ux, cmap='jet')
    fig.colorbar(c, ax=ax)

    plt.show()
    # plt.plot(grid_centers, P)
    # plt.plot(grid_centers, u)
    # plt.show()


if __name__ == "__main__":
    SIMPLE()
