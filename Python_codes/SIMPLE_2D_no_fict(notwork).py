import numpy as np
import matplotlib.pyplot as plt


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
    tau = 0.0000001

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
    # Boundary
    for i in range(Ny):
        ux[i][0] = uxleft
        uy[i][0] = uyleft
        P[i][-1] = pright
    for j in range(Nx):
        uy[0][j] = 0
        uy[-1][j] = 0

    for t in range(1):
        iter = 0
        flag = True
        while (flag and iter < 20):
            new_ux = [[0 for _ in range(Nx)] for _ in range(Ny)]
            new_uy = [[0 for _ in range(Nx)] for _ in range(Ny)]
            for i in range(Ny):
                for j in range(Nx):
                    if i == 0 and j == 0:
                        new_ux[i][j] = ux[i][j] + tau * (-ux[i][j] * (ux[i][j + 1] - ux[i][j]) / dx
                                                         - uy[i][j] * (ux[i + 1][j] - ux[i][j]) / dy
                                                         - (P[i][j + 1] - P[i][j]) / (dx) +
                                                         nu * (ux[i][j + 2] - 2 * ux[i][j + 1] + ux[i][j]) / dx ** 2 +
                                                         nu * (ux[i + 2][j] - 2 * ux[i + 1][j] + ux[i][j]) / dy ** 2)
                        new_uy[i][j] = uy[i][j] + tau * (-ux[i][j] * (uy[i][j + 1] - uy[i][j]) / (dx)
                                                         - uy[i][j] * (uy[i + 1][j] - uy[i][j]) / (dy)
                                                         - (P[i + 1][j] - P[i][j]) / (dy)
                                                         + nu * (uy[i][j + 2] - 2 * uy[i][j + 1] + uy[i][j]) / dx ** 2
                                                         + nu * (uy[i + 2][j] - 2 * uy[i + 1][j] + uy[i][j]) / dy ** 2)
                    elif i == 0 and j == Nx - 1:
                        new_ux[i][j] = ux[i][j] + tau * (-ux[i][j] * (ux[i][j] - ux[i][j - 1]) / (dx)
                                                         - uy[i][j] * (ux[i + 1][j] - ux[i][j]) / (dy)
                                                         - (P[i][j] - P[i][j - 1]) / (dx) +
                                                         nu * (ux[i][j] - 2 * ux[i][j - 1] + ux[i][j - 2]) / dx ** 2 +
                                                         nu * (ux[i + 2][j] - 2 * ux[i + 1][j] + ux[i][j]) / dy ** 2)
                        new_uy[i][j] = uy[i][j] + tau * (-ux[i][j] * (uy[i][j] - uy[i][j - 1]) / (dx)
                                                         - uy[i][j] * (uy[i + 1][j] - uy[i][j]) / (dy)
                                                         - (P[i + 1][j] - P[i][j]) / (dy)
                                                         + nu * (uy[i][j] - 2 * uy[i][j - 1] + uy[i][j - 2]) / dx ** 2
                                                         + nu * (uy[i + 2][j] - 2 * uy[i + 1][j] + uy[i][j]) / dy ** 2)
                    elif i == Ny - 1 and j == 0:
                        new_ux[i][j] = ux[i][j] + tau * (-ux[i][j] * (ux[i][j + 1] - ux[i][j]) / (dx)
                                                         - uy[i][j] * (ux[i][j] - ux[i - 1][j]) / (dy)
                                                         - (P[i][j + 1] - P[i][j]) / (dx) +
                                                         nu * (ux[i][j + 2] - 2 * ux[i][j + 1] + ux[i][j]) / dx ** 2 +
                                                         nu * (ux[i][j] - 2 * ux[i - 1][j] + ux[i - 2][j]) / dy ** 2)
                        new_uy[i][j] = uy[i][j] + tau * (-ux[i][j] * (uy[i][j + 1] - uy[i][j]) / (dx)
                                                         - uy[i][j] * (uy[i][j] - uy[i - 1][j]) / (dy)
                                                         - (P[i][j] - P[i - 1][j]) / (dy)
                                                         + nu * (uy[i][j + 2] - 2 * uy[i][j + 1] + uy[i][j]) / dx ** 2
                                                         + nu * (uy[i][j] - 2 * uy[i - 1][j] + uy[i - 2][j]) / dy ** 2)
                    elif i == Ny - 1 and j == Nx - 1:
                        new_ux[i][j] = ux[i][j] + tau * (-ux[i][j] * (ux[i][j] - ux[i][j - 1]) / (dx)
                                                         - uy[i][j] * (ux[i][j] - ux[i - 1][j]) / (dy)
                                                         - (P[i][j] - P[i][j - 1]) / (dx) +
                                                         nu * (ux[i][j] - 2 * ux[i][j - 1] + ux[i][j - 2]) / dx ** 2 +
                                                         nu * (ux[i][j] - 2 * ux[i - 1][j] + ux[i - 2][j]) / dy ** 2)
                        new_uy[i][j] = uy[i][j] + tau * (-ux[i][j] * (uy[i][j] - uy[i][j - 1]) / (dx)
                                                         - uy[i][j] * (uy[i][j] - uy[i - 1][j]) / (dy)
                                                         - (P[i][j] - P[i - 1][j]) / (dy)
                                                         + nu * (uy[i][j] - 2 * uy[i][j - 1] + uy[i][j - 2]) / dx ** 2
                                                         + nu * (uy[i][j] - 2 * uy[i - 1][j] + uy[i - 2][j]) / dy ** 2)
                    elif i == 0:
                        new_ux[i][j] = ux[i][j] + tau * (-ux[i][j] * (ux[i][j + 1] - ux[i][j - 1]) / (2 * dx)
                                                         - uy[i][j] * (ux[i + 1][j] - ux[i][j]) / (dy)
                                                         - (P[i][j + 1] - P[i][j - 1]) / (2 * dx) +
                                                         nu * (ux[i][j + 1] - 2 * ux[i][j] + ux[i][j - 1]) / dx ** 2 +
                                                         nu * (ux[i + 2][j] - 2 * ux[i + 1][j] + ux[i][j]) / dy ** 2)
                        new_uy[i][j] = uy[i][j] + tau * (-ux[i][j] * (uy[i][j + 1] - uy[i][j - 1]) / (2 * dx)
                                                         - uy[i][j] * (uy[i + 1][j] - uy[i][j]) / (dy)
                                                         - (P[i + 1][j] - P[i][j]) / (dy)
                                                         + nu * (uy[i][j + 1] - 2 * uy[i][j] + uy[i][j - 1]) / dx ** 2
                                                         + nu * (uy[i + 2][j] - 2 * uy[i + 1][j] + uy[i][j]) / dy ** 2)
                    elif i == Ny - 1:
                        new_ux[i][j] = ux[i][j] + tau * (-ux[i][j] * (ux[i][j + 1] - ux[i][j - 1]) / (2 * dx)
                                                         - uy[i][j] * (ux[i][j] - ux[i - 1][j]) / (dy)
                                                         - (P[i][j + 1] - P[i][j - 1]) / (2 * dx) +
                                                         nu * (ux[i][j + 1] - 2 * ux[i][j] + ux[i][j - 1]) / dx ** 2 +
                                                         nu * (ux[i][j] - 2 * ux[i - 1][j] + ux[i - 2][j]) / dy ** 2)
                        new_uy[i][j] = uy[i][j] + tau * (-ux[i][j] * (uy[i][j + 1] - uy[i][j - 1]) / (2 * dx)
                                                         - uy[i][j] * (uy[i][j] - uy[i - 1][j]) / (dy)
                                                         - (P[i][j] - P[i - 1][j]) / (dy)
                                                         + nu * (uy[i][j + 1] - 2 * uy[i][j] + uy[i][j - 1]) / dx ** 2
                                                         + nu * (uy[i][j] - 2 * uy[i - 1][j] + uy[i - 2][j]) / dy ** 2)
                    elif j == 0:
                        new_ux[i][j] = ux[i][j] + tau * (-ux[i][j] * (ux[i][j + 1] - ux[i][j]) / (dx)
                                                         - uy[i][j] * (ux[i + 1][j] - ux[i - 1][j]) / (2 * dy)
                                                         - (P[i][j + 1] - P[i][j]) / (dx) +
                                                         nu * (ux[i][j + 2] - 2 * ux[i][j + 1] + ux[i][j]) / dx ** 2 +
                                                         nu * (ux[i + 1][j] - 2 * ux[i][j] + ux[i - 1][j]) / dy ** 2)
                        new_uy[i][j] = uy[i][j] + tau * (-ux[i][j] * (uy[i][j + 1] - uy[i][j]) / (dx)
                                                         - uy[i][j] * (uy[i + 1][j] - uy[i - 1][j]) / (2 * dy)
                                                         - (P[i + 1][j] - P[i - 1][j]) / (2 * dy)
                                                         + nu * (uy[i][j + 2] - 2 * uy[i][j + 1] + uy[i][j]) / dx ** 2
                                                         + nu * (uy[i + 1][j] - 2 * uy[i][j] + uy[i - 1][j]) / dy ** 2)
                    elif j == Nx - 1:
                        new_ux[i][j] = ux[i][j] + tau * (-ux[i][j] * (ux[i][j] - ux[i][j - 1]) / (dx)
                                                         - uy[i][j] * (ux[i + 1][j] - ux[i - 1][j]) / (2 * dy)
                                                         - (P[i][j] - P[i][j - 1]) / (dx) +
                                                         nu * (ux[i][j] - 2 * ux[i][j - 1] + ux[i][j - 2]) / dx ** 2 +
                                                         nu * (ux[i + 1][j] - 2 * ux[i][j] + ux[i - 1][j]) / dy ** 2)
                        new_uy[i][j] = uy[i][j] + tau * (-ux[i][j] * (uy[i][j] - uy[i][j - 1]) / (dx)
                                                         - uy[i][j] * (uy[i + 1][j] - uy[i - 1][j]) / (2 * dy)
                                                         - (P[i + 1][j] - P[i - 1][j]) / (2 * dy)
                                                         + nu * (uy[i][j] - 2 * uy[i][j - 1] + uy[i][j - 2]) / dx ** 2
                                                         + nu * (uy[i + 1][j] - 2 * uy[i][j] + uy[i - 1][j]) / dy ** 2)
                    else:
                        new_ux[i][j] = ux[i][j] + tau * (-ux[i][j] * (ux[i][j + 1] - ux[i][j - 1]) / (2 * dx)
                                                         - uy[i][j] * (ux[i + 1][j] - ux[i - 1][j]) / (2 * dy)
                                                         - (P[i][j + 1] - P[i][j - 1]) / (2 * dx) +
                                                         nu * (ux[i][j + 1] - 2 * ux[i][j] + ux[i][j - 1]) / dx ** 2 +
                                                         nu * (ux[i + 1][j] - 2 * ux[i][j] + ux[i - 1][j]) / dy ** 2)
                        new_uy[i][j] = uy[i][j] + tau * (-ux[i][j] * (uy[i][j + 1] - uy[i][j - 1]) / (2 * dx)
                                                         - uy[i][j] * (uy[i + 1][j] - uy[i - 1][j]) / (2 * dy)
                                                         - (P[i + 1][j] - P[i - 1][j]) / (2 * dy)
                                                         + nu * (uy[i][j + 1] - 2 * uy[i][j] + uy[i][j - 1]) / dx ** 2
                                                         + nu * (uy[i + 1][j] - 2 * uy[i][j] + uy[i - 1][j]) / dy ** 2)

            coeffs_dp = [[0 for _ in range(Nx * Ny)] for _ in range(Nx * Ny)]
            res_dp = [0 for _ in range(Nx * Ny)]

            for i in range(Ny):
                for j in range(Nx):
                    # уравнение круга: (x - 1.5) ** 2 + (y - 0.5) ** 2 <= 0.5**2
                    if (grid_centers[i][j][0] - 1.5) ** 2 + (grid_centers[i][j][1] - 0.5) ** 2 <= 0.25 ** 2:
                        rho = rhoshar
                    else:
                        rho = rho0
                    if i == 0 and j == 0:
                        coeffs_dp[i * Nx + j][i * Nx + j + 1] = 1 / dx / dx
                        coeffs_dp[i * Nx + j][(i + 1) * Nx + j] = 1 / dy / dy
                        coeffs_dp[i * Nx + j][i * Nx + j] = -(1 / dx / dx + 1 / dy / dy)
                        res_dp[i * Nx + j] = rho / tau * ((ux[i][j + 1] - ux[i][j]) / dx
                                                          + (uy[i + 1][j] - uy[i][j]) / dy)
                    elif i == 0 and j == Nx - 1:
                        coeffs_dp[i * Nx + j][i * Nx + j - 1] = 1 / dx / dx
                        coeffs_dp[i * Nx + j][(i + 1) * Nx + j] = 1 / dy / dy
                        coeffs_dp[i * Nx + j][i * Nx + j] = -(2 / dx / dx + 1 / dy / dy)
                        res_dp[i * Nx + j] = rho / tau * ((ux[i][j] - ux[i][j - 1]) / dx
                                                          + (uy[i + 1][j] - uy[i][j]) / dy) - pright / dx**2
                    elif i == Ny - 1 and j == 0:
                        coeffs_dp[i * Nx + j][i * Nx + j + 1] = 1 / dx / dx
                        coeffs_dp[i * Nx + j][(i - 1) * Nx + j] = 1 / dy / dy
                        coeffs_dp[i * Nx + j][i * Nx + j] = -(1 / dx / dx + 1 / dy / dy)
                        res_dp[i * Nx + j] = rho / tau * ((ux[i][j + 1] - ux[i][j]) / (dx)
                                                          + (uy[i][j] - uy[i - 1][j]) / (dy))
                    elif i == Ny - 1 and j == Nx - 1:
                        coeffs_dp[i * Nx + j][i * Nx + j - 1] = 1 / dx / dx
                        coeffs_dp[i * Nx + j][(i - 1) * Nx + j] = 1 / dy / dy
                        coeffs_dp[i * Nx + j][i * Nx + j] = -(2 / dx / dx + 1 / dy / dy)
                        res_dp[i * Nx + j] = rho / tau * ((ux[i][j] - ux[i][j - 1]) / (dx)
                                                          + (uy[i][j] - uy[i - 1][j]) / (dy)) - pright / dx**2
                    elif i == 0:
                        coeffs_dp[i * Nx + j][i * Nx + j + 1] = 1 / dx / dx
                        coeffs_dp[i * Nx + j][i * Nx + j - 1] = 1 / dx / dx
                        coeffs_dp[i * Nx + j][(i + 1) * Nx + j] = 1 / dy / dy
                        coeffs_dp[i * Nx + j][i * Nx + j] = -(2 / dx / dx + 1 / dy / dy)
                        res_dp[i * Nx + j] = rho / tau * ((ux[i][j + 1] - ux[i][j - 1]) / (2 * dx)
                                                          + (uy[i + 1][j] - uy[i][j]) / (dy))
                    elif i == Ny - 1:
                        coeffs_dp[i * Nx + j][i * Nx + j + 1] = 1 / dx / dx
                        coeffs_dp[i * Nx + j][i * Nx + j - 1] = 1 / dx / dx
                        coeffs_dp[i * Nx + j][(i - 1) * Nx + j] = 1 / dy / dy
                        coeffs_dp[i * Nx + j][i * Nx + j] = -(2 / dx / dx + 1 / dy / dy)
                        res_dp[i * Nx + j] = rho / tau * ((ux[i][j + 1] - ux[i][j - 1]) / (2 * dx)
                                                          + (uy[i][j] - uy[i - 1][j]) / (dy))
                    elif j == 0:
                        coeffs_dp[i * Nx + j][i * Nx + j + 1] = 1 / dx / dx
                        coeffs_dp[i * Nx + j][(i + 1) * Nx + j] = 1 / dy / dy
                        coeffs_dp[i * Nx + j][(i - 1) * Nx + j] = 1 / dy / dy
                        coeffs_dp[i * Nx + j][i * Nx + j] = -(1 / dx / dx + 2 / dy / dy)
                        res_dp[i * Nx + j] = rho / tau * ((ux[i][j + 1] - ux[i][j]) / (dx)
                                                          + (uy[i + 1][j] - uy[i - 1][j]) / (2 * dy))
                    elif j == Nx - 1:
                        coeffs_dp[i * Nx + j][i * Nx + j - 1] = 1 / dx / dx
                        coeffs_dp[i * Nx + j][(i + 1) * Nx + j] = 1 / dy / dy
                        coeffs_dp[i * Nx + j][(i - 1) * Nx + j] = 1 / dy / dy
                        coeffs_dp[i * Nx + j][i * Nx + j] = -2 * (1 / dx / dx + 1 / dy / dy)
                        res_dp[i * Nx + j] = rho / tau * ((ux[i][j] - ux[i][j - 1]) / (dx)
                                                          + (uy[i + 1][j] - uy[i - 1][j]) / (2 * dy)) - pright / dx**2
                    else:
                        coeffs_dp[i * Nx + j][i * Nx + j + 1] = 1 / dx / dx
                        coeffs_dp[i * Nx + j][i * Nx + j - 1] = 1 / dx / dx
                        coeffs_dp[i * Nx + j][(i + 1) * Nx + j] = 1 / dy / dy
                        coeffs_dp[i * Nx + j][(i - 1) * Nx + j] = 1 / dy / dy
                        coeffs_dp[i * Nx + j][i * Nx + j] = -2 * (1 / dx / dx + 1 / dy / dy)
                        res_dp[i * Nx + j] = rho / tau * ((ux[i][j + 1] - ux[i][j - 1]) / (2 * dx)
                                                          + (uy[i + 1][j] - uy[i - 1][j]) / (2 * dy))

            dP = np.linalg.solve(coeffs_dp, res_dp)
            dP_m = [[0 for _ in range(Nx)] for _ in range(Ny)]
            for i in range(Ny):
                for j in range(Nx):
                    dP_m[i][j] = dP[i * Nx + j]

            for i in range(Ny):
                for j in range(Nx):
                    P[i][j] += dP_m[i][j]
                    # уравнение круга: (x - 1.5) ** 2 + (y - 0.5) ** 2 <= 0.5**2
                    if (grid_centers[i][j][0] - 1.5) ** 2 + (grid_centers[i][j][1] - 0.5) ** 2 <= 0.25 ** 2:
                        rho = rhoshar
                    else:
                        rho = rho0

                    if i == 0 and j == 0:
                        ux[i][j] = new_ux[i][j] - tau / rho * (dP_m[i][j + 1] - dP_m[i][j]) / (dx)
                        uy[i][j] = new_uy[i][j] - tau / rho * (dP_m[i + 1][j] - dP_m[i][j]) / (dy)
                    elif i == 0 and j == Nx - 1:
                        ux[i][j] = new_ux[i][j] - tau / rho * (dP_m[i][j] - dP_m[i][j - 1]) / (dx)
                        uy[i][j] = new_uy[i][j] - tau / rho * (dP_m[i + 1][j] - dP_m[i][j]) / (dy)
                    elif i == Ny - 1 and j == 0:
                        ux[i][j] = new_ux[i][j] - tau / rho * (dP_m[i][j + 1] - dP_m[i][j]) / (dx)
                        uy[i][j] = new_uy[i][j] - tau / rho * (dP_m[i][j] - dP_m[i - 1][j]) / (dy)
                    elif i == Ny - 1 and j == Nx - 1:
                        ux[i][j] = new_ux[i][j] - tau / rho * (dP_m[i][j] - dP_m[i][j - 1]) / (dx)
                        uy[i][j] = new_uy[i][j] - tau / rho * (dP_m[i][j] - dP_m[i - 1][j]) / (dy)
                    elif i == 0:
                        ux[i][j] = new_ux[i][j] - tau / rho * (dP_m[i][j + 1] - dP_m[i][j - 1]) / (2 * dx)
                        uy[i][j] = new_uy[i][j] - tau / rho * (dP_m[i + 1][j] - dP_m[i][j]) / (dy)
                    elif i == Ny - 1:
                        ux[i][j] = new_ux[i][j] - tau / rho * (dP_m[i][j + 1] - dP_m[i][j - 1]) / (2 * dx)
                        uy[i][j] = new_uy[i][j] - tau / rho * (dP_m[i][j] - dP_m[i - 1][j]) / (dy)
                    elif j == 0:
                        ux[i][j] = new_ux[i][j] - tau / rho * (dP_m[i][j + 1] - dP_m[i][j]) / (dx)
                        uy[i][j] = new_uy[i][j] - tau / rho * (dP_m[i + 1][j] - dP_m[i - 1][j]) / (2 * dy)
                    elif j == Nx - 1:
                        ux[i][j] = new_ux[i][j] - tau / rho * (dP_m[i][j] - dP_m[i][j - 1]) / (dx)
                        uy[i][j] = new_uy[i][j] - tau / rho * (dP_m[i + 1][j] - dP_m[i - 1][j]) / (2 * dy)
                    else:
                        ux[i][j] = new_ux[i][j] - tau / rho * (dP_m[i][j + 1] - dP_m[i][j - 1]) / (2 * dx)
                        uy[i][j] = new_uy[i][j] - tau / rho * (dP_m[i + 1][j] - dP_m[i - 1][j]) / (2 * dy)

            for i in range(Ny):
                ux[i][0] = uxleft
                uy[i][0] = uyleft
                P[i][-1] = pright
            for j in range(Nx):
                uy[0][j] = 0
                uy[-1][j] = 0

            iter += 1
            R = [[0 for _ in range(Nx)] for _ in range(Ny)]
            for i in range(Ny):
                for j in range(Nx):
                    if i == 0 and j == 0:
                        R[i][j] = abs(
                            (ux[i][j + 1] - ux[i][j]) / (dx) + (uy[i + 1][j] - uy[i][j]) / (dy))
                    elif i == 0 and j == Nx - 1:
                        R[i][j] = abs(
                            (ux[i][j] - ux[i][j - 1]) / (dx) + (uy[i + 1][j] - uy[i][j]) / (dy))
                    elif i == Ny - 1 and j == 0:
                        R[i][j] = abs(
                            (ux[i][j + 1] - ux[i][j]) / (dx) + (uy[i][j] - uy[i - 1][j]) / (dy))
                    elif i == Ny - 1 and j == Nx - 1:
                        R[i][j] = abs(
                            (ux[i][j] - ux[i][j - 1]) / (dx) + (uy[i][j] - uy[i - 1][j]) / (dy))
                    elif i == 0:
                        R[i][j] = abs(
                            (ux[i][j + 1] - ux[i][j - 1]) / (2 * dx) + (uy[i + 1][j] - uy[i][j]) / (dy))
                    elif i == Ny - 1:
                        R[i][j] = abs(
                            (ux[i][j + 1] - ux[i][j - 1]) / (2 * dx) + (uy[i][j] - uy[i - 1][j]) / (dy))
                    elif j == 0:
                        R[i][j] = abs(
                            (ux[i][j + 1] - ux[i][j]) / (dx) + (uy[i + 1][j] - uy[i - 1][j]) / (2 * dy))
                    elif j == Nx - 1:
                        R[i][j] = abs(
                            (ux[i][j] - ux[i][j - 1]) / (dx) + (uy[i + 1][j] - uy[i - 1][j]) / (2 * dy))
                    else:
                        R[i][j] = abs(
                            (ux[i][j + 1] - ux[i][j - 1]) / (2 * dx) + (uy[i + 1][j] - uy[i - 1][j]) / (2 * dy))

            R = max(list(map(lambda x: max(x), R)))
            print("Невязка:", R)
            print()
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
    print(X, Y, ux)
    print(ux)
    c = ax.pcolormesh(X, Y, ux, cmap='jet')
    fig.colorbar(c, ax=ax)

    plt.show()


if __name__ == "__main__":
    SIMPLE()
