import numpy as np
import matplotlib.pyplot as plt

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

grid = [xstart + j * dx for j in range(N + 1)]
P = [P0 for _ in range(N)]
u = [u0 for _ in range(N)]
u[0] = uleft

for t in range(1):
    iter = 0
    while True:
        new_u = [0 for _ in range(N)]
        for j in range(0, N):
            if j == 0:
                new_u[j] = u[j] + tau * (- (P[j + 1] - P[j]) / dx / rho0 - u[j] * (u[j + 1] - u[j]) / dx + nu * (
                            u[j + 2] - 2 * u[j + 1] + u[j]) / dx / dx)
            elif j == N - 1:
                new_u[j] = u[j] + tau * (- (P[j] - P[j - 1]) / dx / rho0 - + nu * (
                            u[j - 2] - u[j - 1]) / dx / dx)
            else:
                new_u[j] = u[j] + tau * (
                            - (P[j + 1] - P[j - 1]) / 2. / dx / rho0 - u[j] * (u[j + 1] - u[j - 1]) / 2 / dx + nu * (
                                u[j + 1] - 2 * u[j] + u[j - 1]) / dx / dx)

        coeffs_dp = [[0 for _ in range(N)] for _ in range(N)]
        res_dp = [0 for _ in range(N)]
        dP = [0 for _ in range(N)]
        for j in range(N):
            if j == 0:
                coeffs_dp[j][j] = -1 / dx**2
                coeffs_dp[j][j + 1] = 1 / dx**2
                res_dp[j] = rho0 / tau * (new_u[j + 1] - new_u[j]) / dx
            elif j == N - 1:
                coeffs_dp[j][j] = -2 / dx**2
                coeffs_dp[j][j - 1] = 1 / dx**2
                res_dp[j] = -pright
            else:
                coeffs_dp[j][j + 1] = 1 / dx**2
                coeffs_dp[j][j] = -2 / dx**2
                coeffs_dp[j][j - 1] = 1 / dx**2
                res_dp[j] = rho0 / tau * (new_u[j + 1] - new_u[j - 1]) / (2 * dx)

        dP = np.linalg.solve(coeffs_dp, res_dp)
        # Boundary condition P[N - 1] = pright and dP / dx [0] = 0

        for j in range(N):
            P[j] += dP[j]
        for j in range(N):
            if j == 0:
                u[j] = new_u[j] - tau / rho0 * (dP[j + 1] - dP[j]) / dx
            elif j == N - 1:
                u[j] = new_u[j] - tau / rho0 * (dP[j] - dP[j - 1]) / dx
            else:
                u[j] = new_u[j] - tau / rho0 * (dP[j + 1] - dP[j - 1]) / (2 * dx)

        # Boundary condition u[0] = uleft and du/dx[N-1] = 0
        # P[-1] = pright, dP/dx[0] = 0
        u[0] = uleft
        P[-1] = pright

        R = [0 for _ in range(N)]
        for j in range(N):
            if j == 0:
                R[j] = abs(u[j + 1] - u[j]) / dx
            elif j == N - 1:
                R[j] = abs(u[j] - u[j - 1]) / dx
            else:
                R[j] = abs(u[j + 1] - u[j - 1]) / (2 * dx)

        R_max = max(R)
        #print("U_starred:", new_u)
        print("Невязка: ", R_max)
        #print("Правая часть давления:", res_dp)
        print("Давление:", P)
        print("U:", u)
        print()
        iter += 1
        if R_max < prec:
            break
    print("Сошлось за", iter, "итераций")

print(*u)
print(*P)
fig = plt.figure()
ax = fig.add_subplot()
fig2 = plt.figure()
ax2 = fig2.add_subplot()
ax.set_title("P")
ax2.set_title("U")

grid_centers = [(grid[i + 1] + grid[i]) / 2 for i in range(N)]

ax.plot(grid_centers[1:-1], P[1:-1])
ax2.plot(grid_centers[1:-1], u[1:-1])
plt.show()
