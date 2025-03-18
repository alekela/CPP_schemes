import numpy as np
import matplotlib.pyplot as plt


def FictBound(A):
    A[0] = A[1]
    A[-1] = A[-2]


xstart = 0
xend = 1
N = 8
fict = 1
N_all = N + 2 * fict
dx = (xend - xstart) / N

rho0 = 1
nu = 0.01
u0 = 0
P0 = 0
uleft = 1
pright = 0
prec = 1e-5
tau = 0.001

grid = [xstart + (j - fict) * dx for j in range(N_all + 1)]
P = [P0 for _ in range(N_all)]
u = [u0 for _ in range(N_all)]
u[1] = uleft
P[-2] = pright
FictBound(u)
FictBound(P)

for t in range(1):
    iter = 0
    while iter < 5:
        new_u = [0 for _ in range(N_all)]
        for j in range(fict, N_all - fict):
            new_u[j] = u[j] + tau * (- (P[j + 1] - P[j - 1]) / 2. / dx / rho0 - u[j] * (u[j] - u[j - 1]) / dx + nu * (
                        u[j + 1] - 2 * u[j] + u[j - 1]) / dx / dx)

        new_u[0] = uleft
        print("U_starred:", new_u)
        coeffs_dp = [[0 for _ in range(N)] for _ in range(N)]
        res_dp = [0 for _ in range(N)]
        dP = [0 for _ in range(N_all)]
        for j in range(N):
            jc = j + 1
            if j == 0:
                coeffs_dp[j][j] = -2
                coeffs_dp[j][j + 1] = 1
                res_dp[j] = rho0 * dx / 2 / tau * (new_u[jc + 1] - new_u[jc - 1])
            elif j == N - 1:
                coeffs_dp[j][j] = -2
                coeffs_dp[j][j - 1] = 1
                res_dp[j] = -pright
            else:
                coeffs_dp[j][j + 1] = 1
                coeffs_dp[j][j] = -2
                coeffs_dp[j][j - 1] = 1
                res_dp[j] = rho0 * dx / 2 / tau * (new_u[jc + 1] - new_u[jc - 1])

        dP = np.linalg.solve(coeffs_dp, res_dp)
        # Boundary condition P[N - 1] = pright and dP / dx [0] = 0
        dP[-1] = 0
        dP[0] = dP[1]

        dP = [0] + list(dP) + [0]
        FictBound(dP)
        for j in range(fict, N_all - fict):
            P[j] += dP[j]
            u[j] = new_u[j] - tau / rho0 / 2 / dx * (dP[j + 1] - dP[j - 1])
        # Boundary condition u[0] = uleft and du/dx[N-1] = 0
        u[1] = uleft
        u[-2] = u[-3]
        P[-2] = pright
        P[1] = P[2]
        FictBound(u)
        FictBound(P)

        R = max([abs(u[j + 1] - u[j - 1]) / 2 / dx for j in range(1, N - 1)])
        print("Невязка: ", R)
        print("Правая часть давления:", res_dp)
        print("U:", u)
        print()
        iter += 1
        if R < prec:
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

grid_centers = [(grid[i + 1] + grid[i]) / 2 for i in range(N_all)]

ax.plot(grid_centers[1:-1], P[1:-1])
ax2.plot(grid_centers[1:-1], u[1:-1])
plt.show()

