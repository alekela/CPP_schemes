import numpy as np
import matplotlib.pyplot as plt


def simple_algorithm(nx, ny, lx, ly, mu, rho, u_inlet, max_iter, tolerance):
    # Grid spacing
    dx = lx / (nx - 1)
    dy = ly / (ny - 1)

    # Initialize fields
    u = np.zeros((ny, nx))  # x-velocity
    v = np.zeros((ny, nx))  # y-velocity
    p = np.zeros((ny, nx))  # Pressure
    u_star = np.zeros((ny, nx))  # Intermediate x-velocity
    v_star = np.zeros((ny, nx))  # Intermediate y-velocity
    p_prime = np.zeros((ny, nx))  # Pressure correction

    # Boundary conditions
    u[:, 0] = u_inlet  # Inlet velocity
    u[:, -1] = u[:, -2]  # Outlet: zero gradient
    u[0, :] = 0  # Bottom wall: no-slip
    u[-1, :] = 0  # Top wall: no-slip

    for iter in range(max_iter):
        # Step 1: Solve momentum equations for intermediate velocities
        for i in range(1, nx - 1):
            for j in range(1, ny - 1):
                u_star[j, i] = u[j, i] + (mu / rho) * (u[j, i + 1] - 2 * u[j, i] + u[j, i - 1]) / dx ** 2 + (
                            mu / rho) * (u[j + 1, i] - 2 * u[j, i] + u[j - 1, i]) / dy ** 2 - (
                                           u[j, i] * (u[j, i + 1] - u[j, i - 1]) / (2 * dx) + v[j, i] * (
                                               u[j + 1, i] - u[j - 1, i]) / (2 * dy))

                v_star[j, i] = v[j, i] + (mu / rho) * (v[j, i + 1] - 2 * v[j, i] + v[j, i - 1]) / dx ** 2 + (
                            mu / rho) * (v[j + 1, i] - 2 * v[j, i] + v[j - 1, i]) / dy ** 2 - (
                                           u[j, i] * (v[j, i + 1] - v[j, i - 1]) / (2 * dx) + v[j, i] * (
                                               v[j + 1, i] - v[j - 1, i]) / (2 * dy))

                # Step 2: Solve pressure correction equation
                for i in range(1, nx - 1):
                    for j in range(1, ny - 1):
                        p_prime[j, i] = (p_prime[j, i + 1] + p_prime[j, i - 1]) * dy ** 2 + (
                                    p_prime[j + 1, i] + p_prime[j - 1, i]) * dx ** 2 - (
                                                    rho * dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2))) * (
                                                    (u_star[j, i + 1] - u_star[j, i - 1]) / (2 * dx) + (
                                                        v_star[j + 1, i] - v_star[j - 1, i]) / (2 * dy))

                # Step 3: Correct velocities and pressure
                for i in range(1, nx - 1):
                    for j in range(1, ny - 1):
                        u[j, i] = u_star[j, i] - (p_prime[j, i + 1] - p_prime[j, i - 1]) * dx / (2 * rho)
                        v[j, i] = v_star[j, i] - (p_prime[j + 1, i] - p_prime[j - 1, i]) * dy / (2 * rho)
                        p[j, i] = p[j, i] + p_prime[j, i]

                # Check for convergence
                residual = np.linalg.norm(p_prime)
                print("Residual:", residual)
                if residual < tolerance:
                    print(f"Converged in {iter} iterations")
                    break
    return u, v, p


# Parameters
nx, ny = 41, 41  # Number of grid points
lx, ly = 1.0, 1.0  # Domain size
mu = 0.01  # Dynamic viscosity
rho = 1.0  # Density
u_inlet = 1.0  # Inlet velocity
max_iter = 1000  # Maximum iterations
tolerance = 1e-5  # Convergence tolerance
dx = lx / (nx - 1)
dy = ly / (ny - 1)
# Run SIMPLE algorithm
u, v, p = simple_algorithm(nx, ny, lx, ly, mu, rho, u_inlet, max_iter, tolerance)

fig, ax = plt.subplots(1, 1)

x = np.array([0 + i * dx for i in range(nx)])
y = np.array([0 + i * dy for i in range(ny)])
X, Y = np.meshgrid(x, y)

c = ax.pcolormesh(X, Y, u, cmap='jet')
fig.colorbar(c, ax=ax)

plt.show()
