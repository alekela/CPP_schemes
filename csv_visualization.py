import matplotlib.pyplot as plt


with open("CPP_codes\\first_step_new.csv") as f:
    labels = f.readline()
    data = f.readlines()

data = list(map(lambda x: x.split(','), data))

t_start = list(map(lambda x: float(x[0]), data))
xs_start = list(map(lambda x: float(x[1]), data))
rhos_start = list(map(lambda x: float(x[2]), data))
Ps_start = list(map(lambda x: float(x[3]), data))
us_start = list(map(lambda x: float(x[4]), data))


with open("CPP_codes\\last_step_new.csv") as f:
    labels = f.readline()
    data2 = f.readlines()

data2 = list(map(lambda x: x.split(','), data2))

t_end = list(map(lambda x: float(x[0]), data2))
xs_end = list(map(lambda x: float(x[1]), data2))
rhos_end = list(map(lambda x: float(x[2]), data2))
Ps_end = list(map(lambda x: float(x[3]), data2))
us_end = list(map(lambda x: float(x[4]), data2))


fig_P = plt.figure()
fig_U = plt.figure()
fig_rho = plt.figure()

ax_P = fig_P.add_subplot()
ax_U = fig_U.add_subplot()
ax_rho = fig_rho.add_subplot()
ax_P.grid()
ax_P.set_xlabel("x")
ax_P.set_ylabel("P")
ax_U.grid()
ax_U.set_xlabel("x")
ax_U.set_ylabel("U")
ax_rho.grid()
ax_rho.set_xlabel("x")
ax_rho.set_ylabel("rho")


ax_P.plot(xs_start, Ps_start)
ax_P.plot(xs_end, Ps_end)

ax_U.plot(xs_start, us_start)
ax_U.plot(xs_end, us_end)

ax_rho.plot(xs_start, rhos_start)
ax_rho.plot(xs_end, rhos_end)

plt.show()
