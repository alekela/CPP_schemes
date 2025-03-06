import matplotlib.pyplot as plt
import os


codename = "KrestSteel"
theory = False

with open(f"CPP_schemes\\CSVs\\{codename}\\Iter=0.csv") as f:
    labels = f.readline()
    data = f.readlines()

data = list(map(lambda x: x.split(','), data))

t_start = list(map(lambda x: float(x[0]), data))
xs_start = list(map(lambda x: float(x[1]), data))
rhos_start = list(map(lambda x: float(x[2]), data))
Ps_start = list(map(lambda x: float(x[3]), data))
us_start = list(map(lambda x: float(x[4]), data))


n = []
for file in os.listdir(f"CPP_schemes\\CSVs\\{codename}"):
    n.append(int(file.split("=")[1].split('.')[0]))

n = max(n)
with open(f"CPP_schemes\\CSVs\\{codename}\\Iter={n}.csv") as f:
    labels = f.readline()
    data2 = f.readlines()

data2 = list(map(lambda x: x.split(','), data2))

t_end = list(map(lambda x: float(x[0]), data2))
xs_end = list(map(lambda x: float(x[1]), data2))
rhos_end = list(map(lambda x: float(x[2]), data2))
Ps_end = list(map(lambda x: float(x[3]), data2))
us_end = list(map(lambda x: float(x[4]), data2))

if theory:
    with open("CPP_schemes\\CSVs\\{codename}\\Theory.csv") as f:
        labels = f.readline()
        data_theory = f.readlines()

    data_theory = list(map(lambda x: x.split(','), data_theory))

    t_theory = list(map(lambda x: float(x[0]), data_theory))
    xs_theory = list(map(lambda x: float(x[1]), data_theory))
    rhos_theory = list(map(lambda x: float(x[2]), data_theory))
    Ps_theory = list(map(lambda x: float(x[3]), data_theory))
    us_theory = list(map(lambda x: float(x[4]), data_theory))



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


if theory:
    ax_P.plot(xs_theory, Ps_theory, color='orange')
    ax_P.plot(xs_end, Ps_end, color='#1f77b4')
    ax_P.scatter(xs_end, Ps_end, s=7, color='#1f77b4')
    ax_P.legend([f"Scheme at t={t_end[0]}", f"Theory at t={t_theory[0]}"])
    
    ax_U.plot(xs_theory, us_theory, color='orange')
    ax_U.plot(xs_end, us_end, color='#1f77b4')
    ax_U.legend([f"Scheme at t={t_end[0]}", f"Theory at t={t_theory[0]}"])
    
    ax_rho.plot(xs_theory, rhos_theory,color='orange')
    ax_rho.plot(xs_end, rhos_end, color='#1f77b4')
    ax_rho.legend([f"Scheme at t={t_end[0]}", f"Theory at t={t_theory[0]}"])

else:
    ax_P.plot(xs_end, Ps_end, color='#1f77b4')
    ax_P.scatter(xs_end, Ps_end, s=7, color='#1f77b4')
    ax_P.legend([f"Scheme at t={t_end[0]}"])
    
    ax_U.plot(xs_end, us_end, color='#1f77b4')
    ax_U.legend([f"Scheme at t={t_end[0]}"])
    
    ax_rho.plot(xs_end, rhos_end, color='#1f77b4')
    ax_rho.legend([f"Scheme at t={t_end[0]}"])

plt.show()
