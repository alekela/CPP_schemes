import numpy as np
import matplotlib.pyplot as plt

ps = list(range(1, 13))
names = ["WENO_2D_p", "Godunov_2D_p", "GK_2D_p", "GKR_2D_p"]
res = {}
for p in ps:
	with open(f"out/out_p_{p}.out") as f:
		data = f.readlines()
	data = list(filter(lambda x: len(x) > 1, map(lambda x: x.split(), data)))
	data_d = {}
	for i in data:
		data_d[i[0]] = float(i[2])
	for name in names:
		if name not in res:
			res[name] = [data_d[name]]
		else:
			res[name].append(data_d[name])

for name in names:
	y = np.array(res[name])
	y = y[0] / y
	
	fig, ax = plt.subplots(1, 1)
	plt.plot(ps, ps)
	plt.plot(ps, y)
	plt.legend(["real", "ideal"])
	plt.grid()
	plt.xlim(0, max(ps))
	plt.ylim(0, max(ps))

	plt.savefig(f"Pics/speed_up_{name}.png")

fig, ax = plt.subplots(1, 1)
plt.plot(ps, ps)

for name in names:
	y = np.array(res[name])
	y = y[0] / y
	
	plt.plot(ps, y)

plt.legend(["ideal"] + names)
plt.grid()
plt.xlim(0, max(ps))
plt.ylim(0, max(ps))

plt.savefig(f"Pics/speed_up_all.png")
