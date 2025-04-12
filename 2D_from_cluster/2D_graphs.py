import numpy as np
import matplotlib.pyplot as plt
import os


names = ["WENO_2D_p", "Godunov_2D_p", "GK_2D_p", "GKR_2D_p"]
for filename in names:
	max_num = -1
	for file in os.listdir(filename):
		max_num = max(max_num, int(file.split("=")[1].split('.')[0]))

	sep = ','
	with open(os.path.join(filename, f"Iter={max_num}.csv")) as f:
		time = float(f.readline().split(":")[1])
		# time = float(f.readline())
		title = f.readline().split(sep)
		data = f.readlines()
	data = list(map(lambda x: list(map(float, x.split(sep))), data))
	x = np.array(list(map(lambda x: x[0], data)))
	y = np.array(list(map(lambda x: x[1], data)))
	P = np.array(list(map(lambda x: x[5], data)))
	Lx = 1
	Ly = len(y)
	for i in range(1, len(x)):
		if y[i] == y[0]:
			Ly = i
			break

	Lx = len(y) // Ly
	x = np.reshape(x, (Lx, Ly)).T
	y = np.reshape(y, (Lx, Ly)).T
	P = np.reshape(P, (Lx, Ly)).T

	fig, ax = plt.subplots(1, 1)
	c = ax.pcolormesh(x, y, P, cmap='jet')
	fig.colorbar(c, ax=ax)

	plt.savefig(f"Pics/2D_{filename}.png")

	fig2, ax2 = plt.subplots(1, 1)
	ax2.plot(x[5], P[5])
	plt.savefig(f"Pics/1D_{filename}.png")
