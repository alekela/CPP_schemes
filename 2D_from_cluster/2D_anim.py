import numpy as np
import matplotlib.pyplot as plt
import os
import imageio
from imageio.v2 import imread


end = "Makh"
# end = "Explosion"
names = [f"WENO_2D_{end}", f"Godunov_2D_{end}", f"GK_2D_{end}", f"GKR_2D_{end}"]
names = [f"Godunov_2D_{end}"]

for filename in names:
	steps = []
	sep = ','
	for file in os.listdir(filename):
		if file[:5] == "Iter=":
			# Получение номера итерации из имени файла
			step = int(file.split("=")[1].split('.')[0])
			steps.append(step)

			with open(os.path.join(filename, f"Iter={step}.csv")) as f :
				time = float(f.readline().split(":")[1])
				# time = float(f.readline())
				title = f.readline().split(sep)
				data = f.readlines()
			data = list(map(lambda x : list(map(float, x.split(sep))), data))
			x = np.array(list(map(lambda x : x[0], data)))
			y = np.array(list(map(lambda x : x[1], data)))
			P = np.array(list(map(lambda x : x[5], data)))
			ux = np.array(list(map(lambda x : x[2], data)))
			uy = np.array(list(map(lambda x : x[3], data)))
			rho = np.array(list(map(lambda x: x[4], data)))


			Lx = 1
			Ly = len(y)
			for i in range(1, len(x)) :
				if y[i] == y[0] :
					Ly = i
					break

			Lx = len(y) // Ly
			x = np.reshape(x, (Lx, Ly))[1:].T
			y = np.reshape(y, (Lx, Ly))[1:].T
			P = np.reshape(P, (Lx, Ly))[1:].T
			ux = np.reshape(ux, (Lx, Ly))[1:].T
			uy = np.reshape(uy, (Lx, Ly))[1:].T
			rho = np.reshape(rho, (Lx, Ly))[1:].T



			fig, ax = plt.subplots(1, 1)
			c = ax.pcolormesh(x, y, rho, cmap = 'jet')
			fig.colorbar(c, ax = ax)

			plt.savefig(f"Pics/2D_{filename}_Iter={step}.png")


	images = []
	steps.sort()
	for step in steps :
		images.append(imread(os.path.join("Pics", f"2D_{filename}_Iter={step}.png")))
	imageio.mimsave(f"Res_{filename}.gif", images)
	print(f"Анимация сохранена в Res_{filename}.gif.")

