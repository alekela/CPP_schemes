import pandas as pd
import matplotlib.pyplot as plt
import os
import imageio
from imageio.v2 import imread


def plot_and_save(x, time, rho, P, u, iteration, im_path):
    '''Функция для построения графиков и сохранения их изображений'''
    plt.figure(figsize=(15, 5))

    # Rho(X)
    plt.subplot(1, 3, 1)
    plt.plot(x, rho, label='Rho')
    plt.xlabel('X')
    plt.ylabel('Rho')

    # P(X)
    plt.subplot(1, 3, 2)
    plt.plot(x, P, label='P')
    plt.title(f'Iteration = {iteration}, time = {time}')
    plt.xlabel('X')
    plt.ylabel('P')

    # U(X)
    plt.subplot(1, 3, 3)
    plt.plot(x, u, label='U')
    plt.xlabel('X')
    plt.ylabel('U')

    # Сохранение изображения
    # plt.tight_layout()
    plt.savefig(f'{im_path}\\Iter={iteration}.png')
    plt.close()


codename = "KrestDetonate"
# Путь к папке с CSV-файлами
data_folder = f'CPP_schemes\\CSVs\\{codename}\\'
video_path = f'CPP_schemes\\Videos\\{codename}_animation.gif'
image_path = f'CPP_schemes\\PNGs\\{codename}'

steps = []
for file in os.listdir(data_folder):
    if file[:5] == "Iter=":
        # Получение номера итерации из имени файла
        step = int(file.split("=")[1].split('.')[0])
        steps.append(step)
        
        data = pd.read_csv(os.path.join(data_folder, file))
        plot_and_save(data["X"], data['Time'][0], data["Rho"], data["P"], data["U"], step, image_path)
print("Графики сохранены.")

images = []
steps.sort()
for step in steps:
    images.append(imread(os.path.join(image_path, f"Iter={step}.png")))
imageio.mimsave(video_path, images)
print(f"Анимация сохранена в файл {video_path}.")
