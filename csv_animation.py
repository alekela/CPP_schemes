import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import subprocess
import imageio


# Функция для построения графиков и сохранения их изображений
def plot_and_save(data, iteration, filename_prefix):
    time = data['Time'][0]
    if int(iteration) < 1000:
        iteration = f'0{iteration}'

    plt.figure(figsize=(15, 5))

    # Rho(X)
    plt.subplot(1, 3, 1)
    plt.plot(data['X'], data['Rho'], label='Rho')
    plt.xlabel('X')
    plt.ylabel('Rho')

    # P(X)
    plt.subplot(1, 3, 2)
    plt.plot(data['X'], data['P'], label='P')
    plt.title(f'Iteration = {iteration}, time = {time}')
    plt.xlabel('X')
    plt.ylabel('P')

    # U(X)
    plt.subplot(1, 3, 3)
    plt.plot(data['X'], data['U'], label='U')
    plt.xlabel('X')
    plt.ylabel('U')

    # Сохранение изображения
    plt.tight_layout()
    plt.savefig(f'CPP_schemes\\PNGs\\{filename_prefix}\\Iter={iteration}.png')
    plt.close()


codename = "KrestDetonate"
# Путь к папке с CSV-файлами
data_folder = f'CPP_schemes\\CSVs\\{codename}\\'

# Поиск всех CSV-файлов по шаблону
csv_files = glob.glob(os.path.join(data_folder, f'Iter=*.csv'))

# Проверка, что найдены файлы
if not csv_files:
    print("CSV-файлы не найдены в папке", os.path.join(data_folder, f'Iter=*.csv'))
else:
    # Проход по каждому найденному файлу
    for csv_file in csv_files:
        # Чтение данных из CSV-файла
        data = pd.read_csv(csv_file)

        # Получение номера итерации из имени файла
        iteration = os.path.basename(csv_file).split('=')[1].split('.')[0]

        # Построение графиков и сохранение их изображений
        plot_and_save(data, iteration, codename)

    print("Графики сохранены.")


    video_path = f'CPP_schemes\\Videos\\{codename}_animation.gif'
    image_path = f'CPP_schemes\\PNGs\\{codename}'
    images = []
    for file in os.listdir(image_path):
        images.append(imageio.v2.imread(os.path.join(image_path, file)))
    imageio.mimsave(video_path, images)
    '''# Команда ffmpeg для создания видео
    cmd = [
        'ffmpeg',
        '-framerate', '1',  # Частота кадров в секунду
        '-i', image_path,
        '-c:v', 'libx264',
        '-pix_fmt', 'yuv420p',
        video_path
    ]

    # Выполнение команды
    subprocess.run(cmd, check=True)
    '''
    print(f"Анимация сохранена в файл {video_path}.")
