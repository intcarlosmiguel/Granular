import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import imageio.v2 as imageio  # Importando a versão 2 especificamente
from tqdm import tqdm

# Carregar dados do arquivo TXT
data = np.loadtxt("./example.txt").T
df = {
    'Tempo': data[0],
    "id": data[1].astype(int),
    "x": data[2],
    "y": data[3],
}
dados = pd.DataFrame(df)
# Configurações iniciais do plot
fig, ax = plt.subplots(figsize= (8,8),dpi = 300)
ax.set_xlim(-304,304)
ax.set_ylim(0, 910)

# Lista para armazenar os frames
filenames = []
dt = 0.25/2
t = 0
count = 0
color = ["red","blue","green"]
m = 150/2
while(t < dados["Tempo"].max()):
    ax.clear()
    ax.set_xlim(-324-m,324+m)
    ax.set_ylim(0, 910)

    ax.plot([0,0],[0,9.8], color='blue')
    ax.plot([150,150],[0,9.8], color='blue')

    ax.plot([0,-154],[9.8,9.8 + 154*np.tan(np.pi*30/180)], color='blue')
    ax.plot([150,304],[9.8,9.8 + 154*np.tan(np.pi*30/180)], color='blue')
    ax.plot([-154,-154],[9.8 + 154*np.tan(np.pi*30/180),910], color='blue')
    ax.plot([304,304],[9.8 + 154*np.tan(np.pi*30/180),910], color='blue')

    data_int_time = dados[dados["Tempo"] == t]
    for particula in data_int_time[["id","x",'y']].values:
        centro_do_circulo = (particula[1], particula[2])
        circulo = plt.Circle(centro_do_circulo, 7.5, color=color[int(particula[0]%3)], fill=False, label='Círculo')
        ax.add_artist(circulo)
    ax.axis('off')
    # Salvar o frame
    ax.text(0.5, 0.95, f'Time: {t:.2f}', transform=ax.transAxes, fontsize=12,verticalalignment='center', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    filename = f'./frames/images/frame_{count}.png'
    plt.gca().set_aspect('equal', adjustable='box')  # Assegurar proporção de aspecto igual
    plt.subplots_adjust(left=0.0, bottom=0.1, right=1.0, top=0.9)
    #plt.tight_layout()
    plt.savefig(filename)
    filenames.append(filename)
    print(t,count)
    if(t >20):
        break
    t += dt
    count += 1

# Criar GIF

with imageio.get_writer('./frames/particula_movimento.gif', mode='I', duration=0.05) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
        # Remover os arquivos de frame para limpeza
        #os.remove(filename)

plt.close()