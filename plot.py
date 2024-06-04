import imageio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from tqdm import tqdm
import os
# Carregar os dados do arquivo example.txt
data = np.loadtxt("./results/15/example_2.50_0.30_0.30.dat").T
df = {
    'Tempo': data[0],
    "id": data[1].astype(int),
    "x": data[2],
    "y": data[3],
}
dados = pd.DataFrame(df)

# Configurações iniciais do plot
dt = 0.001
color = ["red", "blue", "green"]
m = 150 / 2
angulo = 15
fig, ax = plt.subplots(figsize=(8, 8))


# Definir a proporção de aspecto igual e ajustável
ax.set_aspect('equal', adjustable='datalim')
L1 = 3.5*7.5e-3*2
L2 = 154.e-3
# Elementos estáticos
lineas_estaticas = [
    ax.plot([0, 0], [0, 98.], color='blue')[0],
    ax.plot([L1*1000, L1*1000], [0, 98.], color='blue')[0],
    ax.plot([0, -154], [98., 98. + 154 * np.tan(np.pi * angulo / 180)], color='blue')[0],
    ax.plot([L1*1000, (L1+L2)*1000], [98., 98. + 154 * np.tan(np.pi * angulo / 180)], color='blue')[0],
    ax.plot([-154, -154], [98. + 154 * np.tan(np.pi * angulo / 180), 910], color='blue')[0],
    ax.plot([(L1+L2)*1000, (L1+L2)*1000], [98. + 154 * np.tan(np.pi * angulo / 180), 910], color='blue')[0]
]
time = np.unique(df['Tempo'])
# Criar o GIF processando um frame de cada vez
with imageio.get_writer('animation.gif', mode='I', fps=60) as writer:
    for t in tqdm(time):
        # Atualizar o gráfico para o tempo t
        data_int_time = dados[dados["Tempo"] == t]
        ax.clear()
        ax.set_xlim(-324 + m, 324 + m)
        ax.set_ylim(0, 910)
        ax.axis('off')
        # Replotar elementos estáticos
        for linea in lineas_estaticas:
            ax.add_line(linea)

        # Adicionar as partículas
        for _, row in data_int_time.iterrows():
            circ = plt.Circle((row['x'] * 1000, row['y'] * 1000), 7.5, color=color[int(row['id'] % 3)], fill=False)
            ax.add_artist(circ)
        # Salvar o frame atual como imagem
        plt.savefig('temp_frame.png')
        
        # Adicionar a imagem ao GIF
        with Image.open('temp_frame.png') as img:
            writer.append_data(np.array(img))
        
        # Limpar o arquivo temporário
        os.remove('temp_frame.png')

print("GIF criado com sucesso!")
