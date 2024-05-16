import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.animation import FuncAnimation

# Carregar os dados do arquivo example.txt
data = np.loadtxt("./example.txt").T
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
m = 150/2
angulo = 45
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-324 + m, 324 + m)
ax.set_ylim(0, 910)
ax.axis('off')

# Definir a proporção de aspecto igual e ajustável
ax.set_aspect('equal', adjustable='datalim')

# Elementos estáticos
lineas_estaticas = [
    ax.plot([0, 0], [0, 98.], color='blue')[0],
    ax.plot([150, 150], [0, 98.], color='blue')[0],
    ax.plot([0, -154], [98., 98. + 154 * np.tan(np.pi * angulo / 180)], color='blue')[0],
    ax.plot([150, 304], [98., 98. + 154 * np.tan(np.pi * angulo / 180)], color='blue')[0],
    ax.plot([-154, -154], [98. + 154 * np.tan(np.pi * angulo / 180), 910], color='blue')[0],
    ax.plot([304, 304], [98. + 154 * np.tan(np.pi * angulo / 180), 910], color='blue')[0]
]

particulas = []  # Lista para armazenar os círculos das partículas

def init():
    # Inicializa círculos das partículas para garantir que eles existam
    for _ in range(len(dados['id'].unique())):
        circulo = plt.Circle((0, 0), 7.5, color='blue', fill=False)  # Cores e outros detalhes podem ser ajustados
        ax.add_artist(circulo)
        particulas.append(circulo)
    return particulas

def update(t):
    data_int_time = dados[dados["Tempo"] == t]
    for particula, circ in zip(data_int_time[['id', 'x', 'y']].values, particulas):
        circ.center = (particula[1]*1000, particula[2]*1000)
        circ.set_color(color[int(particula[0] % 3)])
    return particulas

ani = FuncAnimation(fig, update, frames=np.arange(dados["Tempo"].min(), dados["Tempo"].max() + dt, dt), init_func=init, blit=True, repeat=False)

# Aplicar tight_layout para garantir que a animação não fique distorcida
plt.tight_layout()

ani.save('animation.gif', writer='imagemagick', fps=200)
