import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- 1. Constantes ---
RT = 6.378160e6   # Radio Tierra (m)
d = 3.844e8       # Distancia Tierra–Luna (m)

RT_norm = RT / d  # Radio Tierra normalizado

# --- 2. Cargar datos ---
cohete_data = np.loadtxt("C:/Users/User/Documents/Fisica_compu/Compu/Obligatorio4 (Voluntario 3)/Obligario-4-Voluntario-3-/coheteaux.txt", delimiter=",")
luna_data = np.loadtxt("C:/Users/User/Documents/Fisica_compu/Compu/Obligatorio4 (Voluntario 3)/Obligario-4-Voluntario-3-/lunaaux.txt", delimiter=",")

x_cohete = cohete_data[:, 0]
y_cohete = cohete_data[:, 1]

x_luna = luna_data[:, 0]
y_luna = luna_data[:, 1]

# --- 3. Crear figura ---
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)

# Dibujar Tierra
earth = plt.Circle((0, 0), RT_norm, color='blue', label='Tierra')
ax.add_artist(earth)

# Inicializar línea y puntos
cohete_line, = ax.plot([], [], 'r-', label='Cohete')
cohete_dot, = ax.plot([], [], 'ro')
luna_dot, = ax.plot([], [], 'go', label='Luna')

# --- 4. Función de animación ---
def animate(i):
    cohete_line.set_data(x_cohete[:i], y_cohete[:i])
    cohete_dot.set_data([x_cohete[i]], [y_cohete[i]])  # CORREGIDO: listas []
    luna_dot.set_data([x_luna[i]], [y_luna[i]])        # CORREGIDO: listas []
    return cohete_line, cohete_dot, luna_dot

# --- 5. Animación ---
frames = min(len(x_cohete), len(x_luna))
ani = FuncAnimation(fig, animate, frames=frames, interval=100, blit=True)

# --- 6. Mostrar ---
ax.set_title("Trayectoria Cohete y Luna (normalizado)")
ax.legend()
plt.show()

ani.save("C:/Users/User/Documents/Fisica_compu/Compu/Obligatorio4 (Voluntario 3)/Obligario-4-Voluntario-3-/trayectoria_cohete_lunarugge_kutta_adaptado.mp4", fps=30, dpi=200)
