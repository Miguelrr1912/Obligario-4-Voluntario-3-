import numpy as np
import matplotlib.pyplot as plt

# 1) Leer los datos de Hprima.txt
# Suponiendo que tienes una columna de números
Hprima = np.loadtxt('C:/Users/User/Documents/Fisica_compu/Compu/Obligatorio4 (Voluntario 3)/Obligario-4-Voluntario-3-/Hprima.txt')

# 2) Crear vector de pasos (índice de cada dato)
pasos = np.arange(len(Hprima))

# 3) Graficar
plt.figure(figsize=(8, 5))
plt.plot(pasos, Hprima, marker='o', linestyle='-', color='blue', label="H'")

plt.xlabel('Paso')
plt.ylabel("Hamiltoniano modificado H'")
plt.title("Evolución de H' con Runge-Kutta Adaptativo")
plt.legend()
plt.grid(True)

# 4) Mostrar o guardar
plt.savefig('C:/Users/User/Documents/Fisica_compu/Compu/Obligatorio4 (Voluntario 3)/Obligario-4-Voluntario-3-/Hprima_grafica.png')  # Guarda la imagen como PNG
plt.show()
