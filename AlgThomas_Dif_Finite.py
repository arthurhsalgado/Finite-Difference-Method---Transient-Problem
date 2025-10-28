import numpy as np
import matplotlib.pyplot as plt

# --- Parâmetros ---
alpha = 5e-6      # difusividade térmica (m²/s)
dx = 0.0014         # passo espacial (m)
dt = 0.01        # passo de tempo (s)
h = 1100          # coeficiente convectivo (W/m².K)
k = 30            # condutividade térmica (W/m.K)
#q_dot = 0    # geração de calor (W/m³)
q_dot = 1e7       # geração de calor (W/m³)
T_inf = 250       # temperatura ambiente (°C)
n = 7             # número de nós
tempo_total = 1  # tempo total de simulação (s)

# --- Parâmetros adimensionais ---
Fo = (alpha * dt) / dx**2
Bi = (h * dx) / k

# --- Solução analítica (regime permanente) ---
A = np.zeros((n, n))
B = np.zeros((n, 1))

for i in range(n):
    if i == 0:  # isolado
        A[i, i] = 1
        A[i, i + 1] = -1
    elif i == n - 1:  # convectivo
        A[i, i - 1] = -1
        A[i, i] = 1 + Bi
        B[i] = Bi * T_inf + (q_dot * dx**2) / k
    else:  # interior
        A[i, i - 1] = -1
        A[i, i] = 2
        A[i, i + 1] = -1
        B[i] = (q_dot * dx**2) / k

T_steady = np.linalg.solve(A, B)  # solução estacionária

# --- Usa a solução estacionária como condição inicial ---
T_prev = T_steady.copy()
n_passos = int(tempo_total / dt)

T_hist = np.zeros((n_passos, n))
T_hist[0, :] = T_prev.flatten()



# --- Loop temporal ---
for passo in range(n_passos):

    # Monta matriz A (tridiagonal) para o método implícito
    A = np.array([
        [(1+2*Fo), -2*Fo, 0, 0, 0, 0, 0],
        [-Fo, (1+2*Fo), -Fo, 0, 0, 0, 0],
        [0, -Fo, (1+2*Fo), -Fo, 0, 0, 0],
        [0, 0, -Fo, (1+2*Fo), -Fo, 0, 0],
        [0, 0, 0, -Fo, (1+2*Fo), -Fo, 0],
        [0, 0, 0, 0, -Fo, (1+2*Fo), -Fo],
        [0, 0, 0, 0, 0, -2*Fo, (1+2*Fo+2*Fo*Bi)]
    ], dtype=float)

    # Monta vetor B
    B = np.zeros((n, 1))
    for i in range(n-1):
        B[i] = T_prev[i] + q_dot * dx**2 / k
    B[n-1] = T_prev[n-1] + Fo*(2*Bi*T_inf) + q_dot * dx**2 / k

    # --- Eliminação progressiva ---
    for i in range(1, n):
        fator = A[i, i-1] / A[i-1, i-1]
        A[i, i-1] = 0
        A[i, i] -= fator * A[i-1, i]
        B[i] -= fator * B[i-1]

    # --- Substituição regressiva ---
    x = np.zeros((n, 1))
    x[-1] = B[-1] / A[-1, -1]
    for i in range(n-2, -1, -1):
        x[i] = (B[i] - A[i, i+1]*x[i+1]) / A[i, i]

    # Atualiza temperatura
    T_prev = x.copy()
    T_hist[passo, :] = T_prev.flatten()

# --- Resultado final ---
print("\n=== Resultado final (transiente iniciado em regime permanente) ===")
for i, Ti in enumerate(T_prev, start=1):
    print(f"T{i} = {Ti[0]:.4f} °C")


# --- Gráfico ---
tempos = np.linspace(0, tempo_total, n_passos)
for i in range(n):
    plt.plot(tempos, T_hist[:, i], label=f"Nó {i+1}")

plt.xlabel("Tempo (s)")
plt.ylabel("Temperatura (°C)")
plt.title("Evolução da temperatura nos nós")
plt.legend()
plt.grid(True)
plt.show()

# --- Gráfico 2: perfil espacial (T(x)) em um instante específico ---
posicoes = np.arange(0, n*dx, dx)
tempo_escolhido = -1  # último instante
plt.plot(posicoes, T_hist[tempo_escolhido, :], 'o-', color='red')
plt.xlabel("Posição ao longo da barra (m)")
plt.ylabel("Temperatura (°C)")
plt.title(f"Distribuição de temperatura na barra (t = {tempo_total:.2f} s)")
plt.grid(True)
plt.show()