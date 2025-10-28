# Método de Diferenças Finitas aplicado a um Problema de Condução de Calor em Regime transiente

O presente trabalho consiste na aplicação do Algoritmo de Thomas para a resolução de um problema de condução de calor em regime transiente.

A figura ilustra o esquemático geral do problema.

![Problema térmico](https://github.com/arthurhsalgado/Finite-Difference-Method---Transient-Problem/blob/main/Problema.png?raw=true)

Com isso, aplicando a formulação de diferenças finitas, tem-se o seguinte problema linear:

![Problema térmico](https://github.com/arthurhsalgado/Finite-Difference-Method---Transient-Problem/blob/main/Formulacao_Matricial.png?raw=true)

Então, considerando os seguintes parâmetros:

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

A distribuição final de temperatura torna-se:

![Problema térmico](https://github.com/arthurhsalgado/Finite-Difference-Method---Transient-Problem/blob/main/Distribuicao_TempBarra.png?raw=true)

A evolução temporal de temperatura torna-se:


![Problema térmico](https://github.com/arthurhsalgado/Finite-Difference-Method---Transient-Problem/blob/main/Temperaturas_Tempo.png?raw=true)

