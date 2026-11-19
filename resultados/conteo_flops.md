# Reporte de Costo Operacional (FLOPs)

## 1. Convención de Conteo Declarada
- **1 FLOP** equivale a una operación aritmética elemental en punto flotante (`+`, `-`, `*`, `/`).
- **Reutilización de multiplicadores:** Cada multiplicador $m_{ik} = a_{ik}^{(k)} / a_{kk}^{(k)}$ se calcula **una única vez por fila** y se reutiliza en toda la fila y en $b$.
- **Barrido triangular:** En el paso $k$, las columnas previas ($j < k$) ya anuladas nunca se tocan. El barrido actualiza la fila desde la columna $k$ hasta $n$, más $b_i$.
- **Omisión algebraica del cero subdiagonal:** Si se asigna `A[i, k] = 0.0` sin evaluar `A[i, k] - m * A[k, k]`, se ahorran exactamente 66 multiplicaciones y 66 restas (pasando de **1486 FLOPs** a **1354 FLOPs**, que es el mínimo teórico).
- **Sobrecosto de estabilidad (pivoteo):** Las comparaciones, cocientes de escala $|a_{ik}|/s_i$ e intercambios de filas se registran en contadores separados.

## 2. Desglose Paso a Paso (escalado)

| Etapa / Paso | Mult. y Div. (`*`, `/`) | Sumas y Restas (`+`, `-`) | Total FLOPs | Comp. Pivoteo | Cocientes Escala | Intercambios |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| Paso 01 (col 1) | 154 | 143 | **297** | 11 | 12 | 1 |
| Paso 02 (col 2) | 130 | 120 | **250** | 10 | 11 | 0 |
| Paso 03 (col 3) | 108 | 99 | **207** | 9 | 10 | 1 |
| Paso 04 (col 4) | 88 | 80 | **168** | 8 | 9 | 1 |
| Paso 05 (col 5) | 70 | 63 | **133** | 7 | 8 | 1 |
| Paso 06 (col 6) | 54 | 48 | **102** | 6 | 7 | 1 |
| Paso 07 (col 7) | 40 | 35 | **75** | 5 | 6 | 1 |
| Paso 08 (col 8) | 28 | 24 | **52** | 4 | 5 | 1 |
| Paso 09 (col 9) | 18 | 15 | **33** | 3 | 4 | 1 |
| Paso 10 (col 10) | 10 | 8 | **18** | 2 | 3 | 0 |
| Paso 11 (col 11) | 4 | 3 | **7** | 1 | 2 | 0 |
| Sust. Regresiva | 78 | 66 | **144** | 0 | 0 | 0 |
| **TOTALES** | **782** | **704** | **1486** | **66** | **77** | **8** |

## 3. Comparación con Cotas de Referencia
- **Término dominante asintótico ($2n^3/3$ para $n=12$):** 1152.00 flops (solo eliminación de $A$).
- **Cota de referencia de la prueba (eliminación densa, barrido desde la columna $k$, incluye $b$ y sustitución regresiva):** 1486 flops; con pivoteo parcial, 1486 flops + 66 comparaciones.
- **Flops contados por el instrumento (esta corrida):** 1486 flops → desviación 0.00 % respecto de la cota (tolerancia ±10 %: [1337, 1635]).
- **Comparaciones de búsqueda de pivote contadas:** 66 (cota: 66); se reportan aparte como sobrecosto de estabilidad.
- **Mínimo estricto si la subdiagonal se asigna a 0 sin operar:** 716 mult./div. + 638 sumas/restas = **1354 flops** (fórmulas $n^3/3+n^2-n/3$ y $n^3/3+n^2/2-5n/6$); el barrido desde la columna $k$ agrega 132 flops (66 mult. + 66 restas). Ambos conteos están dentro de ±10 % de la cota.
- **Comparaciones iniciales para el vector de escalas $s_i$:** 132 comparaciones (aparte; no son flops aritméticos).
