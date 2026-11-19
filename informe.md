# Informe técnico 

## 1. Identificación

- **Estudiante(s):** Jose Herrera, Manuel Soto — Sección: 1
- **Asignatura:** Computación Numérica (INF-321), Ingeniería Civil Informática, UCM
- **Fecha:** 21, noviembre de 2026
- **Repositorio:** https://github.com/ch1mv21/cn-parcial-gauss12--Herrera-Soto-
- **Declaración de uso de IA:** se usó un asistente de IA (Claude, de Anthropic) para revisar el código frente a la rúbrica, auditar los dígitos significativos afirmados y redactar borradores de `informe.md` y `readme.md`. 

## 2. Planteamiento

Se resuelve `Ax = b` con A densa de 12×12 (144 valores no nulos, ingresados con todas sus cifras) y b con 20 dígitos significativos, sin truncar. Tres rasgos lo hacen costoso:

1. **Escalas de fila incompatibles.** El mayor valor absoluto por fila va de ≈ 10⁻⁴ (F1, F5 y F9, entre 8.05·10⁻⁵ y 9.96·10⁻⁵) a ≈ 10⁶ (F6, 9.67·10⁵). Al restar un múltiplo de una fila grande a una pequeña, los bits bajos del operando pequeño salen de los 53 bits de mantisa.
2. **Pivote inicial minúsculo.** a₁₁ = 1.3407·10⁻⁸ y a₂₁/a₁₁ ≈ 8.9·10¹¹.
3. **Casi dependencia lineal.** F7 ≈ 4.173·F3 + 582.7·F5 con ‖δ‖∞ = 2·10⁻⁴ (el ajuste por mínimos cuadrados da α = 4.177, β = 583.25). Esto eleva κ₂(A) ≈ 3.08·10¹³ sin hacer singular la matriz (det ≈ 1.36·10⁹).

## 3. Estrategia de eliminación

- **Orden:** 11 pasos de eliminación hacia adelante (columna k = 1…11) y sustitución regresiva, con **pivoteo parcial escalado**: los factores `sᵢ = maxⱼ|aᵢⱼ|` se calculan una sola vez (132 comparaciones) y en cada paso se elige la fila que maximiza `|aᵢₖ|/sᵢ`.
- **Por qué no excede el mínimo:** cada multiplicador `mᵢₖ = aᵢₖ/aₖₖ` se calcula **una sola vez** y se reutiliza en toda la fila y en b; el barrido del paso k va de la columna k en adelante (nunca desde la 1); no se calcula la inversa; hay **una sola factorización** `PA = LU`, que se reutiliza en el refinamiento iterativo (O(n²), sin refactorizar).
- Las dos corridas (sin pivoteo y escalada) usan el mismo código (`gauss.py`, parámetro `estrategia`).

## 4. Conteo de operaciones

**Convención:** 1 flop = 1 operación aritmética elemental (`+`, `−`, `×`, `÷`). Incluye la actualización de b y la sustitución regresiva. Los conteos salen del contador `ContadorOperaciones` (`conteo.py`), no de una fórmula. Comparaciones, cocientes de escala e intercambios se registran aparte como sobrecosto de estabilidad. Corrida con pivoteo escalado:

| Etapa | Mult./Div. | Sumas/Restas | Total | Comp. pivote | Cocientes escala | Intercambios |
|:--|--:|--:|--:|--:|--:|--:|
| Paso 01 (col 1) | 154 | 143 | 297 | 11 | 12 | 1 |
| Paso 02 (col 2) | 130 | 120 | 250 | 10 | 11 | 0 |
| Paso 03 (col 3) | 108 | 99 | 207 | 9 | 10 | 1 |
| Paso 04 (col 4) | 88 | 80 | 168 | 8 | 9 | 1 |
| Paso 05 (col 5) | 70 | 63 | 133 | 7 | 8 | 1 |
| Paso 06 (col 6) | 54 | 48 | 102 | 6 | 7 | 1 |
| Paso 07 (col 7) | 40 | 35 | 75 | 5 | 6 | 1 |
| Paso 08 (col 8) | 28 | 24 | 52 | 4 | 5 | 1 |
| Paso 09 (col 9) | 18 | 15 | 33 | 3 | 4 | 1 |
| Paso 10 (col 10) | 10 | 8 | 18 | 2 | 3 | 0 |
| Paso 11 (col 11) | 4 | 3 | 7 | 1 | 2 | 0 |
| Sust. Regresiva | 78 | 66 | 144 | 0 | 0 | 0 |
| TOTALES | 782 | 704 | 1486 | 66 | 77 | 8 |

- **Total: 1486 flops** (782 mult./div. + 704 sumas/restas). La cota de la prueba para eliminación densa es ≈ 1486: desviación 0.00 % (dentro de ±10 %).
- **Sobrecosto de estabilidad (aparte):** 66 comparaciones de búsqueda (cota: 66), 132 comparaciones iniciales para los `sᵢ`, 77 cocientes de escala y 8 intercambios de filas.
- El barrido desde la columna k incluye operar el elemento que se anula (132 flops: 66 × + 66 −). Si se asigna 0 directamente, una segunda corrida instrumentada mide el mínimo estricto de 1354 flops (716 + 638); ambos valores están dentro de ±10 % de la cota. El término dominante 2n³/3 vale 1152 para n = 12 y solo cubre la eliminación de A.

## 5. Análisis de mantisa y pivoteo

| Corrida | Multiplicador máx. | Factor de crecimiento | Error relativo ‖·‖₂ | Bits perdidos (de 53) | ‖r‖∞ |
|:--|--:|--:|--:|--:|--:|
| Sin pivoteo | 2.75·10¹³ | 2.64·10³ | 2.38·10⁻¹⁰ | 20.03 | 1.77·10⁻⁷ |
| Pivoteo parcial escalado | 2.47·10⁹ | 1.78 | 2.70·10⁻¹³ | 10.25 | 2.29·10⁻¹⁰ |

(Bits perdidos = log₂(error relativo / ε), con ε = 2.22·10⁻¹⁶.)

**Por qué a(1,1) = 1.34·10⁻⁸ es inaceptable como pivote.** Usarlo obliga a multiplicadores de hasta 10¹³ (m₂₁ ≈ 8.9·10¹¹). Cada fila queda con términos del orden de m·a₁ⱼ que dominan los datos originales, y como dos operandos que difieren en k órdenes de magnitud pierden ≈ 3.32·k bits de la mantisa del menor, la información original queda desplazada fuera de los 53 bits. El factor de crecimiento sube a 2.6·10³ y se pierden 20 bits (≈ 6 dígitos) frente a 10 bits con pivoteo.

**Por qué difieren las corridas.** El pivoteo escalado reduce el multiplicador máximo en cuatro órdenes de magnitud y mantiene el crecimiento en 1.78, de modo que los productos m·aₖⱼ no superan el tamaño de los datos que reciben. Nota: el criterio escalado elige pivotes pequeños en términos absolutos (p. ej. 8.8·10⁻⁵ en el paso 3, con multiplicador máximo 2.5·10⁹) porque su fila también tiene escala pequeña; lo que importa es el cociente |aᵢₖ|/sᵢ, y por eso el crecimiento se mantiene controlado.

## 6. Condicionamiento

- κ₂(A) ≈ 3.08·10¹³, κ∞(A) ≈ 5.37·10¹³ (SVD de NumPy, solo como estimación). Cota clásica: ‖δx‖/‖x‖ ≲ κ·ε ≈ 6.8·10⁻³, que no garantiza ni dos dígitos.
- **Resultado real:** 2.7·10⁻¹³ con pivoteo (≈ 12 dígitos), unas 10¹⁰ veces por debajo de la cota.
- **Por qué la cota es pesimista aquí:** (i) κ(A) está inflado por las escalas de fila: la matriz con cada fila dividida por su máximo tiene κ₂ ≈ 8.1·10³, y el pivoteo escalado trabaja de hecho sobre ese sistema equilibrado; (ii) la cota supone la dirección más desfavorable de la perturbación, que los errores de redondeo no alcanzan; (iii) b se entrega con 20 dígitos, así que el único error de entrada es el redondeo a doble precisión. La conclusión de dígitos recuperables sale del residuo y del error medido contra la solución exacta, no de la cota.

## 7. Resultados

Solución con pivoteo escalado. «Dígitos afirmados» = dígitos cuyo error es ≤ media unidad del último dígito, medido contra la solución exacta obtenida con `fractions.Fraction`; el valor reportado coincide con la solución exacta redondeada a esos dígitos (verificado para las 12 incógnitas).

| Incógnita | Valor reportado | Dígitos afirmados | Error relativo | Dígitos sin pivoteo |
|:--|:--|--:|--:|--:|
| x1 | 0.738491562700 | 12 | 7.25e-14 | 10 |
| x2 | -1.20593741680 | 12 | 5.56e-13 | 9 |
| x3 | 0.004387195200 | 10 | 8.30e-11 | 7 |
| x4 | -0.091726453800 | 11 | 2.55e-12 | 8 |
| x5 | 2.846109372500 | 13 | 1.21e-13 | 10 |
| x6 | -0.00628154090 | 9 | 8.54e-11 | 7 |
| x7 | 1.473820591600 | 13 | 3.10e-13 | 10 |
| x8 | -0.038492716400 | 11 | 1.63e-12 | 8 |
| x9 | 0.915264803700 | 12 | 1.87e-13 | 9 |
| x10 | -2.607139485200 | 13 | 1.40e-13 | 10 |
| x11 | 0.051873629400 | 11 | 1.24e-12 | 8 |
| x12 | -1.084937162500 | 13 | 3.13e-13 | 10 |

Con pivoteo se defienden 9–13 dígitos por incógnita; sin pivoteo, 7–10. Tras el refinamiento iterativo (error 4.9·10⁻¹⁷) las 12 incógnitas son correctas a 10 dígitos o más; con esos valores se arma la tabla resumida de `readme.md`.

## 8. Verificación

- **Residuo en precisión extendida** (80 dígitos decimales, `decimal`): ‖r‖∞ = 2.29·10⁻¹⁰ con pivoteo escalado y 1.77·10⁻⁷ sin pivoteo (vector completo en `resultados/residuo.txt`).
- **Refinamiento iterativo** (1 pasada, reutilizando `PA = LU`): el error relativo baja de 2.70·10⁻¹³ a 4.86·10⁻¹⁷. ‖r‖∞ no baja (2.29·10⁻¹⁰ → 2.60·10⁻¹⁰): tras refinar, x ya está al límite de la doble precisión y el residuo queda acotado por el redondeo de x, |rᵢ| ≲ (ε/2)·Σⱼ|aᵢⱼxⱼ| (≈ 6·10⁻¹⁰ en la fila 6, cuyas entradas son ~10⁶).
- **Control con biblioteca** (`numpy.linalg.solve`, solo como control): error relativo 5.99·10⁻¹³ y ‖r‖∞ = 8.89·10⁻¹¹. La implementación propia con pivoteo escalado tiene menor error relativo (2.70·10⁻¹³).

## 9. Conclusiones

- Se eligió pivoteo parcial escalado con una sola factorización: cuesta 1486 flops (dentro de la cota) más 66 comparaciones de búsqueda, y reduce el error de 2.4·10⁻¹⁰ a 2.7·10⁻¹³. El sobrecosto de estabilidad es O(n²), despreciable frente al O(n³) de la eliminación.
- El refinamiento iterativo, con el residuo en precisión extendida, lleva el error a ≈ 5·10⁻¹⁷ al costo de una pasada O(n²).
- **Con n = 1000:** la eliminación densa costaría ≈ 2n³/3 ≈ 6.7·10⁸ flops, y el almacenamiento sería de 8 MB. Se usaría una LU por bloques de biblioteca (nivel BLAS 3) con pivoteo parcial, equilibrado de filas previo para controlar κ, y refinamiento iterativo con residuo en mayor precisión. Si A fuese rala, se usarían métodos para matrices ralas.

## 10. Reproducibilidad

- **Python:** 3.12.3 (verificado; el código requiere ≥ 3.10) — **NumPy:** 2.5.4 (`requirements.txt`); `decimal` y `fractions` de la biblioteca estándar.
- **Comando:** `pip install -r requirements.txt` y luego `python src/main.py`; regenera los 4 archivos de `resultados/`.
- κ₂ puede variar desde la 6.ª cifra según la versión de LAPACK; el resto es reproducible.
- **Fuentes:** enunciado de la prueba parcial «El pivote» (UCM, noviembre de 2026); documentación de NumPy y de los módulos `fractions` y `decimal` de la biblioteca estándar de Python.
