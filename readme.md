# Gauss 12×12 con pivoteo parcial escalado

Resuelve un sistema lineal denso mal condicionado de 12 incógnitas con eliminación gaussiana propia, con conteo de flops y verificación en precisión extendida.

## Autoría y asignatura

[COMPLETAR nombre(s)] — Computación Numérica (INF-321), Ingeniería Civil Informática, Universidad Católica del Maule. Prueba Parcial «El pivote», noviembre de 2026.

## Requisitos

- Python 3.10 o superior (probado con 3.12.3).
- Instalar dependencias: `pip install -r requirements.txt` (NumPy 2.5.4).

## Ejecución

```bash
python src/main.py
```

Genera en `resultados/`: `solucion.csv`, `traza_eliminacion.txt`, `conteo_flops.md` y `residuo.txt`.

## Estructura de carpetas

```
datos/        A (12×12) y b con todas sus cifras
src/          gauss.py, conteo.py, verificacion.py, main.py
resultados/   salidas generadas por main.py
informe.md    informe técnico y análisis de error
```

## Resultado resumido

| Incógnita | Valor (10 dígitos, tras refinamiento) |
|:--|:--|
| x1 | 0.7384915627 |
| x2 | -1.205937417 |
| x3 | 0.0043871952 |
| x4 | -0.0917264538 |
| x5 | 2.846109372 |
| x6 | -0.0062815409 |
| x7 | 1.473820592 |
| x8 | -0.0384927164 |
| x9 | 0.9152648037 |
| x10 | -2.607139485 |
| x11 | 0.0518736294 |
| x12 | -1.084937162 |

- **Flops (pivoteo escalado):** 1486 (+ 66 comparaciones de búsqueda de pivote)
- **‖r‖∞ (precisión extendida, solución sin refinar):** 2.29·10⁻¹⁰ (tras refinar: 2.60·10⁻¹⁰, ver informe)

## Informe

Análisis completo en [informe.md](informe.md).

