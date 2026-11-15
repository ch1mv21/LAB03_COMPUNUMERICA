"""
Verificación independiente del sistema Ax = b.

Parte 1: carga exacta de los datos, residuo en precisión extendida
(fractions) y control con numpy.linalg.solve.
"""
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import List, Tuple

import numpy as np

Vector = List[Fraction]
Matriz = List[List[Fraction]]


def leer_csv_exacto(ruta: Path) -> Matriz:
    """Lee un CSV como texto y convierte cada celda a Fraccion sin pasar por float.

    Fraccion("1.3407e-08") interpreta el texto tal cual, por lo que no se pierde
    ninguna cifra (los 20 dígitos de b quedan completos).
    """
    filas: Matriz = []
    with open(ruta, encoding="utf-8") as archivo:
        for linea in archivo:
            linea = linea.strip()
            if linea:
                filas.append([Fraction(celda.strip()) for celda in linea.split(",")])
    return filas


def cargar_sistema_exacto(ruta_a: Path, ruta_b: Path) -> Tuple[Matriz, Vector]:
    """Devuelve A y b como racionales exactos y valida sus dimensiones."""
    A = leer_csv_exacto(ruta_a)
    b = [fila[0] for fila in leer_csv_exacto(ruta_b)]
    n = len(A)
    if any(len(fila) != n for fila in A):
        raise ValueError("La matriz A debe ser cuadrada.")
    if len(b) != n:
        raise ValueError("b debe tener tantas filas como A.")
    return A, b


def a_flotante(A: Matriz, b: Vector) -> Tuple[np.ndarray, np.ndarray]:
    """Redondea A y b a doble precisión (único punto donde se pierde información)."""
    return (
        np.array([[float(v) for v in fila] for fila in A], dtype=np.float64),
        np.array([float(v) for v in b], dtype=np.float64),
    )


def residuo_exacto(A: Matriz, b: Vector, x: np.ndarray) -> Vector:
    """Calcula r = b - A x en aritmética racional exacta.

    Cada componente de x (un double) se convierte a Fraccion de forma exacta,
    por lo que el residuo no tiene ningún error de redondeo.
    """
    x_exacto = [Fraction(float(v)) for v in x]
    return [
        b[i] - sum(A[i][j] * x_exacto[j] for j in range(len(x_exacto)))
        for i in range(len(b))
    ]


def norma_inf(v: Vector) -> Fraction:
    """Norma infinito: máximo valor absoluto de las componentes."""
    return max(abs(c) for c in v)


def formatear(q: Fraction, cifras: int = 6) -> str:
    """Escribe un Fraccion en notación científica con 'cifras' decimales."""
    return format(Decimal(q.numerator) / Decimal(q.denominator), f".{cifras}e")


def control_numpy(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Solución de control con numpy.linalg.solve (NO sustituye la implementación propia)."""
    return np.linalg.solve(A, b)


if __name__ == "__main__":
    raiz = Path(__file__).resolve().parent.parent
    A_ex, b_ex = cargar_sistema_exacto(raiz / "datos" / "matriz_a.csv",
                                       raiz / "datos" / "vector_b.csv")
    A_fl, b_fl = a_flotante(A_ex, b_ex)
    x_control = control_numpy(A_fl, b_fl)
    r = residuo_exacto(A_ex, b_ex, x_control)
    print(f"n = {len(b_ex)}")
    print(f"Residuo del control numpy, norma infinito (exacta): {formatear(norma_inf(r))}")
