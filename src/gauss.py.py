from typing import Any, Dict, List, Optional, Tuple
import numpy as np
from conteo import ContadorOperaciones


class EliminacionGaussiana:
    """
    Solver directo para sistemas lineales densos Ax = b mediante eliminación gaussiana.
    Única clase del módulo gauss.py.
    """

    def __init__(self, A: np.ndarray, b: np.ndarray) -> None:
        if A.shape[0] != A.shape[1]:
            raise ValueError("La matriz A debe ser cuadrada.")
        if A.shape[0] != b.shape[0]:
            raise ValueError("Las dimensiones de A y b son incompatibles.")
        self.n: int = A.shape[0]
        self.A_original: np.ndarray = np.array(A, dtype=np.float64, copy=True)
        self.b_original: np.ndarray = np.array(b, dtype=np.float64, copy=True)

    def _calcular_escalas_fila(
        self, A: np.ndarray, contador: ContadorOperaciones
    ) -> np.ndarray:
        """Calcula s_i = max_{1<=j<=n} |a_{ij}| una sola vez antes de eliminar."""
        escalas = np.zeros(self.n, dtype=np.float64)
        comparaciones = 0
        for i in range(self.n):
            max_val = abs(A[i, 0])
            for j in range(1, self.n):
                comparaciones += 1
                val = abs(A[i, j])
                if val > max_val:
                    max_val = val
            if max_val == 0.0:
                raise ZeroDivisionError(f"La fila {i + 1} es completamente nula.")
            escalas[i] = max_val
        contador.registrar_escala_inicial(comparaciones)
        return escalas

    def resolver(
        self,
        estrategia: str = "escalado",
        incluir_cero_subdiagonal: bool = True,
    ) -> Dict[str, Any]:
        """
        Ejecuta eliminación gaussiana hacia adelante y sustitución regresiva.
        Retorna un diccionario estructurado con la solución, factores LU, traza y contador.
        """
        if estrategia not in ("sin_pivoteo", "parcial", "escalado"):
            raise ValueError(f"Estrategia no reconocida: {estrategia}")

        U = self.A_original.copy()
        b_mod = self.b_original.copy()
        L = np.eye(self.n, dtype=np.float64)
        perm = np.arange(self.n, dtype=int)

        contador = ContadorOperaciones(
            estrategia=estrategia,
            incluir_cero_subdiagonal=incluir_cero_subdiagonal,
        )

        if estrategia == "escalado":
            escalas = self._calcular_escalas_fila(U, contador)
        else:
            escalas = np.ones(self.n, dtype=np.float64)

        max_inicial_A = float(np.max(np.abs(U)))
        max_historico_U = max_inicial_A
        multiplicador_max = 0.0
        pivote_minimo = float("inf")
        traza_pasos: List[Tuple[int, np.ndarray, np.ndarray, str]] = []

        # Fase 1: Eliminación hacia adelante (11 pasos para n = 12)
        for k in range(self.n - 1):
            nombre_paso = f"Paso {k + 1:02d} (col {k + 1})"
            fila_pivote = k
            hubo_intercambio = False

            if estrategia == "parcial":
                mejor_val = abs(U[k, k])
                comps = 0
                for i in range(k + 1, self.n):
                    comps += 1
                    candidato = abs(U[i, k])
                    if candidato > mejor_val:
                        mejor_val = candidato
                        fila_pivote = i
                hubo_intercambio = fila_pivote != k
                contador.registrar_busqueda_pivote(
                    nombre_paso, comparaciones=comps, divisiones_escala=0, intercambio=hubo_intercambio
                )

            elif estrategia == "escalado":
                mejor_cociente = abs(U[k, k]) / escalas[k]
                comps = 0
                divs_escala = 1
                for i in range(k + 1, self.n):
                    comps += 1
                    divs_escala += 1
                    cociente = abs(U[i, k]) / escalas[i]
                    if cociente > mejor_cociente:
                        mejor_cociente = cociente
                        fila_pivote = i
                hubo_intercambio = fila_pivote != k
                contador.registrar_busqueda_pivote(
                    nombre_paso,
                    comparaciones=comps,
                    divisiones_escala=divs_escala,
                    intercambio=hubo_intercambio,
                )
            else:
                contador.registrar_busqueda_pivote(
                    nombre_paso, comparaciones=0, divisiones_escala=0, intercambio=False
                )

            detalle_pivote = (
                f"Pivote en F{k + 1} (original F{perm[fila_pivote] + 1}), "
                f"valor = {U[fila_pivote, k]:.9e}"
            )
            if hubo_intercambio:
                detalle_pivote += f" [Intercambio F{k + 1} <-> F{fila_pivote + 1}]"
                U[[k, fila_pivote], :] = U[[fila_pivote, k], :]
                b_mod[[k, fila_pivote]] = b_mod[[fila_pivote, k]]
                escalas[[k, fila_pivote]] = escalas[[fila_pivote, k]]
                perm[[k, fila_pivote]] = perm[[fila_pivote, k]]
                if k > 0:
                    L[[k, fila_pivote], :k] = L[[fila_pivote, k], :k]

            pivote_actual = U[k, k]
            if abs(pivote_actual) == 0.0:
                raise ZeroDivisionError(f"Pivote nulo en el paso {k + 1}.")
            pivote_minimo = min(pivote_minimo, abs(pivote_actual))

            # Eliminación bajo la diagonal en la columna k
            for i in range(k + 1, self.n):
                m_ik = U[i, k] / pivote_actual
                contador.registrar_mul_div(nombre_paso, 1)
                L[i, k] = m_ik

                if abs(m_ik) > multiplicador_max:
                    multiplicador_max = abs(m_ik)

                if incluir_cero_subdiagonal:
                    U[i, k] = U[i, k] - m_ik * U[k, k]
                    contador.registrar_mul_div(nombre_paso, 1)
                    contador.registrar_sum_res(nombre_paso, 1)
                U[i, k] = 0.0

                for j in range(k + 1, self.n):
                    U[i, j] = U[i, j] - m_ik * U[k, j]
                    contador.registrar_mul_div(nombre_paso, 1)
                    contador.registrar_sum_res(nombre_paso, 1)

                b_mod[i] = b_mod[i] - m_ik * b_mod[k]
                contador.registrar_mul_div(nombre_paso, 1)
                contador.registrar_sum_res(nombre_paso, 1)

            max_paso = float(np.max(np.abs(U)))
            if max_paso > max_historico_U:
                max_historico_U = max_paso

            traza_pasos.append((k + 1, U.copy(), b_mod.copy(), detalle_pivote))

        pivote_minimo = min(pivote_minimo, abs(U[self.n - 1, self.n - 1]))
        factor_crecimiento = max_historico_U / max_inicial_A

        # Fase 2: Sustitución regresiva
        x = self._sustitucion_regresiva(U, b_mod, contador)

        return {
            "estrategia": estrategia,
            "x": x,
            "L": L,
            "U": U,
            "perm": perm,
            "factores_escala": escalas,
            "multiplicador_max": multiplicador_max,
            "factor_crecimiento": factor_crecimiento,
            "max_inicial_A": max_inicial_A,
            "max_historico_U": max_historico_U,
            "pivote_minimo": pivote_minimo,
            "traza_pasos": traza_pasos,
            "contador": contador,
        }

    def _sustitucion_regresiva(
        self,
        U: np.ndarray,
        b_mod: np.ndarray,
        contador: Optional[ContadorOperaciones] = None,
    ) -> np.ndarray:
        """Resuelve el sistema triangular superior Ux = b_mod en O(n^2) flops."""
        nombre_etapa = "Sust. Regresiva"
        x = np.zeros(self.n, dtype=np.float64)

        for i in range(self.n - 1, -1, -1):
            acumulador = b_mod[i]
            for j in range(i + 1, self.n):
                acumulador = acumulador - U[i, j] * x[j]
                if contador is not None:
                    contador.registrar_mul_div(nombre_etapa, 1)
                    contador.registrar_sum_res(nombre_etapa, 1)
            x[i] = acumulador / U[i, i]
            if contador is not None:
                contador.registrar_mul_div(nombre_etapa, 1)

        return x

    def resolver_con_lu(
        self, resultado_previo: Dict[str, Any], nuevo_b: np.ndarray
    ) -> np.ndarray:
        """
        Reutiliza la factorización PA = LU ya calculada para resolver A * dx = nuevo_b
        en O(n^2) operaciones (usado en el refinamiento iterativo sin refactorizar A).
        """
        b_permutado = nuevo_b[resultado_previo["perm"]]
        L = resultado_previo["L"]
        U = resultado_previo["U"]
        d = np.zeros(self.n, dtype=np.float64)
        for i in range(self.n):
            suma = b_permutado[i]
            for j in range(i):
                suma -= L[i, j] * d[j]
            d[i] = suma
        return self._sustitucion_regresiva(U, d, contador=None)
