from decimal import Decimal, getcontext
from fractions import Fraction
from typing import Any, Dict, List, Tuple
import math
import numpy as np
from gauss import EliminacionGaussiana

# Configurar precisión extendida de 80 dígitos decimales
getcontext().prec = 80


class VerificadorNumerico:
    """
    Realiza la verificación independiente del sistema lineal en precisión extendida.
    Única clase del módulo verificacion.py.
    """

    def __init__(self, rutas_a_str: List[List[str]], rutas_b_str: List[str]) -> None:
        self.n: int = len(rutas_b_str)
        self.A_dec: List[List[Decimal]] = [
            [Decimal(val.strip()) for val in fila] for fila in rutas_a_str
        ]
        self.b_dec: List[Decimal] = [Decimal(val.strip()) for val in rutas_b_str]
        self.A_np: np.ndarray = np.array(
            [[float(val.strip()) for val in fila] for fila in rutas_a_str],
            dtype=np.float64,
        )
        self.b_np: np.ndarray = np.array(
            [float(val.strip()) for val in rutas_b_str], dtype=np.float64
        )
        self.x_exacto_dec: List[Decimal] = self._resolver_exacto_racional(
            rutas_a_str, rutas_b_str
        )
        self.x_exacto_np: np.ndarray = np.array(
            [float(val) for val in self.x_exacto_dec], dtype=np.float64
        )

    def _resolver_exacto_racional(
        self, a_str: List[List[str]], b_str: List[str]
    ) -> List[Decimal]:
        """Resuelve Ax = b en aritmética racional exacta (fractions.Fraction)."""
        n = len(b_str)
        M = [
            [Fraction(a_str[i][j].strip()) for j in range(n)]
            + [Fraction(b_str[i].strip())]
            for i in range(n)
        ]
        for k in range(n):
            piv = k
            while M[piv][k] == 0:
                piv += 1
            if piv != k:
                M[k], M[piv] = M[piv], M[k]
            piv_val = M[k][k]
            for i in range(k + 1, n):
                factor = M[i][k] / piv_val
                for j in range(k, n + 1):
                    M[i][j] -= factor * M[k][j]
        x_frac = [Fraction(0) for _ in range(n)]
        for i in range(n - 1, -1, -1):
            suma = M[i][n] - sum(M[i][j] * x_frac[j] for j in range(i + 1, n))
            x_frac[i] = suma / M[i][i]
        return [Decimal(f.numerator) / Decimal(f.denominator) for f in x_frac]

    def calcular_residuo_extendido(
        self, x_float: np.ndarray
    ) -> Tuple[List[Decimal], Decimal]:
        """Calcula r = b - A*x en precisión extendida (80 dígitos)."""
        x_dec = [Decimal.from_float(float(val)) for val in x_float]
        residuo: List[Decimal] = []
        for i in range(self.n):
            ax_i = sum(self.A_dec[i][j] * x_dec[j] for j in range(self.n))
            r_i = self.b_dec[i] - ax_i
            residuo.append(r_i)
        norma_inf = max(abs(r) for r in residuo)
        return residuo, norma_inf

    def calcular_error_relativo_norma2(self, x_calc: np.ndarray) -> float:
        """Calcula ||x_calc - x_exacto||_2 / ||x_exacto||_2 en precisión extendida."""
        x_c = [Decimal.from_float(float(v)) for v in x_calc]
        num = sum((x_c[i] - self.x_exacto_dec[i]) ** 2 for i in range(self.n)).sqrt()
        den = sum((self.x_exacto_dec[i]) ** 2 for i in range(self.n)).sqrt()
        return float(num / den)

    def analizar_incognitas(
        self, x_calc: np.ndarray
    ) -> List[Tuple[int, float, Decimal, float, int, float, str]]:
        """Evalúa error relativo, dígitos significativos y bits perdidos por incógnita."""
        eps_maq = np.finfo(np.float64).eps
        resultados = []
        for i in range(self.n):
            val_c = Decimal.from_float(float(x_calc[i]))
            val_e = self.x_exacto_dec[i]
            err_rel = float(abs(val_c - val_e) / abs(val_e))
            if err_rel == 0.0:
                digitos = self.digitos_correctos(0.0)
                bits_perdidos = 0.0
            else:
                digitos = self.digitos_correctos(err_rel)
                bits_perdidos = max(0.0, math.log2(err_rel / eps_maq))
            repr_defendible = f"{float(x_calc[i]):.{max(1, digitos)}g}"
            resultados.append(
                (i + 1, float(x_calc[i]), val_e, err_rel, digitos, bits_perdidos, repr_defendible)
            )
        return resultados

    def ejecutar_refinamiento_iterativo(
        self, solver: EliminacionGaussiana, res_gauss: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Ejecuta una pasada de refinamiento iterativo reutilizando PA = LU."""
        x_ini = res_gauss["x"]
        r_ini, norm_ini = self.calcular_residuo_extendido(x_ini)
        r_float = np.array([float(r) for r in r_ini], dtype=np.float64)

        delta_x = solver.resolver_con_lu(res_gauss, r_float)
        x_ref = x_ini + delta_x
        r_ref, norm_ref = self.calcular_residuo_extendido(x_ref)

        err_ini = self.calcular_error_relativo_norma2(x_ini)
        err_ref = self.calcular_error_relativo_norma2(x_ref)

        return {
            "x_inicial": x_ini.copy(),
            "residuo_inicial": r_ini,
            "norma_inf_inicial": norm_ini,
            "delta_x": delta_x,
            "x_refinado": x_ref,
            "residuo_refinado": r_ref,
            "norma_inf_refinada": norm_ref,
            "error_rel_inicial": err_ini,
            "error_rel_refinado": err_ref,
        }

    @staticmethod
    def digitos_correctos(err_rel: float, maximo: int = 15) -> int:
        """Dígitos significativos correctos, con definición conservadora.

        Se toma la parte entera de -log10(error relativo), sin redondear hacia
        arriba: un error de 1.6e-12 respalda 11 dígitos, no 12.
        """
        if err_rel <= 0.0:
            return maximo
        return max(0, min(maximo, int(math.floor(-math.log10(err_rel)))))

    def ejecutar_refinamiento_multiple(
        self,
        solver: EliminacionGaussiana,
        res_gauss: Dict[str, Any],
        max_iter: int = 3,
    ) -> Dict[str, Any]:
        """Refinamiento iterativo de varias pasadas reutilizando PA = LU.

        En cada pasada: r = b - A x en precisión extendida (80 dígitos), se
        resuelve A d = r con L y U ya calculadas (sin refactorizar) y x <- x + d.
        Devuelve una fila por pasada (la 0 es la solución sin refinar).
        """
        x = res_gauss["x"].copy()
        filas: List[Dict[str, Any]] = []
        for k in range(max_iter + 1):
            r_dec, norma_r = self.calcular_residuo_extendido(x)
            err = self.calcular_error_relativo_norma2(x)
            fila: Dict[str, Any] = {
                "iteracion": k,
                "norma_inf_residuo": norma_r,
                "error_relativo": err,
                "norma_inf_correccion": None,
                "x": x.copy(),
            }
            if k < max_iter:
                r_float = np.array([float(r) for r in r_dec], dtype=np.float64)
                delta = solver.resolver_con_lu(res_gauss, r_float)
                fila["norma_inf_correccion"] = float(np.max(np.abs(delta)))
                x = x + delta
            filas.append(fila)
        return {"iteraciones": filas, "x_final": filas[-1]["x"]}

    def analizar_condicionamiento_y_dependencia(self) -> Dict[str, Any]:
        """Calcula condicionamiento, determinante y la casi dependencia de la fila 7."""
        cond_2 = float(np.linalg.cond(self.A_np, 2))
        cond_inf = float(np.linalg.cond(self.A_np, np.inf))
        det_a = float(np.linalg.det(self.A_np))
        eps = float(np.finfo(np.float64).eps)
        cota_error = cond_2 * eps

        f3 = self.A_np[2, :]
        f5 = self.A_np[4, :]
        f7 = self.A_np[6, :]
        coef, _, _, _ = np.linalg.lstsq(np.column_stack([f3, f5]), f7, rcond=None)
        pert_exacta = f7 - (4.173 * f3 + 582.7 * f5)

        return {
            "det_A": det_a,
            "cond_2": cond_2,
            "cond_inf": cond_inf,
            "eps_maq": eps,
            "cota_relativa": cota_error,
            "alpha_lstsq": float(coef[0]),
            "beta_lstsq": float(coef[1]),
            "alpha_exacto": 4.173,
            "beta_exacto": 582.7,
            "norma_inf_perturbacion": float(np.max(np.abs(pert_exacta))),
            "perturbacion_vector": pert_exacta,
        }