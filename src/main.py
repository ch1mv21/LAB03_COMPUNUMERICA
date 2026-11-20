#Integrantes: Jose Herrera Manuel Soto
#Curso: Computacion Numerica
#Fecha: 21/11/2026

from pathlib import Path
from typing import Any, Dict, List, Tuple
import csv
import math
import numpy as np

from gauss import EliminacionGaussiana
from verificacion import VerificadorNumerico


class EjecutorExperimento:
    """
    Orquestador principal del proyecto. Carga los datos sin truncar, ejecuta las
    corridas del solver, coordina la verificación y exporta resultados/.
    Única clase del módulo main.py.
    """

    def __init__(self, directorio_raiz: Path) -> None:
        self.raiz: Path = directorio_raiz
        self.dir_datos: Path = self.raiz / "datos"
        self.dir_resultados: Path = self.raiz / "resultados"
        self.dir_resultados.mkdir(parents=True, exist_ok=True)

    def cargar_datos_crudos(
        self,
    ) -> Tuple[List[List[str]], List[str], np.ndarray, np.ndarray]:
        """Lee los CSV conservando todas las cifras como cadenas y como float64."""
        ruta_matriz = self.dir_datos / "matriz_a.csv"
        ruta_vector = self.dir_datos / "vector_b.csv"

        filas_a_str: List[List[str]] = []
        with ruta_matriz.open("r", encoding="utf-8") as f_a:
            lector = csv.reader(f_a)
            for fila in lector:
                fila_limpia = [celda.strip() for celda in fila if celda.strip()]
                if fila_limpia:
                    filas_a_str.append(fila_limpia)

        filas_b_str: List[str] = []
        with ruta_vector.open("r", encoding="utf-8") as f_b:
            for linea in f_b:
                texto = linea.strip()
                if texto:
                    filas_b_str.append(texto.split(",")[0].strip())

        if len(filas_a_str) != 12 or any(len(f) != 12 for f in filas_a_str):
            raise ValueError("matriz_a.csv debe tener exactamente 12x12 elementos.")
        if len(filas_b_str) != 12:
            raise ValueError("vector_b.csv debe tener exactamente 12 elementos.")

        A_np = np.array(
            [[float(x) for x in fila] for fila in filas_a_str], dtype=np.float64
        )
        b_np = np.array([float(x) for x in filas_b_str], dtype=np.float64)
        return filas_a_str, filas_b_str, A_np, b_np

    def exportar_solucion_csv(
        self,
        verificador: VerificadorNumerico,
        res_escalado: Dict[str, Any],
        res_sin_piv: Dict[str, Any],
        refinamiento: Dict[str, Any],
    ) -> None:
        """Genera resultados/solucion.csv con la solución refinada y sus dígitos defendibles."""
        ruta_salida = self.dir_resultados / "solucion.csv"
        analisis_ref = verificador.analizar_incognitas(refinamiento["x_refinado"])
        analisis_esc = verificador.analizar_incognitas(res_escalado["x"])
        analisis_sin = verificador.analizar_incognitas(res_sin_piv["x"])

        with ruta_salida.open("w", newline="", encoding="utf-8") as f_out:
            escritor = csv.writer(f_out)
            escritor.writerow([
                "incognita",
                "valor_reportado_refinado",
                "digitos_significativos_afirmados",
                "error_relativo_refinado",
                "digitos_pivoteo_escalado_sin_refinar",
                "digitos_sin_pivoteo",
                "solucion_exacta_referencia",
            ])
            for (idx, _, exacto, err_r, dig_r, _, repr_r), (
                _, _, _, _, dig_e, _, _
            ), (_, _, _, _, dig_s, _, _) in zip(
                analisis_ref, analisis_esc, analisis_sin
            ):
                escritor.writerow([
                    f"x_{idx}",
                    repr_r,
                    dig_r,
                    f"{err_r:.4e}",
                    dig_e,
                    dig_s,
                    f"{exacto:.10f}",
                ])
    def exportar_traza_eliminacion(
        self,
        res_escalado: Dict[str, Any],
        res_sin_piv: Dict[str, Any],
    ) -> None:
        """Genera resultados/traza_eliminacion.txt con la matriz tras cada paso."""
        ruta_salida = self.dir_resultados / "traza_eliminacion.txt"
        lineas: List[str] = [
            "================================================================================",
            "TRAZA COMPLETA DE ELIMINACIÓN GAUSSIANA (11 PASOS)",
            "================================================================================",
            "",
        ]

        for titulo, resultado in [
            ("CORRIDA 1: CON PIVOTEO PARCIAL ESCALADO (ÓPTIMA EN ESTABILIDAD)", res_escalado),
            ("CORRIDA 2: SIN PIVOTEO (DEGRADACIÓN DE MANTISA POR PIVOTE a(1,1))", res_sin_piv),
        ]:
            lineas.append(f"### {titulo} ###")
            lineas.append(f"Multiplicador máximo: {resultado['multiplicador_max']:.6e}")
            lineas.append(f"Factor de crecimiento: {resultado['factor_crecimiento']:.6e}")
            lineas.append(
                f"Permutación final de filas (1-indexado): {[int(p + 1) for p in resultado['perm']]}"
            )
            lineas.append("")

            for idx_paso, (paso_num, mat_u, vec_b, info_piv) in enumerate(
                resultado["traza_pasos"]
            ):
                ceros_bajo_diag = sum(
                 1
                 for col in range(paso_num)
                 for fila in range(col + 1, mat_u.shape[0])
                 if mat_u[fila, col] == 0.0
                )
                lineas.append(
                    f"--- Paso {paso_num:02d}/11 | {info_piv} | Ceros acumulados bajo diag: {ceros_bajo_diag} ---"
                )
                lineas.append(
                    "Multiplicadores m(i,k) de este paso (filas ya intercambiadas, 17 cifras):"
                )
                for desplaz, m_ik in enumerate(resultado["multiplicadores_pasos"][idx_paso]):
                    lineas.append(f"   m({paso_num + 1 + desplaz:02d},{paso_num:02d}) = {m_ik:+.16e}")
                lineas.append("Matriz [A | b] tras el paso (17 cifras; los elementos anulados se fijan en 0):")
                for i in range(mat_u.shape[0]):
                    fila_str = " ".join(
                        f"{mat_u[i, j]:+.16e}" for j in range(mat_u.shape[1])
                    )
                    lineas.append(f"F{i + 1:02d}: [{fila_str} | {vec_b[i]:+.16e}]")
                lineas.append("")

        ruta_salida.write_text("\n".join(lineas), encoding="utf-8")

    def exportar_residuo_txt(
        self,
        verificador: VerificadorNumerico,
        res_escalado: Dict[str, Any],
        res_sin_piv: Dict[str, Any],
        refinamiento: Dict[str, Any],
        info_cond: Dict[str, Any],
        x_numpy: np.ndarray,
    ) -> None:
        """Genera resultados/residuo.txt con el residuo en precisión extendida."""
        ruta_salida = self.dir_resultados / "residuo.txt"
        r_sin, norm_sin = verificador.calcular_residuo_extendido(res_sin_piv["x"])
        r_np, norm_np = verificador.calcular_residuo_extendido(x_numpy)
        err_np = verificador.calcular_error_relativo_norma2(x_numpy)
        err_sin = verificador.calcular_error_relativo_norma2(res_sin_piv["x"])

        eps_maq = info_cond["eps_maq"]
        bits_perdidos_esc = max(
            0.0, math.log2(refinamiento["error_rel_inicial"] / eps_maq)
        )
        bits_perdidos_sin = max(0.0, math.log2(err_sin / eps_maq))

        lineas: List[str] = [
            "================================================================================",
            "VERIFICACIÓN INDEPENDIENTE EN PRECISIÓN EXTENDIDA (80 DÍGITOS DECIMALES)",
            "================================================================================",
            "",
            "1. RESUMEN DE ESTABILIDAD, ERROR RELATIVO Y MANTISA",
            f"   - Condicionamiento L2 kappa_2(A)         : {info_cond['cond_2']:.6e}",
            f"   - Condicionamiento Linf kappa_inf(A)     : {info_cond['cond_inf']:.6e}",
            f"   - Determinante det(A)                    : {info_cond['det_A']:.6e}",
            f"   - Cota clásica pesimista kappa_2(A)*eps  : {info_cond['cota_relativa']:.6e}",
            f"   - Casi dependencia lineal detectada      : F7 ≈ {info_cond['alpha_exacto']}*F3 + {info_cond['beta_exacto']}*F5 (||delta||_inf = {info_cond['norma_inf_perturbacion']:.6e}; ajuste MCO: alpha={info_cond['alpha_lstsq']:.6f}, beta={info_cond['beta_lstsq']:.6f})",
            "",
            "2. COMPARACIÓN DE CORRIDAS",
            "   a) Sin Pivoteo:",
            f"      * Multiplicador máximo     : {res_sin_piv['multiplicador_max']:.6e}",
            f"      * Factor de crecimiento    : {res_sin_piv['factor_crecimiento']:.6e}",
            f"      * Error relativo ||dx||/||x||: {err_sin:.6e}",
            f"      * Bits de mantisa perdidos : {bits_perdidos_sin:.2f} de 53 bits",
            f"      * Norma infinito ||r||_inf : {norm_sin:.6e}",
            "",
            "   b) Con Pivoteo Parcial Escalado:",
            f"      * Multiplicador máximo     : {res_escalado['multiplicador_max']:.6e}",
            f"      * Factor de crecimiento    : {res_escalado['factor_crecimiento']:.6e}",
            f"      * Error relativo ||dx||/||x||: {refinamiento['error_rel_inicial']:.6e}",
            f"      * Bits de mantisa perdidos : {bits_perdidos_esc:.2f} de 53 bits",
            f"      * Norma infinito ||r||_inf : {refinamiento['norma_inf_inicial']:.6e}",
            "",
            "   c) Control con Biblioteca (numpy.linalg.solve):",
            f"      * Error relativo ||dx||/||x||: {err_np:.6e}",
            f"      * Norma infinito ||r||_inf : {norm_np:.6e}",
            "",
            "3. REFINAMIENTO ITERATIVO (1 PASADA REUTILIZANDO PA = LU)",
            f"   - Norma infinito del residuo antes   ||r^(0)||_inf : {refinamiento['norma_inf_inicial']:.6e}",
            f"   - Norma infinito del residuo después ||r^(1)||_inf : {refinamiento['norma_inf_refinada']:.6e}",
            f"   - Error relativo antes               ||e^(0)||_2   : {refinamiento['error_rel_inicial']:.6e}",
            f"   - Error relativo después             ||e^(1)||_2   : {refinamiento['error_rel_refinado']:.6e}",
            "",
            "4. VECTOR RESIDUO DETALLADO r = b - Ax (PRECISIÓN EXTENDIDA)",
            "Fila | r_i (Pivoteo Escalado)               | r_i (Refinado 1 pasada)              | r_i (Sin Pivoteo)",
            "-----+--------------------------------------+--------------------------------------+--------------------------------------",
        ]

        for i in range(12):
            lineas.append(
                f"F{i + 1:02d}  | {refinamiento['residuo_inicial'][i]:+36.28e} | "
                f"{refinamiento['residuo_refinado'][i]:+36.28e} | {r_sin[i]:+36.28e}"
            )

        ruta_salida.write_text("\n".join(lineas) + "\n", encoding="utf-8")

    def ejecutar(self) -> None:
        """Ejecuta el flujo completo y genera los 4 archivos en resultados/."""
        a_str, b_str, A_np, b_np = self.cargar_datos_crudos()

        solver = EliminacionGaussiana(A_np, b_np)
        res_sin_piv = solver.resolver(
            estrategia="sin_pivoteo", incluir_cero_subdiagonal=True
        )
        res_escalado = solver.resolver(
            estrategia="escalado", incluir_cero_subdiagonal=True
        )

        # Segunda corrida instrumentada: subdiagonal asignada a 0 sin operar (mínimo medido)
        res_minimo = solver.resolver(
            estrategia="escalado", incluir_cero_subdiagonal=False
        )

        verificador = VerificadorNumerico(a_str, b_str)
        refinamiento = verificador.ejecutar_refinamiento_iterativo(
            solver, res_escalado
        )
        info_cond = verificador.analizar_condicionamiento_y_dependencia()
        x_numpy = np.linalg.solve(A_np, b_np)

        self.exportar_solucion_csv(verificador, res_escalado, res_sin_piv, refinamiento)
        self.exportar_traza_eliminacion(res_escalado, res_sin_piv)
        (self.dir_resultados / "conteo_flops.md").write_text(
            res_escalado["contador"].generar_tabla_markdown(
                n=12, alternativo=res_minimo["contador"]
            ), encoding="utf-8"
        )
        self.exportar_residuo_txt(
            verificador,
            res_escalado,
            res_sin_piv,
            refinamiento,
            info_cond,
            x_numpy,
        )

        print("Ejecución completada exitosamente. Archivos generados en 'resultados/':")
        print(f" - FLOPs contabilizados: {res_escalado['contador'].total_flops}")
        print(
            f" - Error relativo (Sin pivoteo)     : {verificador.calcular_error_relativo_norma2(res_sin_piv['x']):.4e}"
        )
        print(
            f" - Error relativo (Pivoteo escalado): {refinamiento['error_rel_inicial']:.6e}"
        )
        print(
            f" - Error relativo (Tras 1 ref. iter): {refinamiento['error_rel_refinado']:.6e}"
        )


if __name__ == "__main__":
    raiz_proyecto = Path(__file__).resolve().parent.parent
    app = EjecutorExperimento(raiz_proyecto)
    app.ejecutar()  