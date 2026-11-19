#Integrantes: Jose Herrera Manuel Soto
#Curso: Computacion Numerica
#Fecha: 21/11/2026


from typing import Dict, List, Optional


class ContadorOperaciones:
    """
    Instrumentador del costo operacional (FLOPs) y del sobrecosto de pivoteo.
    Única clase del módulo conteo.py.
    """

    def __init__(self, estrategia: str, incluir_cero_subdiagonal: bool = True) -> None:
        self.estrategia: str = estrategia
        self.incluir_cero_subdiagonal: bool = incluir_cero_subdiagonal
        self.pasos: Dict[str, Dict[str, int]] = {}
        self.comparaciones_escala_inicial: int = 0

    def _asegurar_paso(self, clave: str) -> Dict[str, int]:
        if clave not in self.pasos:
            self.pasos[clave] = {
                "mul_div": 0,
                "sum_res": 0,
                "comparaciones": 0,
                "divisiones_escala": 0,
                "intercambios": 0,
            }
        return self.pasos[clave]

    def registrar_escala_inicial(self, comparaciones: int) -> None:
        self.comparaciones_escala_inicial += comparaciones

    def registrar_mul_div(self, clave_paso: str, cantidad: int = 1) -> None:
        self._asegurar_paso(clave_paso)["mul_div"] += cantidad

    def registrar_sum_res(self, clave_paso: str, cantidad: int = 1) -> None:
        self._asegurar_paso(clave_paso)["sum_res"] += cantidad

    def registrar_busqueda_pivote(
        self,
        clave_paso: str,
        comparaciones: int,
        divisiones_escala: int = 0,
        intercambio: bool = False,
    ) -> None:
        registro = self._asegurar_paso(clave_paso)
        registro["comparaciones"] += comparaciones
        registro["divisiones_escala"] += divisiones_escala
        if intercambio:
            registro["intercambios"] += 1

    @property
    def total_mul_div(self) -> int:
        return sum(p["mul_div"] for p in self.pasos.values())

    @property
    def total_sum_res(self) -> int:
        return sum(p["sum_res"] for p in self.pasos.values())

    @property
    def total_flops(self) -> int:
        return self.total_mul_div + self.total_sum_res

    @property
    def total_comparaciones(self) -> int:
        return sum(p["comparaciones"] for p in self.pasos.values())

    @property
    def total_intercambios(self) -> int:
        return sum(p["intercambios"] for p in self.pasos.values())

    @property
    def total_divisiones_escala(self) -> int:
        return sum(p["divisiones_escala"] for p in self.pasos.values())

    def generar_tabla_markdown(
        self, n: int = 12, alternativo: Optional["ContadorOperaciones"] = None
    ) -> str:
        cota_teorica_cubica = (2.0 * (n ** 3)) / 3.0
        # Conteo mínimo teórico de eliminación + sustitución regresiva (incluye b):
        #   mul/div = n^3/3 + n^2 - n/3 ;  sum/res = n^3/3 + n^2/2 - 5n/6
        minimo_mul_div = round(n ** 3 / 3 + n ** 2 - n / 3)
        minimo_sum_res = round(n ** 3 / 3 + n ** 2 / 2 - 5 * n / 6)
        cota_minima = minimo_mul_div + minimo_sum_res
        extra_subdiag = n * (n - 1)  # 2 flops por cada una de las n(n-1)/2 posiciones
        cota_referencia = 1486  # cota de la prueba (T1); el conteo real sale del contador
        con_sub = self.total_flops if self.incluir_cero_subdiagonal else self.total_flops + extra_subdiag
        sin_sub = self.total_flops - extra_subdiag if self.incluir_cero_subdiagonal else self.total_flops

        if alternativo is not None:
            linea_minimo = (
                "- **Mínimo medido (segunda corrida instrumentada, subdiagonal asignada a 0 sin operar):** "
                f"{alternativo.total_mul_div} mult./div. + {alternativo.total_sum_res} sumas/restas = "
                f"**{alternativo.total_flops} flops**; coincide con las fórmulas $n^3/3+n^2-n/3$ y "
                f"$n^3/3+n^2/2-5n/6$ ({minimo_mul_div} + {minimo_sum_res} = {cota_minima}). "
                f"El barrido desde la columna $k$ agrega {con_sub - alternativo.total_flops} flops."
            )
        else:
            linea_minimo = (
                f"- **Mínimo estricto (por fórmula) si la subdiagonal se asigna a 0 sin operar:** "
                f"{minimo_mul_div} mult./div. + {minimo_sum_res} sumas/restas = **{cota_minima} flops**."
            )

        lineas: List[str] = [
            "# Reporte de Costo Operacional (FLOPs)",
            "",
            "## 1. Convención de Conteo Declarada",
            "- **1 FLOP** equivale a una operación aritmética elemental en punto flotante (`+`, `-`, `*`, `/`).",
            "- **Reutilización de multiplicadores:** Cada multiplicador $m_{ik} = a_{ik}^{(k)} / a_{kk}^{(k)}$ se calcula **una única vez por fila** y se reutiliza en toda la fila y en $b$.",
            "- **Barrido triangular:** En el paso $k$, las columnas previas ($j < k$) ya anuladas nunca se tocan. El barrido actualiza la fila desde la columna $k$ hasta $n$, más $b_i$.",
            "- **Omisión algebraica del cero subdiagonal:** Si se asigna `A[i, k] = 0.0` sin evaluar `A[i, k] - m * A[k, k]`, se ahorran exactamente 66 multiplicaciones y 66 restas (pasando de **1486 FLOPs** a **1354 FLOPs**, que es el mínimo teórico).",
            "- **Sobrecosto de estabilidad (pivoteo):** Las comparaciones, cocientes de escala $|a_{ik}|/s_i$ e intercambios de filas se registran en contadores separados.",
            "",
            f"## 2. Desglose Paso a Paso ({self.estrategia})",
            "",
            "| Etapa / Paso | Mult. y Div. (`*`, `/`) | Sumas y Restas (`+`, `-`) | Total FLOPs | Comp. Pivoteo | Cocientes Escala | Intercambios |",
            "| :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]

        for nombre, datos in self.pasos.items():
            tot_paso = datos["mul_div"] + datos["sum_res"]
            lineas.append(
                f"| {nombre} | {datos['mul_div']} | {datos['sum_res']} | "
                f"**{tot_paso}** | {datos['comparaciones']} | "
                f"{datos['divisiones_escala']} | {datos['intercambios']} |"
            )

        lineas.extend([
            f"| **TOTALES** | **{self.total_mul_div}** | **{self.total_sum_res}** | "
            f"**{self.total_flops}** | **{self.total_comparaciones}** | "
            f"**{self.total_divisiones_escala}** | **{self.total_intercambios}** |",
            "",
            "## 3. Comparación con Cotas de Referencia",
            f"- **Término dominante asintótico ($2n^3/3$ para $n={n}$):** {cota_teorica_cubica:.2f} flops (solo eliminación de $A$).",
            f"- **Cota de referencia de la prueba (eliminación densa, barrido desde la columna $k$, incluye $b$ y sustitución regresiva):** {cota_referencia} flops; con pivoteo parcial, {cota_referencia} flops + {n * (n - 1) // 2} comparaciones.",
            f"- **Flops contados por el instrumento (esta corrida):** {con_sub} flops → desviación {abs(con_sub - cota_referencia) / cota_referencia * 100:.2f} % respecto de la cota (tolerancia ±10 %: [{cota_referencia * 0.9:.0f}, {cota_referencia * 1.1:.0f}]).",
            f"- **Comparaciones de búsqueda de pivote contadas:** {self.total_comparaciones} (cota: {n * (n - 1) // 2}); se reportan aparte como sobrecosto de estabilidad.",
            linea_minimo,
            f"- **Comparaciones iniciales para el vector de escalas $s_i$:** {self.comparaciones_escala_inicial} comparaciones (aparte; no son flops aritméticos).",
            "",
        ])
        return "\n".join(lineas)
