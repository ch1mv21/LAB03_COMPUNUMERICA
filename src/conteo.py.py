from typing import Dict, List


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

    def generar_tabla_markdown(self, n: int = 12) -> str:
        cota_teorica_cubica = (2.0 * (n ** 3)) / 3.0
        cota_referencia = 1486

        lineas: List[str] = [
            "# Reporte de Costo Operacional (FLOPs)",
            "",
            "## 1. Convención de Conteo Declarada",
            "- **1 FLOP** equivale a una operación aritmética elemental en punto flotante (`+`, `-`, `*`, `/`).",
            "- **Reutilización de multiplicadores:** Cada multiplicador $m_{ik} = a_{ik}^{(k)} / a_{kk}^{(k)}$ se calcula **una única vez por fila** y se reutiliza en toda la fila y en $b$.",
            "- **Barrido triangular:** En el paso $k$, las columnas previas ($j < k$) ya anuladas nunca se tocan. El barrido actualiza la fila desde la columna $k$ hasta $n$, más $b_i$.",
            "- **Omisión algebraica del cero subdiagonal:** Si se asigna `A[i, k] = 0.0` sin evaluar `A[i, k] - m * A[k, k]`, se ahorran exactamente 66 multiplicaciones y 66 restas (pasando de **1486 FLOPs** a **1354 FLOPs**).",
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
            f"- **Cota mínima alcanzable de referencia (incluyendo $b$ y sustitución regresiva):** {cota_referencia} flops.",
            f"- **FLOPs contabilizados (con cálculo explícito de subdiagonal):** {self.total_flops} flops (desviación: {abs(self.total_flops - cota_referencia) / cota_referencia * 100:.2f} %).",
            f"- **FLOPs netos sin operar la columna subdiagonal ($1486 - 132$):** {self.total_flops - 132 if self.incluir_cero_subdiagonal else self.total_flops} flops.",
            f"- **Comparaciones iniciales para vector de escalas $s_i$:** {self.comparaciones_escala_inicial} comparaciones.",
            "",
        ])
        return "\n".join(lineas)
