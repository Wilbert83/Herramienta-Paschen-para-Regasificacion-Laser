# 🔬 Laser Re-Gassing & Paschen Law Tool

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![Flet](https://img.shields.io/badge/Flet-GUI-00BCD4?style=flat-square)](https://flet.dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](LICENSE)
[![Status](https://img.shields.io/badge/status-active-brightgreen?style=flat-square)]()
[![Physics](https://img.shields.io/badge/domain-plasma%20physics-blueviolet?style=flat-square)]()

Herramienta de ingeniería para calcular el **voltaje de ruptura dieléctrica** de mezclas de gases mediante la Ley de Paschen, y simular el proceso completo de **re-gassing** de tubos láser CO₂ tipo RECI. Incluye una versión de terminal (CLI) y una aplicación de escritorio con panel SCADA en tiempo real.

---

## 🧠 Contexto Técnico

### Ley de Paschen

El voltaje de ruptura eléctrica en un gas depende del producto presión × distancia entre electrodos ($p \cdot d$). La fórmula de Paschen modela la ionización por colisiones en avalancha de Townsend:

$$V_b = \frac{B \cdot p \cdot d}{\ln(A \cdot p \cdot d) - \ln\!\left(\ln\!\left(1 + \dfrac{1}{\gamma}\right)\right)}$$

Donde:
| Símbolo | Descripción | Unidades |
|---------|-------------|----------|
| $V_b$ | Voltaje de ruptura (breakdown) | V |
| $p$ | Presión del gas | Torr |
| $d$ | Distancia entre electrodos | cm |
| $A$ | Coeficiente de ionización primaria | (Torr·cm)⁻¹ |
| $B$ | Coeficiente de excitación | V/(Torr·cm) |
| $\gamma$ | Coeficiente de Townsend secundario | adimensional |

El **mínimo analítico** de la curva (punto óptimo de operación):

$$\left(p \cdot d\right)_{\min} = \frac{e \cdot \ln\!\left(1 + \frac{1}{\gamma}\right)}{A} \qquad V_{b,\min} = B \cdot (p \cdot d)_{\min}$$

**Condición física del denominador** (verificada en código):

$$A \cdot p \cdot d > \ln\!\left(1 + \frac{1}{\gamma}\right)$$

### Mezcla de Gases — Ponderación Molar

Para una mezcla de gases, los coeficientes se obtienen por fracción molar:

$$A_{\text{mezcla}} = \sum_i x_i \cdot A_i \qquad B_{\text{mezcla}} = \sum_i x_i \cdot B_i$$

### Masa de Gas — Ley de Gas Ideal

Cantidad de gas para alcanzar la presión objetivo en el volumen total del sistema:

$$n = \frac{P_{\text{obj}} \cdot V_{\text{total}}}{R \cdot T} \qquad m = n \cdot PM_{\text{mezcla}}$$

Con $R = 62.364\,\text{L·Torr/mol·K}$ y $T = 298\,\text{K}$.

### Constantes de Ionización Incorporadas

| Gas | A [1/(Torr·cm)] | B [V/(Torr·cm)] | PM [g/mol] |
|-----|:-:|:-:|:-:|
| H₂  | 5  | 130 | 2.016  |
| N₂  | 12 | 342 | 28.014 |
| CO₂ | 20 | 466 | 44.01  |
| He  | 3  | 34  | 4.003  |
| Hg  | 20 | 370 | 200.59 |
| Air | 15 | 365 | 28.97  |

---

## 🗂️ Estructura del Repositorio

```
laser-regassing-paschen-tool/
│
├── paschen_calculator.py      # CLI: calculadora de terminal (standalone)
├── paschen_app_v12.py         # GUI: aplicación de escritorio con Flet + SCADA
├── manifold.png               # Diagrama del sistema de manifold (requerido por la GUI)
│
├── requirements.txt           # Dependencias Python
├── .gitignore
└── README.md
```

---

## ⚙️ Instalación

```bash
# 1. Clonar el repositorio
git clone https://github.com/WilbertMiguel/laser-regassing-paschen-tool.git
cd laser-regassing-paschen-tool

# 2. Crear entorno virtual (recomendado)
python -m venv venv
source venv/bin/activate        # Linux / macOS
venv\Scripts\activate           # Windows

# 3. Instalar dependencias
pip install -r requirements.txt
```

---

## 🚀 Uso

### Versión CLI — Terminal

```bash
python paschen_calculator.py
```

Menú interactivo en terminal con 4 opciones:

```
=======================================================
    CALCULADORA DE LA LEY DE PASCHEN
=======================================================
  1. Conversión de unidades de presión
  2. Calcular A y B de una mezcla de gases
  3. Calcular presión de operación a partir de Vb
  4. Agregar gas personalizado
  0. Salir
=======================================================
Selecciona una opción:
```

### Versión GUI — Aplicación de Escritorio

> ⚠️ `manifold.png` debe estar en la misma carpeta que `paschen_app_v12.py`.

```bash
python paschen_app_v12.py
```

---

## 🖥️ Funcionalidades — GUI (v12)

La aplicación tiene 6 módulos accesibles desde la barra de navegación:

### 🔄 Conversión de Presión
Convierte entre Pa, Torr, atm, mbar, psi, mmHg, kg/cm², inHg.

### ⚗️ Mezcla de Gases
- Agregar gases por porcentaje
- Calcula $A$ y $B$ ponderados por fracción molar
- Exporta resultados a `.txt`
- Botón **→ Paschen** y **→ Simulador** para transferir valores

### ⚡ Paschen
- Calcula $V_b$ dado $p$, $d$, $A$, $B$, $\gamma$
- Calcula $p$ dado $V_b$ (bisección numérica, error < 10⁻⁸ V)
- Grafica la curva completa $V_b$ vs $p \cdot d$ con tema oscuro/claro
- Botón **→ Simulador** para transferir $A$, $B$, $\gamma$

### 🧪 Gases Personalizados
Agregar, editar y eliminar gases con constantes $A$, $B$ y peso molecular propios.

### 🔬 Simulador Re-Gassing (SCADA)
Simulación física completa del proceso de re-gassing con **4 fases**:

| Fase | Descripción | Duración aprox. |
|------|-------------|:---:|
| **A** | Evacuación exponencial + prueba de fugas (60 s) | ~5 min |
| **B** | Flushing: 3 ciclos de inyección a 50 Torr + evacuación | ~18 min |
| **C** | Llenado final controlado hasta presión objetivo | ~5 min |
| **D** | Estabilización y verificación de micro-fugas | 10 min |

Incluye:
- **Panel SCADA** con indicadores en tiempo real (presión, fase, tiempo, masa de gas, $V_b$)
- **Animación matplotlib** a velocidad ×1/×3/×5/×10
- **Control manual de válvulas** V1, V2, V3 con validación de secuencia
- **4 gráficas simultáneas:** Presión vs Tiempo, Curva de Paschen, Característica V-I, Composición de mezcla
- **Curvas suavizadas** con interpolación spline en escala logarítmica
- Exporta gráfica y reporte a escritorio

#### Parámetros configurables del simulador

| Parámetro | Default | Descripción |
|-----------|---------|-------------|
| Longitud tubo | 140 cm | Longitud del tubo RECI |
| Diámetro tubo | 8 cm | Diámetro interno |
| Volumen manual | 5500 cm³ | Volumen medido real |
| d electrodos | 125 cm | Distancia entre electrodos |
| P objetivo | 15 Torr | Presión de operación |
| P atmosférica | 760 Torr | Presión ambiente |
| P vacío | 0.001 Torr | Vacío objetivo (10⁻³ Torr) |
| V ignición | 24 kV | Voltaje de arranque |
| V operación | 18 kV | Voltaje de plasma estable |
| I máxima | 28 mA | Corriente máxima |
| γ (Townsend 2°) | 0.01 | Coef. emisión secundaria |
| Tasa evacuación | 25 Torr/s | Velocidad de la bomba |
| Tasa llenado | 0.15 Torr/s | Apertura de válvula de aguja |

### 📖 Guía de Uso
Documentación embebida con ecuaciones LaTeX renderizadas, descripción de válvulas y procedimiento paso a paso.

---

## 🎨 Interfaz

- Tema **oscuro / claro** con paleta coherente en toda la aplicación
- Cambio de tema en tiempo real (textos, gráficas y campos se actualizan automáticamente)
- Ecuaciones LaTeX renderizadas con Matplotlib
- Exportación de gráficas `.png` y reportes `.txt` al escritorio

---

## 📦 Tecnologías

| Herramienta | Uso |
|-------------|-----|
| `Python 3.10+` | Lenguaje principal |
| `Flet` | Framework GUI multiplataforma |
| `Matplotlib` | Gráficas científicas + renderizado LaTeX |
| `NumPy` | Cálculo numérico y arrays |
| `math` (stdlib) | Funciones matemáticas (bisección, exponencial) |
| `threading` | Animación SCADA en hilo separado |
| `tempfile` | Caché de imágenes LaTeX |

---

## 🔧 Compatibilidad

| SO | Probado |
|----|---------|
| Windows 10/11 | ✅ |
| Linux (Ubuntu) | ✅ |
| macOS | ⚠️ No probado |

Python 3.10 o superior. Flet 0.21+.

---

## 📄 Licencia

MIT © 2025 Wilbert Miguel Nahuatlato

---

## 👤 Autor

**Wilbert Miguel Nahuatlato**  
Ingeniero Mecatrónico  
[![GitHub](https://img.shields.io/badge/GitHub-WilbertMiguel-181717?style=flat-square&logo=github)](https://github.com/WilbertMiguel)

---

> *Proyecto desarrollado como herramienta de laboratorio para el proceso de rehabilitación y re-gassing de tubos láser CO₂ sellados (Reci W4, 100–130 W).*
