## CAPÍTULO 1 - SISTEMA DE VACÍO Y CINÉTICA DE GASES EN RESONADORES RECI

### La Situación y el Objetivo

La rehabilitación de un resonador molecular de $\text{CO}_2$ (arquitectura RECI W4) mediante el proceso de *re-gassing* constituye un desafío de ingeniería que integra la cinética de gases, la mecánica y termodinámica del vacío, y la física de descargas luminiscentes (*glow discharges*). La construcción de una estación de rehabilitación para modelos de alta exigencia requiere una integración interdisciplinaria de la ciencia de materiales, la metrología de vacío ultra-alto (UHV) y la dinámica de fluidos en regímenes de transición.

El objetivo central de la estación es recrear la atmósfera interna original, eliminando contaminantes críticos como el vapor de agua y el oxígeno residual que envenenan los sitios activos del catalizador interno. Este proceso no es una simple recarga, sino una operación de descontaminación profunda y acondicionamiento químico y térmico de la cavidad de borosilicato. El objetivo primordial es restaurar la inversión de población mediante el control preciso del medio activo, minimizando las pérdidas térmicas y garantizando que la transición dieléctrica del gas no derive en un arco eléctrico destructivo.

El éxito del diseño de la estación radica en no concebirla como una simple línea de llenado, sino como un puente de conmutación cinético que permita alternar entre la evacuación profunda y la inyección controlada sin romper la estanqueidad del sistema.

---

### Dinámica del Medio Activo: Cinética Molecular y Sinergia de la Mezcla Gaseosa

En la ingeniería de plasmas, el medio activo no es un fluido estático, sino un sistema cinético fuera del equilibrio donde la transferencia de energía vibracional dicta la eficiencia cuántica del dispositivo. El gas actúa como un transductor de energía cinética electrónica a energía vibracional molecular, requiriendo una mezcla estequiométrica precisa donde cada componente desempeña una función crítica en el ciclo de bombeo y relajación:

* **Dióxido de Carbono ($\text{CO}_2$):** Actúa como el medio de ganancia. Representa típicamente entre el **5%** y el **10%** de la mezcla total. Su función es la emisión de fotones a $10.6\,\mu\text{m}$ mediante la transición $(001) \to (100)$.
* **Nitrógeno ($\text{N}_2$):** Es el reservorio de energía. Debido a que el primer nivel vibratorio del nitrógeno ($v=1$) está en resonancia casi exacta con el nivel (001) del $\text{CO}_2$ (discrepancia de solo $18\,\text{cm}^{-1}$), la transferencia de energía por colisión es extremadamente eficiente. Se excita fácilmente por impacto electrónico debido a su gran sección eficaz colisional.
* **Helio ($\text{He}$):** Es fundamental para la conductividad térmica. Constituye el componente mayoritario, aproximadamente el **75%** a **80%**. Su alta conductividad disipa el calor hacia las paredes del tubo, y aumenta la tasa de de-excitación del nivel inferior del láser (estado 010), eliminando el embotellamiento energético que cesaría la emisión si la temperatura superara los **200°C**.
* **Hidrógeno ($\text{H}_2$) o Vapor de Agua ($\text{H}_2\text{O}$):** Presentes en trazas ($\sim$ **1%**), actúan como catalizadores para la recombinación de los productos de la disociación. Bajo impacto electrónico, el $\text{CO}_2$ se disocia en $\text{CO}$ y $\text{O}$; el $\text{H}_2$ facilita el retorno a $\text{CO}_2$.
* **Xenón ($\text{Xe}$):** Se añade en trazas para reducir el potencial de ionización de la mezcla, disminuyendo la temperatura electrónica y la disociación molecular.

La reacción de transferencia resonante fundamental entre los componentes principales se describe mediante la ecuación:

$$\text{N}_2(v=1) + \text{CO}_2(000) \to \text{N}_2(v=0) + \text{CO}_2(001) - 18\,\text{cm}^{-1}$$

Para comprender la termodinámica de este ecosistema sinérgico, es imperativo utilizar el Modelo de las Tres Temperaturas, donde se definen reservorios energéticos que interactúan mediante colisiones de primera y segunda especie:

| Reservorio de Energía | Símbolo | Descripción Cinética |
| --- | --- | --- |
| Traslacional-Rotacional | $T$ | Energía cinética de las moléculas neutras; define la temperatura del gas. |
| Vibracional (Modo $v_3$) | $T_3$ | Asociado al nivel superior del láser (001); excitado por colisión con $\text{N}_2$. |
| Vibracional (Modos $v_1, v_2$) | $T_{12}$ | Niveles inferiores del láser; acoplados por resonancia de Fermi. |
| Vibracional del Nitrógeno | $T_4$ | Almacén metaestable de energía que bombea al $\text{CO}_2$ por colisión resonante. |

La función de partición vibracional ($Q$) para el sistema en equilibrio local se define como:

$$Q(q,r,s) = [(1-q)(1-r)^2(1-s)]^{-1}$$

Donde $q, r, s$ son los factores de Boltzmann asociados a las frecuencias fundamentales de los modos vibracionales del $\text{CO}_2$.

---

### Impacto de la Alteración de Proporciones en la Eficiencia y Cinética

Experimentar con proporciones no convencionales altera drásticamente la función de distribución de energía electrónica (EEDF) y la ganancia del láser. Aunque la mezcla estándar suele situarse en rangos de 1:1:3 a 1:1:8 ($\text{CO}_2$:$\text{N}_2$:$\text{He}$), los cambios producen los siguientes efectos cinéticos:

* **Exceso de $\text{CO}_2$:** Incrementa el valor de $E/N$ cuasi-estacionario. Los electrones pierden más energía en excitaciones vibracionales, elevando la temperatura del gas, incrementando la razón de decaimiento del nivel superior y colapsando la ganancia.
* **Exceso de $\text{N}_2$:** Favorece el bombeo del nivel superior, pero una concentración desproporcionada inestabiliza el plasma, provocando tonalidades incorrectas (púrpura/blanco) y aumentando el riesgo de transición a arco destructivo.
* **Reducción de $\text{He}$:** Provoca un embotellamiento en el nivel inferior (010). La inversión de población se reduce significativamente debido al aumento del factor de Boltzmann inducido por el calor no disipado.
* **Variación de $\text{H}_2$:** Su ausencia total acelerará la disociación del $\text{CO}_2$ en un sistema sellado, reduciendo la potencia tras pocas horas. Un exceso puede actuar como supresor de electrones.

---

### Análisis de la Mezcla Industrial y Seguridad Termodinámica

Para aplicaciones de alto rendimiento, se emplean mezclas pre-certificadas como la línea *SparkLaser® Premix* de Grupo Infra, garantizando una relación optimizada (cercana a 1:1.7:8.4). Desde la perspectiva de la seguridad industrial, el manejo de estos cilindros en una estación de re-gassing exige protocolos rigurosos.

* **Clasificación de Riesgo:** Gas comprimido, no inflamable, asfixiante simple. Desplaza el oxígeno atmosférico por debajo del umbral de supervivencia (**19.5%** $\text{O}_2$), provocando asfixia sin advertencia.
* **Peligro Termodinámico:** El calentamiento del cilindro puede provocar un incremento exponencial de la presión interna. Se deben evitar temperaturas superiores a **50°C** en el área de almacenamiento para prevenir explosiones mecánicas.

---

### Niveles de Vacío Críticos, Dinámica de Fluidos y Descontaminación

En la ingeniería de láseres gaseosos, el límite entre una rehabilitación exitosa y un fallo operativo está definido por el vacío final alcanzado. Si bien una bomba mecánica rotatoria puede alcanzar presiones de $10^{-3}\,\text{Torr}$ (limpieza gruesa), es estrictamente necesario llegar al régimen de vacío ultra-alto (UHV), situándose en el rango de $10^{-6}$ a $10^{-7}\,\text{Torr}$. A estos niveles, la trayectoria libre media ($\lambda$) de las moléculas es lo suficientemente larga como para asegurar que la probabilidad de colisión con impurezas residuales de oxígeno ($\text{O}_2$) sea despreciable durante la ignición.

En este régimen de transición, la dinámica de los gases está dictada por el número de Knudsen ($Kn = \lambda/d$). Cuando el sistema entra en flujo molecular libre ($Kn > 0.5$), el transporte ya no es viscoso, sino que depende de las colisiones con las paredes de borosilicato, limitando la conductancia y requiriendo tiempos prolongados de evacuación para vencer la desgasificación (*outgassing*) de las superficies.

#### Envenenamiento del Catalizador Interno

La infiltración de vapor de agua y $\text{O}_2$ representa una barrera termodinámica significativa. Estos contaminantes fisisorbidos bloquean los sitios activos del catalizador interno de los tubos RECI (basado en nanopartículas de oro dopadas o depósitos de platino sobre dióxido de estaño, $\text{SnO}_2$), impidiendo la reacción exergónica:

$$2\text{CO} + \text{O}_2 \to 2\text{CO}_2$$

Esto deriva en una acumulación de monóxido de carbono y el colapso de la inversión de población. Para regenerar activamente el sistema se requiere dopaje con hidrógeno (añadiendo una fracción molar cercana al **1%** de $\text{H}_2$) durante la primera descarga de limpieza para reducir los óxidos superficiales.

---

### Arquitectura del Módulo de Re-gassing y Gestión de Fluidos

Para procesar un tubo de entrada única, el diseño debe centrarse en un *manifold* (colector) en forma de "puente de conmutación" de tres vías, minimizando el volumen muerto para evitar atrapar aire atmosférico.

* **Cuerpo del Manifold:** Se recomienda topología en "T" fabricada en acero inoxidable 316L o vidrio borosilicato (Pyrex) de pared gruesa por sus bajas tasas de *outgassing*. El vidrio ofrece la ventaja de ser aislante eléctrico natural.
* **Válvulas de Control:** Uso imperativo de válvulas micrométricas de aguja (transición suave de presión) o de fuelle de alta estanqueidad.
* **Interfaces Metal-Vidrio:** La conexión entre el sistema de acero y el tubo debe realizarse mediante un rabillo de transición de *vidrio de uranio*. Posee un coeficiente de expansión térmica intermedio entre el borosilicato y el tungsteno/Kovar de los electrodos, minimizando fracturas por fatiga térmica.
* **Sellado y Lubricación:** Juntas de Viton o polímeros fluorocarbonados, epóxicos de baja presión de vapor (Torr-Seal). Se prohíben las grasas de hidrocarburos; usar exclusivamente grasa de silicona tipo Apiezon.

| Componente | Material Recomendado | Razón de Ingeniería |
| --- | --- | --- |
| Cuerpo del Manifold | Acero Inox 316L / Pyrex | Baja porosidad y resistencia química. |
| Válvulas de Control | Aguja Micrométrica (Whitley) | Control de presiones parciales a nivel de mTorr. |
| Sellos | O-rings de Viton / Torr-Seal | Estanqueidad en vacío profundo. |
| Mirillas | ZnSe (Selenuro de Zinc) | Alta transmitancia a $10.6\,\mu\text{m}$ y visibilidad. |

---

### Recomendación de Instrumentación y Metrología de Vacío

La precisión en la medición dicta la ganancia del medio activo. La presión final del sistema debe situarse por debajo de la presión atmosférica para aprovechar la ventaja del "autosellado" mecánico sobre los espejos de ZnSe. La presión absoluta se rige por:

$$P_{abs} = P_{atm} - P_{vacío}$$

Se requieren dos niveles de medición:

1. **Monitoreo de Limpieza:** Sensor Pirani ($10^2$ a $10^{-3}\,\text{Torr}$) para vacío primario, y un Ion Gauge (Bayer-Alpert) para regímenes UHV ($10^{-6}\,\text{Torr}$).
2. **Control de Mezcla:** Manómetro de capacitancia o diafragma con resolución de $0.5\,\text{Torr}$. Para tubos RECI, el punto óptimo se documenta consistentemente entre **18 y 20 Torr**.

| Instrumento | Rango de Aplicación | Función Crítica |
| --- | --- | --- |
| Pirani Gauge | $10^{-3}$ a $100\,\text{Torr}$ | Monitoreo de limpieza gruesa y llenado final. |
| Ionization Gauge | $10^{-3}$ a $10^{-10}\,\text{Torr}$ | Verificación de vacío ultra-alto y desgasificación. |
| Manómetro Absoluto | $0$ a $760\,\text{Torr}$ | Lectura visual y control fino de la mezcla. |

---

### Protocolo Operativo: Control, Lavado Gaseoso y Gestión Termodinámica

El procedimiento debe ser secuencial para mitigar choques térmicos o implosiones mecánicas (es obligatorio el uso de pantallas de policarbonato por el estrés de compresión del vidrio).

1. **I: Evacuación y Desgasificación:** Descenso asintótico de presión hasta $10^{-6}\,\text{Torr}$ por al menos 30 minutos.
2. **II: Purga y Lavado Gaseoso por Nitrógeno:** Inyectar $\text{N}_2$ de ultra-alta pureza hasta $100\,\text{Torr}$ y re-evacuar (3 ciclos). La fracción de impurezas remanentes se rige por: $f_n = (P_{residual} / P_{purga})^n$.
3. **III: Plasma Cleaning (Pasivación):** Establecer una descarga luminiscente de baja corriente ($I < 5\,\text{mA}$) con resistencia de balasto. El bombardeo iónico desprende óxidos del cátodo.
4. **IV: Llenado Controlado y Acondicionamiento (Aging):** Abrir la válvula micrométrica extremadamente lento ("ladle") considerando el Efecto Joule-Thomson (el enfriamiento adiabático de la válvula altera la densidad local). Si el gas entra muy frío, al calentarse en el resonador elevará su presión según la Ley de Gay-Lussac: $\frac{P_1}{T_1} = \frac{P_2}{T_2}$. El llenado debe detenerse en $\sim$ **18 Torr**.

Para asegurar la estabilidad durante estas pruebas de ignición y envejecimiento, el sistema debe integrar refrigeración por agua de circuito cerrado (**2 a 3 L/min** a **15-18°C**), ya que el calor altera el balance energético isocórico ($pV = NkT$), reduciendo la densidad local ($N$) y elevando el voltaje de sostenimiento.

---

### Fundamentos Cinéticos de la Ruptura (Ignición) y Criterio de Paschen

Desde la perspectiva de la ingeniería de potencia, el estado inicial del tubo es un dieléctrico gaseoso. Al aplicar un campo eléctrico, los electrones libres residuales inician una avalancha dictada por el criterio de Townsend:

$$\gamma (e^{\alpha d} - 1) = 1$$

Donde $\alpha$ es el coeficiente de ionización primaria, $d$ es la distancia interelectródica y $\gamma$ es la emisión secundaria del cátodo. La transición hacia el plasma (ruptura) exige un voltaje masivo ($V_r$) modelado por la **Ley de Paschen**:

$$V_r = \frac{B pd}{\ln(A pd) - \ln[\ln(1 + \frac{1}{\gamma})]}$$

Para determinar la presión de llenado ideal ($P_{ideal}$) en mezclas cuaternarias, los coeficientes empíricos $A$ y $B$ no son constantes universales, sino promedios ponderados por la fracción molar ($x_i$) de los constituyentes para calcular las colisiones ionizantes por unidad de longitud ($A_{mix}$) y la pérdida de energía no ionizante ($B_{mix}$):

$$A_{mix} = \sum_{i=1}^{n} x_i A_i \quad ; \quad B_{mix} = \sum_{i=1}^{n} x_i B_i$$

| Gas | $A$ ($\text{cm}^{-1}\text{Torr}^{-1}$) | $B$ ($\text{V}\cdot\text{cm}^{-1}\text{Torr}^{-1}$) | Rol Eléctrico en la Mezcla |
| --- | --- | --- | --- |
| Dióxido de Carbono ($\text{CO}_2$) | 20 | 466 | Incrementa el $E/N$ operativo. |
| Nitrógeno ($\text{N}_2$) | 12 | 342 | Facilita descarga estable a voltajes menores. |
| Helio ($\text{He}$) | 3 | 34 | Reduce la impedancia (diluyente eléctrico). |
| Hidrógeno ($\text{H}_2$) | 5 | 130 | Trazas catalíticas; mínimo impacto en Paschen. |

A presiones operativas de $\sim$ **20 Torr**, el tubo RECI se encuentra a la derecha del "mínimo de Paschen" ($pd \approx 2.5$ a $4\,\text{Torr}\cdot\text{cm}$). Cualquier incremento de presión elevará drásticamente el umbral de ignición (frecuentemente superando los **15 kV** a **20 kV**), definiendo el límite de cumplimiento ("compliance") absoluto que deberá satisfacer el diseño de tu fuente conmutada (SMPS).
