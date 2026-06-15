# INGENIERÍA DE SISTEMAS LÁSER DE TRÁNSITO MOLECULAR: DISEÑO OPERATIVO DE ESTACIONES DE RE-GASSING Y ELECTRÓNICA DE POTENCIA EN CAVIDADES DE BOROSILICATO

---

## CAPÍTULO 1: DINÁMICA DEL PLASMA Y CINÉTICA GASODINÁMICA EN EL PROCESO DE RE-GASSING

### 1.1. Fenomenología de la Descarga Luminescente y Ruptura Dieléctrica No Convencional

La rehabilitación de un resonador molecular de $\text{CO}_2$ con arquitectura RECI W4 mediante técnicas de *re-gassing* exige la transición controlada de un gas neutro multiespecie a un estado de plasma parcialmente ionizado, clasificado termodinámicamente como un plasma frío o no térmico ($T_e \gg T_g$). En este régimen, la energía suministrada por el campo eléctrico acelera preferencialmente a los electrones libres debido a su alta movilidad en comparación con los iones pesados. Bajo condiciones experimentales no convencionales—caracterizadas típicamente por desviaciones en las presiones parciales nominales o por el uso de subsistemas de vacío de presupuesto limitado que operan en las fronteras del régimen de transición (flujo de Knudsen)—, la predicción del voltaje de ruptura dieléctrica ($V_b$) se aparta del comportamiento clásico lineal. La ley de Paschen, que gobierna matemáticamente este fenómeno, establece que el potencial de ignición es una función no lineal del producto de la presión hidrostática total ($p$) y la distancia interelectródica ($d$), corregida por los coeficientes de ionización de Townsend:

$$V_b = \frac{B \cdot p \cdot d}{\ln(A \cdot p \cdot d) - \ln\left[\ln\left(1 + \frac{1}{\gamma_{se}}\right)\right]}$$

Donde $A$ y $B$ son constantes empíricas dependientes de la composición molecular específica del gas, y $\gamma_{se}$ representa el segundo coeficiente de Townsend, el cual cuantifica la emisión electrónica secundaria en el cátodo por impacto iónico. Cuando se trabaja con restricciones de instrumentación, entender este balance es crítico: una estación de vacío que deje trazas microscópicas de vapor de agua u oxígeno alterará drásticamente las secciones eficaces de colisión electrónica ($\sigma_j$), elevando el umbral de $V_b$ y sometiendo a la fuente de alimentación conmutada (SMPS) a un estrés dieléctrico destructivo antes de consolidar el canal conductor.

Para clarificar el impacto de los contaminantes y la variación del gas de aporte durante el proceso de desalojo y recarga en el laboratorio, se deben identificar los siguientes vectores fenomenológicos:

* **Envenenamiento de Sitios Activos:** La presencia de $\text{H}_2\text{O}$ residual actúa como un sumidero energético de alta eficiencia que desexcita por colisión los niveles vibracionales del nitrógeno ($\text{N}_2$), bloqueando la transferencia resonante hacia el nivel superior del láser de $\text{CO}_2$ ($00^01$).
* **Electronegatividad del Oxígeno Molecular:** El $\text{O}_2$ generado por la disociación inherente del $\text{CO}_2$ ($2\text{CO}_2 \rightleftharpoons 2\text{CO} + \text{O}_2$) captura electrones libres para formar iones negativos estables, reduciendo la densidad electrónica efectiva ($n_e$) y alterando de forma abrupta la conductividad eléctrica del medio.
* **Efecto de Confinamiento Térmico:** El helio ($\text{He}$), gracias a su elevada conductividad térmica, desaloja el calor del eje óptico de la descarga hacia las paredes refrigeradas por agua de la chaqueta de borosilicato; una deficiencia de helio incrementa la temperatura rotacional del gas ($T_g$), ensanchando las líneas espectrales por efecto Doppler y disminuyendo drásticamente la ganancia óptica del medio activo.

---

### 1.2. Impedancia Dinámica del Plasma y Modelado Equivalente

Una vez superado el transitorio de ruptura dieléctrica, el gas experimenta un colapso de impedancia en la escala de los nanosegundos, transitando de una resistencia cuasi-infinita a un estado de conducción con una impedancia caracterizada por una pendiente fuertemente negativa. Este fenómeno de resistencia dinámica negativa ($\mathbb{R}_{dyn} = \frac{dV}{dI} < 0$) define la región de descarga luminescente normal y anormal donde opera el tubo RECI W4. Matemáticamente, el plasma no se comporta como una carga óhmica lineal, sino como una fuente de contratensión no lineal en serie con una resistencia remanente de carácter fuertemente dinámico.

La ecuación analítica que describe el voltaje operativo del tubo ($V_{tube}$) en función de la corriente de descarga ($I_d$) se parametriza mediante la relación empírica de descarga en gases estables:

$$V_{tube} = V_c + \frac{K \cdot p \cdot L}{\sqrt{I_d}}$$

Donde $V_c$ representa la caída de potencial catódica (una constante ligada al material del electrodo y al tipo de gas), $L$ es la longitud efectiva de la columna de plasma positiva, y $K$ es un factor de escala geométrico y composicional.

```
   Ignición (V_b) 
        ▲
        │      ═══ Balastrado Requerido
        │     /
        │    /     Punto de Operación Estable
Voltaje │   /     /
  (V)   │  ▼     ▼
        │  * .
        │      . * . _ _ _ _ _ _ _ _ _ _ _ (Carga del Plasma)
        │                                  R_dyn < 0
        └────────────────────────────────────────► Corriente (I)

```

Para propósitos de diseño electrónico con herramientas accesibles de simulación, la carga plasmática se modela en régimen permanente mediante un circuito equivalente no lineal, cuyos parámetros se sintetizan en la siguiente tabla de diseño operacional:

| Parámetro Eléctrico del Plasma | Expresión Matemática / Valor Típico | Impacto en la Topología de la Fuente | Método de Validación Económico |
| --- | --- | --- | --- |
| **Resistencia Estática ($R_{st}$)** | $R_{st} = \frac{V_{op}}{I_{op}} \approx 1.5 \text{ M}\Omega \text{ a } 2 \text{ M}\Omega$ | Define el punto de polarización inicial y el dimensionamiento de aislamiento. | Voltímetro electrostático artesanal o divisor resistivo calibrado de alta disipación. |
| **Resistencia Dinámica ($R_{dyn}$)** | $R_{dyn} = \frac{\partial V_{tube}}{\partial I_d} $  | ${I_{op}} \approx -20 \text{ k}\Omega \text{ a } -50 \text{ k}\Omega$ | Provoca oscilaciones destructivas en lazo abierto si la fuente actúa en modo voltaje. |
| **Voltaje de Sostenimiento ($V_{sus}$)** | $V_{sus} \approx 15 \text{ kV} - 22 \text{ kV}$ | Determina la relación de transformación mínima del secundario en régimen estacionario. | Monitoreo directo con sonda de alto voltaje comercial o construida con resistores de película de carbón. |
| **Capacitancia de Envoltura ($C_{sheath}$)** | $C_{sheath} = \frac{\varepsilon_0 A_{sh}}{d_{sh}} \approx 5 \text{ pF} - 15 \text{ pF}$ | Introduce desfases reactivos en esquemas de modulación por ancho de pulso a alta frecuencia. | Cálculo analítico mediante la geometría del cátodo y medición indirecta por tiempo de subida. |

---

## CAPÍTULO 2: INGENIERÍA DE POTENCIA AVANZADA: DISEÑO DE TOPOLOGÍAS DE CONMUTACIÓN DE BAJO COSTO Y ALTA ESTABILIDAD

### 2.1. Topologías Conmutadas (SMPS) y Control por Lazo Cerrado con TL494

El núcleo de la electrónica de potencia requerida para estabilizar una descarga con resistencia dinámica negativa consiste en una fuente conmutada configurada intrínsecamente como una fuente de corriente constante. Para un proyecto universitario con restricciones severas de costo, el uso de arquitecturas de puente completo (Full-Bridge) o medio puente (Half-Bridge) controladas por el circuito integrado clásico TL494 representa la solución óptima en términos de coste, robustez y comprensibilidad pedagógica. La topología Half-Bridge aprovecha la división de tensión capacitiva del bus de continua, reduciendo a la mitad el estrés de voltaje de bloqueo en los MOSFETs de potencia en comparación con una topología Push-Pull. El transformador elevador de alta frecuencia, que opera en el rango de los $20 \text{ kHz}$ a $50 \text{ kHz}$, transfiere energía hacia un bloque multiplicador de voltaje encargado de alcanzar los $>35 \text{ kV}$ necesarios para la ignición.

El lazo de control debe configurarse estrictamente en Modo de Corriente de Pico (Peak Current Mode Control) o, en su defecto mediante el TL494, un Modo de Control por Corriente Promedio implementado a través de sus amplificadores de error internos. La función de transferencia del lazo de control se ve severamente afectada por la resistencia negativa de la descarga, lo que exige una compensación de tipo II (polo en el origen y un par cero-polo) para estabilizar el sistema. La corriente de descarga se sensa en el retorno del cátodo a través de una resistencia no inductiva ($R_{sense}$), generando un voltaje de realimentación que se compara con una referencia fija, según se define en la siguiente ecuación cinético-eléctrica de control:

$$V_{feedback} = I_d \cdot R_{sense} \implies \Delta D = A_{v(ol)} \cdot \left( V_{ref} - I_d \cdot R_{sense} \right)$$

Donde $\Delta D$ es la variación del ciclo de trabajo (*duty cycle*) de los pulsos de compuerta y $A_{v(ol)}$ es la ganancia en lazo abierto del amplificador de error configurado. Si el ciclo de trabajo no responde de forma instantánea a los incrementos microscópicos de corriente, la descarga luminescente transitará exponencialmente hacia un arco eléctrico concentrado, destruyendo el cátodo por evaporación térmica localizada.

* **Configuración del TL494 para Control de Corriente:** Los dos amplificadores de error del circuito integrado deben usarse de forma segregada: uno dedicado exclusivamente al lazo de regulación estricta de corriente constante (pines 1 y 2), y el segundo configurado como un lazo de protección contra sobrevoltaje en circuito abierto (pines 15 y 16) para evitar el colapso del aislamiento del transformador si el tubo no enciende.
* **Aislamiento del Driver de Compuerta:** Dadas las limitaciones de presupuesto, la implementación de transformadores de pulso artesanales bobinados sobre toroides de ferrita recuperados de fuentes ATX viejas es preferible frente a los circuitos integrados *drivers* flotantes dedicados (como la serie IR2110), los cuales son propensos a destruirse por enganche de baja impedancia (*latch-up*) ante los ruidos electromagnéticos (EMI) de la descarga.
* **Diseño del Transformador Elevador:** El devanado secundario debe seccionarse en múltiples galletas físicas para minimizar la capacitancia parásita distribuida ($C_{para}$); de lo contrario, la corriente de desplazamiento resonará con la inductancia de fuga ($L_{leak}$), disipando la energía útil en forma de calor sobre los interruptores de silicio del primario.

---

### 2.2. Mitigación de Transitorios mediante Redes de Amortiguamiento y Filtrado de Rizado

El instante exacto en que ocurre la ruptura dieléctrica genera un transitorio electromagnético severo que se propaga tanto de forma radiada como conducida hacia las etapas de baja tensión de la electrónica de control. La inductancia de fuga del transformador de potencia ($L_{leak}$), al no encontrar una descarga óhmica simétrica durante los primeros nanosegundos del colapso del gas, induce picos de sobrevoltaje destructivos en los drenadores de los MOSFETs primarios debido al fenómeno de almacenamiento de energía magnética residual ($E_m = \frac{1}{2} L_{leak} I_{peak}^2$). Para proteger el silicio sin recurrir a costosos componentes activos de grado militar, se requiere el diseño riguroso de redes de amortiguamiento pasivas RCD (Snubbers). La constante de tiempo de la red de amortiguamiento ($\tau = R_{snub} \cdot C_{snub}$) debe calcularse para ser significativamente menor que el tiempo de apagado ($t_{off}$) del interruptor, absorbiendo el transitorio inductivo de conmutación:

$$V_{snub(max)} = V_{bus} + I_{peak} \sqrt{\frac{L_{leak}}{C_{snub}}}$$

A la salida de la etapa multiplicadora de alta tensión, la estabilidad espacial y temporal del haz láser de $\text{CO}_2$ depende críticamente del factor de rizado de la corriente ($\Delta I_d / I_d$). Un rizado excesivo modula la densidad de población de los niveles energéticos excitados, degradando la calidad del haz óptico (parámetro $M^2$) y generando inestabilidades acústicas en la columna de plasma. La inclusión de un filtro LC de alta tensión en combinación con una resistencia de balasto física ($R_{balast}$) actúa amortiguando el remanente de impedancia negativa. El balasto proporciona un acoplamiento de carga balanceado, desplazando el punto de operación nulo del sistema hacia una región de estabilidad pasiva absoluta donde la pendiente total vista por la fuente conmutada sea estrictamente positiva.

```
       MOSFET Drenador
            │
            ├───┐
            │   ├──┐
           / \  │  │
    Diodo  ▲   █  █ Resistor (R_snub)
   Ultrafast│  █  █
            └───┼──┘
               --- Capacitor (C_snub)
                │
            ────┴──── GND Bus

```

Para implementar este filtrado y amortiguamiento bajo un esquema de optimización de recursos, se sugieren los siguientes lineamientos prácticos basados en fundamentos de diseño electrónico:

* **Selección de Capacitores en el Multiplicador:** El uso de capacitores cerámicos de disco de bajo costo tipo Y5V debe evitarse debido a su drástica pérdida de capacitancia bajo polarización de DC elevada; en su lugar, se deben recuperar o adquirir capacitores de polipropileno o cerámicos de grado industrial NPO/COG que mantengan su permitividad eléctrica estable.
* **Construcción de la Resistencia de Balasto:** Ante la imposibilidad de adquirir resistencias cerámicas especializadas de alta tensión y gran disipación, es técnicamente viable construir un arreglo en serie-paralelo de múltiples resistencias de película de óxido metálico estándar de $2\text{ W}$, montadas sobre una placa de circuito impreso espaciada y sumergidas en aceite mineral dieléctrico o resina epóxica para prevenir el arco de efecto corona superficial.
* **Blindaje Electromagnético de la Electrónica de Control:** La sección lógica del circuito (TL494 y amplificadores operacionales de acondicionamiento) debe confinarse dentro de un chasis metálico ferroso conectado a una tierra física independiente (tierra de protección); todas las señales de realimentación analógica provenientes del lado de alta tensión deben filtrarse localmente mediante redes pasivas pasabajas RC con frecuencias de corte situadas una década por debajo de la frecuencia de conmutación de la SMPS.








