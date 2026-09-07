<p align="center">
  <img src="logo/Macc.png" alt="Logo Universidad del Rosario - MACC" width="160"/>
</p>

<h1 align="center">Modelado de la Propagación de Epidemias en Redes Dinámicas Usando Teoría de Grafos</h1>

<p align="center">
  <em>Dinámica Epidémica SIS Estocástica, Caminos Causalmente Respetuosos del Tiempo (TBFS), Umbrales Epidémicos Espectrales e Identificación de Supercontagiadores basada en Centralidad sobre Redes de Contacto SocioPatterns.</em>
</p>

<p align="center">
  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/Python-3.9%2B-blue.svg?style=flat&logo=python" alt="Versión de Python"/></a>
  <a href="https://networkx.org/"><img src="https://img.shields.io/badge/NetworkX-3.0%2B-blueviolet.svg?style=flat" alt="NetworkX"/></a>
  <a href="https://jupyter.org/"><img src="https://img.shields.io/badge/Jupyter-Notebook-orange.svg?style=flat&logo=jupyter" alt="Jupyter"/></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/Licencia-MIT-green.svg?style=flat" alt="Licencia: MIT"/></a>
  <a href="https://urosario.edu.co/"><img src="https://img.shields.io/badge/Instituci%C3%B3n-Universidad%20del%20Rosario-red.svg?style=flat" alt="Universidad del Rosario"/></a>
</p>

<p align="center">
  <strong><a href="README.md">English Version</a></strong> | <strong><a href="README.es.md">Versión en Español</a></strong>
</p>

---

## 🔬 Resumen del Proyecto y Motivación

En la epidemiología real, las interacciones físicas entre individuos son intrínsecamente **dinámicas, variables en el tiempo y no estacionarias**. Los enfoques tradicionales que agregan contactos temporales en grafos estáticos introducen severas distorsiones estructurales: generan caminos topológicos que violan la flecha del tiempo y sobreestiman artificialmente la velocidad y el alcance de la propagación epidémica.

En este proyecto de investigación se modela la dinámica de propagación del **COVID-19** implementando un **proceso estocástico Susceptible-Infectado-Susceptible (SIS) sobre redes temporales**, utilizando registros reales de telemetría de proximidad física del conjunto de datos **SocioPatterns High-School** (Netzschleuder / SocioPatterns 2012).

### Principales Contribuciones:
1. **Discretización Temporal de Contactos:** Procesamiento de eventos continuos en ventanas discretas estandarizadas ($\Delta t = 20\text{ s}$), preservando la duración y el peso de las interacciones $w_t(u,v)$.
2. **Grafo Causal de Contagios ($G_{\text{inf}}$):** Reconstrucción del árbol dirigido exacto de transmisión, registrando pares infector-infectado y la estampa temporal $t_{\text{infect}}$.
3. **Algoritmo TBFS (Time-respecting Breadth-First Search):** Implementación de una búsqueda en anchura temporal que respeta estrictamente la causalidad para calcular tiempos mínimos de llegada ($\text{arrival}(v)$) y alcanzabilidad real desde el Paciente Cero.
4. **Análisis del Umbral Epidémico Espectral:** Comparación entre el umbral teórico estático ($\tau_c \approx 1/\lambda_{\max}(A_{\text{agg}})$) y la tasa efectiva simulada ($\tau = \beta / \delta$).
5. **Identificación de Supercontagiadores por Centralidad:** Aplicación de **PageRank**, **Betweenness Centrality** y **Grado de Salida (Out-Degree)** sobre $G_{\text{inf}}$, demostrando que el grado estático en $G_{\text{agg}}$ es insuficiente para predecir a los verdaderos difusores del virus.
6. **Evaluación de Coarse-Graining Temporal:** Cuantificación del sesgo en la alcanzabilidad al variar el tamaño del bin de agregación temporal ($\Delta t \in [1, 1000]$ ventanas).

---

## 📐 Formulación Matemática y Teórica

### 1. Representación de la Red Temporal
La red de contactos se modela como una secuencia ordenada de grafos instantáneos (snapshots):
$$\mathcal{G} = \{G_t = (V, E_t, w_t)\}_{t=1}^T$$
donde $V$ representa el conjunto de individuos, $E_t \subseteq V \times V$ son los contactos no dirigidos activos en el intervalo $[(t-1)\Delta t, t\Delta t]$, y $w_t(u,v)$ denota la duración acumulada del contacto.

### 2. Modelo Epidémico SIS Estocástico en Tiempo Discreto
Cada individuo $v \in V$ se encuentra en un estado $\sigma_v(t) \in \{S, I\}$ en el instante $t$:
- **Infección ($S \to I$):** Si el nodo $v$ es susceptible en $t$ ($\sigma_v(t) = S$), la probabilidad de contraer la infección en $t+1$ depende de sus vecinos infectados $N_t(v) \cap I_t$ y el tiempo de exposición:
  $$P(v \in I_{t+1} \mid v \in S_t) = 1 - \prod_{u \in N_t(v) \cap I_t} (1 - \beta)^{\max\left(1, \left\lfloor \frac{w_t(u,v)}{20} \right\rceil\right)}$$
  donde $\beta \in (0,1)$ es la probabilidad de transmisión por bloque temporal de 20s.
- **Recuperación ($I \to S$):** Un individuo infectado se recupera y vuelve a ser susceptible (sin inmunidad duradera) con probabilidad $\delta \in (0,1)$ en cada ventana temporal:
  $$P(v \in S_{t+1} \mid v \in I_t) = \delta$$

```
   +-------------------+       Tasa de Infección: P(infección | w_t, beta)       +------------------+
   |  Susceptible (S)  | ------------------------------------------------------> |   Infectado (I)  |
   +-------------------+                                                         +------------------+
             ^                                                                             |
             |                                Tasa de Recuperación: delta                  |
             +-----------------------------------------------------------------------------+
```

### 3. Grafo Causal Dirigido de Contagios ($G_{\text{inf}}$)
Cuando un nodo susceptible $v$ es infectado por el vecino $u$ en la ventana $t$, se crea una arista dirigida causal:
$$e = (u \to v, t_{\text{infect}}) \in E_{\text{inf}}$$
Se obtiene así el subgrafo acíclico de transmisión $G_{\text{inf}} = (V, E_{\text{inf}})$, preservando la trazabilidad epidemiológica completa.

### 4. Algoritmo TBFS (Temporal Breadth-First Search)
Para determinar si un nodo $v$ es alcanzable desde la raíz $v_0$ (Paciente Cero), las transiciones deben cumplir la condición de orden temporal no decreciente $t_1 \le t_2 \le \dots \le t_k$:
$$\text{arrival}(v) = \min \left\{ t \in [1, T] \mid \exists \text{ un camino temporalmente respetuoso de } v_0 \text{ a } v \text{ en o antes de } t \right\}$$
A diferencia de la distancia geodésica estática $d_{\text{estatica}}(v_0, v)$, TBFS garantiza estricta causalidad física.

### 5. Umbral Epidémico Espectral
En la matriz de adyacencia del grafo agregado estático $A_{\text{agg}}$, el umbral crítico está determinado por el autovalor principal $\lambda_{\max}(A_{\text{agg}})$:
$$\tau_c = \frac{1}{\lambda_{\max}(A_{\text{agg}})}$$
Al contrastarlo con $\tau = \frac{\beta}{\delta}$, se revela la discrepancia entre la teoría espectral estática y la evolución real sobre la red temporal no homogénea.

---

## 📊 Resultados y Visualizaciones

<p align="center">
  <img src="figs/sis_I_curve.png" alt="Curva de Infección SIS" width="700"/>
  <br/>
  <em>Figura 1: Evolución temporal de individuos infectados activos I(t), evidenciando ondas de rebrote y estabilización dinámica.</em>
</p>

### Snapshots de la Expansión de Contagios
Evolución de la red de transmisión en diferentes cortes temporales ($t = 40, 255, 471, 902$). El color de los nodos indica el tiempo de primera infección (escala $\text{plasma}$: amarillo/naranja $\to$ temprano, morado $\to$ tardío, gris $\to$ susceptible):

<p align="center">
  <img src="figs/snapshot_t40.png" width="45%" alt="Snapshot t=40"/>
  <img src="figs/snapshot_t255.png" width="45%" alt="Snapshot t=255"/>
  <br/>
  <img src="figs/snapshot_t471.png" width="45%" alt="Snapshot t=471"/>
  <img src="figs/snapshot_t902.png" width="45%" alt="Snapshot t=902"/>
  <br/>
  <em>Figura 2: Snapshots dirigidos de contagio. El grosor de las aristas y el color del nodo reflejan la secuencia cronológica.</em>
</p>

### Supercontagiadores y Centralidad PageRank ($G_{\text{inf}}$)
<p align="center">
  <img src="figs/infections_ALL_overlay_pagerank.png" alt="Superposición de PageRank en Grafo de Contagio" width="750"/>
  <br/>
  <em>Figura 3: Identificación del Top-10 de supercontagiadores mediante PageRank sobre el grafo dirigido de contagios G_inf, con halos dorados y etiquetas destacadas.</em>
</p>

---

## 📁 Estructura del Repositorio

```text
MODELADO-DE-LA-PROPAGACION-DE-EPIDEMIAS-EN-REDES-DINAMICAS/
├── .gitignore                      # Reglas de exclusión de Git (.venv, pycache, zips)
├── LICENSE                         # Licencia de código abierto MIT
├── README.md                       # Documentación científica en Inglés (MITACS-ready)
├── README.es.md                    # Documentación científica en Español
├── requirements.txt                # Dependencias de Python reproducibles
│
├── data/
│   ├── raw/                        # Telemetría cruda SocioPatterns High School
│   │   ├── edges.csv               # Interacciones temporales con estampas de tiempo
│   │   ├── gprops.csv              # Metadatos del grafo
│   │   └── nodes.csv               # Identificadores de nodos
│   └── processed/                  # Datos preprocesados (ventanas de 20s)
│       ├── contacts.csv            # (source, target, t_start, t_end, weight, directed)
│       └── nodes.csv               # (id, label, group)
│
├── docs/                           # Documentación académica del proyecto
│   ├── Paper_Modelacion_Epidemias_Redes_Dinamicas.pdf  # Artículo científico final (11 páginas)
│   ├── Poster_Proyecto_Teoria_de_Grafos.pdf            # Póster académico de investigación
│   └── Presentacion_Proyecto_Teoria_de_Grafos.pdf      # Diapositivas de sustentación
│
├── figs/                           # Figuras generadas en alta resolución
│   ├── sis_I_curve.png             # Curva de casos activos I(t)
│   ├── snapshot_t*.png             # Snapshots secuenciales de contagio
│   ├── infections_ALL_overlay_pagerank.png # Overlay de centralidad PageRank
│   └── tbfs_coarse_reach_vs_bin.png # Gráfica de coarse-graining temporal
│
├── logo/
│   └── Macc.png                    # Escudo institucional (UR - MACC)
│
├── media/
│   └── Grabacion_sis_animacion.mp4 # Video MP4 completo con la animación de la propagación
│
├── results/                        # Salidas tabulares del modelo
│   ├── infection_edges.csv         # Aristas dirigidas de contagio (source, target, t_infect)
│   ├── states_by_t.csv             # Conteo de estados (t, S, I)
│   ├── rank_contact_degree.csv     # Ranking por grado en red estática
│   ├── rank_infection_centrality.csv # Métricas de centralidad en G_inf (Out-degree, PageRank, Betweenness)
│   └── rank_pagerank.csv           # Ranking de PageRank
│
└── src/                            # Código fuente y Notebook Interactivo
    ├── centrality_rankings.py      # Cálculo de rankings de centralidad
    ├── prepare_robusto_2012.py     # Limpieza y discretización de contactos temporales
    ├── simulate_sis_and_infection_graph.py # Motor de simulación estocástica SIS
    ├── sis_animacion.ipynb         # Jupyter Notebook interactivo (Pipeline + TBFS + Animación)
    ├── viz_expansion.py            # Gráficas básicas de expansión de casos
    ├── viz_expansion_colored.py    # Generador de snapshots con mapa de color plasma
    ├── viz_expansion_directed.py   # Visualizador de aristas dirigidas ponderadas
    └── viz_overlay_centrality.py   # Renderizador de superposición de PageRank
```

---

## 🚀 Guía de Instalación y Reproducibilidad

### 1. Clonar el Repositorio
```bash
git clone https://github.com/JuanFeletes24/MODELADO-DE-LA-PROPAGACION-DE-EPIDEMIAS-EN-REDES-DINAMICAS.git
cd MODELADO-DE-LA-PROPAGACION-DE-EPIDEMIAS-EN-REDES-DINAMICAS
```

### 2. Configurar el Entorno Virtual
```bash
# Crear entorno virtual
python -m venv .venv

# Activar en Linux/macOS
source .venv/bin/activate

# Activar en Windows PowerShell
.venv\Scripts\Activate.ps1

# Instalar dependencias
pip install -r requirements.txt
```

### 3. Pipeline de Ejecución

Ejecuta el flujo científico completo en consola:

```bash
# Paso 1: Preprocesar la telemetría en ventanas discretas de 20s
python src/prepare_robusto_2012.py

# Paso 2: Ejecutar la simulación estocástica SIS en la red dinámica
python src/simulate_sis_and_infection_graph.py

# Paso 3: Calcular los rankings de centralidad sobre el grafo causal
python src/centrality_rankings.py

# Paso 4: Generar visualizaciones, snapshots y superposiciones de PageRank
python src/viz_expansion_colored.py
python src/viz_overlay_centrality.py
```

### 4. Jupyter Notebook Interactivo
Inicia el entorno Jupyter para explorar el código paso a paso y reproducir las animaciones interactivas:
```bash
jupyter notebook src/sis_animacion.ipynb
```

---

## 💻 Contenido del Notebook (`sis_animacion.ipynb`)

El notebook `src/sis_animacion.ipynb` constituye el núcleo demostrativo e interactivo del proyecto:
1. **Encabezado y Metadatos Académicos:** Identificación institucional y del grupo.
2. **Carga y Simulación SIS:** Ajuste de parámetros epidemiológicos ($\beta, \delta$, siembras).
3. **Definición de Layouts y Colores:** Layout de resortes basado en contactos con separación mínima forzada.
4. **Funciones de Dibujo de Snapshots:** Renderizado dinámico de la red de transmisión.
5. **Animación en Tiempo Real:** Reproductor cuadro a cuadro de la propagación viral.
6. **Búsqueda de Nodos Centrales:** Tablas de ranking por PageRank, Betweenness y Out-Degree.
7. **Algoritmo TBFS:** Implementación de BFS temporal para trazabilidad desde el Paciente Cero.
8. **Comparativa Distancia Estática vs. Temporal:** Demostración cuantitativa de la distorsión del grafo estático.
9. **Cálculo del Umbral Epidémico:** Determinación analítica de $\lambda_{\max}(A_{\text{agg}})$ y $\tau_c$.
10. **Análisis de Agregación Temporal ($\Delta t$):** Estudio del impacto de coarse-graining en alcanzabilidad.

---

## 👥 Autores y Afiliación Académica

**Matemáticas Aplicadas y Ciencias de la Computación (MACC)**  
*Escuela de Ingeniería, Ciencia y Tecnología*  
**Universidad del Rosario, Bogotá, Colombia**

- **Juan Felipe Amaya Montaña**
- **Juan Felipe Rojas Manjarres**
- **Sebastian Alejandro Sanchez Urrego**

*Proyecto desarrollado para la asignatura de Teoría de Grafos (2025-II) y presentado como portafolio académico para postulaciones a pasantías internacionales de investigación (MITACS Globalink Research Internship).*

---

## 📚 Referencias Bibliográficas

1. **Pastor-Satorras, R., Castellano, C., Van Mieghem, P., & Vespignani, A.** (2015). *Epidemic processes in complex networks*. Reviews of Modern Physics, 87(3), 925.
2. **Holme, P., & Saramäki, J.** (2012). *Temporal networks*. Physics Reports, 519(3), 97-125.
3. **Fournet, J., & Barrat, A.** (2014). *Contact patterns among high school students*. PLoS ONE, 9(9), e107878. SocioPatterns Collaboration.
4. **Kiss, I. Z., Miller, J. C., & Simon, P. L.** (2017). *Mathematics of Epidemics on Networks: From Stochastic to Deterministic Models*. Springer.
5. **Newman, M. E. J.** (2018). *Networks: An Introduction*. Oxford University Press.

---

## 📄 Licencia
Este proyecto se distribuye bajo los términos de la [Licencia MIT](LICENSE).
