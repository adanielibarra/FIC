# FIC Shoreline Change Analysis

[English](#english) · [Español](#español)

---

## English

QGIS plugin for shoreline change analysis, developed by
**Daniel Ibarra Marinas**, Facultad de Ingeniería y Ciencias, Universidad
Autónoma de Tamaulipas (Mexico).

### What it does

It casts transects perpendicular to a user-digitized baseline, every X
metres. Each transect is intersected with a multitemporal shoreline layer
(lines with a date field) and the following DSAS-style metrics are
computed:

- **NSM (Net Shoreline Movement, m):** distance between the oldest and the
  most recent shoreline crossing the transect.
- **EPR (End Point Rate, m/yr):** NSM divided by the years elapsed between
  the first and the last shoreline.
- **LRR (Linear Regression Rate, m/yr):** slope of the linear regression of
  all shoreline positions against time. With 3 or more dates, also its R²,
  the standard error of the estimate (LSE) and the confidence interval of
  the slope (LCI) at the chosen level (90, 95 or 99 %).
- **WLR (Weighted Linear Regression, m/yr)**, optional: like LRR, but dates
  with less positional error weigh more. Requires positional uncertainty.

Also:

- **Merges by date:** if the shoreline of a date is digitized in several
  pieces, they are merged before intersecting, so each date gives a single
  position per transect.
- **Multiple intersections:** if a transect crosses the shoreline of a date
  more than once, the intersection closest to the baseline is used (like
  the "closest" option in DSAS).
- **Reprojection:** shorelines in a different CRS are reprojected to the
  baseline CRS automatically.
- **DSAS-style sign:** you state whether the baseline is onshore or
  offshore, and each transect is oriented seaward. Positive rates mean
  accretion and negative rates mean erosion, whatever the digitizing
  direction of the baseline.

### Installation

1. Install it from the QGIS plugin manager (Plugins → Manage and Install
   Plugins), searching for "FIC Shoreline Change Analysis". Since it is
   experimental, first tick "Show also experimental plugins" in the
   Settings tab.
2. Or download the ZIP of the latest version and use "Install from ZIP"
   in the same manager.
3. It appears in the Processing Toolbox, under
   "FIC Shoreline Change Analysis" → "Coastal".

### Usage

You need two line layers:

- A **baseline**: a single continuous line, digitized by hand, in a
  **projected CRS in metres** (for example, the corresponding UTM zone).
  If it is in degrees or another unit, the algorithm stops and says so.
- A **multitemporal shoreline** layer: shorelines of different dates, with
  a Date field (or ISO text `YYYY-MM-DD`). Features with an invalid date
  are discarded and reported in the log.

The algorithm also asks for:

- Transect spacing and total transect length, both in metres. Transects
  are symmetric, half on each side of the baseline.
- The **baseline position**: onshore or offshore. The baseline must lie
  entirely on one side of all shorelines.
- The **smoothing distance** (m): the direction of each transect is
  measured over a stretch of baseline of this length, so transects follow
  the general shape of the coast rather than digitizing kinks. A
  reasonable starting point is 2 to 5 times the spacing. 0 = no smoothing.
- The **confidence level** of the intervals: 90 % (default), 95 % or 99 %.
- The **positional uncertainty** of the shorelines (optional), in metres.
  It is taken, in order of priority, from:
  1. A field with the **total error** already computed.
  2. Fields with its **components**, combined in quadrature (square root of
     the sum of squares): source image **resolution** (counted as 1 pixel),
     **georeferencing RMSE** and **digitizing error**. Missing ones count
     as 0.
  3. A **single value** for all shorelines.

  This way each date can have its own error: an old aerial photograph and a
  recent satellite image do not have the same accuracy. If a date is split
  into pieces with different errors, the largest is used. Resolution alone
  underestimates the error (the plugin warns about it): add at least the
  RMSE. Without uncertainty, the plugin works the same but does not
  compute EPRunc or WLR.

#### How the seaward side is decided

For each transect, the plugin looks at which side of the baseline the
closest intersection of each date falls on, and takes the majority side.
If the baseline is onshore, the sea is on that side; if offshore, on the
opposite one. Then, for each date, the closest intersection within that
side is used, so a shore on the other side of the baseline (for example,
a lagoon shore) does not get mixed into the calculation.

### Output

A line layer (the transects) with these fields:

| Field | Content |
|---|---|
| `id` | Transect identifier |
| `n_pts` | Number of dates crossing the transect |
| `yr_min`, `yr_max` | Decimal year of the first and the last date |
| `dist_min`, `dist_max` | Distance to the baseline (m) of the first and the last date, positive seaward |
| `NSM` | Net Shoreline Movement (m) |
| `EPR` | End Point Rate (m/yr) |
| `EPRunc` | EPR uncertainty (m/yr): square root of the sum of the squared uncertainties of the first and last date, divided by the years. Empty without positional uncertainty |
| `LRR` | Linear Regression Rate (m/yr) |
| `LRR_R2` | R² of the regression (empty with fewer than 3 dates) |
| `LSE` | Standard error of the estimate (m): how far, on average, positions lie from the trend line. Empty with fewer than 3 dates |
| `LCI90` (or `LCI95`, `LCI99`) | Half-width of the LRR confidence interval (m/yr). If the absolute LRR is smaller than LCI, the rate is not different from zero at that level. Empty with fewer than 3 dates |
| `WLR` | Weighted Linear Regression (m/yr): like LRR, but dates with less error weigh more (weights 1/U²). Empty with fewer than 3 dates or without uncertainty |
| `WR2`, `WSE` | R² and standard error of the estimate of the WLR |
| `WCI90` (or `WCI95`, `WCI99`) | Half-width of the WLR confidence interval (m/yr) |
| `flag` | 0 = OK; 1 = the baseline crosses a shoreline at that transect; 2 = the seaward side could not be determined |
| `tr_cross` | 1 = the transect crosses another transect (usually at sharp bends) |

Positive NSM, EPR and LRR = accretion; negative = erosion. Only transects
crossing at least 2 dates are kept. Transects are drawn from land to sea:
with an arrow symbol you can quickly check the orientation.

### Limitations

- At least 2 dates are needed. With only 2 there is no R², LSE or
  confidence interval: 3 or more are recommended.
- **The baseline must lie entirely on one side of all shorelines.** Where
  it does not, the transect gets `flag` 1 or 2. Filter those transects out
  before interpreting results.
- Positional uncertainty must be supplied by the user: the plugin does not
  estimate it. EPRunc and WLR are only as good as those values.
- At very sharp bends, transects may cross even with smoothing. They are
  marked in `tr_cross`: filter them, increase smoothing or shorten the
  transects.

### Status

Experimental. Not yet tested extensively on real case studies. If you find
a bug, please open an [issue](https://github.com/adanielibarra/FIC/issues).

---

## Español

Plugin de QGIS para el análisis del cambio de línea de costa, desarrollado
por **Daniel Ibarra Marinas**, Facultad de Ingeniería y Ciencias,
Universidad Autónoma de Tamaulipas (México).

### Qué hace

Genera transectos perpendiculares a una línea base (baseline) digitalizada
por el usuario, cada X metros. Cada transecto se corta con una capa
multitemporal de líneas de costa (shorelines) con un campo de fecha y se
calculan estas métricas, al estilo DSAS:

- **NSM (Net Shoreline Movement, m):** distancia entre la línea de costa
  más antigua y la más reciente que cruzan el transecto.
- **EPR (End Point Rate, m/año):** el NSM dividido entre los años
  transcurridos entre la primera y la última línea de costa.
- **LRR (Linear Regression Rate, m/año):** pendiente de la regresión
  lineal de la posición de todas las líneas de costa frente al tiempo. Con
  3 o más fechas, también su R², el error estándar de la estimación (LSE)
  y el intervalo de confianza de la pendiente (LCI) al nivel que elijas
  (90, 95 o 99 %).
- **WLR (Weighted Linear Regression, m/año)**, opcional: como el LRR, pero
  las fechas con menos error de posición pesan más. Necesita la
  incertidumbre de posición.

Además:

- **Agrupa por fecha:** si la línea de costa de una fecha está
  digitalizada en varios trozos, se juntan antes de cortar, de modo que
  cada fecha aporta una sola posición por transecto.
- **Cortes múltiples:** si un transecto corta la línea de costa de una
  fecha más de una vez, se usa el corte más cercano a la línea base (como
  la opción "closest" de DSAS).
- **Reproyección:** si las líneas de costa están en otro SRC, se
  reproyectan automáticamente al SRC de la línea base.
- **Signo al estilo DSAS:** indicas si la línea base está en tierra o en
  el mar, y cada transecto se orienta hacia el mar. Las tasas positivas son
  acreción y las negativas, erosión, sin importar el sentido en que
  digitalizaste la línea base.

### Instalación

1. Instálalo desde el gestor de complementos de QGIS (Complementos →
   Administrar e instalar complementos), buscando "FIC Shoreline Change
   Analysis". Al ser experimental, activa antes "Mostrar también los
   complementos experimentales" en la pestaña de configuración.
2. O descarga el ZIP de la última versión e instálalo con "Instalar
   desde ZIP" en el mismo gestor.
3. Aparece en la Caja de herramientas de Processing, en
   "FIC Shoreline Change Analysis" → "Coastal".

### Uso

Necesitas dos capas de líneas:

- Una **línea base**: una única línea continua, digitalizada a mano, en un
  **SRC proyectado en metros** (por ejemplo, la zona UTM que corresponda).
  Si está en grados u otra unidad, el algoritmo se detiene y lo avisa.
- Una capa de **líneas de costa multitemporales**: líneas de distintas
  fechas, con un campo de tipo fecha (o texto ISO `AAAA-MM-DD`). Las
  features con fecha no válida se descartan y se avisa en el registro.

El algoritmo pide también:

- El espaciado entre transectos y la longitud total de cada transecto,
  ambos en metros. Los transectos se generan simétricos, mitad a cada
  lado de la línea base.
- La **posición de la línea base**: en tierra o en el mar. La línea base
  debe quedar entera a un lado de todas las líneas de costa.
- La **distancia de suavizado** (m): la dirección de cada transecto se
  mide sobre un tramo de línea base de esa longitud, así los transectos
  siguen la forma general de la costa y no los quiebros de la
  digitalización. Un punto de partida razonable es entre 2 y 5 veces el
  espaciado. Con 0 no se suaviza.
- El **nivel de confianza** de los intervalos: 90 % (por defecto),
  95 % o 99 %.
- La **incertidumbre de posición** de las líneas de costa (opcional), en
  metros. Se toma, por orden de prioridad:
  1. Un campo con el **error total** ya calculado.
  2. Campos con sus **componentes**, que se combinan en cuadratura
     (raíz de la suma de cuadrados): **resolución** de la imagen de
     origen (cuenta como 1 píxel), **RMSE de georreferenciación** y
     **error de digitalización**. Los que falten cuentan como 0.
  3. Un **valor único** para todas las líneas de costa.

  Así cada fecha puede tener su propio error: una foto aérea antigua y
  una imagen de satélite reciente no tienen la misma precisión. Si una
  fecha está en varios trozos con errores distintos, se toma el mayor.
  La resolución sola subestima el error (el plugin lo avisa): conviene
  añadir al menos el RMSE. Sin incertidumbre, el plugin funciona igual
  pero no calcula EPRunc ni WLR.

#### Cómo se decide el lado del mar

En cada transecto se mira a qué lado de la línea base cae el corte más
cercano de cada fecha, y se toma el lado de la mayoría. Si la línea base
está en tierra, el mar está en ese lado; si está en el mar, en el
contrario. Después, para cada fecha se usa el corte más cercano dentro de
ese lado, de modo que una orilla que quede al otro lado de la línea base
(por ejemplo, la de una laguna) no se mezcla en el cálculo.

### Resultado

Una capa de líneas (los transectos) con estos campos:

| Campo | Contenido |
|---|---|
| `id` | Identificador del transecto |
| `n_pts` | Número de fechas que cortan el transecto |
| `yr_min`, `yr_max` | Año decimal de la primera y la última fecha |
| `dist_min`, `dist_max` | Distancia a la línea base (m) de la primera y la última fecha, positiva hacia el mar |
| `NSM` | Net Shoreline Movement (m) |
| `EPR` | End Point Rate (m/año) |
| `EPRunc` | Incertidumbre del EPR (m/año): raíz de la suma de los cuadrados de la incertidumbre de la primera y la última fecha, dividida entre los años. Vacío sin incertidumbre de posición |
| `LRR` | Linear Regression Rate (m/año) |
| `LRR_R2` | R² de la regresión (vacío con menos de 3 fechas) |
| `LSE` | Error estándar de la estimación (m): cuánto se separan, de media, las posiciones de la recta. Vacío con menos de 3 fechas |
| `LCI90` (o `LCI95`, `LCI99`) | Semiamplitud del intervalo de confianza del LRR (m/año). Si el valor absoluto del LRR es menor que el LCI, la tasa no es distinta de cero a ese nivel. Vacío con menos de 3 fechas |
| `WLR` | Weighted Linear Regression (m/año): como el LRR, pero las fechas con menos error pesan más (pesos 1/U²). Vacío con menos de 3 fechas o sin incertidumbre |
| `WR2`, `WSE` | R² y error estándar de la estimación de la WLR |
| `WCI90` (o `WCI95`, `WCI99`) | Semiamplitud del intervalo de confianza de la WLR (m/año) |
| `flag` | 0 = correcto; 1 = la línea base cruza alguna línea de costa en ese transecto; 2 = no se pudo saber hacia qué lado está el mar |
| `tr_cross` | 1 = el transecto se cruza con otro transecto (suele pasar en curvas cerradas) |

NSM, EPR y LRR positivos = acreción; negativos = erosión. Solo se guardan
los transectos que cortan al menos 2 fechas. Los transectos se dibujan de
tierra a mar: con un símbolo de flecha se ve enseguida si la orientación
es correcta.

### Limitaciones

- Hacen falta al menos 2 fechas. Con solo 2 no hay R², LSE ni intervalo
  de confianza: se recomiendan 3 o más.
- **La línea base tiene que quedar entera a un lado de todas las líneas
  de costa.** Donde no es así, el transecto sale con `flag` 1 o 2.
  Filtra esos transectos antes de interpretar resultados.
- La incertidumbre de posición la tiene que aportar el usuario: el
  plugin no la estima. Los resultados de EPRunc y WLR son tan buenos como
  esos valores.
- En curvas muy cerradas, los transectos pueden cruzarse aunque se
  suavice. Se marcan en `tr_cross`: fíltralos, sube el suavizado o acorta
  los transectos.

### Estado

Experimental. Aún no se ha probado a fondo en casos de uso reales. Si
encuentras un error, abre un [issue](https://github.com/adanielibarra/FIC/issues).

---

## Changelog / Cambios

- **0.1.2:**
  - EN: positional uncertainty per shoreline (total error, components in
    quadrature: resolution, RMSE and digitizing, or single value) with
    EPRunc and weighted regression (WLR, WR2, WSE, WCI); R², LSE and LRR
    confidence interval (90, 95 or 99 %), empty with fewer than 3 dates;
    baseline smoothing distance and `tr_cross` field; DSAS-style sign
    (onshore/offshore baseline, transects oriented seaward, positive =
    accretion) with a `flag` quality field; baseline CRS check (projected,
    metres) and automatic reprojection of shorelines; features with the
    same date are merged; inputs restricted to line layers; bilingual
    interface (English and Spanish).
  - ES: incertidumbre de posición por línea de costa (error total,
    componentes en cuadratura: resolución, RMSE y digitalización, o valor
    único) con EPRunc y regresión ponderada (WLR, WR2, WSE, WCI); R², LSE e
    intervalo de confianza del LRR (90, 95 o 99 %), vacíos con menos de 3
    fechas; distancia de suavizado y campo `tr_cross`; signo al estilo
    DSAS (línea base en tierra o en el mar, transectos orientados hacia el
    mar, positivo = acreción) con campo de calidad `flag`; comprobación del
    SRC de la línea base (proyectado, en metros) y reproyección automática
    de las líneas de costa; las features con la misma fecha se agrupan;
    las entradas solo admiten capas de líneas; interfaz bilingüe (inglés y
    español).
- **0.1.1:** first version / primera versión.

## License / Licencia

GPL v3. See / Ver [LICENSE](LICENSE).
