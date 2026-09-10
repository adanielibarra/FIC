"""
FIC Coastal Change Analysis
===========================
Flujo básico al estilo DSAS: transectos perpendiculares a una baseline
digitalizada por el usuario y estadísticas de cambio de línea de costa
por transecto (NSM, EPR, LRR).

Supuestos (cambian el resultado):

1. SRC: la baseline debe estar en un SRC proyectado con unidades en metros.
   Si las shorelines están en otro SRC, se reproyectan al de la baseline.

2. La baseline debe ser una única línea (una feature, una parte). Si es
   multiparte se intenta fusionar con mergeLines(); si no se puede, el
   algoritmo se detiene.

3. Signo al estilo DSAS: el usuario indica si la baseline está en tierra
   o en el mar. En cada transecto se mira a qué lado caen las shorelines
   (por mayoría de fechas) y el eje se orienta hacia el mar. Resultado:
   distancias positivas hacia el mar, tasas > 0 = acreción y < 0 =
   erosión, sea cual sea el sentido de digitalización. Requiere que la
   baseline quede entera a un lado de todas las shorelines; donde no se
   cumple, el transecto se marca en el campo flag.

4. Los transectos son simétricos (mitad a cada lado de la baseline).

5. Las features de shorelines con la misma fecha se agrupan en una sola
   geometría antes de cortar. Si un transecto corta varias veces la
   línea de una fecha, se usa el corte más cercano a la baseline
   (equivalente a la opción "closest" de DSAS).

6. La capa de shorelines necesita un campo de fecha (tipo Date o texto
   ISO YYYY-MM-DD). Las features con fecha no interpretable se descartan
   y se avisa en el registro.

7. Mínimo 2 fechas para NSM/EPR/LRR. R2, LSE y LCI solo se calculan con
   3 o más fechas (con 2 quedan vacíos).

8. LCI: semiamplitud del intervalo de confianza de la pendiente,
   t(1 - alfa/2, n - 2) * SE(LRR), con el cuantil de la t calculado sin
   dependencias externas. LSE: sqrt(SSres / (n - 2)).

9. Suavizado: la dirección de la baseline en cada transecto se mide entre
   dos puntos separados por la distancia de suavizado (mínimo +-1 m).
   Los transectos que se cruzan con otro se marcan en tr_cross.

10. Incertidumbre de posición: por feature. Prioridad: campo de error
    total > componentes en cuadratura, sqrt(res^2 + rmse^2 + digit^2), con
    la resolución contada como 1 píxel y los que falten a 0 > valor por
    defecto. Por fecha se toma la mayor de sus trozos; si algún
    trozo no tiene, la fecha queda sin incertidumbre. EPRunc =
    sqrt(U1^2 + U2^2) / años. WLR: mínimos cuadrados ponderados con
    w = 1 / U^2; WSE = sqrt(sum(w * res^2) / (n - 2)).

11. Tendencia (campo trend): erosion / accretion / stable / unclassified.
    Criterio estadístico: stable si |tasa| <= intervalo, con la primera
    disponible de WLR/WCI, LRR/LCI, EPR/EPRunc. Criterio umbral: stable si
    |tasa| <= umbral, con WLR, LRR (n >= 3) o EPR. flag != 0 siempre
    unclassified. trend_src guarda la tasa usada.
"""

import math
import os

from qgis.PyQt.QtCore import QUrl, QVariant
from qgis.PyQt.QtGui import QIcon
from qgis.core import (
    Qgis,
    QgsCoordinateTransform,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsUnitTypes,
    QgsProcessingParameterVectorLayer,
    QgsProcessingParameterField,
    QgsProcessingParameterNumber,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSink,
    QgsProcessingException,
    QgsFeature,
    QgsFeatureSink,
    QgsSpatialIndex,
    QgsFields,
    QgsField,
    QgsGeometry,
    QgsPointXY,
    QgsWkbTypes,
)

# Unidad "metros": en QGIS >= 3.30 esta en Qgis.DistanceUnit; en versiones
# anteriores en QgsUnitTypes.
try:
    METERS = Qgis.DistanceUnit.Meters
except AttributeError:
    METERS = QgsUnitTypes.DistanceMeters


# ---------------------------------------------------------------------------
# Estadística sin dependencias externas (scipy no viene en todas las
# instalaciones de QGIS). Cuantil de la t de Student a partir de la función
# beta incompleta regularizada (fracción continua de Lentz).
# ---------------------------------------------------------------------------

def _betacf(a, b, x, max_iter=300, eps=3e-14):
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _betainc_reg(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbt = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
           + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return math.exp(lbt) * _betacf(a, b, x) / a
    return 1.0 - math.exp(lbt) * _betacf(b, a, 1.0 - x) / b


def t_cdf(t, df):
    """Función de distribución de la t de Student."""
    x = df / (df + t * t)
    tail = 0.5 * _betainc_reg(df / 2.0, 0.5, x)
    return 1.0 - tail if t >= 0 else tail


def t_ppf(p, df):
    """Cuantil de la t de Student (p > 0.5), por bisección."""
    lo, hi = 0.0, 1.0
    while t_cdf(hi, df) < p:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if t_cdf(mid, df) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def linreg_stats(xs, ys, conf):
    """Regresión lineal de la posición frente al tiempo, al estilo DSAS.

    Devuelve (LRR, R2, LSE, LCI):
      LRR: pendiente (m/año).
      R2:  coeficiente de determinación. None con menos de 3 fechas.
      LSE: error estándar de la estimación (m), sqrt(SSres / (n - 2)).
           None con menos de 3 fechas.
      LCI: semiamplitud del intervalo de confianza de la pendiente (m/año)
           al nivel conf (0.90, 0.95...): t(1 - alfa/2, n - 2) * SE(LRR).
           None con menos de 3 fechas.
    """
    n = len(xs)
    if n < 2:
        return None, None, None, None
    mx, my = sum(xs) / n, sum(ys) / n
    ss_xx = sum((x - mx) ** 2 for x in xs)
    if ss_xx == 0:
        return None, None, None, None
    ss_xy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    lrr = ss_xy / ss_xx
    if n < 3:
        return lrr, None, None, None
    ss_res = sum((y - (my + lrr * (x - mx))) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - my) ** 2 for y in ys)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else None
    lse = math.sqrt(ss_res / (n - 2))
    se_slope = lse / math.sqrt(ss_xx)
    lci = t_ppf(1.0 - (1.0 - conf) / 2.0, n - 2) * se_slope
    return lrr, r2, lse, lci


def _positive(v):
    """Devuelve v como float si es un número > 0; si no, None."""
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if v > 0 and not math.isnan(v) else None


def position_uncertainty(total, res, rmse, digit, default):
    """Incertidumbre de posición de una shoreline (m).

    Orden de prioridad:
      1. total: error total ya calculado por el usuario (> 0).
      2. Componentes combinados en cuadratura:
         sqrt(res^2 + rmse^2 + digit^2), donde res es la resolución de la
         imagen (1 píxel), rmse el error de georreferenciación y digit el
         error de digitalización. Los que falten cuentan como 0.
      3. default: valor único por defecto (> 0).
    Devuelve None si no hay nada.
    """
    t = _positive(total)
    if t is not None:
        return t
    comps = [c for c in (_positive(res), _positive(rmse), _positive(digit))
             if c is not None]
    if comps:
        return math.sqrt(sum(c * c for c in comps))
    return _positive(default)


def epr_uncertainty(u_first, u_last, dyears):
    """Incertidumbre del EPR (m/año) al estilo DSAS:
    sqrt(U1^2 + U2^2) / años entre la primera y la última fecha."""
    if u_first is None or u_last is None or dyears <= 0:
        return None
    return math.sqrt(u_first ** 2 + u_last ** 2) / dyears


def wlr_stats(xs, ys, us, conf):
    """Regresión lineal ponderada (WLR), pesos w = 1 / U^2.

    Devuelve (WLR, WR2, WSE, WCI):
      WLR: pendiente ponderada (m/año).
      WR2: R2 ponderado.
      WSE: error estándar ponderado de la estimación,
           sqrt(sum(w * res^2) / (n - 2)).
      WCI: semiamplitud del intervalo de confianza de la pendiente,
           t(1 - alfa/2, n - 2) * WSE / sqrt(sum(w * (x - xw)^2)).
    Todo None con menos de 3 fechas o si falta alguna incertidumbre.
    """
    n = len(xs)
    if n < 3 or any(u is None or u <= 0 for u in us):
        return None, None, None, None
    ws = [1.0 / (u * u) for u in us]
    sw = sum(ws)
    xw = sum(w * x for w, x in zip(ws, xs)) / sw
    yw = sum(w * y for w, y in zip(ws, ys)) / sw
    sxx = sum(w * (x - xw) ** 2 for w, x in zip(ws, xs))
    if sxx == 0:
        return None, None, None, None
    sxy = sum(w * (x - xw) * (y - yw) for w, x, y in zip(ws, xs, ys))
    b = sxy / sxx
    ss_res = sum(w * (y - (yw + b * (x - xw))) ** 2
                 for w, x, y in zip(ws, xs, ys))
    ss_tot = sum(w * (y - yw) ** 2 for w, y in zip(ws, ys))
    wr2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else None
    wse = math.sqrt(ss_res / (n - 2))
    wci = t_ppf(1.0 - (1.0 - conf) / 2.0, n - 2) * wse / math.sqrt(sxx)
    return b, wr2, wse, wci


TREND_STATISTICAL = 0
TREND_THRESHOLD = 1


def classify_trend(flag, method, threshold, n,
                   epr, eprunc, lrr, lci, wlr, wci):
    """Clasifica un transecto en erosion / accretion / stable / unclassified.

    Devuelve (clase, tasa_usada) con tasa_usada = 'WLR', 'LRR', 'EPR' o
    None.

    - Transectos con flag distinto de 0: 'unclassified' (signo no fiable).
    - Estadístico: se usa la mejor tasa con intervalo disponible, por este
      orden: WLR con WCI, LRR con LCI, EPR con EPRunc. 'stable' si el valor
      absoluto de la tasa es menor o igual que su intervalo (no distinta de
      cero). Si no hay ninguna con intervalo: 'unclassified'.
    - Umbral: se usa WLR si existe; si no, LRR con 3 o más fechas; si no,
      EPR. 'stable' si el valor absoluto de la tasa es menor o igual que el
      umbral (m/año).
    - En los demás casos: 'accretion' si la tasa es positiva, 'erosion' si
      es negativa.
    """
    if flag != 0:
        return 'unclassified', None

    if method == TREND_STATISTICAL:
        if wlr is not None and wci is not None:
            rate, band, src = wlr, wci, 'WLR'
        elif lrr is not None and lci is not None:
            rate, band, src = lrr, lci, 'LRR'
        elif epr is not None and eprunc is not None:
            rate, band, src = epr, eprunc, 'EPR'
        else:
            return 'unclassified', None
    else:
        if wlr is not None:
            rate, src = wlr, 'WLR'
        elif lrr is not None and n >= 3:
            rate, src = lrr, 'LRR'
        elif epr is not None:
            rate, src = epr, 'EPR'
        else:
            return 'unclassified', None
        band = threshold

    if abs(rate) <= band:
        return 'stable', src
    return ('accretion' if rate > 0 else 'erosion'), src


def lbl(en, es):
    """Etiqueta bilingüe en una línea: 'English / Español'."""
    return f'{en} / {es}'


def msg(en, es):
    """Mensaje bilingüe en dos líneas: primero inglés, luego español."""
    return f'{en}\n{es}'


FLAG_OK = 0         # todas las fechas al mismo lado de la baseline
FLAG_CROSSES = 1    # alguna fecha solo aparece al otro lado: la baseline
                    # cruza las shorelines en este transecto
FLAG_AMBIGUOUS = 2  # empate: tantas fechas a un lado como al otro, la
                    # orientación de este transecto no es fiable


def orient_seaward(per_date, onshore):
    """Orienta un transecto hacia el mar al estilo DSAS.

    per_date: lista de (año, [distancias con signo de todos los cortes de
              esa fecha]), medidas sobre la perpendicular sin orientar.
    onshore:  True si la baseline está en tierra, False si está en el mar.

    Devuelve (pts, f, flag):
      pts:  lista de (año, distancia) con distancia positiva hacia el mar.
      f:    +1 o -1, factor que convierte la perpendicular original en el
            eje hacia el mar.
      flag: FLAG_OK, FLAG_CROSSES o FLAG_AMBIGUOUS.

    Lado de las shorelines: se decide por mayoría, con el corte más
    cercano a la baseline de cada fecha. Con la baseline en tierra, el mar
    está en ese lado; con la baseline en el mar, en el lado contrario.
    Para cada fecha se usa el corte más cercano DENTRO del lado
    mayoritario, así un corte con otra orilla al otro lado de la baseline
    (por ejemplo, una laguna) no se cuela.
    """
    votes = [min(sd, key=abs) for _, sd in per_date]
    n_pos = sum(1 for v in votes if v > 0)
    n_neg = sum(1 for v in votes if v < 0)
    side = 1 if n_pos >= n_neg else -1

    pts = []
    crosses = False
    for year, sd in per_date:
        same_side = [v for v in sd if v * side > 0]
        if same_side:
            best = min(same_side, key=abs)
        else:
            best = min(sd, key=abs)
            if best != 0:
                crosses = True
        pts.append((year, best))

    f = side if onshore else -side
    pts = [(yr, f * d) for yr, d in pts]

    if n_pos == n_neg:
        flag = FLAG_AMBIGUOUS
    elif crosses:
        flag = FLAG_CROSSES
    else:
        flag = FLAG_OK
    return pts, f, flag


class ShorelineChangeAlgorithm(QgsProcessingAlgorithm):

    BASELINE = 'BASELINE'
    SHORELINES = 'SHORELINES'
    DATE_FIELD = 'DATE_FIELD'
    SPACING = 'SPACING'
    LENGTH = 'LENGTH'
    BASELINE_POS = 'BASELINE_POS'
    SMOOTH = 'SMOOTH'
    CONF = 'CONF'
    CONF_LEVELS = [90, 95, 99]
    UNC_FIELD = 'UNC_FIELD'
    UNC_DEFAULT = 'UNC_DEFAULT'
    TREND_METHOD = 'TREND_METHOD'
    TREND_THRESH = 'TREND_THRESH'
    RES_FIELD = 'RES_FIELD'
    GEOREF_FIELD = 'GEOREF_FIELD'
    DIGIT_FIELD = 'DIGIT_FIELD'
    OUTPUT = 'OUTPUT'

    def createInstance(self):
        return ShorelineChangeAlgorithm()

    def name(self):
        return 'shoreline_change_stats'

    def displayName(self):
        return 'FIC Coastal Change Analysis'

    def group(self):
        return 'Coastal'

    def groupId(self):
        return 'coastal'

    def icon(self):
        import os
        icon_path = os.path.join(os.path.dirname(__file__), 'icon.png')
        return QIcon(icon_path)

    def shortHelpString(self):
        logo = QUrl.fromLocalFile(
            os.path.join(os.path.dirname(__file__), 'icon.png')).toString()
        credit = (
            '<table><tr>'
            f'<td valign="middle"><img src="{logo}" width="64" height="64"></td>'
            '<td valign="middle" style="padding-left:10px">'
            'Created by / Creado por Daniel Ibarra Marinas<br>'
            '<span style="color:#2e8a3d"><b>Facultad de Ingeniería '
            'y Ciencias</b></span><br>'
            '<span style="color:#d97706"><b>Universidad Autónoma de '
            'Tamaulipas</b></span>'
            '</td></tr></table>'
        )
        en = (
            '<p>Casts <b>transects</b> perpendicular to a baseline (a single '
            'line, in a projected CRS in metres) every X metres. Each '
            'transect is intersected with the multitemporal shorelines '
            '(features with the same date are merged) and the following are '
            'computed:</p>'
            '<p><b>NSM (Net Shoreline Movement, m):</b> distance between the '
            'oldest and the most recent shoreline. Total change over the '
            'period, ignoring intermediate dates.</p>'
            '<p><b>EPR (End Point Rate, m/yr):</b> NSM divided by the years '
            'elapsed. It only uses the two extreme dates: if one of them was '
            'an unusual year, the EPR is biased.</p>'
            '<p><b>LRR (Linear Regression Rate, m/yr):</b> slope of the '
            'linear regression of all shoreline positions against time. With '
            '3 or more dates it adds <b>LRR_R2</b> (fit to a straight line), '
            '<b>LSE</b> (standard error of the estimate, m) and <b>LCI</b> '
            '(half-width of the LRR confidence interval, m/yr, at the chosen '
            'level). If the absolute LRR is smaller than LCI, the rate is not '
            'different from zero at that level. With only 2 dates these three '
            'fields are empty.</p>'
            '<p><b>Positional uncertainty (optional):</b> error in metres of '
            'each shoreline. In order of priority: (1) a field with the total '
            'error; (2) fields with its components, combined in quadrature: '
            'source image resolution (counted as 1 pixel), georeferencing '
            'RMSE and digitizing error (missing ones count as 0); (3) a '
            'single value for all. Resolution alone underestimates the error: '
            'add at least the RMSE. If a date is split into several pieces, '
            'the largest error is used. With it the plugin computes '
            '<b>EPRunc</b> (EPR uncertainty, m/yr) and the weighted '
            'regression <b>WLR</b> (m/yr, weights 1/U²), with <b>WR2</b>, '
            '<b>WSE</b> and <b>WCI</b>. Without it, these fields are '
            'empty.</p>'
            '<p><b>Sign:</b> state whether the baseline is onshore or '
            'offshore. Positive rates mean accretion and negative rates mean '
            'erosion, whatever the digitizing direction of the baseline. '
            'Transects are drawn from land to sea.</p>'
            '<p><b>Smoothing:</b> the direction of each transect is measured '
            'over a stretch of baseline of this length. Increase it if the '
            'baseline has kinks and transects come out skewed or crossing. '
            'A reasonable starting point is 2 to 5 times the spacing; '
            '0 = no smoothing.</p>'
            '<p><b>Trend:</b> the <b>trend</b> field classifies each '
            'transect as erosion, accretion, stable or unclassified. '
            '<i>Statistical</i> criterion (default): stable if the rate is '
            'not different from zero, using the best rate with an interval: '
            'WLR with WCI, else LRR with LCI, else EPR with EPRunc; with 2 '
            'dates and no positional uncertainty there is no interval and '
            'the transect is unclassified. <i>Threshold</i> criterion: stable '
            'if the absolute rate is at most the threshold (m/yr), using WLR, '
            'else LRR (3 or more dates), else EPR. The 0.5 m/yr default is '
            'only an example: choose a value that makes sense for your '
            'coast. <b>trend_src</b> says which rate was used. Transects with '
            'flag 1 or 2 are always unclassified.</p>'
            '<p style="color:#6b6b6b"><b>IMPORTANT:</b> check two fields '
            'before interpreting. <b>flag</b>: 0 = OK; 1 = the baseline '
            'crosses a shoreline at that transect (unreliable rates); 2 = the '
            'seaward side could not be determined (do not use its rates). '
            '<b>tr_cross</b>: 1 = the transect crosses another one.</p>'
        )
        es = (
            '<p>Genera <b>transectos</b> perpendiculares a una línea base '
            '(una sola línea, en un SRC proyectado en metros) cada X metros. '
            'Cada transecto se corta con las líneas de costa multitemporales '
            '(las de la misma fecha se agrupan) y se calcula:</p>'
            '<p><b>NSM (Net Shoreline Movement, m):</b> distancia entre la '
            'línea de costa más antigua y la más reciente. Es el cambio total '
            'en el periodo, sin contar las fechas intermedias.</p>'
            '<p><b>EPR (End Point Rate, m/año):</b> el NSM dividido entre '
            'los años transcurridos. Solo usa las dos fechas extremas: si '
            'una de ellas fue un año anómalo, el EPR sale sesgado.</p>'
            '<p><b>LRR (Linear Regression Rate, m/año):</b> pendiente de la '
            'regresión lineal de la posición de todas las líneas de costa '
            'frente al tiempo. Con 3 o más fechas se añaden: <b>LRR_R2</b> '
            '(ajuste a una recta), <b>LSE</b> (error estándar de la '
            'estimación, m) y <b>LCI</b> (semiamplitud del intervalo de '
            'confianza del LRR, m/año, al nivel elegido). Si el valor '
            'absoluto del LRR es menor que el LCI, la tasa no es distinta de '
            'cero a ese nivel. Con solo 2 fechas estos tres campos quedan '
            'vacíos.</p>'
            '<p><b>Incertidumbre de posición (opcional):</b> error en metros '
            'de cada línea de costa. Por orden de prioridad: (1) un campo con '
            'el error total; (2) campos con sus componentes, que se combinan '
            'en cuadratura: resolución de la imagen de origen (cuenta como '
            '1 píxel), RMSE de georreferenciación y error de digitalización '
            '(los que falten cuentan como 0); (3) un valor único para todas. '
            'La resolución sola subestima el error: conviene añadir al menos '
            'el RMSE. Si una fecha está en varios trozos, se toma el error '
            'mayor. Con ella se calculan <b>EPRunc</b> (incertidumbre del '
            'EPR, m/año) y la regresión ponderada <b>WLR</b> (m/año, pesos '
            '1/U²), con <b>WR2</b>, <b>WSE</b> y <b>WCI</b>. Sin ella, estos '
            'campos quedan vacíos.</p>'
            '<p><b>Signo:</b> indica si la línea base está en tierra o en el '
            'mar. Las tasas positivas son acreción y las negativas, erosión, '
            'sin importar el sentido en que se digitalizó la línea base. Los '
            'transectos se dibujan de tierra a mar.</p>'
            '<p><b>Suavizado:</b> la dirección de cada transecto se mide '
            'sobre un tramo de línea base de esa longitud. Súbelo si la línea '
            'base tiene quiebros y los transectos salen torcidos o se cruzan. '
            'Un punto de partida razonable es entre 2 y 5 veces el '
            'espaciado; 0 = sin suavizar.</p>'
            '<p><b>Tendencia:</b> el campo <b>trend</b> clasifica cada '
            'transecto como erosion (erosión), accretion (acreción), stable '
            '(estable) o unclassified (sin clasificar). Criterio '
            '<i>estadístico</i> (por defecto): estable si la tasa no es '
            'distinta de cero, con la mejor tasa que tenga intervalo: WLR con '
            'WCI; si no, LRR con LCI; si no, EPR con EPRunc. Con 2 fechas y '
            'sin incertidumbre de posición no hay intervalo y el transecto '
            'queda sin clasificar. Criterio <i>umbral</i>: estable si el '
            'valor absoluto de la tasa no pasa del umbral (m/año), usando '
            'WLR; si no, LRR (3 o más fechas); si no, EPR. El valor por '
            'defecto de 0,5 m/año es solo un ejemplo: elige uno que tenga '
            'sentido en tu costa. <b>trend_src</b> dice qué tasa se usó. Los '
            'transectos con flag 1 o 2 quedan siempre sin clasificar.</p>'
            '<p style="color:#6b6b6b"><b>IMPORTANTE:</b> revisa dos campos '
            'antes de interpretar. <b>flag</b>: 0 = correcto; 1 = la línea '
            'base cruza alguna línea de costa en ese transecto (tasas poco '
            'fiables); 2 = no se pudo saber hacia qué lado está el mar (no '
            'uses sus tasas). <b>tr_cross</b>: 1 = el transecto se cruza con '
            'otro.</p>'
        )
        return credit + en + '<hr>' + es


    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterVectorLayer(
            self.BASELINE, lbl('Baseline (single line)', 'Línea base (una sola línea)'),
            types=[QgsProcessing.TypeVectorLine]))
        self.addParameter(QgsProcessingParameterVectorLayer(
            self.SHORELINES, lbl('Multitemporal shorelines (lines)',
                                'Líneas de costa multitemporales (líneas)'),
            types=[QgsProcessing.TypeVectorLine]))
        self.addParameter(QgsProcessingParameterField(
            self.DATE_FIELD, lbl('Date field', 'Campo de fecha'),
            parentLayerParameterName=self.SHORELINES))
        self.addParameter(QgsProcessingParameterField(
            self.UNC_FIELD,
            lbl('Total positional uncertainty field (m)',
                'Campo de incertidumbre total de posición (m)'),
            parentLayerParameterName=self.SHORELINES,
            type=QgsProcessingParameterField.Numeric,
            optional=True))
        self.addParameter(QgsProcessingParameterField(
            self.RES_FIELD,
            lbl('Source image resolution field (m/pixel)',
                'Campo de resolución de la imagen de origen (m/píxel)'),
            parentLayerParameterName=self.SHORELINES,
            type=QgsProcessingParameterField.Numeric,
            optional=True))
        self.addParameter(QgsProcessingParameterField(
            self.GEOREF_FIELD,
            lbl('Georeferencing RMSE field (m)',
                'Campo de RMSE de georreferenciación (m)'),
            parentLayerParameterName=self.SHORELINES,
            type=QgsProcessingParameterField.Numeric,
            optional=True))
        self.addParameter(QgsProcessingParameterField(
            self.DIGIT_FIELD,
            lbl('Digitizing error field (m)',
                'Campo de error de digitalización (m)'),
            parentLayerParameterName=self.SHORELINES,
            type=QgsProcessingParameterField.Numeric,
            optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            self.UNC_DEFAULT,
            lbl('Default positional uncertainty (m, 0 = no EPRunc or WLR)',
                'Incertidumbre de posición por defecto (m, 0 = sin EPRunc ni WLR)'),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(
            self.SPACING, lbl('Transect spacing (m)', 'Espaciado entre transectos (m)'),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=50.0, minValue=0.1))
        self.addParameter(QgsProcessingParameterNumber(
            self.LENGTH, lbl('Total transect length (m)', 'Longitud total del transecto (m)'),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=500.0, minValue=1.0))
        self.addParameter(QgsProcessingParameterEnum(
            self.BASELINE_POS, lbl('Baseline position', 'Posición de la línea base'),
            options=[lbl('Onshore', 'En tierra'), lbl('Offshore', 'En el mar')],
            defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            self.SMOOTH,
            lbl('Baseline smoothing distance (m, 0 = no smoothing)',
                'Distancia de suavizado de la línea base (m, 0 = sin suavizar)'),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterEnum(
            self.CONF, lbl('Confidence level', 'Nivel de confianza'),
            options=[f'{c} %' for c in self.CONF_LEVELS],
            defaultValue=0))
        self.addParameter(QgsProcessingParameterEnum(
            self.TREND_METHOD,
            lbl('Stability criterion (trend field)',
                'Criterio de estabilidad (campo trend)'),
            options=[lbl('Statistical: rate not different from zero',
                         'Estadístico: tasa no distinta de cero'),
                     lbl('Threshold', 'Umbral')],
            defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            self.TREND_THRESH,
            lbl('Stability threshold (m/yr, only with Threshold criterion)',
                'Umbral de estabilidad (m/año, solo con criterio Umbral)'),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=0.5, minValue=0.0))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT, lbl('Transects with shoreline change rates',
                             'Transectos con tasas de cambio')))

    def _output_fields(self, conf_pct):
        fields = QgsFields()
        fields.append(QgsField('id', QVariant.Int))
        fields.append(QgsField('n_pts', QVariant.Int))
        fields.append(QgsField('yr_min', QVariant.Double))
        fields.append(QgsField('yr_max', QVariant.Double))
        fields.append(QgsField('dist_min', QVariant.Double))
        fields.append(QgsField('dist_max', QVariant.Double))
        fields.append(QgsField('NSM', QVariant.Double))
        fields.append(QgsField('EPR', QVariant.Double))
        fields.append(QgsField('EPRunc', QVariant.Double))
        fields.append(QgsField('LRR', QVariant.Double))
        fields.append(QgsField('LRR_R2', QVariant.Double))
        fields.append(QgsField('LSE', QVariant.Double))
        fields.append(QgsField(f'LCI{conf_pct}', QVariant.Double))
        fields.append(QgsField('WLR', QVariant.Double))
        fields.append(QgsField('WR2', QVariant.Double))
        fields.append(QgsField('WSE', QVariant.Double))
        fields.append(QgsField(f'WCI{conf_pct}', QVariant.Double))
        fields.append(QgsField('flag', QVariant.Int))
        fields.append(QgsField('tr_cross', QVariant.Int))
        fields.append(QgsField('trend', QVariant.String, len=12))
        fields.append(QgsField('trend_src', QVariant.String, len=3))
        return fields

    def _decimal_year(self, value, feedback):
        """Convierte QDate o texto ISO a año decimal. Devuelve None si no
        se puede interpretar (y avisa por el log)."""
        if value is None:
            return None
        try:
            if hasattr(value, 'year') and hasattr(value, 'dayOfYear'):
                return value.year() + (value.dayOfYear() - 1) / 365.25
            from datetime import datetime
            dt = datetime.fromisoformat(str(value)[:10])
            return dt.year + (dt.timetuple().tm_yday - 1) / 365.25
        except Exception:
            feedback.pushWarning(
                msg(f'Could not read {value!r} as a date: shoreline discarded.',
                    f'No se pudo leer {value!r} como fecha: se descarta esa '
                    'línea de costa.')
            )
            return None

    def processAlgorithm(self, parameters, context, feedback):
        baseline_layer = self.parameterAsVectorLayer(parameters, self.BASELINE, context)
        shorelines_layer = self.parameterAsVectorLayer(parameters, self.SHORELINES, context)
        date_field = self.parameterAsString(parameters, self.DATE_FIELD, context)
        spacing = self.parameterAsDouble(parameters, self.SPACING, context)
        length = self.parameterAsDouble(parameters, self.LENGTH, context)
        onshore = self.parameterAsEnum(parameters, self.BASELINE_POS, context) == 0
        smooth = self.parameterAsDouble(parameters, self.SMOOTH, context)
        conf_pct = self.CONF_LEVELS[
            self.parameterAsEnum(parameters, self.CONF, context)]
        conf = conf_pct / 100.0
        unc_field = self.parameterAsString(parameters, self.UNC_FIELD, context)
        unc_default = self.parameterAsDouble(parameters, self.UNC_DEFAULT, context)
        trend_method = self.parameterAsEnum(parameters, self.TREND_METHOD, context)
        trend_thresh = self.parameterAsDouble(parameters, self.TREND_THRESH, context)
        res_field = self.parameterAsString(parameters, self.RES_FIELD, context)
        georef_field = self.parameterAsString(parameters, self.GEOREF_FIELD, context)
        digit_field = self.parameterAsString(parameters, self.DIGIT_FIELD, context)
        unc_fields_used = [x for x in (unc_field, res_field, georef_field,
                                       digit_field) if x]

        def attr(feat, name):
            return feat.attribute(name) if name else None

        if baseline_layer is None or shorelines_layer is None:
            raise QgsProcessingException(msg(
                'The baseline or the shoreline layer is missing.',
                'Falta la línea base o la capa de líneas de costa.'))

        # --- SRC: proyectado y en metros; shorelines al SRC de la baseline ---
        base_crs = baseline_layer.crs()
        if not base_crs.isValid():
            raise QgsProcessingException(msg(
                'The baseline has no valid CRS.',
                'La línea base no tiene un SRC válido.'))
        if base_crs.isGeographic():
            raise QgsProcessingException(msg(
                f'The baseline is in a geographic CRS ({base_crs.authid()}), '
                'in degrees. Reproject it to a projected CRS in metres (for '
                'example, the corresponding UTM zone) and run again.',
                f'La línea base está en un SRC geográfico ({base_crs.authid()}), '
                'en grados. Reproyéctala a un SRC proyectado en metros (por '
                'ejemplo, la zona UTM que corresponda) y vuelve a lanzarlo.'))
        if base_crs.mapUnits() != METERS:
            raise QgsProcessingException(msg(
                f'The baseline CRS ({base_crs.authid()}) is not in metres. '
                'Spacing, length and rates are expressed in metres: reproject '
                'to a CRS in metres.',
                f'El SRC de la línea base ({base_crs.authid()}) no está en '
                'metros. El espaciado, la longitud y las tasas se expresan en '
                'metros: reproyecta a un SRC en metros.'))

        xform = None
        sl_crs = shorelines_layer.crs()
        if sl_crs.isValid() and sl_crs != base_crs:
            xform = QgsCoordinateTransform(sl_crs, base_crs,
                                           context.transformContext())
            feedback.pushInfo(msg(
                f'Shorelines ({sl_crs.authid()}) are reprojected to the '
                f'baseline CRS ({base_crs.authid()}).',
                f'Las líneas de costa ({sl_crs.authid()}) se reproyectan al '
                f'SRC de la línea base ({base_crs.authid()}).'))
        elif not sl_crs.isValid():
            feedback.pushWarning(msg(
                'The shoreline layer has no valid CRS: the baseline CRS is '
                'assumed.',
                'La capa de líneas de costa no tiene SRC válido: se asume el '
                'de la línea base.'))

        # --- baseline: una sola geometria continua ---
        feats = list(baseline_layer.getFeatures())
        if not feats:
            raise QgsProcessingException(msg(
                'The baseline layer has no features.',
                'La capa de línea base no tiene features.'))
        if len(feats) > 1:
            feedback.pushWarning(msg(
                'The baseline has more than one feature; only the first one '
                'is used. Merge the lines first if you want to use them all.',
                'La línea base tiene más de una feature; solo se usa la '
                'primera. Fusiona las líneas antes si quieres usarlas todas.'))
        baseline_geom = feats[0].geometry()
        if baseline_geom.isMultipart():
            merged = baseline_geom.mergeLines()
            if merged and not merged.isEmpty():
                baseline_geom = merged
            else:
                raise QgsProcessingException(msg(
                    'The baseline is multipart and could not be merged into '
                    'a single continuous line.',
                    'La línea base es multiparte y no se pudo fusionar en una '
                    'sola línea continua.'))

        total_length = baseline_geom.length()
        if total_length <= 0:
            raise QgsProcessingException(msg(
                'The baseline has zero length.',
                'La línea base tiene longitud 0.'))

        # --- shorelines con fecha válida, agrupadas por fecha ---
        # Si una fecha está digitalizada en varios trozos, se juntan en una
        # sola geometría para que cada fecha aporte un único punto por
        # transecto.
        by_date = {}
        unc_by_date = {}  # incertidumbre por fecha (la mayor de sus trozos)
        n_feats = 0
        n_unc_missing = 0
        for f in shorelines_layer.getFeatures():
            year = self._decimal_year(f.attribute(date_field), feedback)
            if year is None:
                continue
            geom = f.geometry()
            if geom is None or geom.isEmpty():
                continue
            geom = QgsGeometry(geom)
            if xform is not None:
                geom.transform(xform)
            key = round(year, 6)
            by_date.setdefault(key, []).append(geom)
            n_feats += 1

            # Incertidumbre: total > componentes en cuadratura > por defecto.
            u_fields = position_uncertainty(
                attr(f, unc_field), attr(f, res_field),
                attr(f, georef_field), attr(f, digit_field), None)
            if unc_fields_used and u_fields is None:
                n_unc_missing += 1
            u = u_fields if u_fields is not None else _positive(unc_default)
            prev = unc_by_date.get(key, 'nada')
            if prev == 'nada':
                unc_by_date[key] = u
            elif prev is None or u is None:
                unc_by_date[key] = None
            else:
                unc_by_date[key] = max(prev, u)

        if res_field and not (unc_field or georef_field or digit_field):
            feedback.pushWarning(msg(
                'Only image resolution was given. It is just one part of the '
                'positional error (georeferencing and digitizing are '
                'missing): EPRunc and WLR will underestimate the error.',
                'Solo se ha indicado la resolución de la imagen. Es solo una '
                'parte del error de posición (faltan georreferenciación y '
                'digitalización): EPRunc y WLR saldrán con el error '
                'subestimado.'))
        if n_unc_missing:
            if unc_default > 0:
                feedback.pushWarning(msg(
                    f'{n_unc_missing} features without a valid uncertainty in '
                    f'the given fields: the default value ({unc_default} m) '
                    'is used.',
                    f'{n_unc_missing} features sin incertidumbre válida en los '
                    f'campos indicados: se usa el valor por defecto '
                    f'({unc_default} m).'))
            else:
                feedback.pushWarning(msg(
                    f'{n_unc_missing} features without a valid uncertainty in '
                    'the given fields and no default value: transects using '
                    'those dates will have no EPRunc or WLR.',
                    f'{n_unc_missing} features sin incertidumbre válida en los '
                    'campos indicados y sin valor por defecto: los transectos '
                    'que usen esas fechas no tendrán EPRunc ni WLR.'))
        if not unc_fields_used and unc_default <= 0:
            feedback.pushInfo(msg(
                'No positional uncertainty: EPRunc and WLR are left empty.',
                'Sin incertidumbre de posición: EPRunc y WLR quedan vacíos.'))

        shorelines = []
        for year, geoms in by_date.items():
            geom = geoms[0] if len(geoms) == 1 else QgsGeometry.collectGeometry(geoms)
            shorelines.append((year, geom))

        if len(shorelines) < 2:
            raise QgsProcessingException(msg(
                'At least 2 different dates with valid shorelines are needed '
                'to compute NSM/EPR/LRR.',
                'Hacen falta al menos 2 fechas distintas con líneas de costa '
                'válidas para calcular NSM/EPR/LRR.'))
        shorelines.sort(key=lambda t: t[0])
        feedback.pushInfo(msg(
            f'{n_feats} shoreline features grouped into '
            f'{len(shorelines)} different dates.',
            f'{n_feats} features de líneas de costa agrupadas en '
            f'{len(shorelines)} fechas distintas.'))

        fields = self._output_fields(conf_pct)
        (sink, dest_id) = self.parameterAsSink(
            parameters, self.OUTPUT, context, fields,
            QgsWkbTypes.LineString, base_crs
        )
        if sink is None:
            raise QgsProcessingException(msg(
                'The output layer could not be created.',
                'No se pudo crear la capa de salida.'))

        half_len = length / 2.0
        n_transects = int(total_length // spacing) + 1
        tid = 0
        n_flag = {FLAG_OK: 0, FLAG_CROSSES: 0, FLAG_AMBIGUOUS: 0}
        n_trend = {}
        out_feats = []

        for i in range(n_transects):
            if feedback.isCanceled():
                break
            dist = i * spacing
            if dist > total_length:
                break
            feedback.setProgress(int(100 * i / max(n_transects, 1)))

            # Dirección de la baseline medida entre dos puntos separados
            # por la distancia de suavizado (mínimo +-1 m). Con suavizado,
            # los transectos siguen la forma general de la costa y no los
            # quiebros de la digitalización.
            half_smooth = max(smooth / 2.0, 1.0)
            d0 = max(0.0, dist - half_smooth)
            d1 = min(total_length, dist + half_smooth)
            if d1 == d0:
                continue

            p0 = baseline_geom.interpolate(d0).asPoint()
            p1pt = baseline_geom.interpolate(d1).asPoint()
            dx, dy = p1pt.x() - p0.x(), p1pt.y() - p0.y()
            norm = (dx ** 2 + dy ** 2) ** 0.5
            if norm == 0:
                continue
            ux, uy = dx / norm, dy / norm
            perp = (-uy, ux)  # rotacion +90 grados del vector tangente

            p_center = baseline_geom.interpolate(dist).asPoint()
            pa = QgsPointXY(p_center.x() + perp[0] * half_len, p_center.y() + perp[1] * half_len)
            pb = QgsPointXY(p_center.x() - perp[0] * half_len, p_center.y() - perp[1] * half_len)
            transect_geom = QgsGeometry.fromPolylineXY([pa, pb])

            per_date = []
            for year, sl_geom in shorelines:
                inter = transect_geom.intersection(sl_geom)
                if inter is None or inter.isEmpty():
                    continue

                candidates = []
                if inter.type() == QgsWkbTypes.PointGeometry:
                    if inter.isMultipart():
                        candidates = inter.asMultiPoint()
                    else:
                        candidates = [inter.asPoint()]
                if not candidates:
                    # corte en linea (solapamiento) u otro tipo raro: se omite
                    continue

                signed = [(c.x() - p_center.x()) * perp[0]
                          + (c.y() - p_center.y()) * perp[1]
                          for c in candidates]
                per_date.append((year, signed))

            if len(per_date) < 2:
                continue

            # Orientación hacia el mar: distancias positivas hacia el mar,
            # así NSM/EPR/LRR > 0 = acreción y < 0 = erosión (como DSAS).
            pts, f, flag = orient_seaward(per_date, onshore)
            n_flag[flag] += 1

            # Transecto de salida dibujado de tierra a mar.
            sx, sy = f * perp[0] * half_len, f * perp[1] * half_len
            out_geom = QgsGeometry.fromPolylineXY([
                QgsPointXY(p_center.x() - sx, p_center.y() - sy),
                QgsPointXY(p_center.x() + sx, p_center.y() + sy),
            ])

            pts.sort(key=lambda t: t[0])
            yr_min, dist_min_v = pts[0]
            yr_max, dist_max_v = pts[-1]
            nsm = dist_max_v - dist_min_v
            dyears = yr_max - yr_min
            epr = nsm / dyears if dyears > 0 else None

            n = len(pts)
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            us = [unc_by_date.get(p[0]) for p in pts]
            eprunc = epr_uncertainty(us[0], us[-1], dyears)
            lrr, r2, lse, lci = linreg_stats(xs, ys, conf)
            wlr, wr2, wse, wci = wlr_stats(xs, ys, us, conf)

            trend, trend_src = classify_trend(
                flag, trend_method, trend_thresh, n,
                epr, eprunc, lrr, lci, wlr, wci)
            n_trend[trend] = n_trend.get(trend, 0) + 1

            feat = QgsFeature(fields)
            feat.setId(tid)
            feat.setGeometry(out_geom)
            feat.setAttributes([tid, n, yr_min, yr_max, dist_min_v, dist_max_v,
                                 nsm, epr, eprunc, lrr, r2, lse, lci,
                                 wlr, wr2, wse, wci, flag, 0,
                                 trend, trend_src])
            out_feats.append(feat)
            tid += 1

        # --- transectos que se cruzan entre sí (curvas cerradas) ---
        index = QgsSpatialIndex()
        for feat in out_feats:
            index.addFeature(feat)
        cross_col = fields.indexOf('tr_cross')
        n_cross = 0
        for k, feat in enumerate(out_feats):
            geom = feat.geometry()
            for other_id in index.intersects(geom.boundingBox()):
                if other_id == feat.id():
                    continue
                if geom.intersects(out_feats[other_id].geometry()):
                    feat.setAttribute(cross_col, 1)
                    n_cross += 1
                    break
        for feat in out_feats:
            sink.addFeature(feat, QgsFeatureSink.FastInsert)

        if tid:
            parts = ', '.join(f'{k}: {n_trend.get(k, 0)}' for k in
                              ('erosion', 'accretion', 'stable', 'unclassified'))
            crit_en = ('statistical' if trend_method == TREND_STATISTICAL
                       else f'threshold {trend_thresh} m/yr')
            crit_es = ('estadístico' if trend_method == TREND_STATISTICAL
                       else f'umbral {trend_thresh} m/año')
            feedback.pushInfo(msg(
                f'trend ({crit_en}), {tid} transects: {parts}.',
                f'trend ({crit_es}), {tid} transectos: {parts}.'))
            if (trend_method == TREND_STATISTICAL
                    and n_trend.get('unclassified', 0) > n_flag[FLAG_CROSSES]
                    + n_flag[FLAG_AMBIGUOUS]):
                feedback.pushInfo(msg(
                    'Some transects are unclassified because they have no '
                    'confidence interval (2 dates without positional '
                    'uncertainty). Add uncertainty or use the Threshold '
                    'criterion.',
                    'Algunos transectos quedan sin clasificar porque no '
                    'tienen intervalo de confianza (2 fechas sin '
                    'incertidumbre de posición). Añade la incertidumbre o '
                    'usa el criterio Umbral.'))

        if n_cross:
            feedback.pushWarning(msg(
                f'{n_cross} transects cross another transect (tr_cross = 1). '
                'This usually happens at sharp bends of the baseline: '
                'increase the smoothing distance or shorten the transects.',
                f'{n_cross} transectos se cruzan con otro transecto '
                '(tr_cross = 1). Suele pasar en curvas cerradas de la línea '
                'base: sube la distancia de suavizado o acorta los '
                'transectos.'))
        if n_flag[FLAG_CROSSES]:
            feedback.pushWarning(msg(
                f'{n_flag[FLAG_CROSSES]} transects with flag = 1: the '
                'baseline crosses a shoreline. Their rates may be unreliable; '
                'check the baseline in those stretches.',
                f'{n_flag[FLAG_CROSSES]} transectos con flag = 1: la línea '
                'base cruza alguna línea de costa. Sus tasas pueden no ser '
                'fiables; revisa la línea base en esos tramos.'))
        if n_flag[FLAG_AMBIGUOUS]:
            feedback.pushWarning(msg(
                f'{n_flag[FLAG_AMBIGUOUS]} transects with flag = 2: the '
                'seaward side could not be determined. Do not use their rates.',
                f'{n_flag[FLAG_AMBIGUOUS]} transectos con flag = 2: no se '
                'pudo saber hacia qué lado está el mar. No uses sus tasas.'))

        if tid == 0:
            feedback.pushWarning(msg(
                'No transect crossed at least 2 valid shorelines. Check the '
                'spacing, the transect length and that the shorelines '
                'actually lie within reach of the baseline.',
                'Ningún transecto cortó al menos 2 líneas de costa válidas. '
                'Revisa el espaciado, la longitud del transecto y que las '
                'líneas de costa queden al alcance de la línea base.'))

        return {self.OUTPUT: dest_id}
