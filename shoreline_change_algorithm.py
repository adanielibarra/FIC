"""
Shoreline Change Stats
=======================
Replica el flujo basico de DSAS: transectos perpendiculares a una baseline
digitalizada a mano + estadisticas de cambio de linea de costa por transecto.

SUPUESTOS QUE ASUMO (léelos antes de usarlo, cambian el resultado):

1. La baseline debe ser una unica linea simple (una sola feature, una sola
   parte). Si tiene varias partes intento fusionarlas con mergeLines(), si
   no puedo, el algoritmo falla y te lo dice.

2. El SENTIDO del signo de NSM/EPR/LRR depende de en qué dirección
   digitalizaste la baseline. La perpendicular se calcula rotando +90°
   (sentido antihorario) el vector tangente de la baseline en cada punto.
   Si digitalizaste la baseline siempre "de izquierda a derecha mirando
   hacia el mar", el signo positivo = el mar se aleja = acreción; si la
   digitalizaste al revés, el signo se invierte. ESO NO LO CALCULO YO,
   lo tienes que comprobar tú mirando un par de transectos en el mapa.

3. Los transectos se generan simétricos (mitad hacia cada lado del punto
   de la baseline), no solo "hacia el mar" como en algunas config. de DSAS.

4. Para cada shoreline, si el transecto la corta más de una vez, me quedo
   con el punto de corte más cercano al punto de la baseline (igual que
   hace DSAS con la opción "closest").

5. La capa de shorelines debe tener un campo de FECHA (tipo Date, o texto
   en formato ISO YYYY-MM-DD). Si el campo no se puede interpretar como
   fecha en alguna feature, esa feature se descarta (te avisa en el log).

6. Hacen falta minimo 2 shorelines con fecha válida para poder calcular
   NSM/EPR, y para que el R2 de la regresión LRR tenga algún sentido
   necesitas 3 o más (con 2 puntos el R2 siempre da 1.0 y no dice nada).

No lo he podido probar dentro de QGIS real (aquí no tengo QGIS instalado),
así que verifica con un caso pequeño antes de correrlo sobre todo tu área
de estudio.
"""

from qgis.PyQt.QtCore import QVariant
from qgis.PyQt.QtGui import QIcon
from qgis.core import (
    QgsProcessingAlgorithm,
    QgsProcessingParameterVectorLayer,
    QgsProcessingParameterField,
    QgsProcessingParameterNumber,
    QgsProcessingParameterFeatureSink,
    QgsProcessingException,
    QgsFeature,
    QgsFeatureSink,
    QgsFields,
    QgsField,
    QgsGeometry,
    QgsPointXY,
    QgsWkbTypes,
)


class ShorelineChangeAlgorithm(QgsProcessingAlgorithm):

    BASELINE = 'BASELINE'
    SHORELINES = 'SHORELINES'
    DATE_FIELD = 'DATE_FIELD'
    SPACING = 'SPACING'
    LENGTH = 'LENGTH'
    OUTPUT = 'OUTPUT'

    def createInstance(self):
        return ShorelineChangeAlgorithm()

    def name(self):
        return 'shoreline_change_stats'

    def displayName(self):
        return 'FIC Shoreline Change Analisys'

    def group(self):
        return 'Coastal'

    def groupId(self):
        return 'coastal'

    def icon(self):
        import os
        icon_path = os.path.join(os.path.dirname(__file__), 'icon.png')
        return QIcon(icon_path)

    def shortHelpString(self):
        return (
            '<p>Plugin para el calculo de la erosion costera.<br>'
            'Creado por Daniel Ibarra Marinas, '
            '<span style="color:#2e8a3d"><b>Facultad de Ingenieria '
            'y Ciencias</b></span>, '
            '<span style="color:#d97706"><b>Universidad Autonoma de '
            'Tamaulipas</b></span>.</p>'

            '<p>Genera <b>transectos</b> perpendiculares a una baseline '
            '(una sola linea, digitalizada por ti) espaciados cada X '
            'metros. Para cada transecto, corta cada shoreline de la capa '
            'multitemporal y calcula tres metricas de cambio:</p>'

            '<p><b>NSM (Net Shoreline Movement):</b> distancia entre la '
            'shoreline mas antigua y la mas reciente que cruzan ese '
            'transecto, en metros. Es el cambio total en el periodo '
            'estudiado, sin tener en cuenta cuantas shorelines intermedias '
            'haya.</p>'

            '<p><b>EPR (End Point Rate):</b> el NSM dividido entre los '
            'años transcurridos entre la primera y la ultima shoreline. '
            'Solo usa esos dos extremos, ignora las shorelines '
            'intermedias aunque existan. Es rapido de calcular pero '
            'sensible a que esas dos fechas concretas sean '
            'representativas (si una de las dos fue un año anomalo, '
            'el EPR sale sesgado).</p>'

            '<p><b>LRR (Linear Regression Rate):</b> pendiente de la '
            'recta de regresion lineal de la distancia de todas las '
            'shorelines frente al año, no solo las dos extremas. Viene '
            'acompañada del <b>R2</b> de esa regresion (que tan bien se '
            'ajustan los puntos a una linea recta). Con solo 2 '
            'shorelines el R2 siempre da 1.0 y no significa nada, hacen '
            'falta 3 o mas fechas para que el LRR y su R2 tengan sentido '
            'real.</p>'

            '<p style="color:#6b6b6b">'
            '<b>IMPORTANTE:</b> el signo de las tasas (erosion o acrecion) '
            'depende del sentido en que digitalizaste la baseline. '
            'Revisalo en un par de transectos antes de interpretar el '
            'mapa completo.</p>'
        )


    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterVectorLayer(
            self.BASELINE, 'Baseline (una sola linea)'))
        self.addParameter(QgsProcessingParameterVectorLayer(
            self.SHORELINES, 'Shorelines multitemporales (lineas)'))
        self.addParameter(QgsProcessingParameterField(
            self.DATE_FIELD, 'Campo de fecha en la capa de shorelines',
            parentLayerParameterName=self.SHORELINES))
        self.addParameter(QgsProcessingParameterNumber(
            self.SPACING, 'Espaciado entre transectos (m)',
            type=QgsProcessingParameterNumber.Double,
            defaultValue=50.0, minValue=0.1))
        self.addParameter(QgsProcessingParameterNumber(
            self.LENGTH, 'Longitud total del transecto (m)',
            type=QgsProcessingParameterNumber.Double,
            defaultValue=500.0, minValue=1.0))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT, 'Transectos con NSM/EPR/LRR'))

    def _output_fields(self):
        fields = QgsFields()
        fields.append(QgsField('id', QVariant.Int))
        fields.append(QgsField('n_pts', QVariant.Int))
        fields.append(QgsField('yr_min', QVariant.Double))
        fields.append(QgsField('yr_max', QVariant.Double))
        fields.append(QgsField('dist_min', QVariant.Double))
        fields.append(QgsField('dist_max', QVariant.Double))
        fields.append(QgsField('NSM', QVariant.Double))
        fields.append(QgsField('EPR', QVariant.Double))
        fields.append(QgsField('LRR', QVariant.Double))
        fields.append(QgsField('LRR_R2', QVariant.Double))
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
                f"No pude interpretar como fecha el valor {value!r}, "
                "se descarta esa shoreline."
            )
            return None

    def processAlgorithm(self, parameters, context, feedback):
        baseline_layer = self.parameterAsVectorLayer(parameters, self.BASELINE, context)
        shorelines_layer = self.parameterAsVectorLayer(parameters, self.SHORELINES, context)
        date_field = self.parameterAsString(parameters, self.DATE_FIELD, context)
        spacing = self.parameterAsDouble(parameters, self.SPACING, context)
        length = self.parameterAsDouble(parameters, self.LENGTH, context)

        if baseline_layer is None or shorelines_layer is None:
            raise QgsProcessingException('Falta la baseline o la capa de shorelines.')

        # --- baseline: una sola geometria continua ---
        feats = list(baseline_layer.getFeatures())
        if not feats:
            raise QgsProcessingException('La capa de baseline no tiene features.')
        if len(feats) > 1:
            feedback.pushWarning(
                'La baseline tiene mas de una feature, solo se usa la primera. '
                'Fusiona las lineas antes si querias usarlas todas.'
            )
        baseline_geom = feats[0].geometry()
        if baseline_geom.isMultipart():
            merged = baseline_geom.mergeLines()
            if merged and not merged.isEmpty():
                baseline_geom = merged
            else:
                raise QgsProcessingException(
                    'La baseline es multiparte y no se pudo fusionar en una '
                    'sola linea continua.'
                )

        total_length = baseline_geom.length()
        if total_length <= 0:
            raise QgsProcessingException('La baseline tiene longitud 0.')

        # --- shorelines con fecha valida ---
        shorelines = []
        for f in shorelines_layer.getFeatures():
            year = self._decimal_year(f.attribute(date_field), feedback)
            if year is None:
                continue
            geom = f.geometry()
            if geom is None or geom.isEmpty():
                continue
            shorelines.append((year, geom))

        if len(shorelines) < 2:
            raise QgsProcessingException(
                'Hacen falta al menos 2 shorelines con fecha valida para '
                'calcular NSM/EPR/LRR.'
            )
        shorelines.sort(key=lambda t: t[0])

        fields = self._output_fields()
        (sink, dest_id) = self.parameterAsSink(
            parameters, self.OUTPUT, context, fields,
            QgsWkbTypes.LineString, baseline_layer.crs()
        )
        if sink is None:
            raise QgsProcessingException('No se pudo crear la capa de salida.')

        half_len = length / 2.0
        n_transects = int(total_length // spacing) + 1
        tid = 0

        for i in range(n_transects):
            if feedback.isCanceled():
                break
            dist = i * spacing
            if dist > total_length:
                break
            feedback.setProgress(int(100 * i / max(n_transects, 1)))

            d0 = max(0.0, dist - 1.0)
            d1 = min(total_length, dist + 1.0)
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

            pts = []
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

                best_signed = None
                best_abs = None
                for c in candidates:
                    vx, vy = c.x() - p_center.x(), c.y() - p_center.y()
                    signed = vx * perp[0] + vy * perp[1]
                    if best_abs is None or abs(signed) < best_abs:
                        best_abs = abs(signed)
                        best_signed = signed
                pts.append((year, best_signed))

            if len(pts) < 2:
                continue

            pts.sort(key=lambda t: t[0])
            yr_min, dist_min_v = pts[0]
            yr_max, dist_max_v = pts[-1]
            nsm = dist_max_v - dist_min_v
            dyears = yr_max - yr_min
            epr = nsm / dyears if dyears > 0 else None

            n = len(pts)
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            mean_x = sum(xs) / n
            mean_y = sum(ys) / n
            ss_xy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
            ss_xx = sum((x - mean_x) ** 2 for x in xs)
            if ss_xx == 0:
                lrr, r2 = None, None
            else:
                lrr = ss_xy / ss_xx
                y_pred = [mean_y + lrr * (x - mean_x) for x in xs]
                ss_res = sum((y - yp) ** 2 for y, yp in zip(ys, y_pred))
                ss_tot = sum((y - mean_y) ** 2 for y in ys)
                r2 = 1 - ss_res / ss_tot if ss_tot > 0 else None

            feat = QgsFeature(fields)
            feat.setGeometry(transect_geom)
            feat.setAttributes([tid, n, yr_min, yr_max, dist_min_v, dist_max_v,
                                 nsm, epr, lrr, r2])
            sink.addFeature(feat, QgsFeatureSink.FastInsert)
            tid += 1

        if tid == 0:
            feedback.pushWarning(
                'No se genero ningun transecto con al menos 2 shorelines '
                'validas. Revisa el espaciado, la longitud del transecto y '
                'que las shorelines realmente crucen la baseline.'
            )

        return {self.OUTPUT: dest_id}
