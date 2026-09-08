import os
from qgis.PyQt.QtGui import QIcon
from qgis.core import QgsProcessingProvider
from .shoreline_change_algorithm import ShorelineChangeAlgorithm


class ShorelineChangeProvider(QgsProcessingProvider):

    def loadAlgorithms(self):
        self.addAlgorithm(ShorelineChangeAlgorithm())

    def id(self):
        return 'shoreline_change_stats'

    def name(self):
        return 'FIC Shoreline Change Analisys'

    def longName(self):
        return 'FIC Shoreline Change Analisys'

    def icon(self):
        icon_path = os.path.join(os.path.dirname(__file__), 'icon.png')
        return QIcon(icon_path)
