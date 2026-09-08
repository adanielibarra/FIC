from qgis.core import QgsApplication
from .shoreline_provider import ShorelineChangeProvider


class ShorelineChangePlugin:
    """Clase minima: solo registra/desregistra el provider de Processing.
    No añade dialogos propios, todo se usa desde la Caja de herramientas
    de Processing (busca 'Shoreline change')."""

    def __init__(self, iface):
        self.iface = iface
        self.provider = None

    def initGui(self):
        self.provider = ShorelineChangeProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def unload(self):
        QgsApplication.processingRegistry().removeProvider(self.provider)
