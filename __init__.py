def classFactory(iface):
    from .shoreline_change_plugin import ShorelineChangePlugin
    return ShorelineChangePlugin(iface)
