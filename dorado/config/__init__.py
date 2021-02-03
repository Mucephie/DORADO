'''
This is the configuration class for Dorado which utilizes the Astropy configuration 
system.
'''

import astropy.config as astropyconfig


class ConfigNamespace(astropyconfig.ConfigNamespace):
    rootname = 'dorado'


class ConfigItem(astropyconfig.ConfigItem):
    rootname = 'dorado'