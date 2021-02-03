'''
This is the configuration class for Dorado which utilizes the Astropy configuration 
system.
'''

from astropy import config as _config

class ConfigNamespace(_config.ConfigNamespace):
    rootname = 'dorado'
    api_key = _config.ConfigItem(
        '',
        "The Astrometry.net API key."
    )

class ConfigItem(_config.ConfigItem):
    rootname = 'dorado'




conf = ConfigNamespace()


print('Dorado directory created: ', _config.get_config_dir('dorado'))

