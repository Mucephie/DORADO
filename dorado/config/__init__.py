'''
This is the configuration class for Dorado which utilizes the Astropy configuration 
system.
'''

from astropy import config as _config
from configparser import ConfigParser
import os

rootname = 'dorado'

class ConfigItem(_config.ConfigItem):
    rootname = 'dorado'

class ConfigNamespace(_config.ConfigNamespace):
    rootname = 'dorado'
    get_config = _config.get_config(packageormod = rootname, rootname = rootname)
    config_dir = _config.get_config_dir(rootname)
    _config.generate_config(pkgname = 'dorado')
    # config = ConfigParser()
    # api_key = config.read('api_key')
    api_key = ConfigItem(
        '',
        "The Astrometry.net API key."
    )
    use_logger = ConfigItem(
        True, 'Whether to use the Dorado logger system.')

    # _config.get_config.write()






conf = ConfigNamespace()

print('Dorado directory created: ', conf.config_dir)
for items in conf.values():
    print(items)
conf.use_logger = False
print(conf.use_logger)

__all__ = ['conf']