'''
This is the configuration class for Dorado which utilizes the Astropy configuration 
system.
'''

from astropy import config as _config

rootname = 'dorado'

class ConfigNamespace(_config.ConfigNamespace):
    rootname = 'dorado'
    get_config = _config.get_config(packageormod = rootname, rootname = rootname)
    config_dir = _config.get_config_dir(rootname)
    api_key = _config.ConfigItem(
        '',
        "The Astrometry.net API key."
    )
    use_logger = _config.ConfigItem(
        True, 'Whether to use the Dorado logger system.')
    # _config.generate_config(pkgname = rootname)
    _config.get_config.write()

class ConfigItem(_config.ConfigItem):
    rootname = 'dorado'




conf = ConfigNamespace()
print('Dorado directory created: ', conf.config_dir)
# conf.get_config.write()

__all__ = ['conf']