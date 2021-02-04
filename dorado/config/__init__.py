'''
This is the configuration class for Dorado which utilizes the Astropy configuration 
system.
'''

from astropy import config as _config

class ConfigNamespace(_config.ConfigNamespace):
    rootname = 'dorado'
    config_dir = _config.get_config_dir(rootname)
    _config.get_config(packageormod = rootname, rootname = rootname)
    api_key = _config.ConfigItem(
        '',
        "The Astrometry.net API key."
    )
    use_logger = _config.ConfigItem(
        True, 'Whether to use the Dorado logger system.')
    # _config.generate_config(pkgname = rootname)


class ConfigItem(_config.ConfigItem):
    rootname = 'dorado'




conf = ConfigNamespace()


# print('Dorado directory created: ', _config.get_config_dir('dorado'))

__all__ = ['conf']