import numpy as np
from astropy import config as _config
import os
import datetime
'''
Clippy is the handler of the Dorado system,
'''

# __all__ = []
rootname = 'dorado'
config_dir = _config.get_config_dir(rootname)


# os.environ['HOME'] is also the home directory if needed
# os.uname() for returning operating system information

# os.open os.read os.write
# os.get_terminal_size (gets columns, lines)

# os.chdir for choosing current directory path
# os.getcwd() get the current working directory
# os.listdir(path='.') list directory contents by name

# os.mkdir(path) make a directory
# os.mkdirs(name, exist_ok = False) make directory and intermediate directories

# os.remove(path) remove/delete a file path (not for directories)
# os.rmdir(path) remove/delete a directory path (if dir is empty)
# os.removedirs(name) recursively removes directories until it finds a non empty directory

# os.rename(src, dst) rename a file or directory from src to dst
# os.scandir(path='.') better than listdir, you can specify a file descriptor
# with os.scandir(path) as it:
#     for entry in it:
#         if not entry.name.startswith('.') and entry.is_file():
#             print(entry.name)

# os.startfile(path[, operation]) acts like double clicking file in explorer to open with default
# os.system(command) execute a system command string in subshell

# getpass.getuser() similar to os.getlogin() gets username or login name. 


class clippy:

    def __init__(self):
        # open and use logger
        # make function to create data class from hardware, processed, or raw data folder
        # needs function to import data into raw
        # fix get_night()
        # add new calibration frames to folders
        # find most recetn calibration frame if none
        # -> figure out if none
        # clear cache function
        self.init_dir()

    def init_dir(self):
        stardir = os.getcwd()
        os.chdir(config_dir)
        os.makedirs('../data/wrk', exist_ok = True)
        os.makedirs('../data/flats', exist_ok = True)
        os.makedirs('../data/bias', exist_ok = True)
        os.makedirs('../data/darks', exist_ok = True)
        os.makedirs('../data/raw', exist_ok = True)
        os.makedirs('../data/graphical', exist_ok = True)
        os.makedirs('../data/projects', exist_ok = True)
        os.makedirs('../targets', exist_ok = True)
        os.makedirs('../logs', exist_ok = True)
        os.makedirs('../cache', exist_ok = True)
        os.chdir(stardir)

    def newdat(self):
        # find data that hasn't been processed yet
        print('searching for unprocessed data...')
    
    def make_ceres(self, input):
        # grab data from hardware storage and construct instance of ceres class
        output = input ()
        return output

    def get_night(self):
        """
        get_night obtains a timestring for the most recent(previous) night based on local/hardware
        time provided by datetime. The format follows yyyy-mm-(dd-1)+dd where (dd-1) is last nights 
        day of the month. This currently does not support dates at which the start of the observing
        night was the last day of the month.

        Parameters
        ----------
        None

        Returns
        -------
        night: str
                Timestring for the most recent night.
        """
        # currently does not support first/last of the month
        date = datetime.date.today()
        year = date.year
        month = date.month
        date2 = date.day
        date1 = date2 - 1
        night = str(year) + '-' + str(month) + '-' + str(date1) + '+' + str(date2)

        return night



clip = clippy()
# try:
#     fp = open("myfile")
# except PermissionError:
#     return "some default data"
# else:
#     with fp:
#         return fp.read()