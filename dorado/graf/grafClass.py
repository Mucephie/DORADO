import numpy as np

import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
from matplotlib.transforms import Affine2D
from matplotlib.textpath import TextPath

from astropy.time import Time

DATA_SYMBOLS = [
    TextPath((0, 0), "☹"),
    TextPath((0, 0), "😒"),
    TextPath((0, 0), "☺"),
]

__all__ = ['TSPlot']
class TSPlot:
    def __init__(self, target, cmap_str = 'plasma'):
        cmap = plt.colormaps[cmap_str]

        self.target = target
        self.ts = self.target.ts
        if self.ts.table == None:
            self.ts.toTable()
    def plot(self, fi):
        fig, ax = plt.subplots(1, 2, subplot_kw={'aspect': 'equal'})
        ax[0].set_xlabel('Time')
        ax[0].set_ylabel('Flux')
        time = self.ts[self.ts.filters[fi]].times.mjd
        flux = self.ts[self.ts.filters[fi]].flux
        time_fit = self.ts[self.ts.filters[fi]].fit_flux
        flux_fit = self.ts[self.ts.filters[fi]].fit_times
        m = []
        t = []
        for i in range(len(flux)):
            t.append(Affine2D().scale(10 * np.cos(time[i])+1))
            m.append(MarkerStyle('☺', transform=t)) # DATA_SYMBOLS[]
        ax[0].plot(time, flux, marker=m, color=cmap(flux), label = 'raw data')
        ax[0].plot(time_fit, flux_fit, label = 'fit data')
        fig.colorbar(plt.cm.ScalarMappable(cmap=cmap), ax=ax[0], label="Amplitude")
        ax[0].grid()

        ax[1].set_xlabel('Frequency')
        ax[1].set_ylabel('Power')
        freq = self.target.freq
        power = self.target.power_vec

        ## add peak markers and text


        plt.tight_layout()
        plt.show()