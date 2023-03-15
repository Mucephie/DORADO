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

__all__ = ['TSPlot', 'star_chart']
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


from astropy.visualization import astropy_mpl_style
from astropy.wcs import WCS
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
import astropy.units as u
from astropy.visualization.wcsaxes import add_scalebar, Quadrangle

font = {'family': 'serif',
        'color':  'k',
        'weight': 'normal',
        'size': 20,
        }
marker_colours = ['orangered', 'gold', 'darkorange', 'orange', 'yellow']

class star_chart:
    def __init__(self, im, title = '', wcs = None, cm = 'Greys', grid_c = 'firebrick'):
        self.fig = plt.figure(layout='constrained')
        self.im = im
        self.wcs = wcs
        if wcs:
            # ax.projection=wcs
            self.ax = self.fig.add_subplot(projection=self.wcs)
            self.ax.coords.grid(True, color=grid_c, ls='dotted')
            self.ax.coords[0].set_axislabel('Right Accension')
            self.ax.coords[1].set_axislabel('Declination')
        else:
            self.ax = self.fig.add_subplot()
        plt.title(title + ' Star Chart', fontdict = font)
        up, down = plt_eye(self.im.data)
        self.ax.imshow(self.im.data, cmap = cm, vmin = down, vmax = 3 * up)
        self.legend_labels = []
        self.legend_patchs = []

    def plt_stars(self, stars, label = 'stars'):
        c = marker_colours[len(self.legend_labels)]
        apertures = [Circle((x, y), r, fill = False, facecolor = None, 
            edgecolor = c) for x, y, r in zip(stars['x'], stars['y'], stars['r'])]
        # Create patch collection with specified colour/alpha
        pc = PatchCollection(apertures, facecolor = 'none', edgecolor = c, 
            alpha = 0.75, label = label)
        # Add collection to axes
        self.ax.add_collection(pc)
        self.legend_patchs.append(Circle((0,0), 10, fill = False, 
            facecolor = None, edgecolor = c))
        self.legend_labels.append(label + ' : ' + str(len(stars)))

    def plot(self, legend = True):
        if legend:
            self.fig.legend(self.legend_patchs, self.legend_labels, 
                loc='center right')
        plt.show()
    
    def add_compass(self, length = 0.7 * u.arcmin, loc = 'br'):
        y, x = self.im.data.shape
        co = self.wcs.pixel_to_world(x- 200, y - 300)
        self.barH = Quadrangle((co.ra, co.dec), 1.25 * u.arcmin, 0.09 * u.arcmin , transform=self.ax.get_transform('fk5'), color = 'k')
        self.barW = Quadrangle((co.ra, co.dec), 0.2 * u.arcmin, 0.55 * u.arcmin , transform=self.ax.get_transform('fk5'), color = 'k')
        self.ax.add_patch(self.barH)
        self.ax.add_patch(self.barW)
        # self.ax.arrow(co.ra.value, co.dec.value, 0, 0.01, 
        #     head_width=0, head_length=0, 
        #     fc='firebrick', ec='k', width=0.003, 
        #     transform=self.ax.get_transform('icrs'))
        # self.ax.arrow(co.ra.value, co.dec.value, 0.02, 0, 
        #     head_width=0, head_length=0, 
        #     fc='firebrick', ec='k', width=0.003, 
        #     transform=self.ax.get_transform('icrs'))
        fontComp = {'family': 'serif',
            'color':  'k',
            'weight': 'normal',
            'size': 10,
        }
        plt.text(co.ra.value - 0.003, co.dec.value + 0.012, 'N', 
            color='k', rotation=67, 
            transform=self.ax.get_transform('icrs'), fontdict = fontComp)
        plt.text(co.ra.value + 0.022, co.dec.value + 0.002, 'E', 
            color='k', rotation=67, 
            transform=self.ax.get_transform('icrs'), fontdict = fontComp)

    def add_scale(self):
        add_scalebar(self.ax, 1 * u.arcmin, label="1\"", color="k")


def plt_eye(data):
    mean = np.mean(data)
    std = np.std(data)
    mean_up = mean + 1.5 * std
    mean_down = mean - 1.5 * std
    return mean_up, mean_down

