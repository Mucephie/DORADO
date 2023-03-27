import numpy as np

import matplotlib as mpl
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

# this will most likely be renamed visualization 
__all__ = ['TimeSeries_Plot', 'star_chart', 'dorado_mpl_style_1', 'dorado_mpl_style', 'mlty_map']

import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This module contains dictionaries that can be used to set a matplotlib
# plotting style. It is no longer documented/recommended as of Astropy v3.0
# but is kept here for backward-compatibility.

# Version 1 dorado plotting style for matplotlib
dorado_mpl_style_1 = {
    # Lines
    "lines.linewidth": 1.7,
    "lines.antialiased": True,
    # Patches
    "patch.linewidth": 1.0,
    "patch.facecolor": "#348ABD",
    "patch.edgecolor": "#CCCCCC",
    "patch.antialiased": True,
    # Images
    "image.cmap": "mlty",
    "image.origin": "upper",
    # Font
    "font.size": 12.0,
    "font.family": 'serif',
    # Axes
    "axes.facecolor": 'whitesmoke',
    "axes.edgecolor": "#AAAAAA",
    "axes.linewidth": 1.0,
    "axes.grid": True,
    "grid.color": 'gold',
    "grid.linestyle": ':',
    "axes.titlesize": "x-large",
    "axes.labelsize": "large",
    "axes.labelcolor": "k",
    "axes.axisbelow": True,
    # Ticks
    "xtick.major.size": 0,
    "xtick.minor.size": 0,
    "xtick.major.pad": 6,
    "xtick.minor.pad": 6,
    "xtick.color": "#565656",
    "xtick.direction": "in",
    "ytick.major.size": 0,
    "ytick.minor.size": 0,
    "ytick.major.pad": 6,
    "ytick.minor.pad": 6,
    "ytick.color": "#565656",
    "ytick.direction": "in",
    # Legend
    "legend.fancybox": True,
    "legend.loc": "best",
    # Figure
    "figure.figsize": [8, 6],
    "figure.facecolor": '1.0',
    "figure.edgecolor": "0.50",
    "figure.subplot.hspace": 0.5,
    # Other
    "savefig.dpi": 72,
}
color_cycle = [
    "#FB9062",  # pink peach
    "#FF9248",  # peach
    "#BDFD23",  # limish yellow green
    "#F8E604",  # lemon
    "#FDD023",  # bumble bee yellow
    "#FF6700",  # tangerine
    "#FF4E50",  # watermelon red
    "#6A0D83",  # porple
    "#602320",
]  # burgandy

try:
    # This is a dependency of matplotlib, so should be present if matplotlib
    # is installed.
    from cycler import cycler

    dorado_mpl_style_1["axes.prop_cycle"] = cycler("color", color_cycle)
except ImportError:
    dorado_mpl_style_1["axes.color_cycle"] = color_cycle


dorado_mpl_style = dorado_mpl_style_1
"""This is a modified version of the astropy plotting style for use in dorado."""


# map_colours : name MLTY comes from low mass spectral type M, L, T, and Y. 
# it is pronounced 'melty'
mlty_map = ListedColormap(color_cycle, name="mlty")
mpl.colormaps.register(cmap=mlty_map)

# might not need font dict anymore since dorado_mpl_style
font = {'family': 'serif',
        'color':  'k',
        'weight': 'normal',
        'size': 20,
        }
plt.style.use(dorado_mpl_style)

class TimeSeries_Plot:
    def __init__(self, target, fi, cmap_str = 'mlty'):
        self.cmap = plt.colormaps[cmap_str]
        self.target = target
        self.ts = self.target.ts
        if self.ts.table == None:
            self.ts.toTable()
        self.plot(fi)
        # Comparison star(s) add function?

    def plot(self, fi):
        self.fig, self.ax = plt.subplots(1, 2, subplot_kw={'aspect': 'equal'})
        # Produce title string : plt.title('BL-Cam (2023-03-09+10) $\ ^{J.\ Parsons\ @\ Allan\ I.\ Carswell\ Observatory}$', fontdict=font)
        title_str = 'Timeseries' # This needs target name, datestr, observer name, observatory etc.
        plt.title(title_str) 
        self.ax[0].set_xlabel('Time (mjd)')
        self.ax[0].set_ylabel('Flux')
        time = self.ts[self.ts.filters[fi]].times.mjd
        flux = self.ts[self.ts.filters[fi]].flux
        time_fit = self.ts[self.ts.filters[fi]].fit_flux
        flux_fit = self.ts[self.ts.filters[fi]].fit_times
        m = [] # marker styles
        t = [] # marker size transforms
        for i in range(len(flux)):
            # This may be better as optional, maybe size should be inversly proportional to FWHM of image
            # or error on comparison average.
            t.append(Affine2D().scale(10 * np.cos(time[i])+1))
            # This could be star markers ?? or saturn icons
            m.append(MarkerStyle('☺', transform=t)) # DATA_SYMBOLS[]
        self.ax[0].plot(time, flux, marker=m, color=self.cmap(flux), label = 'raw data')
        self.ax[0].plot(time_fit, flux_fit, label = 'fit data')
        self.fig.colorbar(plt.cm.ScalarMappable(cmap=self.cmap), ax=self.ax[0], label="Amplitude")
        self.ax[0].grid() # will this really be needed with the new dorado style?
        # set axis limits
        

        # Fourier plot
        self.ax[1].set_xlabel('Frequency')
        self.ax[1].set_ylabel('Power')
        freq = self.target.freq
        power = self.target.power_vec

        ## add peak markers and text

        plt.tight_layout()
        plt.show()


from astropy.wcs import WCS
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
import astropy.units as u
from astropy.visualization.wcsaxes import add_scalebar, Quadrangle


# marker_colours = ['orangered', 'gold', 'darkorange', 'orange', 'yellow'] # old colour library

class star_chart:
    def __init__(self, im, title = '', wcs = None, cm = 'mlty', grid_c = 'w'):
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
        c = color_cycle[len(color_cycle) % len(self.legend_labels)]
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






