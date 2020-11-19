## Processing imports
import numpy as np

## plotting imports
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from astropy.visualization import astropy_mpl_style

__all__ = ['plt_eye', 'plt_stars', 'plt_fits', 'plt_flat', 'plt_bias']

def plt_eye(data, num_sigma = 3):
    """
    plt_eye takes image data for plotting and computes the optimal colourbar range.

    Parameters
    ----------
    data: image/array
            image data in 2-dimensional array format.
    num_sigma: float
            Positive float value that defines how many standard deviations are to be considered when computing the range.
    Returns
    -------
    mean_up, mean_down: float
           Upper and lower colourbar limits for plotting.
    """
    mean = np.mean(data)
    std = np.std(data)
    mean_up = mean + num_sigma * std
    mean_down = mean - num_sigma * std

    return mean_up, mean_down

def plt_stars(data, x, y, r, num_sigma = 3):
        plt.style.use(astropy_mpl_style)
        plt.figure()
        # plt.imshow(data, cmap='viridis', vmin=0)
        up, down = plt_eye(data, num_sigma=num_sigma)
        plt.imshow(data, cmap='viridis', vmin=down, vmax=up)
        cbar = plt.colorbar()
        plt.grid(False)
        
        # for i in range(len(x)):
        #         circlei=plt.Circle((x[i],y[i]), r[i], edgecolor=, alpha = 0.75, linewidth = 1)
        #         plt.gcf().gca().add_artist(circlei)
        colours = ['r', 'g', 'b', 'y', 'cyan', 'w', 'm']
        plt.scatter(x, y, s = r ** 2, edgecolors = colours, alpha = 0.2) 
        for p in range(len(x)):
                plt.text(x[p], y[p], str(p))

        plt.show()

def plt_fits(image_data, cmap_str = 'viridis', num_sigma = 3):
    plt.style.use(astropy_mpl_style)
    #print(image_data.shape)

    plt.figure()
    #plt.imshow(np.transpose(np.fliplr(image_data)), cmap='viridis')
    #plt.imshow(image_data, cmap='viridis', norm=LogNorm())
    # plt.imshow(image_data, cmap=cmap_str, vmin=0)
    # edited
    up, down = plt_eye(image_data, num_sigma = num_sigma)
    plt.imshow(image_data, cmap=cmap_str, vmin=down, vmax=up)
    cbar = plt.colorbar() #.ColorbarBase(ax, cmap=cmap_str,norm=mpl.colors.Normalize(vmin=down, vmax=up))
    #cbar.set_clim(down, up)

    # .ColorbarBase(ax, cmap=cm,norm=mpl.colors.Normalize(vmin=-0.5, vmax=1.5))
    #
    plt.grid(False)

    plt.show()
    
def plt_flat(image_data, cmap_str = 'viridis'):
    plt.style.use(astropy_mpl_style)
    #print(image_data.shape)

    plt.figure()
    #plt.imshow(np.transpose(np.fliplr(image_data)), cmap='viridis')
    #plt.imshow(image_data, cmap='viridis', norm=LogNorm())
    plt.imshow(image_data, cmap=cmap_str, norm=LogNorm(vmin=image_data.min(), vmax=image_data.max()))
    plt.colorbar()
    plt.grid(False)

    plt.show()
    
def plt_bias(image_data, cmap_str = 'viridis'):

    plt.style.use(astropy_mpl_style)
    #print(image_data.shape)

    plt.figure()
    #plt.imshow(np.transpose(np.fliplr(image_data)), cmap='viridis')
    #plt.imshow(image_data, cmap='viridis', norm=LogNorm())
    plt.imshow(image_data, cmap=cmap_str, vmin=0)
    plt.colorbar()
    plt.grid(False)

    plt.show()


# def canyon (curve analysis @ York )
# remove uneeded functions