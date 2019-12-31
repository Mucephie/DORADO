## Processing imports
import numpy as np

## plotting imports
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.visualization import astropy_mpl_style

def plt_eye(data):
    mean = np.mean(data)
    std = np.std(data)
    mean_up = mean + 2 * std
    mean_down = mean - 2 * std

    return mean_up, mean_down

def plt_stars(data, x, y, r):
        plt.style.use(astropy_mpl_style)
        plt.figure()
        # plt.imshow(data, cmap='viridis', vmin=0)
        up, down = plt_eye(data)
        plt.imshow(data, cmap='viridis', vmin=down, vmax=up)
        cbar = plt.colorbar()
        plt.grid(False)
        
        # for i in range(len(x)):
        #         circlei=plt.Circle((x[i],y[i]), r[i], edgecolor=, alpha = 0.75, linewidth = 1)
        #         plt.gcf().gca().add_artist(circlei)
        colours = ['r', 'g', 'b', 'y', 'cyan', 'w', 'm']
        plt.scatter(x, y, s = r ** 2, edgecolors = colours, alpha = 0.5) 
        for p in range(len(x)):
                plt.text(x[p], y[p], str(p))

        plt.show()

def plt_fits(image_data, cmap_str):
    plt.style.use(astropy_mpl_style)
    #print(image_data.shape)

    plt.figure()
    #plt.imshow(np.transpose(np.fliplr(image_data)), cmap='viridis')
    #plt.imshow(image_data, cmap='viridis', norm=LogNorm())
    # plt.imshow(image_data, cmap=cmap_str, vmin=0)
    # edited
    up, down = plt_eye(image_data)
    plt.imshow(image_data, cmap=cmap_str, vmin=down, vmax=up)
    cbar = plt.colorbar() #.ColorbarBase(ax, cmap=cmap_str,norm=mpl.colors.Normalize(vmin=down, vmax=up))
    #cbar.set_clim(down, up)

    # .ColorbarBase(ax, cmap=cm,norm=mpl.colors.Normalize(vmin=-0.5, vmax=1.5))
    #
    plt.grid(False)

    plt.show()
    
def plt_flat(image_data, cmap_str):
    plt.style.use(astropy_mpl_style)
    #print(image_data.shape)

    plt.figure()
    #plt.imshow(np.transpose(np.fliplr(image_data)), cmap='viridis')
    #plt.imshow(image_data, cmap='viridis', norm=LogNorm())
    plt.imshow(image_data, cmap=cmap_str, norm=LogNorm(vmin=image_data.min(), vmax=image_data.max()))
    plt.colorbar()
    plt.grid(False)

    plt.show()
    
def plt_bias(image_data, cmap_str):

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