import numpy as np
from astropy.utils.misc import isiterable


def integrate(f, a, b, N):
        # integrate takes a function, an upper and lower bound and a number of
        # panels and then performs a basic numerical integration technique 
        # consisting of summing the area of panels created from the function 
        # height and a panel width defined by a, b, and N.
        #
        # variables:
        # variable1 = description of variable (data type)
        # return_var = description of what is returned (data type)
        #
        # Invoke example: integrate(E_inv, 0, 2, 200)
        #
        # Author: June Parsons
        # Date: 20190602
        # Version: 1.1.0
        # Last update: 20190602
        # 
        # Comments: 
        if isiterable(b):
                b = np.asarray(b)
                area = []
                for n in range(len(b)):
                        aa = (a+(b[n]-a)/(2*N))
                        bb = (b[n]-(b[n]-a)/(2*N))
                        x = np.linspace(aa, bb, N)
                        fx = f(x)
                        area.append(np.sum(fx)*(b[n]-a)/N)
                area = np.array(area)
            
        else:
                aa = (a+(b-a)/(2*N))
                bb = (b-(b-a)/(2*N))
                
                x = np.linspace(aa, bb, N)
                fx = f(x)
                area = np.sum(fx)*(b-a)/N
        return area

def gauss_func(x):
        a, b, c, d, e, F_0 = [coeffs[1][1], coeffs[1][2], coeffs[1][3], coeffs[1][4], coeffs[1][5], F_0]
        func = 1 - (a*np.exp((-np.power((x-b), 2))/2 * np.power(c,2)) + d * x + e)/F_0

        return func