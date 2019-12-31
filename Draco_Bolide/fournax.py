import numpy as np

def fournax(x, y, terms):
    # Fournax is an abbreviation of Fourier numerical astronomy extension, its name is a backronym styled to match the constellation 'fornax'.
    
    # Fourier series domain
    tau = (np.max(x) - np.min(x))
    # Compute real valued Fourier transform
    f = np.fft.rfft(y)
    # Null or zero coefficients above ammount of series "terms"
    # This corresponds to undesired high-frequency terms
    f[terms+1:] = 0
    # Collapse back into function space, result is smoothed Fourier curve
    Fcurve = np.fft.irfft(f)
    
    return Fcurve


# def powerranger
# def maxlight
# def properties

