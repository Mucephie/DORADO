import warnings
warnings.filterwarnings('ignore')

## file imports
import os
from astropy.io import fits 
#from astropy.nddata import CCDData
from astropy.table import Table

## Processing imports
import numpy as np

## photometry imports
#from astropy import units as u

## other imports
import datetime
#import requests
from astroquery.simbad import Simbad
import dracoOP2


# ## testing
START_DATE_TIME = datetime.datetime.now()
print('\nStarting time: ', START_DATE_TIME)

home_dir = os.getcwd()
print('\nHome dir: ', home_dir)
target = input('Enter the target name (e.g dypeg): ')
usetdy = input('Reduce todays data (y/n)? ')
if (usetdy == 'y'):
        night = dracoOP2.get_night()
else:
        night = input('Enter night code (YYYY-MM-DD+DD): ')

directory = 'C:/Users/mucep/Offline/Draco-test/'
## get VSTAR data from Simbad

def get_field(target_list):
        ra = []
        dec = []
        for g in range(len(target_list)):
                query = Simbad.query_object(target_list[g])
                ra.append(str(query[0]['RA']))
                dec.append(str(query[0]['DEC']).replace('+','').replace('-',''))

        results = Table([np.transpose(np.linspace(0, len(ra)-1, len(ra)))])
        results.rename_column('col0', 'field#') 
        results['ra'] = ra
        results['dec'] = dec
        # results['field#'] = np.linspace(0, len(ra)-1, len(ra))
        results['name'] = target_list
        isV = np.zeros(len(ra))
        isC = np.ones(len(ra))
        isV[0] = 1
        isC[0] = 0
        results['isV'] = isV
        results['isC'] = isC
        # print(results.colnames)

        return results


dypeg = ['dy peg', 'ic 1472', 'GSC 01712-00542']
xxcyg = ['xx cyg', 'GSC 03948-02563', 'TYC 3948-2105-1', 'TYC 3948-2018-1', 'GSC 03948-02266']
blcam = []
yzboo = []
belyn = []

if (target == 'dypeg'):
        inp = dypeg
        x = [1147.43, 613.48, 1067]
        y = [581.12, 100.7, 900]
elif (target == 'xxcyg'):
        inp = xxcyg
        x = []
        y = []
elif (target == 'blcam'):
        inp = blcam
        x = []
        y = []
elif (target == 'yzboo'):
        inp = yzboo
        x = []
        y = []
elif (target == 'belyn'):
        inp = belyn
        x = []
        y = []
else:
        print('\nNot a recognized star!\nexiting...')
        exit

results = get_field(inp)
print('\n{} comp stars found'.format(len(results)))

print('\n')
print(results)

# write results

# get distances

# x = [1147.43, 613.48, 1097.32]
# y = [581.12, 100.7, 912.89]
distances = np.zeros((len(x), len(y)))
for u in range(len(x)):
        for v in range(len(x)):
                separation = np.sqrt((x[u] - x[v])**2 + (y[u] - y[v])**2)
                distances[u, v] = separation


disTable = Table(distances)
for r in range(len(x)):
        disTable.rename_column('col' + str(r), 'star' + str(r+1))

print(disTable)
stardir = directory + night + '/wrk/stars/'
os.chdir(stardir)
disTable.write(str(target) + '-separations.fit')
results.write(str(target) + '-fieldinfo.fit')

# x: 1147.43 y: 581.12 DYPEG xs: 1124.96 ys: 562.11
# x: 1097.32 y: 912.89 GSC xs: 1080 ys: 900
# x: 613.48 y: 100.7 IC 1472 xs: 604.84 ys: 92.09

END_DATE_TIME = datetime.datetime.now()

print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME), '\n')