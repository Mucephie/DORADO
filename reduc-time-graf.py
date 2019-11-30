import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

lights = np.sort([222, 186, 104, 170, 155])
time = np.sort([((8*60)+41) / 60, ((7*60)+15) / 60, ((4*60)+43) / 60, ((6*60)+46) / 60, ((6*60)+1) / 60])

plt.figure()
plt.plot(time, lights, 'r--')
plt.xlabel('time (min)')
plt.ylabel('lights (#)')
plt.title('Computing time vs. number of images')

plt.show()