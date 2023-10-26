import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
T = [[[11, 12, 5, 2], [15, 6,10, 3], [10, 8, 12, 5], [12,15,8,6]], [[11, 12, 5, 2], [15, 6,10, 3], [10, 8, 12, 5], [12,15,8,6]]]
shape = (4,4)
Intensity = []
for i in range(0, len(T)):
    Target = T[i]
    Intensity = Intensity + [Target[1][1]]
print(Intensity)