import numpy as np
import matplotlib.pyplot as plt 

x1 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
y1 = np.array([44.82, 22.729, 15.836, 12.434, 19.706, 16.107, 11.154, 9.940, 16.768, 12.009])

y2 = np.array([2.608, 1.584, 1.304, 1.257, 1.305, 1.051, 1.032, 1.126, 1.157, 1.160])

plt.plot(x1, y2)
plt.xlabel("Number of threads")
plt.ylabel("Time [s]")
plt.title("Timer versus number of threads")
plt.show()