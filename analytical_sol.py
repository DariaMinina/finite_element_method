import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 16,
    "text.usetex": True
})

#fig = plt.figure(figsize=(20, 6))
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.xlabel(r"$x$")
plt.ylabel(r"$u$")
plt.grid(True)
major_ticks = np.arange(-70, 21, 5)
minor_ticks = np.arange(3, 5.5, 0.5)

ax.set_xticks(minor_ticks)
ax.set_yticks(major_ticks)
x = np.arange(3., 5.05, 0.05)
t1 = np.arange(9., 10., 0.5)
c2 = (5 - 12/37)/((37/22)*np.exp(37*3/22))
c1 = 10 - c2 * np.exp(37 * 5 / 22) - 12 * 5 / 37
y = c1 + c2*np.exp(37*x/22) + 12*x/37
plt.plot(x, y, color="red")
#y = (-3806*np.exp(1665/11) - 26270*np.exp(111/22))/(1369*np.exp(111/222)) + (3806/(1369*np.exp(111/22)))*np.exp((37/22)*x) + (12/37)*x
#y = 2*(37*(6*x - 355) + 1903*np.exp(37*(x-3)/22) - 1903*np.exp(3219/22))/1369

plt.savefig('image.png')
