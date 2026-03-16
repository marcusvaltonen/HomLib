import numpy as np
import matplotlib.pyplot as plt
import cv2




img1_name = 'DSC_2519'
img1 = cv2.cvtColor(cv2.imread(f'/home/elramav/Downloads/drive-download-20250816T102202Z-1-001/fisheye_grossmunster/images/{img1_name}.JPG'), cv2.COLOR_BGR2RGB)

kappa = -1.1302514855113225e-07

x = np.arange(0,img1.shape[1]) - img1.shape[1]/2
y = np.arange(0,img1.shape[0]) - img1.shape[0]/2
xv, yv = np.meshgrid(x, y)

for i in range(img1.shape[1]):
    for j in range(img1.shape[0]):
        r2 = xv[j,i]**2 + yv[j,i]**2
        xv[j, i] /= 1 + kappa*r2
        yv[j, i] /= 1 + kappa * r2

fig, ax = plt.subplots()
ax.pcolormesh(xv, yv, img1.astype(float) / 255.0)
ax.yaxis.set_inverted(True)

plt.axis('off')
plt.savefig(f"{img1_name}.png", bbox_inches='tight')

plt.show()