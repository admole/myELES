import numpy as np
import matplotlib.pyplot as plt

Nx = 20
Ny = 10


x = np.arange(0, 20, 20/Nx)
y = np.arange(0, 10, 10/Ny)

positionx, positiony = np.meshgrid(x, y, indexing='ij')
position = np.array([positionx, positiony])
print(np.shape(position))

bmin = np.array([np.min(x), np.min(y)])
bmax = np.array([np.max(x), np.max(y)])
bsize = bmax - bmin
# max_dist = np.max(bsize)
max_dist = np.sqrt(bsize[0]**2+bsize[1]**2)
print(max_dist)

vel = np.array([position[1]*2-10, position[0]*1-5])
# vel = np.array([(position[0]-5)*-1, (position[1]-5)*-1])
# vel = np.array([position[1]*0+10, position[0]*0])
# vel = np.array([position[1]*1-5, position[0]*0])

flowDir = vel / np.maximum(np.sqrt(vel[0]**2 + vel[1]**2), 0.0001) #

dist = np.empty(shape=np.shape(position[0]))
# max_dists = np.empty(shape=np.shape(position[0]))

for i in range(Nx):
    for j in range(Ny):
        dists = np.empty(2)
        # print(flowDir[:, i, j])

        # if np.all(flowDir[:, i, j] == 0):
        #     print('fixing 0 flowDir')
        #     flowDir[:, i, j] = 0.5*np.array([np.sqrt(1), np.sqrt(1)])

        for k in [0, 1]:
            if flowDir[k, i, j] > 0.0:
                dist_normal = bmax[k] - position[k, i, j]
                dists[k] = dist_normal / (flowDir[k, i, j])
            elif flowDir[k, i, j] < 0.0:
                dist_normal = position[k, i, j] - bmin[k]
                dists[k] = dist_normal / (flowDir[k, i, j])
            elif np.all(flowDir[:, i, j] == 0):
                # when vel = 0 set flow dir straight towards nearest wall
                dist_normal = min(position[k, i, j] - bmin[k], bmax[k] - position[k, i, j])# , key=abs)
                dists[k] = dist_normal
                # flowDir[k, i, j] = 0.0001
            else:
                dists[k] = max_dist


            # dists[k] = dist_normal/(flowDir[k, i, j]) #
            # dists[k] = dist_normal
        dist[i, j] = np.min(abs(dists))
        # dist[i, j] = np.dot(dists, flowDir[:, i, j])
        # dist[i, j] = np.sqrt(dists[0]**2 + dists[1]**2)
        # max_dists[i, j] = bsize[np.argmin(abs(dists))]
norm_dist = abs(dist) / max_dist
wi = (1 - norm_dist)
wiexp_orig = 1 - (np.exp(norm_dist**3.5)-1) / (np.exp(1)-1)
wiexp = 1 - (np.exp(norm_dist**2)-1) / (np.exp(1)-1)


fig, axes = plt.subplots(1, 1)

im = axes.imshow(wi.T, cmap='viridis', origin='lower', alpha=1,
                 extent=(bmin[0]-10/(2*Nx), bmax[0]+10/(2*Nx), bmin[1]-10/(2*Nx), bmax[1]+10/(2*Nx)),
                 vmin=0.0, vmax=1.0)
fig.colorbar(im, ax=axes)
plt.quiver(position[0], position[1], vel[0], vel[1])
# plt.quiver(position[0], position[1], flowDir[0], flowDir[1], color='r')

plt.show()
