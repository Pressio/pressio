import numpy as np
import matplotlib.pyplot as plt


dt = 0.1
t0 = 0.0
L = 128.0

data = np.loadtxt("primal.dat")

nnode = data.shape[0]
nt = data.shape[1]

# append zeros
u = np.vstack([np.zeros(nt),data,np.zeros(nt)])

print(u.shape)

# create grid
x = np.linspace(0.0,L,nnode+2)

# create time
t = np.linspace(t0,t0 + nt*dt,nt)

x_grid,t_grid = np.meshgrid(x,t)

plt.contourf(x_grid,t_grid,u.T,81)

plt.xlabel('x')
plt.ylabel('t')

plt.show()

