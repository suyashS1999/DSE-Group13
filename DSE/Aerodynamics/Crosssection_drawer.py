import numpy as np
import matplotlib.pyplot as plt

#theta_upper = np.linspace(0,106.142,106143)
#theta_upper = theta_upper-16.142
theta_upper = np.linspace(0,90,91)
theta_upper = np.radians(theta_upper)
#theta_lower = np.linspace(0,89.882,89883)
theta_lower = np.linspace(0,90,91)
theta_lower = theta_lower-90
theta_lower = np.radians(theta_lower)


r_upper = 1909.271 + 50
r_lower = 1834.004 + 50
k_upper = 430.811
k_lower = -96.223

x_upper = r_upper*np.cos(theta_upper)
x_lower = r_lower*np.cos(theta_lower)
y_upper = k_upper + r_upper*np.sin(theta_upper)
y_lower = k_lower + r_lower*np.sin(theta_lower)

plt.plot(x_upper,y_upper)
plt.plot(x_lower,y_lower)
plt.show()
print(x_upper[0])

x_upper = np.flip(x_upper, 0)
x_upper[0] = 0
y_upper = np.flip(y_upper, 0)
x_lower = np.flip(x_lower, 0)
x_lower[-1] = 0
y_lower = np.flip(y_lower, 0)

x = np.concatenate([x_upper,x_lower])
y = np.concatenate([y_upper,y_lower])
x = np.array([x])
y = np.array([y])
x = x/np.max(x)
y = y/np.max(y)

# func = np.polyfit(x[0],y[0],12)
# nodes = np.linspace(0,1,10001)
# def reconstruct_f(coeff, plot_nodes):
# 	F = 0;
# 	coeff = coeff[::-1];
# 	for i in range(len(coeff)):
# 		F += coeff[i]*plot_nodes**i;
# 	return F;

# funcdraw = reconstruct_f(func,nodes)

xy = np.vstack([x,y])
xy = np.transpose(xy)
print(len(xy[:]))
np.savetxt("Crosssection",xy)

