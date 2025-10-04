import numpy as np
import matplotlib.pyplot as plt

# Load data from a tab-separated file
data = np.loadtxt("charge.csv", delimiter="\t")  # replace with your filename

x = data[:, 0]
y = data[:, 1]

# Compute derivative
dy_dx = np.gradient(y, x)  # works for uneven spacing

# Plot y vs x
plt.figure(figsize=(10,5))
plt.plot(x, y, label='y(x)')
plt.plot(x, dy_dx, label='dy/dx', linestyle='--')
plt.xlabel('x')
plt.ylabel('y and dy/dx')
plt.title('Data and its Derivative')
plt.legend()
plt.grid(True)
plt.show()
