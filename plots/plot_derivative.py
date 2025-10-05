import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Directory containing the TSV files
directory = os.path.dirname(os.path.abspath(__file__))

plt.figure(figsize=(10, 6))

for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        filepath = os.path.join(directory, filename)

        # Read tab-separated file
        df = pd.read_csv(filepath, sep="\t", header=None)
        
        # Assume first col = x, second col = y
        x = df.iloc[:, 0].values
        y = df.iloc[:, 1].values

        # Compute derivative dy/dx
        dydx = np.gradient(y, x)

        # Plot derivative
        plt.plot(x, dydx, label=filename)

plt.xlabel("x")
plt.ylabel("dy/dx")
plt.title("Derivatives of TSV files")
plt.legend()
plt.grid(True)
plt.show()
