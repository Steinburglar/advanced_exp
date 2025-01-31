import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the function to fit
def gaussian(x, a, mu, sigma):
    return a * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))



# Generate synthetic data
np.random.seed(42)
x_data = np.linspace(-5, 5, 100)
y_data = gaussian(x_data, 2.5, 0, 1) + 0.2 * np.random.normal(size=len(x_data))

# Perform the curve fitting
popt, pcov = curve_fit(gaussian, x_data, y_data, p0=[2, 0, 1])

# Extract fitted parameters
a_fit, mu_fit, sigma_fit = popt
print(f"Fitted parameters: a={a_fit:.3f}, mu={mu_fit:.3f}, sigma={sigma_fit:.3f}")

# Plot the data and the fitted function
plt.scatter(x_data, y_data, label='Data', color='blue', alpha=0.5)
plt.plot(x_data, gaussian(x_data, *popt), label='Fitted function', color='red')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Curve Fitting Example')
plt.show()
