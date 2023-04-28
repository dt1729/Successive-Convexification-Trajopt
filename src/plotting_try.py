import matplotlib.pyplot as plt
import numpy as np

Controls = np.loadtxt("Controls_final.txt",delimiter=',')
v, ω, α = [], [], []

for i in range(np.shape(Controls)[0]):
    v.append(Controls[i,:][0])
    ω.append(Controls[i,:][1])
    α.append(Controls[i,:][2])

# Define the time array
t = np.arange(0.1, 10, 0.1)

# Define the x and y arrays
x = np.sin(t)
y = np.cos(t)

# Create a new figure
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

# Plot the x and y arrays against time
ax1.plot(t, v)
ax1.set(title='Subplot 1', ylabel='Velocity')
ax1.grid(visible=True)
ax2.plot(t, ω)
ax2.set(title='Subplot 2', ylabel='omega')
ax2.grid(visible=True)
ax3.plot(t, α)
ax3.set(title='Subplot 3', xlabel='Time (s)', ylabel='Alpha')


# Add a legend and axis labels
fig.suptitle('Control Signals')
ax3.grid(visible=True)
# Show the plot
plt.show()