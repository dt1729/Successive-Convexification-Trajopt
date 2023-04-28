import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.transform import Rotation as R

states = np.loadtxt("States_landing.txt",delimiter=',')
points, angles_rad_o = [], []

for i in range(np.shape(states)[0]):
    points.append(states[i,:][0:3].tolist())
    angles_rad_o.append(states[i,:][3:6].tolist())

def rot_xyz(theta_x,theta_y,theta_z):
    rot_x = np.array([
        [1, 0, 0, 0],
        [0, np.cos(theta_x), -np.sin(theta_x), 0],
        [0, np.sin(theta_x), np.cos(theta_x), 0],
        [0, 0, 0, 1]
    ])

    rot_y = np.array([
        [np.cos(theta_y), 0, np.sin(theta_y), 0],
        [0, 1, 0, 0],
        [-np.sin(theta_y), 0, np.cos(theta_y), 0],
        [0, 0, 0, 1]
    ])

    rot_z = np.array([
        [np.cos(theta_z), -np.sin(theta_z), 0, 0],
        [np.sin(theta_z), np.cos(theta_z), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])

    rot_xyz = rot_z.dot(rot_y).dot(rot_x)
    return rot_xyz


points = np.array(points)
print(points)
breakpoint()
angles_deg = np.array([[0, 0, 0], [45, 60, 30], [-30, 90, 0]])

# Convert the Euler angles to Cartesian coordinates
angles_rad = np.radians(angles_deg)
rotations = R.from_euler('zyx', angles_rad_o, degrees=True)
vectors = []
for i in range(np.shape(points)[0]):
    r1 = rot_xyz(angles_rad_o[i][0], angles_rad_o[i][1],angles_rad_o[i][2])
    vectors.append(r1.dot(np.append(points[i],1))[0:3]/(np.sum(r1.dot(np.append(points[i],1))**2)))
vectors = np.array(vectors)
# vectors = rotations.apply(np.eye(3))

# Create a 3D figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the points
ax.scatter(points[:,0], points[:,1], points[:,2], color='b')

# Plot the direction arrows
for i in range(len(points)):
    ax.quiver(points[i,0], points[i,1], points[i,2], vectors[i,0], vectors[i,1], vectors[i,2], length=1.0, color='r')

# Set the axis limits
ax.set_xlim([np.min(points[:,0])-1, np.max(points[:,0])+1])
ax.set_ylim([np.min(points[:,1])-1, np.max(points[:,1])+1])
ax.set_zlim([np.min(points[:,2])-1, np.max(points[:,2])+1])

# Add labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Positions with Direction Arrows')

# Show the plot
plt.show()