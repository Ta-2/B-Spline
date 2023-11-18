from B_spline3D import B_spline3D
import numpy as np
import matplotlib.pyplot as plt

class B_spline_surface3D:
    pass

if __name__ == "__main__":
    control_points_dim = 2
    control_points_num = 12
    control_points = np.random.rand(control_points_num, control_points_num)
    print(control_points)
    bs = B_spline3D(control_points_dim, control_points_dim)
    bs.set_control_points(control_points, close=True)

    curve_num = 50
    t = np.linspace(0.0, 1.0, curve_num)
    x, y = np.meshgrid(t, t)
    z = bs.calc_face(x, y)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_wireframe(x, y, z)
    #ax.plot_surface(x, y, z)
    plt.show()