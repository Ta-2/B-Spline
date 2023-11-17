import B_spline_module as BS
import numpy as np
import matplotlib.pyplot as plt

class B_spline2D:
    def __init__(self, dim):
        self.B_sp_x = BS.B_spline(dim)
        self.B_sp_y = BS.B_spline(dim)
    
    def set_control_points(self, cp_x, cp_y, close):
        self.B_sp_x.set_control_points(cp_x, close=close)
        self.B_sp_y.set_control_points(cp_y, close=close)

    def calc_curve(self, t):
        x = self.B_sp_x.calc_curve(t)
        y = self.B_sp_y.calc_curve(t)
        return x, y

if __name__ == "__main__":
    control_points_dim = 2
    control_points_num = 8
    control_points_x = np.random.rand(control_points_num)
    control_points_y = np.random.rand(control_points_num)
    
    bs2 = B_spline2D(control_points_dim)
    bs2.set_control_points(control_points_x, control_points_y, True)
    
    curve_num = 200
    t = np.linspace(0.0, 1.0, curve_num)
    x, y = bs2.calc_curve(t)

    ax = plt.gcf().gca()
    ax.set_xticks(control_points_x)
    ax.set_xlim(-0.1, 1.1)
    ax.set_yticks(control_points_y)
    ax.set_xlim(-0.1, 1.1)

    ax.grid()
    plt.plot(x, y)
    plt.scatter(control_points_x, control_points_y)
    plt.show()