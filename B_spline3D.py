import B_spline_module as BS
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class B_spline3D:
    def __init__(self, dim_u, dim_v=None):
        if type(dim_v) != int:
            dim_v = dim_u
        print(str(dim_u) + ", " + str(dim_v))
        self.B_sp_u = BS.B_spline(dim_u)
        self.B_sp_v = BS.B_spline(dim_v)
    
    def set_control_points(self, cp_uv, close):
        self.cp_uv = cp_uv
        print(cp_uv.shape)
        self.cp_u, self.cp_v = cp_uv.shape
        self.B_sp_u.set_control_points(np.ones(self.cp_u), close=close)
        self.B_sp_v.set_control_points(np.ones(self.cp_v), close=close)
        
    def calc_curve_uv(self, u_list, v_list):
        i_u = u_list.shape[0]
        for i in range(i_u):
            u_list[i] = self.B_sp_u.calc_curve(u_list[i])

        i_v = v_list.shape[0]
        for i in range(i_v):
            v_list[i] = self.B_sp_v.calc_curve(v_list[i])
        return u_list, v_list

    def calc_face(self, t_u, t_v=None):
        if type(t_v) == type(None):
            t_v = np.array(t_u)
        u, v = np.meshgrid(t_u, t_v)
        w_x, w_y = self.calc_curve_uv(u, v)
        w_xy = (w_x * w_y)
        return w_xy

if __name__ == "__main__":
    control_points_dim = 2
    control_points_num = 5
    #control_points_x = np.random.rand(control_points_num)
    #control_points_y = np.random.rand(control_points_num)
    control_points_xy = np.array(
        [
            np.random.rand(control_points_num),
            np.random.rand(control_points_num),
            np.random.rand(control_points_num),
            np.random.rand(control_points_num),
            np.random.rand(control_points_num)
        ]
    )
    
    bs2 = B_spline3D(control_points_dim)
    bs2.set_control_points(control_points_xy, True)
    
    curve_num = 9
    t = np.linspace(0.0, 1.0, curve_num)
    x, y = np.meshgrid(t, t)
    z = bs2.calc_face(t)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_wireframe(x, y, z)
    #ax.set_xticks(control_points_x)
    #ax.set_xlim(-0.1, 1.1)
    #ax.set_yticks(control_points_y)
    #ax.set_xlim(-0.1, 1.1)

    #ax.grid()
    #plt.plot(x, y)
    #plt.scatter(control_points_x, control_points_y)
    plt.show()