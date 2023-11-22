import numpy as np
import matplotlib.pyplot as plt

class B_spline3D:
    def __init__(self, dim_u, dim_v):
        self.dim_u = dim_u
        self.dim_v = dim_v

    def set_control_points(self, control_points, close=True):
        self.control_size_u = control_points.shape[0]
        self.control_size_v = control_points.shape[1]
        dim_u = self.dim_u
        dim_v = self.dim_v
        self.knot_size_u = self.control_size_u+dim_u-1
        self.knot_size_v = self.control_size_v+dim_v-1
        #端を閉じているか
        close = True

        if close:
            self.control_vec_uv = control_points
            self.knot_vec_u = np.linspace(0.0, 1.0, self.knot_size_u)
            self.knot_vec_v = np.linspace(0.0, 1.0, self.knot_size_v)

            r_add = self.control_vec_uv[:, -1].reshape(-1, 1)
            l_add = self.control_vec_uv[:, 0].reshape(-1, 1)
            for _ in range(dim_u-1):
                self.control_vec_uv = np.hstack([
                    l_add,
                    self.control_vec_uv
                    ])
                self.control_vec_uv = np.hstack([
                    self.control_vec_uv,
                    r_add
                    ])
                
            u_add = self.control_vec_uv[0, :].reshape(1, -1)
            b_add = self.control_vec_uv[-1, :].reshape(1, -1)
            for _ in range(dim_u-1):
                self.control_vec_uv = np.vstack([
                    u_add,
                    self.control_vec_uv
                    ])
                self.control_vec_uv = np.vstack([
                    self.control_vec_uv,
                    b_add
                    ])

            self.knot_vec_u = np.append([self.knot_vec_u[0]]*dim_u, self.knot_vec_u)
            self.knot_vec_u = np.append(self.knot_vec_u, [self.knot_vec_u[-1]]*dim_u)
            self.knot_vec_v = np.append([self.knot_vec_v[0]]*dim_v, self.knot_vec_v)
            self.knot_vec_v = np.append(self.knot_vec_v, [self.knot_vec_v[-1]]*dim_v)
        else:
            self.control_vec = control_points
            self.knot_vec = np.linspace(0.0, 1.0, self.knot_size)

        self.control_size_u = self.control_vec_uv.shape[0]
        self.control_size_v = self.control_vec_uv.shape[1]
        self.knot_size_u = self.knot_vec_u.size
        self.knot_size_v = self.knot_vec_v.size
        """
        print(self.control_vec_uv)
        
        print (self.control_size_u)
        print (self.control_size_v)
        print (self.knot_size_u)
        print (self.knot_size_v)
        """
    def k_index_u(self, k):
        if k < self.knot_vec_u[0] or self.knot_vec_u[-1] < k:
            return None
        
        for idx, knot in enumerate(self.knot_vec_u):
            if k < knot:
                return idx-1
        if self.knot_vec_u[-1] == k:
            return self.knot_size_u-1
        
    def k_index_v(self, k):
        if k < self.knot_vec_v[0] or self.knot_vec_v[-1] < k:
            return None
        
        for idx, knot in enumerate(self.knot_vec_v):
            if k < knot:
                return idx-1
        if self.knot_vec_v[-1] == k:
            return self.knot_size_v-1

    def base_func_u(self, j, t, k):
        t_idx = self.k_index_u(t)
        if k == 0:
            if j <= t_idx and t_idx < j+1:
                return 1.0
            else:
                return 0.0
        else:
            t_j   = self.knot_vec_u[j]
            t_j1  = self.knot_vec_u[j+1]
            t_jk  = self.knot_vec_u[j+k]
            t_jk1 = self.knot_vec_u[j+k+1]

            if (t_jk  - t_j ) != 0.0:
                left  = (t     - t_j)/(t_jk  - t_j )*self.base_func_u(j,   t, k-1)
            else:
                left = 0.0
            if (t_jk1 - t_j1) != 0.0:
                right = (t_jk1 -   t)/(t_jk1 - t_j1)*self.base_func_u(j+1, t, k-1)
            else:
                right = 0.0

            return left + right

    def base_func_v(self, j, t, k):
        t_idx = self.k_index_v(t)
        if k == 0:
            if j <= t_idx and t_idx < j+1:
                return 1.0
            else:
                return 0.0
        else:
            t_j   = self.knot_vec_v[j]
            t_j1  = self.knot_vec_v[j+1]
            t_jk  = self.knot_vec_v[j+k]
            t_jk1 = self.knot_vec_v[j+k+1]

            if (t_jk  - t_j ) != 0.0:
                left  = (t     - t_j)/(t_jk  - t_j )*self.base_func_v(j,   t, k-1)
            else:
                left = 0.0
            if (t_jk1 - t_j1) != 0.0:
                right = (t_jk1 -   t)/(t_jk1 - t_j1)*self.base_func_v(j+1, t, k-1)
            else:
                right = 0.0

            return left + right

    def B_spline(self, t_u, t_v):
        curve = 0.0
        for u in range(self.control_size_u):
            for v in range(self.control_size_v):
                if t_u == 1.0 and u==self.control_size_u-1:
                    weight_u = 1.0
                else:
                    weight_u = self.base_func_u(u, t_u, self.dim_u)
                if t_v == 1.0 and v==self.control_size_v-1:
                    weight_v = 1.0
                else:
                    weight_v = self.base_func_v(v, t_v, self.dim_v)
                cp_uv = self.control_vec_uv[u][v]
                curve += weight_u * weight_v * cp_uv
                #print("u: " + str(u) + ", v: " + str(v) + ", t_u: " + str(t_u) + ", t_v: " + str(t_v)+ ", weight_u: " + str(weight_u) + ", weight_v: " + str(weight_v))
        #print()
        return curve
    
    def calc_face(self, t_x, t_y):
        z = np.zeros(t_x.shape)
        for i, (x_list, y_list) in enumerate(zip(t_x, t_y)):
            z[i] = [self.B_spline(x, y) for x, y in zip(x_list, y_list)]
        return z
    
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