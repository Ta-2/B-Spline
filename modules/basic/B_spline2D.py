import numpy as np
import matplotlib.pyplot as plt
import time

class B_spline:
    def __init__(self, dim):
        self.dim = dim

    def set_control_points(self, control_points, knot_vector=None, close=True):
        self.control_size = control_points.size
        dim = self.dim
        if type(knot_vector) is np.ndarray:
            if knot_vector.size != self.control_size+dim-1:
                print("error")
        self.knot_size = self.control_size+dim-1
        #端を閉じているか
        close = True

        if close:
            self.control_vec = control_points
            if type(knot_vector) is not np.ndarray:
                self.knot_vec = np.linspace(0.0, 1.0, self.knot_size)
            else:
                self.knot_vec = knot_vector

            self.control_vec = np.append([self.control_vec[0]]*(dim-1), self.control_vec)
            self.control_vec = np.append(self.control_vec, [self.control_vec[-1]]*(dim-1))
            self.knot_vec = np.append([self.knot_vec[0]]*dim, self.knot_vec)
            self.knot_vec = np.append(self.knot_vec, [self.knot_vec[-1]]*dim)
        else:
            self.control_vec = control_points
            if type(knot_vector) is not np.ndarray:
                self.knot_vec = np.linspace(0.0, 1.0, self.knot_size)
            else:
                self.knot_vec = knot_vector

        self.control_size = self.control_vec.size
        self.knot_size = self.knot_vec.size

        print (self.knot_vec)
        print (self.control_vec)

    def k_index(self, k):
        if k < self.knot_vec[0] or self.knot_vec[-1] < k:
            return None
        
        for idx, knot in enumerate(self.knot_vec):
            if k < knot:
                return idx-1
        if self.knot_vec[-1] == k:
            return self.knot_size-1

    def base_func(self, j, t, k):
        t_idx = self.k_index(t)
        if k == 0:
            if j <= t_idx and t_idx < j+1:
                return 1.0
            else:
                return 0.0
        else:
            t_j   = self.knot_vec[j]
            t_j1  = self.knot_vec[j+1]
            t_jk  = self.knot_vec[j+k]
            t_jk1 = self.knot_vec[j+k+1]

            if (t_jk  - t_j ) != 0.0:
                left  = (t     - t_j)/(t_jk  - t_j )*self.base_func(j,   t, k-1)
            else:
                left = 0.0
            if (t_jk1 - t_j1) != 0.0:
                right = (t_jk1 -   t)/(t_jk1 - t_j1)*self.base_func(j+1, t, k-1)
            else:
                right = 0.0

            return left + right

    def B_spline(self, t, base_num=-1):
        curve = 0.0
        if base_num == -1:
            for i in range(self.control_size):
                if t == 1.0 and i==self.control_size-1:
                    curve += self.control_vec[i]
                else:
                    curve += self.control_vec[i]*self.base_func(i, t, self.dim)
        else:
            if t == 1.0 and base_num==self.control_size-1:
                curve = 1.0
            else:
                curve = self.base_func(base_num, t, self.dim)
        return curve
    
    def calc_curve(self, t, base_num=-1):
        return np.array([self.B_spline(i, base_num=base_num) for i in t])
    
if __name__ == "__main__":
    control_points_dim = 2
    control_points_num = 500
    control_points = np.arange(control_points_num)
    #control_points = np.array([  2,    4,   3,    2,   4,    2])
    #knot =           np.array([0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 1.0])
    #control_points = np.random.rand(control_points_num)
    bs = B_spline(control_points_dim)
    #bs.set_control_points(control_points, knot_vector=knot, close=True)
    bs.set_control_points(control_points, close=True)

    curve_num = 2000
    x = np.linspace(0.0, 1.0, curve_num)

    start = time.time()
    y = bs.calc_curve(x)
    end = time.time()
    print("time: " + str(end-start))
    """
    bases = []
    for i in range(bs.control_size):
        bases.append([bs.B_spline(t, base_num=i) for t in x])
    """
    """ax = plt.gcf().gca()
    ax.set_xticks(bs.knot_vec)
    ax.set_xlim(0.0, 1.0)
    ax.set_yticks(bs.control_vec)

    ax.grid()
    """
    plt.plot(x, y)
    #plt.scatter(np.linspace(0.0, 1.0, control_points_num), control_points)
    #plt.scatter(np.linspace(0.0, 1.0, bs.control_size), bs.control_vec)
    """
    for b in bases:
        plt.plot(x, b)
    """
    plt.show()