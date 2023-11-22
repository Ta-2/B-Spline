import numpy as np
import matplotlib.pyplot as plt
import time

class B_spline_ad:
    def __init__(self, dim):
        self.dim = dim

    def set_control_points(self, control_points, knot_vector=None, close=True):
        self.close = close
        self.control_size = control_points.size
        dim = self.dim
        if type(knot_vector) is np.ndarray:
            if knot_vector.size != self.control_size+dim-1:
                print("error")
        self.knot_size = self.control_size+dim-1

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
            if self.close:
                return self.knot_size-1-(self.dim+1)
            else:
                return self.knot_size-1
    
    def calc_weight_of_base(self, t, t_idx, Up_cnt, Down_cnt, L=0, R=0, space=""):
        if Up_cnt==0 and Down_cnt==0:
            return 1.0
        up_weight, down_weight = 0.0, 0.0
        t_l = self.knot_vec[t_idx    -L]
        t_r = self.knot_vec[t_idx +1 +R]
        #space += "  "
        #print(space + "t_ixd: " + str(t_idx) + ", up_cnt: " + str(Up_cnt) + ", down_cnt: " + str(Down_cnt) + ", L: " + str(L) + ", R: " + str(R))

        t_dif = t_r-t_l
        if t_dif == 0.0:
            return 0.0
        else:
            if Up_cnt   > 0:
                up_weight   = self.calc_weight_of_base(t, t_idx, Up_cnt-1, Down_cnt, L, R+1, space)
                up_weight   *= (t - t_l)/t_dif
            if Down_cnt > 0:
                down_weight = self.calc_weight_of_base(t, t_idx, Up_cnt, Down_cnt-1, L+1, R, space)
                down_weight *= (t_r -t)/t_dif
        
        return up_weight + down_weight
    
    def B_spline_ad(self, t):
        t_idx = self.k_index(t)
        curve = 0.0
        for i in range(self.dim+1):
            cp_val = self.control_vec[t_idx-(self.dim-i)]
            #print("t_idx: " + str(t_idx) + ", cp_idx: " + str(t_idx-i))
            curve += cp_val*self.calc_weight_of_base(t, t_idx, i, self.dim-i)

        return curve

if __name__ == "__main__":
    control_points_dim = 2
    control_points_num = 500
    control_points = np.arange(control_points_num)
    #control_points = np.random.randint(2, 5, (control_points_num))
    bs = B_spline_ad(control_points_dim)
    bs.set_control_points(control_points, close=True)

    curve_num = 2000
    x = np.linspace(0.0, 1.0, curve_num)
    start = time.time()
    y = [bs.B_spline_ad(i) for i in x]
    end = time.time()
    print("time: " + str(end-start))

    plt.plot(x, y)
    """
    sc_x = np.linspace(0.0, 1.0, bs.control_size - 2*(bs.dim-1))
    sc_y = bs.control_vec[bs.dim-1 : -bs.dim+1]
    print(sc_x)
    print(sc_y)
    plt.scatter(sc_x, sc_y)
    
    ax = plt.gcf().gca()
    ax.grid()
    ax.set_xticks(sc_x)
    ax.set_xlim(0.0, 1.0)
    ax.set_yticks(np.append(bs.control_vec, [0.5, 1.0]))
    """
    plt.show()