import numpy as np
import matplotlib.pyplot as plt
import time

class B_spline_ad:
    def __init__(self, dim):
        self.dim = dim

    def generate_dim_fitted_knot(self, cp_num=2):
        knot_vec = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]).astype("float64")/4.0
        knot_vec = []
        knot_num = cp_num+self.dim+1
        if (knot_num-2*(self.dim+1)) % self.dim != 0:
            print("error")
        middle_num = (knot_num-2*(self.dim+1)) // self.dim
        num = 0
        knot_vec += [num]*(self.dim+1)
        for _ in range(middle_num):
            num += 1
            knot_vec += [num]*self.dim
        num += 1
        knot_vec += [num]*(self.dim+1)
        knot_vec = np.array(knot_vec)/knot_vec[-1]

        return knot_vec

    def set_control_points(self, control_points, knot_vec=None, weight_vec=None, close=True, gradient=True):
        self.close = close
        self.control_size = control_points.size
        dim = self.dim

        #両端を閉じないときの要素数チェック
        if type(knot_vec) is np.ndarray and not close:
            if knot_vec.size != self.control_size+dim+1:
                print("error")
            else:
                close = False

        #両端を閉じるときの要素数チェック
        if type(knot_vec) is np.ndarray and close:
            if knot_vec.size != self.control_size+dim-1:
                print("error")

        #重みベクトルがある時の要素数チェック
        if weight_vec is np.ndarray:
            if weight_vec.size != control_points.size:
                print("error")
            else:
                close = False
                
        #両端を閉じさせる
        if close:
            self.knot_size = self.control_size+dim-1
            self.control_vec = control_points

            #knotベクトルがない時
            if type(knot_vec) is not np.ndarray:
                self.knot_vec = np.linspace(0.0, 1.0, self.knot_size)
            else:
                self.knot_vec = knot_vec
            
            cp_front = self.control_vec[0]
            cp_back  = self.control_vec[-1]

            #勾配を再現する
            if gradient:
                self.control_vec[0]  = (self.control_vec[1]/3.0  + cp_front*2.0/3.0)
                self.control_vec[-1] = (self.control_vec[-2]/3.0 + cp_back*2.0/3.0 )
            #制御点を閉じる
            self.control_vec = np.append([cp_front]*(dim-1), self.control_vec)
            self.control_vec = np.append(self.control_vec, [cp_back]*(dim-1))

            #knotベクトルを閉じる
            knot_front = self.knot_vec[0]
            knot_back  = self.knot_vec[-1]
            self.knot_vec = np.append([knot_front]*dim, self.knot_vec)
            self.knot_vec = np.append(self.knot_vec, [knot_back]*dim)

        #そのまま使う
        else:
            self.control_vec = control_points
            if type(knot_vec) is not np.ndarray:
                self.knot_vec = np.linspace(0.0, 1.0, self.knot_size)
            else:
                self.knot_vec = knot_vec

        self.control_size = self.control_vec.size
        self.knot_size = self.knot_vec.size

        print(self.knot_vec)
        print(self.control_vec)
        #print("cp size: " + str(self.control_size) + ", knot size: " + str(self.knot_size) + ", dim: " + str(self.dim))

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
    
    def k_index_ad(self, k, close=True):
        knot_size = self.knot_vec.size
        if k < self.knot_vec[0] or self.knot_vec[-1] < k:
            return None
        
        L, R = self.dim, knot_size-1-self.dim

        while( 1 < R-L ):
            C = (L+R)//2
            pivot = self.knot_vec[C]
            if pivot == k:
                return C
            elif pivot < k:
                L = C
            else:
                R = C

        return (L+R)//2

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
        #t_idx = self.k_index(t)
        t_idx = self.k_index_ad(t)
        curve = 0.0
        for i in range(self.dim+1):
            cp_val = self.control_vec[t_idx-(self.dim-i)]
            #print("t_idx: " + str(t_idx) + ", cp_idx: " + str(t_idx-i))
            curve += cp_val*self.calc_weight_of_base(t, t_idx, i, self.dim-i)

        return curve
    
    def base(self, t, base_num, cp=False):
        t_idx = self.k_index_ad(t)
        curve = 0.0
        for i in range(self.dim+1):
            if base_num == t_idx-(self.dim-i):
                if cp:
                    cp_val = self.control_vec[t_idx-(self.dim-i)]
                else:
                    cp_val = 1.0
                #print("t_idx: " + str(t_idx) + ", cp_idx: " + str(t_idx-i))
                curve += cp_val*self.calc_weight_of_base(t, t_idx, i, self.dim-i)

        return curve

if __name__ == "__main__":
    control_points_dim = 2
    control_points_num = 9
    bs = B_spline_ad(control_points_dim)
    close = False
    grad = False
    control_points = np.random.randint(2, 5, (control_points_num))
    control_points_not_close = np.random.randint(2, 5, (control_points_num-2*(control_points_dim-1)))
    knot_vec = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]).astype("float64")/4.0
    knot_vec_not_close = np.array([0, 1, 1, 2, 2, 3, 3, 4]).astype("float64")/4.0

    #case 1
    print("case 1")
    bs.set_control_points(control_points, close=True, gradient=False)
    #case 2
    print("case 2")
    bs.set_control_points(control_points, close=True, gradient=True)
    #case 3
    print("case 3")
    bs.set_control_points(control_points, knot_vec=knot_vec, close=False, gradient=False)
    #case 4
    print("case 4")
    bs.set_control_points(control_points_not_close, knot_vec=knot_vec_not_close, close=True, gradient=False)
    #case 5
    print("case 5")
    control_points = np.array([3, 3, 2, 1, 1, 1, 2, 3, 3]).astype("float64")
    weight_vec     = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1]).astype("float64")
    knot_vec = bs.generate_dim_fitted_knot(cp_num=control_points.size)
    bs.set_control_points(control_points, knot_vec=knot_vec, close=close, gradient=grad)

    curve_num = 200
    x = np.linspace(0.0, 1.0, curve_num)
    start = time.time()
    y = [bs.B_spline_ad(i) for i in x]
    end = time.time()
    print("time: " + str(end-start))
    
    plt.plot(x, y)
    
    bases = []
    for i in range(control_points_num + 2*(control_points_dim-1)):
        base = [bs.base(t, i, cp=True) for t in x]
        bases.append(base)
    for b in bases:
        plt.plot(x, b)
        
    if close:
        sc_x = np.linspace(0.0, 1.0, bs.knot_vec[bs.dim : -bs.dim].size)
        sc_x += (sc_x[1]-sc_x[0])/2
        sc_x = sc_x[:-1]
        sc_x = np.append([0.0]*(bs.dim-1), np.append(sc_x, [1.0]*(bs.dim-1)))
        sc_y = bs.control_vec[bs.dim-1 : -(bs.dim-1)]
    else:
        sc_y = bs.control_vec
        sc_x = np.linspace(0.0, 1.0, sc_y.size)

    #print(sc_x)
    #print(sc_y)
    plt.scatter(sc_x, sc_y)
    ax = plt.gcf().gca()
    ax.grid()
    ax.set_xticks(sc_x)
    ax.set_xlim(0.0, 1.0)
    ax.set_yticks(np.append(bs.control_vec, [0.5, 1.0]))
    
    plt.show()