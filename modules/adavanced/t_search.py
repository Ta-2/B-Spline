import numpy as np
import time

def k_index(knot_vec, k, close=True, dim=2):
        knot_size = knot_vec.size
        if k < knot_vec[0] or knot_vec[-1] < k:
            return None
        
        for idx, knot in enumerate(knot_vec):
            if k < knot:
                return idx-1
        if knot_vec[-1] == k:
            if close:
                return knot_size-1-(dim+1)
            else:
                return knot_size-1

def k_index_ad(knot_vec, k, close=True, dim=2, comm = False):
    knot_size = knot_vec.size
    if k < knot_vec[0] or knot_vec[-1] < k:
        return None
    
    L, R = dim, knot_size-1-dim

    while( 1 < R-L ):
        C = (L+R)//2
        pivot = knot_vec[C]
        if comm:
            print("L: " + str(L) + ", C: " + str(C) + ", R: " + str(R) + ", pivot: " + str(pivot))
        if pivot == k:
            return C
        elif pivot < k:
            L = C
        else:
            R = C

    if comm:
        print("just before return L: " + str(L) + ", C: " + str(C) + ", R: " + str(R) + ", pivot: " + str(pivot))
    return (L+R)//2

def accuracy_test():
    dim = 2
    for num in np.arange(5, 3000):
        print("accuracy test: dim:" + str(dim))
        knot_vec = np.linspace(0.0, 1.0, num)
        knot_vec = np.append([knot_vec[0]]*dim, knot_vec)
        knot_vec = np.append(knot_vec, [knot_vec[-1]]*dim)
        if num % 20 == 0:
            print("processing num: " + str(num))

        test_num = 7000
        for t in np.linspace(0.0, 1.0, test_num):
            old_var = k_index(knot_vec, t)
            advenced_var = k_index_ad(knot_vec, t)

            t_prev = knot_vec[advenced_var]
            t_next = knot_vec[advenced_var+1]
            if(old_var != advenced_var and (t_prev <= t and t < t_next)):
                print("dim: " + str(dim) +  ", k_num: " + str(num) +  ", knot_vec: ", knot_vec)
                print("t: " + str(t) + ", old: " + str(old_var) + ", ad: " + str(advenced_var))
                advenced_var = k_index_ad(knot_vec, t, dim=dim, comm=True)

def speed_test(dim, knot_num, test_num):
    print("speed_test: dim:" + str(dim) + ", knot_num: " + str(knot_num) + ", test_num: " + str(test_num))

    knot_vec = np.linspace(0.0, 1.0, knot_num)
    knot_vec = np.append([knot_vec[0]]*dim, knot_vec)
    knot_vec = np.append(knot_vec, [knot_vec[-1]]*dim)

    start = time.time()
    for t in np.linspace(0.0, 1.0, test_num):
        k_index(knot_vec, t)
    end = time.time()
    print("time old: " + str(end-start))

    
    start = time.time()
    for t in np.linspace(0.0, 1.0, test_num):
        k_index_ad(knot_vec, t)
    end = time.time()
    print("time adv: " + str(end-start))

if __name__ == "__main__":
    speed_test(2, 2000, 5000)
    speed_test(2, 5000, 10000)