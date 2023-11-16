import numpy as np
import matplotlib.pyplot as plt

dim = 2
control_size = 5
control_size += 2*(dim-1)
knot_size = control_size+dim+1
#端を閉じているか
close = True

if close:
    control_size = control_size-2*(dim-1)
    knot_size = knot_size-2*dim
    control_vec = np.linspace(1.0, 1.0, control_size)
    control_vec = np.random.rand(control_size)
    knot_vec = np.linspace(0.0, 1.0, knot_size)

    control_vec = np.append([control_vec[0]]*(dim-1), control_vec)
    control_vec = np.append(control_vec, [control_vec[-1]]*(dim-1))
    knot_vec = np.append([knot_vec[0]]*dim, knot_vec)
    knot_vec = np.append(knot_vec, [knot_vec[-1]]*dim)
else:
    control_vec = np.linspace(1.0, 1.0, control_size)
    knot_vec = np.linspace(0.0, 1.0, knot_size)

control_size = control_vec.size
knot_size = knot_vec.size

print (knot_vec)
print (control_vec)

def k_index(k):
    if k < knot_vec[0] or knot_vec[-1] < k:
        return None
    
    for idx, knot in enumerate(knot_vec):
        if k < knot:
            return idx-1
    if knot_vec[-1] == k:
        return knot_size-1

def base_func(j, t, k, close=False):
    t_idx = k_index(t)
    if k == 0:
        if j <= t_idx and t_idx < j+1:
            return 1.0
        else:
            return 0.0
    else:
        t_j   = knot_vec[j]
        t_j1  = knot_vec[j+1]
        t_jk  = knot_vec[j+k]
        t_jk1 = knot_vec[j+k+1]

        if (t_jk  - t_j ) != 0.0:
            left  = (t     - t_j)/(t_jk  - t_j )*base_func(j,   t, k-1, close)
        else:
            left = 0.0
        if (t_jk1 - t_j1) != 0.0:
            right = (t_jk1 -   t)/(t_jk1 - t_j1)*base_func(j+1, t, k-1, close)
        else:
            right = 0.0

        return left + right

def B_spline(t, base_num=-1, close=False):
    curve = 0.0
    if base_num == -1:
        for i in range(control_size):
            if t == 1.0 and i==control_size-1:
                curve += control_vec[i]
            else:
                curve += control_vec[i]*base_func(i, t, dim, close)
    else:
        if t == 1.0 and base_num==control_size-1:
            curve = 1.0
        else:
            curve = base_func(base_num, t, dim, close)
    return curve

curve_num = 200
x = np.linspace(0.0, 1.0, curve_num)
y = np.array(
    [B_spline(t) for t in x]
    )

bases = []
for i in range(control_size):
    bases.append([B_spline(t, base_num=i, close=close) for t in x])

ax = plt.gcf().gca()
ax.set_xticks(knot_vec)
ax.set_xlim(0.0, 1.0)
ax.set_yticks(control_vec)

ax.grid()
plt.plot(x, y)
for b in bases:
    plt.plot(x, b)
plt.show()