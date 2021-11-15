import numpy as np
import matplotlib.pyplot as plt

A0 = 1.0
B0 = 0.0
C0 = 0.0
k1 = 1
k2 = 1

def euler_expl(h, stop=10):
    # u_i+1 = u_i + h f(t_i, u_i)
    xs = np.arange(stop, step=h)
    a_ys = np.array([A0])
    b_ys = np.array([B0])
    c_ys = np.array([C0])
    for x in xs:
        new_a = a_ys[-1] + h * (-k1) * a_ys[-1]
        new_b = b_ys[-1] + h * (k1 * a_ys[-1] - k2 * b_ys[-1])
        new_c = A0 - new_a - new_b

        a_ys = np.append(a_ys, new_a)
        b_ys = np.append(b_ys, new_b)
        c_ys = np.append(c_ys, new_c)

    return np.append(xs, [stop]), a_ys, b_ys, c_ys

def euler_impl(h, stop=10):
    # u_i+1 = u_i + h f(t_i+1, u_i+1)
    xs = np.arange(stop, step=h)
    a_ys = np.array([A0])
    b_ys = np.array([B0])
    c_ys = np.array([C0])
    for x in xs:
        new_a = a_ys[-1] + h * (-k1) * a_ys[-1]
        new_b = (b_ys[-1] + h * k1 * new_a) / (1 + h * k2)
        new_c = A0 - new_a - new_b

        a_ys = np.append(a_ys, new_a)
        b_ys = np.append(b_ys, new_b)
        c_ys = np.append(c_ys, new_c)

    return np.append(xs, [stop]), a_ys, b_ys, c_ys

def heun(h, stop=10):
    # u_i+1 = u_i + h/2 (u_i + f(t_i+1, u_ip))
    # u_ip = u_i + h f(t_i, u_i)
    xs = np.arange(stop, step=h)
    a_ys = np.array([A0])
    b_ys = np.array([B0])
    c_ys = np.array([C0])
    for x in xs:
        a_ip = a_ys[-1] + h * (-k1) * a_ys[-1]
        new_a = a_ys[-1] + h/2 * (-k1) * (a_ys[-1] + a_ip)
        
        b_ip = b_ys[-1] + h * (k1 * a_ys[-1] - k2 * b_ys[-1])
        new_b = b_ys[-1] + h/2 * (k1 * a_ys[-1] - k2 * b_ys[-1] + k1 * new_a - k2 * b_ip)

        new_c = A0 - new_a - new_b

        a_ys = np.append(a_ys, new_a)
        b_ys = np.append(b_ys, new_b)
        c_ys = np.append(c_ys, new_c)

    return np.append(xs, [stop]), a_ys, b_ys, c_ys

def C(t):
    if k2 == k1:
        return A0 * (1 - np.exp(-k1 * t) * (1 + k1 * t))
    return A0 * (1 - (k1 * np.exp(-k2 * t) - k2 * np.exp(-k1 * t)) / (k2 - k1))

def analytical(h, stop=10):
    xs = np.arange(stop, step=h)
    return xs, C(xs)

if __name__ == "__main__":
    fig, ax = plt.subplots(3)

    # ee_xs, ee_a_ys, ee_b_ys, ee_c_ys = euler_expl(0.1)
    # ax[0].plot(ee_xs, ee_a_ys)
    # ax[1].plot(ee_xs, ee_b_ys)
    # ax[2].plot(ee_xs, ee_c_ys)

    # ei_xs, ei_a_ys, ei_b_ys, ei_c_ys = euler_impl(0.01)
    # ax[0].plot(ei_xs, ei_a_ys)
    # ax[1].plot(ei_xs, ei_b_ys)
    # ax[2].plot(ei_xs, ei_c_ys)

    he_xs, he_a_ys, he_b_ys, he_c_ys = heun(0.001)
    ax[0].plot(he_xs, he_a_ys)
    ax[1].plot(he_xs, he_b_ys)
    ax[2].plot(he_xs, he_c_ys)

    ana_xs, ana_c_ys = analytical(0.001)
    # ax[2].plot(ana_xs, ana_c_ys, label="Analytic solution")

    plt.show()
