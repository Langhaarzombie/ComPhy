import numpy as np
import matplotlib.pyplot as plt

A0 = 1.0
B0 = 0.0
C0 = 0.0
k1 = 1
k2 = 1

def euler_forward(h, stop=10):
    # u_i+1 = u_i + h f(t_i, u_i)
    xs = np.arange(stop, step=h)
    a_ys = np.array([A0])
    b_ys = np.array([B0])
    c_ys = np.array([C0])
    for x in xs:
        new_a = a_ys[-1] + h * (-k1) * A0 * np.exp(-k1 * x)
        new_b = b_ys[-1] + h * (-k2) * (b_ys[-1] + new_a) * np.exp(-k2 * h)
        new_c = A0 - new_a - new_b

        a_ys = np.append(a_ys, new_a)
        b_ys = np.append(b_ys, np.abs(new_b))
        c_ys = np.append(c_ys, new_c)

    return np.append(xs, [stop]), a_ys, b_ys, c_ys

def euler_backward(h, stop=10):
    # u_i+1 = u_i + h f(t_i+1, u_i+1)
    xs = np.arange(stop, step=h)
    a_ys = np.array([A0])
    b_ys = np.array([B0])
    c_ys = np.array([C0])
    for x in xs:
        new_a = a_ys[-1] + h * (-k1) * A0 * np.exp(-k1 * (x + h))
        new_b = b_ys[-1] + h * (-k2) * (b_ys[-1] + new_a) * np.exp(-k2 * 2 * h)
        new_c = A0 - new_a - new_b

        a_ys = np.append(a_ys, new_a)
        b_ys = np.append(b_ys, np.abs(new_b))
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
        new_a = a_ys[-1] + h/2 * ( a_ys[-1] + (-k1) * A0 * np.exp(-k1 * (x + h)))
        new_b = b_ys[-1] + h/2 * ( b_ys[-1] + (-k2) * (b_ys[-1] + new_a) * np.exp(-k2 * 2 * h))
        new_c = A0 - new_a - new_b

        a_ys = np.append(a_ys, new_a)
        b_ys = np.append(b_ys, np.abs(new_b))
        c_ys = np.append(c_ys, new_c)

    return np.append(xs, [stop]), a_ys, b_ys, c_ys

if __name__ == "__main__":
    fig, ax = plt.subplots(3)

    ef_xs, ef_a_ys, ef_b_ys, ef_c_ys = euler_forward(0.01)
    # ax[0].plot(ef_xs, ef_a_ys)
    # ax[1].plot(ef_xs, ef_b_ys)
    # ax[2].plot(ef_xs, ef_c_ys)

    eb_xs, eb_a_ys, eb_b_ys, eb_c_ys = euler_backward(0.01)
    # ax[0].plot(eb_xs, eb_a_ys)
    # ax[1].plot(eb_xs, eb_b_ys)
    # ax[2].plot(eb_xs, eb_c_ys)

    he_xs, he_a_ys, he_b_ys, he_c_ys = heun(0.01)
    ax[0].plot(he_xs, he_a_ys)
    ax[1].plot(he_xs, he_b_ys)
    ax[2].plot(he_xs, he_c_ys)

    plt.show()
