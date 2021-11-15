import numpy as np
import matplotlib.pyplot as plt

A0 = 1.0
B0 = 0.0
C0 = 0.0

def euler_expl(h, k2, k1=1, stop=10):
    """
    Solves the decay ode using explicit euler method:
        u_i+1 = u_i + h f(t_i, u_i)
    h    ... timestep
    k2   ... decay ratio for B
    k1   ... decay ratio for A
    stop ... time to which to simulate

    Returns:
    xs   ... timesteps of simulation
    a_ys ... amount of A over time
    b_ys ... amount of B over time
    c_ys ... amount of C over time

    Decay is happening from A - B - C.

    """
    xs = np.arange(stop, step=h)
    a_ys = np.array([A0])
    b_ys = np.array([B0])
    c_ys = np.array([C0])
    for _ in xs:
        new_a = a_ys[-1] + h * (-k1) * a_ys[-1]
        new_b = b_ys[-1] + h * (k1 * a_ys[-1] - k2 * b_ys[-1])
        new_c = A0 - new_a - new_b

        a_ys = np.append(a_ys, new_a)
        b_ys = np.append(b_ys, new_b)
        c_ys = np.append(c_ys, new_c)

    return np.append(xs, [stop]), a_ys, b_ys, c_ys

def euler_impl(h, k2, k1=1, stop=10):
    """
    Solves the decay ode using implicit euler method:
        u_i+1 = u_i + h f(t_i+1, u_i+1)
    h    ... timestep
    k2   ... decay ratio for B
    k1   ... decay ratio for A
    stop ... time to which to simulate

    Returns:
    xs   ... timesteps of simulation
    a_ys ... amount of A over time
    b_ys ... amount of B over time
    c_ys ... amount of C over time

    Decay is happening from A - B - C.
    The equation to get an explicit form for u_i+1 was transformed
    on paper.
    """
    xs = np.arange(stop, step=h)
    a_ys = np.array([A0])
    b_ys = np.array([B0])
    c_ys = np.array([C0])
    for _ in xs:
        new_a = a_ys[-1] + h * (-k1) * a_ys[-1]
        new_b = (b_ys[-1] + h * k1 * new_a) / (1 + h * k2)
        new_c = A0 - new_a - new_b

        a_ys = np.append(a_ys, new_a)
        b_ys = np.append(b_ys, new_b)
        c_ys = np.append(c_ys, new_c)

    return np.append(xs, [stop]), a_ys, b_ys, c_ys

def heun(h, k2, k1=1, stop=10):
    """
    Solves the decay ode using implicit euler method:
        u_i+1 = u_i + h/2 (u_i + f(t_i+1, u_ip))
        u_ip = u_i + h f(t_i, u_i)
    h    ... timestep
    k2   ... decay ratio for B
    k1   ... decay ratio for A
    stop ... time to which to simulate

    Returns:
    xs   ... timesteps of simulation
    a_ys ... amount of A over time
    b_ys ... amount of B over time
    c_ys ... amount of C over time

    Decay is happening from A - B - C.
    """
    xs = np.arange(stop, step=h)
    a_ys = np.array([A0])
    b_ys = np.array([B0])
    c_ys = np.array([C0])
    for _ in xs:
        a_ip = a_ys[-1] + h * (-k1) * a_ys[-1]
        new_a = a_ys[-1] + h/2 * (-k1) * (a_ys[-1] + a_ip)
        
        b_ip = b_ys[-1] + h * (k1 * a_ys[-1] - k2 * b_ys[-1])
        new_b = b_ys[-1] + h/2 * (k1 * a_ys[-1] - k2 * b_ys[-1] + k1 * new_a - k2 * b_ip)

        new_c = A0 - new_a - new_b

        a_ys = np.append(a_ys, new_a)
        b_ys = np.append(b_ys, new_b)
        c_ys = np.append(c_ys, new_c)

    return np.append(xs, [stop]), a_ys, b_ys, c_ys

def C(t, k2, k1=1):
    """
    Calculates the analytic solution for the amount of C over time.

    t  ... timesteps
    k2 ... decay ratio for B
    k1 ... decay ratio for A

    Returns:
    amount of C over time
    """
    if k2 == k1:
        return A0 * (1 - np.exp(-k1 * t) * (1 + k1 * t))
    return A0 * (1 - (k2 * np.exp(-k1 * t) - k1 * np.exp(-k2 * t)) / (k2 - k1))

def analytical(h, k2, k1=1, stop=10):
    xs = np.arange(stop+h, step=h)
    return xs, C(xs, k2, k1)

if __name__ == "__main__":
    fig, ax = plt.subplots(5)

    # ee_xs, ee_a_ys, ee_b_ys, ee_c_ys = euler_expl(0.1, 1)
    # ax[0].plot(ee_xs, ee_a_ys)
    # ax[1].plot(ee_xs, ee_b_ys)
    # ax[2].plot(ee_xs, ee_c_ys)

    # ei_xs, ei_a_ys, ei_b_ys, ei_c_ys = euler_impl(0.01, 1)
    # ax[0].plot(ei_xs, ei_a_ys)
    # ax[1].plot(ei_xs, ei_b_ys)
    # ax[2].plot(ei_xs, ei_c_ys)

    for k2 in [1000, 100, 10, 1]:
        he_xs, he_a_ys, he_b_ys, he_c_ys = heun(0.0005, k2)
        ax[0].plot(he_xs, he_a_ys, label=f"k={k2}")
        ax[1].plot(he_xs, he_b_ys, label=f"k={k2}")
        ax[2].plot(he_xs, he_c_ys, label=f"k={k2}")

        ana_xs, ana_c_ys = analytical(0.01, k2)
        ax[3].plot(ana_xs, ana_c_ys, label=f"Analytic k={k2}")

    ax[4].plot(ana_xs, ana_c_ys - he_c_ys[::20], label=f"Difference k=1")

    ax[0].set_xlabel("Time [s]")
    ax[1].set_xlabel("Time [s]")
    ax[2].set_xlabel("Time [s]")
    ax[3].set_xlabel("Time [s]")
    ax[4].set_xlabel("Time [s]")
    ax[0].set_ylabel("A(t)")
    ax[1].set_ylabel("B(t)")
    ax[2].set_ylabel("C(t)")
    ax[3].set_ylabel("C_ana(t)")
    ax[4].set_ylabel("C_ana(t) - C(t)")

    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    ax[3].legend()
    ax[4].legend()

    plt.show()
