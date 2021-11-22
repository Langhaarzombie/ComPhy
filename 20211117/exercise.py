import numpy as np
import numpy.linalg as linalg
import scipy.constants as const
import matplotlib.pyplot as plt

me = 5.972e24 # mass of the earth
r0 = 42157 # the initial distance from earth

def euler_expl(dt, stop=15):
    """
    Args:
        dt   ... timestep for expl. Euler method
        stop ... number of secods to calculate (stop / dt is number of y values)
    Returns:
        ts   ... timestamps of calculated values
        rs   ... calculated positions of satellite
    
    Uses expl. Euler method to calculate the position of the satellite.
    """
    # timesteps the satellite will go through
    ts = np.arange(stop, step=dt)
    # initial values for the satellite
    rs = np.array([[r0, 0]])

    # Calculated values so that the satellite stays in orbit
    # for any kind of orbit
    T = 2 * const.pi * np.sqrt(r0**3 / (const.G * me))
    # vn = [0, 2 * const.pi * r0 / T]
    # for perfectly circular orbit
    vn = np.array([0, np.sqrt((const.G * me) / r0)])

    for t in ts:
        r_abs = linalg.norm(rs[-1])
        rn = rs[-1] + dt * vn
        vn = vn + dt * (- const.G * me / r_abs**3 ) * rs[-1]

        rs = np.concatenate((rs, [rn]), axis=0)

        if t/T > 2 and t/T < 3:
            v_abs = linalg.norm(vn)
            vn = vn + (vn/v_abs) * dt * 5000

    return ts, rs

def euler_sympl(dt, stop=15):
    """
    Args:
        dt   ... timestep for expl. Euler method
        stop ... number of secods to calculate (stop / dt is number of y values)
    Returns:
        ts   ... timestamps of calculated values
        rs   ... calculated positions of satellite
    
    Uses semi implicit Euler method to calculate the position of the satellite.
    """
    # timesteps the satellite will go through
    ts = np.arange(stop, step=dt)
    # initial values for the satellite
    rs = np.array([[r0, 0]])

    # Calculated values so that the satellite stays in orbit
    # for any kind of orbit
    T = 2 * const.pi * np.sqrt(r0**3 / (const.G * me))
    # vn = [0, 2 * const.pi * r0 / T]
    # for perfectly circular orbit
    vn = np.array([0, np.sqrt((const.G * me) / r0)])

    for t in ts:
        r_abs = linalg.norm(rs[-1])
        rn = rs[-1] + dt * vn
        vn = vn + dt * (- const.G * me / r_abs**3 ) * rn

        rs = np.concatenate((rs, [rn]), axis=0)

        if t/T > 2 and t/T < 3:
            v_abs = linalg.norm(vn)
            vn = vn + (vn/v_abs) * dt * 5000

    return ts, rs

if __name__ == "__main__":
    _, ee_rs = euler_expl(0.01)
    _, es_rs = euler_sympl(0.01)

    _, ax = plt.subplots(2)

    ax[0].plot(ee_rs[:,0], ee_rs[:,1])
    ax[0].plot([0], [0], "o")
    ax[0].plot(ee_rs[-1,0], ee_rs[-1,1], "o")

    ax[1].plot(es_rs[:,0], es_rs[:,1])
    ax[1].plot([0], [0], "o")
    ax[1].plot(es_rs[-1,0], es_rs[-1,1], "o")

    plt.show()
