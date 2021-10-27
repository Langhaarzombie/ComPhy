import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

time = np.array([])
m_free = np.array([])


def init(file='20211019/oscillator.dat'):
    global time
    global m_free
    data = np.loadtxt(file)
    time = data.T[0]
    m_free = data.T[1:4].T


def calculate_averages():
    avg = np.ndarray((3, 1))
    length = np.size(m_free)
    avg[0] = np.sum(m_free.T[0]) / length
    avg[1] = np.sum(m_free.T[1]) / length
    avg[2] = np.sum(m_free.T[2]) / length
    print("The average for the free layer:")
    print(avg)
    print("Magintude of that vector:")
    print(np.linalg.norm(avg))


def interpolate():
    f = interp1d(time, m_free.T[0])
    plt.plot(time, f(time), label="Interpolated data")


def fit():
    fit_points = 12500
    fit, _ = curve_fit(sinus, time[-fit_points:], m_free.T[0]
                       [-fit_points:])
    plt.plot(time[-fit_points:], sinus(time[-fit_points:], *fit))


def sinus(x, a, b, c, d):
    return a*np.sin(b*x + c) + d


if __name__ == "__main__":
    init()
    calculate_averages()
    interpolate()
    fit()
    plt.show()
