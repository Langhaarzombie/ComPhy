import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
from scipy.misc import derivative
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt

time = np.array([])
m_free = np.array([])
interpolation = None

_, ax = plt.subplots(3)


def init(file='20211019/oscillator.dat'):
    global time
    global m_free
    data = np.loadtxt(file)
    time = data.T[0]
    m_free = data.T[1:4].T


def calculate_averages():
    avg = np.ndarray((3, 1))
    avg[0] = np.mean(m_free.T[0])
    avg[1] = np.mean(m_free.T[1])
    avg[2] = np.mean(m_free.T[2])
    print("The average for the free layer:")
    print(avg)
    print("Magintude of that vector:")
    print(np.linalg.norm(avg))


def interpolate():
    global interpolation
    interpolation = UnivariateSpline(time, m_free.T[0], s=0)
    xs = np.linspace(time[0], time[-1], num=10000, endpoint=False)
    ax[0].plot(xs, interpolation(xs), "x", label="equidistant data")


def fit():
    fit_points = 12500
    fit, _ = curve_fit(sinus, time[-fit_points:], m_free.T[0]
                       [-fit_points:], p0=[0, 22e9, np.pi/2, 0])
    ax[0].plot(time[-fit_points:], sinus(time[-fit_points:], *fit), label="fit curve")

def ddx():
    print(interpolation.derivative)
    # ax[1].plot(time, interpolation.derivative)

def fourier():
    # timestep defines how fine the ftt resolution is
    # N is the number of sample points
    # yf contains the frequencies embedded in the interpolation of the original data
    # xf contains the frequencies included in ftt
    # The usage of [:N//2] is necessary as the lower half of the array
    # contains the positiv frequencies and the upper half the negative.
    # As the spectrum is symmetric we can reduce the dataset by half.
    # When plotting, we multiply by 2 to represent this combination of positive
    # and negative frequencies and divide by N to get the relative occurance.

    deltaT = 5e-11
    N = int((time[-1] - time[0])//deltaT)
    yf = fft(interpolation(np.linspace(time[0], time[-1], num=N, endpoint=False)))
    xf = fftfreq(N, d=deltaT)[:N//2]
    ax[2].plot(xf, 2/N * np.abs(yf[:N//2]))

def sinus(x, a, b, c, d):
    return a*np.sin(b*x + c) + d


if __name__ == "__main__":
    init()
    calculate_averages()
    interpolate()
    fit()
    ddx()
    fourier()
    plt.show()
