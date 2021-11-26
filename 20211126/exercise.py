import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt

def with_dense_matrix(N):
    D = N**2 #size of matrix A
    y0 = 1 #boundary condition for y = 0

    A = (-4) * np.eye(D) + np.eye(D, k=-1) + np.eye(D, k=1) + np.eye(D, k=N) + np.eye(D, k=-N)
    for n in np.arange(N, D, step=N):
        A[n, n-1] = 0
        A[n-1, n] = 0

    b = np.zeros(D)
    b[:N] = -y0

    np.savetxt('dense_matrix_solved.dat', linalg.solve(A, b))

    # One problem I see here is that the boundary is not included in the plot.
    # We could solve tha by simply appending to the results array, but it would
    # make the solution so much uglier. That is why I omitted it.


if __name__ == "__main__":
    N = 50 #size of square grid

    with_dense_matrix(N)

    DS = np.array(np.loadtxt('dense_matrix_solved.dat'))

    # makes a heat map from the values in the results matrix
    # origin lower sets the start to be from the lower left (else it would
    # start at the top)
    # extent tells the axis to go from 0 to 1
    # cmap chooses a nice looking color map for the heat map
    plt.imshow(np.reshape(DS, (N, N)), extent=[0,1,0,1], cmap='jet', origin='lower')
    plt.colorbar()
    plt.show()

