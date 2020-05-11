import math
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import curve_fit
import statistics as stat

# welke vragen runnen?
runVraag2a = False
runVraag2b = False
runVraag2c = True

# animate chain and plot ending config of chain after runs?
video = False
endplot = False

# for average end_to_end distances
average_end_to_end_squared = []

# define a dot product function used for the rotate operation
def v_dot(a): return lambda b: np.dot(a, b)

class lattice_SAW:
    def __init__(self, N):
        self.N = N
        # initial configuration. Usually we just use a straight chain as inital configuration
        self.init_state = np.dstack((np.arange(N), np.zeros(N)))[0]
        self.state = self.init_state.copy()
        # 3 possible rotations: rotate angles(90,180,270)
        self.rotate_matrix = np.array([[[0, -1], [1, 0]],
                                       [[-1, 0], [0, -1]],
                                       [[0, 1], [1, 0]]])

    # define pivot algorithm process where t is the number of successful steps
    def walk(self, t, c):
        afstand_squared = []
        timestep = 0
        x = []
        y = []
        while timestep < t:
            pick_pivot = np.random.randint(1, self.N - 1)  # pick a pivot site

            # divide the chain up into two parts
            chain_after_pivot = self.state[pick_pivot:]
            chain_before_pivot = self.state[0:pick_pivot]

            # pick a random turn
            symtry_oprtr = self.rotate_matrix[np.random.randint(len(self.rotate_matrix))]

            # turn the chain (always the chain before the pivot, wouldn't matter)
            turned_chain_half = np.apply_along_axis(v_dot(symtry_oprtr), 1, chain_before_pivot - self.state[pick_pivot]) \
                        + self.state[pick_pivot]

            # check if chain overlaps with itself with cdist
            distances_between_points = cdist(turned_chain_half, chain_after_pivot).flatten()

            # determine whether the new state is accepted or rejected
            # aka if no zero distances between points exist
            if len(np.nonzero(distances_between_points)[0]) != len(distances_between_points):
                continue
            else:
                self.state = np.concatenate((turned_chain_half, chain_after_pivot), axis=0)

                # Warm-up steps before statistics usage; dont append distance before c*N timesteps
                if timestep < c * self.N:
                    afstand_squared.append(((self.state[N - 1][0] - self.state[0][0]) ** 2)
                                    + ((self.state[N - 1][1] - self.state[0][1]) ** 2))

                # animate the steps
                if video:
                    x = []
                    y = []
                    for i in range(N):
                        x.append(self.state[i][0])
                        y.append(self.state[i][1])
                    plt.plot(x, y, "-", markersize=3)

                    plt.ylabel('y-as')
                    plt.xlabel('x-as')

                    plt.xlim(0, 2*N)
                    plt.ylim(0.5 * -N, 0.5 * N)
                    plt.grid(b=True, axis='both')

                    plt.draw()
                    plt.pause(0.2)
                    plt.clf()

                timestep += 1

        # plot the ending config of the chain
        if endplot:
            for i in range(N):
                x.append(self.state[i][0])
                y.append(self.state[i][1])
            plt.plot(x, y, '-')
            plt.xlabel("x-as")
            plt.ylabel("y-as")
            plt.grid(b=True, axis='both')
            plt.show()
            
        # append average end_to_end after these amount of steps
        average_end_to_end_squared.append(stat.mean(afstand_squared))
        # reset distance list for next iteration
        afstand_squared = []

        print(f"N = {N} done.")

# Weight for warm-up rotations
c = 1

# amount of steps per total walk
t = 1000

################################ Vraag 2a

if runVraag2a:
    N_tests = [10, 50, 100, 500]

    # run alle N tests
    for N in N_tests:
        chain = lattice_SAW(N)
        chain.walk(t, c)
    

    # plot x = N, y = avg end-to-end distance
    def endtoend2a():
        average_distances = [np.sqrt(value) for value in average_end_to_end_squared]

        plt.plot(N_tests, average_distances, 'o-')
        plt.legend(('Simulatie'))
        plt.xlabel("N")
        plt.ylabel("Average End to end distance")
        plt.grid(b=True, axis='both')
        plt.show()

    endtoend2a()

################################ Vraag 2b

if runVraag2b:
    repeat = 200
    N = 10

    # run chains repeat amount of times
    for _ in range(repeat):
        chain = lattice_SAW(N)
        chain.walk(t, c)

    # plot probability density function of the end-to-end
    def histogram():
        # plot the histogram
        sqrt_end_to_end = [np.sqrt(i) for i in average_end_to_end_squared]
        plt.hist(sqrt_end_to_end, bins=20, density=True)
        plt.xlim(0, 10)
        plt.ylabel("Probability density")
        plt.xlabel("Average End-To-End Distance")
        plt.show()

    histogram()

################################ Vraag 2c

if runVraag2c:
    N_values = []

    N_begin = 10
    N_end = 300
    steps = 40

    # Run different N's
    for N in range(N_begin, N_end, int((N_end - N_begin) / steps)):
        chain = lattice_SAW(N)
        chain.walk(t, c)
        N_values.append(N)

    # function for fitting
    def func(x, a, b):
        return a + (b * np.sqrt(x))

    # fit data with this function
    def fit_data(xdata, ydata):
        popt, pcov = curve_fit(func, xdata, ydata)
        return popt, pcov

    # mean value of end_to_end distances
    def endtoend2c():
        average_distances = [np.sqrt(value) for value in average_end_to_end_squared]

        # fit data with asymptotic function and plot
        popt, pcov = fit_data(N_values, average_distances)
        plt.plot(N_values, func(N_values, *popt), 'r--', label='Fit: y = %5.3f + (%5.3f * sqrt(x))' % tuple(popt))

        # plot simulation data
        plt.plot(N_values, average_distances, 'o-', label="Simulatie")
        plt.xlabel("Nodes N")
        plt.ylabel("Average end to end Distance")
        plt.grid(b=True, axis='both')
        plt.legend()
        plt.show()

    endtoend2c()
