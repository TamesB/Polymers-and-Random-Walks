import math
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
import statistics as stat

runVraag1a = False
runVraag1b = False
runVraag1c = True

class chain:
    def __init__(self, N):
        # length of chain
        self.N = N

        # first node at (0,0)
        self.current_state = [[0, 0]]

        # all steps a new node can take
        self.steps = [[1,0], [0,1], [-1,0], [0,-1]]
        
    def generate_random_chain(self, animate):
        # place nodes
        while len(self.current_state) < self.N:
            # choose random step
            random_step = self.steps[np.random.randint(len(self.steps))]

            last_node_coords = self.current_state[-1]

            # get the new node coordinates to place a new node
            new_node_x = last_node_coords[0] + random_step[0]
            new_node_y = last_node_coords[1] + random_step[1]
            # append to chain
            self.current_state.append([new_node_x, new_node_y])
            
            # animate drawing
            if animate:
                x = [coord[0] for coord in self.current_state]
                y = [coord[1] for coord in self.current_state]
                plt.plot(x, y, "-", markersize=3)
                plt.ylabel('y-as')
                plt.xlabel('x-as')
                plt.grid(b=True, axis='both')
                plt.draw()
                plt.pause(0.001)
                plt.clf()

    # after generation, calculate end-to-end
    def calculate_end_to_end(self):
        afstand_squared = ((self.current_state[0][1] - self.current_state[-1][1])**2 + (self.current_state[0][0] - self.current_state[-1][0])**2)
        return afstand_squared

    # generate repeat amount of chains for average end to end
    def repeat_generation(self, repeat):
        distances = []

        # generate repeat amount of chains and calculate end-to-end
        for _ in range(repeat):
            self.current_state = [[0,0]]
            self.generate_random_chain()
            self.calculate_end_to_end()
            distances.append(self.afstand_squared)

        # return the average of all of these distances
        return stat.mean(distances)

    # plot the chain
    def plot_chain(self):
        x = [coord[0] for coord in self.current_state]
        y = [coord[1] for coord in self.current_state]
        plt.plot(x, y, '-')
        plt.xlabel("x-as")
        plt.ylabel("y-as")
        plt.grid(b=True, axis='both')
        plt.show()

############# Tests
#N = 300
#
#ketting = chain(N)
#ketting.generate_random_chain(True)
#ketting.plot_chain()

##############


############ Vraag 1a
if runVraag1a:
    N_tests = [10, 50, 100, 500, 1000, 5000, 10000]

    # run alle N tests
    for N in N_tests:
        ketting = chain(N)
        ketting.generate_random_chain(False)
        ketting.plot_chain()

    # amount of chains to be generated
    runs_per_chain = 1000

    avg_distance_list = []

    for N in N_tests:
        ketting = chain(N)
        avg_distance = ketting.repeat_generation(runs_per_chain)
        avg_distance_list.append(avg_distance)
        print(f"N = {N} Done.")

    # plot all the chains with (x = N, y = endtoend distance)
    def endtoend1a():
        average_distances = [np.sqrt(value) for value in avg_distance_list]

        plt.plot(N_tests, average_distances, 'o-')
        plt.legend('Simulatie')
        plt.xlabel("N")
        plt.ylabel("Average End to end distance")
        plt.grid(b=True, axis='both')
        plt.show()

    endtoend1a()


############### Vraag 1b

if runVraag1b:
    chain_crossed_list = []
    N_list = []
    # go from N=1 to N=10
    for N in range(1, 10):
        chain_crossed = False
        N_list.append(N)

        ketting = chain(N)
        ketting.generate_random_chain(False)

        # check of de chain zichzelf heeft overlapt
        seen = []
        for number in ketting.current_state:
            if number in seen:
                chain_crossed = True
            else:
                seen.append(number)

        if chain_crossed:
            chain_crossed_list.append(1)
        else:
            chain_crossed_list.append(0)
    
    plt.scatter(N_list, chain_crossed_list)
    plt.xlabel("Chain length")
    plt.ylabel("Crossed (1 = Yes, 0 = No)")
    plt.title("Chain crossed itself vs chain length")
    plt.show()

############## Vraag 1c
if runVraag1c:
    N_start = 10
    N_end = 10000
    steps = 100

    gyration_radius_list = []
    N_list = []
    # run alle N tests
    for N in range(N_start, N_end, int((N_end - N_start) / steps)):
        N_list.append(N)
        ketting = chain(N)
        ketting.generate_random_chain(False)

        gyration_sum = 0

        # bereken de centre of mass
        x_sum = 0
        y_sum = 0
        for point in ketting.current_state:
            x_sum += point[0]
            y_sum += point[1]
        
        centre_of_mass_x = x_sum / N
        centre_of_mass_y = y_sum / N

        # bereken de som
        for point in ketting.current_state:
            dist_x = point[0] - centre_of_mass_x
            dist_y = point[1] - centre_of_mass_y

            dist_squared = (dist_x**2 + dist_y**2)
            gyration_sum += dist_squared

        # bereken de gyration radius squared
        gyration_radius = math.sqrt((1 / (N - 1)) * gyration_sum)
        gyration_radius_list.append(gyration_radius)

    # analytical sqrt(N)
    analytical_N = [0.4*math.sqrt(N) for N in N_list]

    plt.plot(N_list, gyration_radius_list, 'o-', label="Simulation")
    plt.plot(N_list, analytical_N, 'r--', label="y=0.4*sqrt(N)")
    plt.xlabel("Chain length N")
    plt.ylabel("Gyration radius")
    plt.title("Gyration radius vs chain length")
    plt.legend()
    plt.show()
