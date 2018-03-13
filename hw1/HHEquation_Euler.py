# Copyright 2018 dufengtong
#
# a neuron class using HH equation as neuron model
# solving HH equation using Euler method
import matplotlib.pyplot as plt
import numpy as np
from Neurons import HHNeuron


if __name__ == '__main__':
    neuron = HHNeuron(dt=0.001, inj_type='slope', inj_I=1.0)

    time = 100 # s
    steps = int(time / neuron.dt) + 1
    x = neuron.dt * np.arange(steps)
    # numerical solution
    y = neuron.v * np.ones(x.shape)

    fig = plt.figure()
    ax = fig.add_subplot(111)


    for i in range(0, steps):
        y[i] = neuron.v
        print('time: %d v: %f ' % (neuron.t, neuron.v))
        neuron.update_all_params()

    ax.set_xlim([0, time])
    ax.set_ylim([-80, 80])
    ax.plot(x, y, color='k')
    plt.show()