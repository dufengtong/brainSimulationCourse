import matplotlib.pyplot as plt
import numpy as np
from Neurons import HHNeuron_RK, HHNeuron



if __name__ == '__main__':
    time_step = 0.05
    type = 'slope'
    inject = 0.4
    rk_neuron = HHNeuron_RK(dt=time_step, inj_type=type, inj_I=inject)
    euler_neuron = HHNeuron(dt=time_step, inj_type=type, inj_I=inject)
    time = 100 # ms
    steps = int(time / rk_neuron.dt) + 1
    x = rk_neuron.dt * np.arange(steps)
    # numerical solution
    rk_y = rk_neuron.v * np.ones(x.shape)
    euler_y = euler_neuron.v * np.ones(x.shape)
    inject_I = np.zeros(x.shape)
    fig = plt.figure()
    ax_inj = fig.add_subplot(211)
    ax = fig.add_subplot(212)


    for i in range(0, steps):
        rk_y[i] = rk_neuron.v
        euler_y[i] = euler_neuron.v
        inject_I[i] = rk_neuron.I
        print('time: %d rk_v: %f euler_v: %f diff: %f' % (rk_neuron.t, rk_neuron.v, euler_neuron.v,
        rk_neuron.v-euler_neuron.v))
        rk_neuron.update_all_params()
        euler_neuron.update_all_params()

    ax.set_xlim([0, time])
    ax.set_ylim([-80, 80])
    ax.set_xlabel('time(ms)')
    ax.set_ylabel('v(mV)')
    ax.set_title('neuron response dt:%.2f current type:%s'%(rk_neuron.dt, rk_neuron.inj_type))
    rk_handle, = ax.plot(x, rk_y, color='r', label='Runge Kutta')
    euler_handle, = ax.plot(x, euler_y, color='y', label='Euler')
    ax.legend(handles=[euler_handle, rk_handle], loc=1)
    ax_inj.plot(x, inject_I)
    ax_inj.set_xlim([0, time])
    ax_inj.set_ylim([0, 50])
    ax_inj.set_xlabel('time(ms)')
    ax_inj.set_ylabel('I(mA)')
    ax_inj.set_title('Inject Current')
    plt.show()