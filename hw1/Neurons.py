import numpy as np

class HHNeuron:
    '''Neuron class using HH equation with Euler method'''
    def __init__(self, dt=0.05, inj_type=None, inj_I=0.0):
        '''
        :param dt: time step length used for estimating potential
        :param inj_type: type of inject current, support 'direct', 'noise', 'sin' and 'slope'
        :param inj_I: value of inject direct current, or direct component of sin,
        or slope of inject current proportion to time, noise is adding white noise to inject direct current.
        '''
        # set parameters
        self.v = -65.0
        self.Cm = 1.0
        self.gl = 0.0
        self.gNa = 120.0
        self.gK = 36.0
        self.ENa = 50.0
        self.EK = -77.0
        self.El = -54.387

        self.dt = dt
        self.inj_type = inj_type
        self.inj_I = inj_I
        self.I = 0

        self.t = 0
        self.m = 0
        self.n = 0
        self.h = 0
        self.alpha_m = 0
        self.alpha_n = 0
        self.alpha_h = 0
        self.beta_m = 0
        self.beta_n = 0
        self.beta_h = 0
        self.INa = 0
        self.IK = 0

    def print_params(self):
        print('INa:%f, IK:%f' % (self.INa, self.IK))
        print('m:%f, am:%f, bm:%f' % (self.m, self.alpha_m, self.beta_m))

    def update_all_params(self):
        self.t += self.dt

        self.update_alpha_m()
        self.update_alpha_n()
        self.update_alpha_h()
        self.update_beta_m()
        self.update_beta_n()
        self.update_beta_h()

        self.update_m()
        self.update_n()
        self.update_h()

        self.update_INa()
        self.update_IK()

        # inject current
        self.update_I()

        self.update_v()

    def update_I(self):
        if self.inj_type == 'direct':
            self.I = self.inj_I
        elif self.inj_type == 'sin':
            self.I = (np.sin(self.t) + 1) * self.inj_I
        elif self.inj_type == 'noise':
            self.I = self.inj_I * (1 + 0.1*np.random.normal(0, 1))
            if self.I < 0 :
                self.I = 0
        elif self.inj_type == 'slope':
            self.I = self.inj_I * self.t
        else:
            raise ('Unsupported inject current type!')

    def update_v(self):
        delta_v = (self.gl*(self.El-self.v) + self.INa + self.IK + self.I) / self.Cm
        self.v += delta_v * self.dt

    def update_INa(self):
        self.INa = self.gNa * (self.m**3) * self.h * (self.ENa-self.v)
    def update_IK(self):
        self.IK = self.gK * (self.n**4) * (self.EK - self.v)

    def update_m(self):
        delta_m = self.alpha_m*(1-self.m) - self.beta_m*self.m
        self.m += delta_m * self.dt
    def update_n(self):
        delta_n = self.alpha_n * (1 - self.n) - self.beta_n * self.n
        self.n += delta_n * self.dt
    def update_h(self):
        delta_h = self.alpha_h * (1 - self.h) - self.beta_h * self.h
        self.h += delta_h * self.dt

    def update_alpha_m(self):
        self.alpha_m = 0.1 * (self.v+40) / \
                       (1 - np.exp(-(self.v+40)/10))
    def update_beta_m(self):
        self.beta_m = 4 * np.exp(-(self.v+65)/18)
    def update_alpha_n(self):
        self.alpha_n = 0.01 * (self.v+55) / \
                       (1 - np.exp(-(self.v+55)/10))
    def update_beta_n(self):
        self.beta_n = 0.125 * np.exp(-(self.v+65)/80)
    def update_alpha_h(self):
        self.alpha_h = 0.07 * np.exp(-(self.v+65)/20)
    def update_beta_h(self):
        self.beta_h = 1 / (np.exp(-(self.v+35)/10) + 1)

def runge_kutta(f, y, dt):
    y1 = f(y) * dt
    y2 = f(y + y1 / 2) * dt
    y3 = f(y + y2 / 2) * dt
    y4 = f(y + y3) * dt
    return y + (y1 + 2 * y2 + 2 * y3 + y4) / 6

class HHNeuron_RK(HHNeuron):
    '''Neuron class using HH equation with Runge Kutta method'''
    def __init__(self, dt=0.05, inj_type=None, inj_I=0.0):
        HHNeuron.__init__(self, dt, inj_type, inj_I)
    def update_v(self):
        self.v = runge_kutta(self.d_v, self.v, self.dt)
    def update_m(self):
        self.m = runge_kutta(self.d_m, self.m, self.dt)
    def update_n(self):
        self.n = runge_kutta(self.d_n, self.n, self.dt)
    def update_h(self):
        self.h = runge_kutta(self.d_h, self.h, self.dt)
    def d_v(self, v):
        return (self.gl * (self.El - v) + self.INa + self.IK + self.I) / self.Cm
    def d_m(self, m):
        return self.alpha_m * (1 - m) - self.beta_m * m
    def d_n(self, n):
        return self.alpha_n * (1 - n) - self.beta_n * n
    def d_h(self, h):
        return self.alpha_h * (1 - h) - self.beta_h * h
