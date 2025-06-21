import numpy as np

# Initializes the SEIR Model with its relevant parameters.
# Alpha: Determines the transition rate from exposed to infectious compartment. Its reciprocal 1/alpha is the incubation time in days.
# Beta: Determines the transition rate from susceptible to exposed compartment. It is usually connected to the reproduction number R0 and the recovery rate gamma.
# Gamma: Determines the transition rate from infectious/hospitalized to recovered compartment. Its reciprocal 1/gamma is the recovery time in days.
# Delta : Determines the transition rate from recovered to susceptible compartment. Its reciprocal 1/delta is the immunity time in days.
# Zeta: Determines the transition rate from infectious to hospitalized compartment. Its reciprocal 1/zeta is the transition time from mild to severe desease in days.
# Eta: Determines the transition rate from hospitalized to recovered compartment. Its reciprocal 1/eta is the recovery time for a hospitalized patient.
# Epsilon: Determines the transition rate from hospitalized to dead compartment. Its reciprocal 1/epsilon is the decease time in days.
# Ny: Determines the natural birth rate for the total population.
# My: Determines the natural death rate for the total population.
def initialize_seir_model(alpha, beta, gamma, delta=0, zeta=0, eta=0, epsilon=0, ny=0, my=0):
    def f(t, x):
        n, s, e, i, ih, r, _, _, _, _, _, _, _ = x

        dn = (ny - my) * n - epsilon(ih) * ih

        ds = ny * n + delta * r - 1 / n * beta(t) * s * i - my * s
        de = 1 / n  * beta(t) * s * i - alpha * e - my * e
        di = alpha * e - gamma * i - zeta * i - my * i
        dih = zeta * i - eta * ih - epsilon(ih) * ih - my * ih
        dr = gamma * i + eta * ih - delta * r - my * r
        dd = epsilon(ih) * ih

        dc_s = ny * n + delta * r
        dc_e = 1 / n  * beta(t) * s * i
        dc_i = alpha * e
        dc_ih = zeta * i
        dc_r = eta * ih + gamma * i
        dc_d = epsilon(ih) * ih

        dx = np.array([dn, ds, de, di, dih, dr, dd, dc_s, dc_e, dc_i, dc_ih, dc_r, dc_d])

        return dx

    return f