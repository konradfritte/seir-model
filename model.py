import numpy as np

# Initializes the SEIHRS-D Model.
# Alpha: Determines the transition rate from exposed to infectious compartment. Its reciprocal 1/alpha is the incubation time in days.
# Beta: Determines the transition rate from susceptible to exposed compartment. It is usually connected to the reproduction number R0 and the recovery rate gamma.
# Gamma: Determines the transition rate from infectious to recovered compartment. Its reciprocal 1/gamma is the recovery time in days.
# Delta : Determines the transition rate from recovered to susceptible compartment. Its reciprocal 1/delta is the immunity time in days.
# Epsilon: Determines the transition rate from hospitalized to dead compartment. Its reciprocal 1/epsilon is the decease time in days.
# Epsilon H: Determines the transition rate from hospitalized to dead compartment. Its reciprocal 1/epsilon is the decease time in days.
# Zeta: Determines the transition rate from infectious to hospitalized compartment. Its reciprocal 1/zeta is the transition time from mild to severe desease in days.
# Eta: Determines the transition rate from hospitalized to recovered compartment. Its reciprocal 1/eta is the recovery time for a hospitalized patient.
# Ny: Determines the natural birth rate for the total population.
# My: Determines the natural death rate for the total population.
def seihrsd_model(alpha, beta, gamma, delta, epsilon, epsilon_h, zeta, eta, ny, my):
    def f(t, x):
        n, s, e, i, ih, r, _, _, _, _, _, _, _ = x

        dn = (ny - my) * n - epsilon(t,ih) * ih - epsilon() * i

        ds = ny * n + delta() * r - 1 / n * beta(t) * s * i - my * s
        de = 1 / n  * beta(t) * s * i - alpha() * e - my * e
        di = alpha() * e - gamma() * i - zeta() * i - epsilon() * i - my * i
        dih = zeta() * i - eta() * ih - epsilon_h(t, ih) * ih - my * ih
        dr = gamma() * i + eta() * ih - delta() * r - my * r
        dd = epsilon_h(ih) * ih + epsilon() * i

        dc_s = ny * n + delta() * r
        dc_e = 1 / n  * beta(t) * s * i
        dc_i = alpha() * e
        dc_ih = zeta() * i
        dc_r = eta() * ih + gamma() * i
        dc_d = epsilon_h(t, ih) * ih + epsilon() * i

        dx = np.array([dn, ds, de, di, dih, dr, dd, dc_s, dc_e, dc_i, dc_ih, dc_r, dc_d])

        return dx

    return f

def alpha_moderator(parameter):
    return basic_moderator(parameter)

def beta_moderator(parameter, amplitude=0, phi=0, period=365):
    def f(t=0,x=0):
        return parameter * (1 + amplitude * np.sin(2 * np.pi * (t + phi) / period))

    return f

def gamma_moderator(parameter):
    return basic_moderator(parameter)

def delta_moderator(parameter):
    return basic_moderator(parameter)

def epsilon_moderator(parameter, limit=0, magnitude=1):
    def f(t=0,x=0):
        return parameter if x <= limit else magnitude * parameter

    return f

def zeta_moderator(parameter):
    return basic_moderator(parameter)

def eta_moderator(parameter):
    return basic_moderator(parameter)

def basic_moderator(parameter):
    def f(t=0,x=0):
        return parameter

    return f