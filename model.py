import numpy as np

# Initializes the SEhRS-D Model.
# Alpha: Determines the transition rate from exposed to infectious compartment. Its reciprocal 1/alpha is the incubation time in days.
# Beta: Determines the transition rate from susceptible to exposed compartment. It is usually connected to the reproduction number R0 and the recovery rate gamma.
# Gamma: Determines the transition rate from infectious to recovered compartment. Its reciprocal 1/gamma is the recovery time in days.
# Delta : Determines the transition rate from recovered to susceptible compartment. Its reciprocal 1/delta is the immunity time in days.
# Epsilon I: Determines the transition rate from infectious to dead compartment. Its reciprocal 1/epsilon_i is the decease time in days.
# Epsilon H: Determines the transition rate from hospitalized to dead compartment. Its reciprocal 1/epsilon_h is the decease time in days.
# Zeta: Determines the transition rate from infectious to hospitalized compartment. Its reciprocal 1/zeta is the transition time from mild to severe desease in days.
# Eta: Determines the transition rate from hospitalized to recovered compartment. Its reciprocal 1/eta is the recovery time for a hospitalized patient.
# Ny: Determines the natural birth rate for the total population.
# My: Determines the natural death rate for the total population.
def seihrsd_model(alpha, beta, gamma, delta, epsilon_i, epsilon_h, zeta, eta, ny, my):
    def f(t, x):
        n, s, e, i, h, r, _, _, _, _, _, _, _ = x

        dn = (ny - my) * n - epsilon_i(t,h) * h - epsilon_i() * i

        ds = ny * n + delta() * r - (1 / n * beta(t) * i + my) * s
        de = 1 / n  * beta(t) * s * i - (alpha() + my) * e
        di = alpha() * e - (gamma() + zeta() + epsilon_i() + my) * i
        dh = zeta() * i - (eta() + epsilon_h(t, h) + my) * h
        dr = gamma() * i + eta() * h - (delta() + my) * r
        dd = epsilon_h(h) * h + epsilon_i() * i

        # Rates for cumulated cases
        dc_s = ny * n + delta() * r
        dc_e = 1 / n  * beta(t) * s * i
        dc_i = alpha() * e
        dc_h = zeta() * i
        dc_r = eta() * h + gamma() * i
        dc_d = epsilon_h(t, h) * h + epsilon_i() * i

        dx = np.array([dn, ds, de, di, dh, dr, dd, dc_s, dc_e, dc_i, dc_h, dc_r, dc_d])

        return dx

    return f

# Moderator function for the alpha model parameter
# Alpha: The transition rate from exposed to infectious compartment
def alpha_moderator(alpha):
    return basic_moderator(alpha)

# Moderator function for the beta model parameter
# Beta: The transition rate from susceptible to exposed compartment
# Amplitude: The amplitude by which beta should be modulated
# Phi: The phase shift for the modulator
# Period: The period of the modulator
def beta_moderator(beta, amplitude=0, phi=0, period=365):
    def f(t=0,x=0):
        return beta * (1 + amplitude * np.sin(2 * np.pi * (t + phi) / period))

    return f

# Moderator function for the gamma model parameter
# Gamma: The transition rate from infectious to recovered compartment
def gamma_moderator(parameter):
    return basic_moderator(parameter)

# Moderator function for the delta model parameter
# Delta: The transition rate from recovered to susceptible compartment
def delta_moderator(parameter):
    return basic_moderator(parameter)

# Moderator function for the epsilon model parameter
# Epsilon: The transition rate from infectious/hospitalized to dead compartment
# Limit: The threshold value at which the rate should be assigned a new value
def epsilon_moderator(parameter, limit=0, magnitude=1):
    def f(t=0,x=0):
        return parameter if x <= limit else magnitude * parameter

    return f

# Moderator function for the zeta model parameter
# Zeta: The transition rate from infectious to hospitalized compartment
def zeta_moderator(parameter):
    return basic_moderator(parameter)

# Moderator function for the eta model parameter
# Eta: The transition rate from hospitalized to recovered compartment
def eta_moderator(parameter):
    return basic_moderator(parameter)

def basic_moderator(parameter):
    def f(t=0,x=0):
        return parameter

    return f