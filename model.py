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
# def seihrsd_model(alpha, beta, gamma, delta, epsilon, epsilon_h, zeta, eta, ny, my):
#     def f(t, x):
#         n, s, e, i, ih, r, _, _, _, _, _, _, _ = x

#         dn = (ny - my) * n - epsilon(t,ih) * ih - epsilon() * i

#         ds = ny * n + delta() * r - 1 / n * beta(t) * s * i - my * s
#         de = 1 / n  * beta(t) * s * i - alpha() * e - my * e
#         di = alpha() * e - gamma() * i - zeta() * i - epsilon() * i - my * i
#         dih = zeta() * i - eta() * ih - epsilon_h(t, ih) * ih - my * ih
#         dr = gamma() * i + eta() * ih - delta() * r - my * r
#         dd = epsilon_h(ih) * ih + epsilon() * i

#         # Rates for cumulated cases
#         dc_s = ny * n + delta() * r
#         dc_e = 1 / n  * beta(t) * s * i
#         dc_i = alpha() * e
#         dc_ih = zeta() * i
#         dc_r = eta() * ih + gamma() * i
#         dc_d = epsilon_h(t, ih) * ih + epsilon() * i

#         dx = np.array([dn, ds, de, di, dih, dr, dd, dc_s, dc_e, dc_i, dc_ih, dc_r, dc_d])

#         return dx

#     return f

def seihrsd_model(alpha, beta, gamma, delta, epsilon, epsilon_h, zeta, eta, ny, my):
    def f(t, x):
        n, s, e, i, ih, r, d, _, _, _, _, _, _ = x

        p = [s, e, i, ih, r, d]

        w_ss = - 1 / n * beta(t) * i
        w_se = 0
        w_si = 0
        w_sh = 0
        w_sr = delta()
        w_sd = 0

        w_es = 1 / n  * beta(t) * i
        w_ee = - alpha()
        w_ei = 0
        w_eh = 0
        w_er = 0
        w_ed = 0

        w_is = 0
        w_ie = alpha()
        w_ii = - gamma() - zeta() - epsilon()
        w_ih = 0
        w_ir = 0
        w_id = 0

        w_hs = 0
        w_he = 0
        w_hi = zeta()
        w_hh = - eta() - epsilon_h(t, ih)
        w_hr = 0
        w_hd = 0

        w_rs = 0
        w_re = 0
        w_ri = gamma()
        w_rh = eta()
        w_rr = - delta()
        w_rd = 0

        w_ds = 0
        w_de = 0
        w_di = epsilon()
        w_dh = epsilon_h(t, ih)
        w_dr = 0
        w_dd = 0

        q = np.array([
            [w_ss, w_se, w_si, w_sh, w_sr, w_sd],
            [w_es, w_ee, w_ei, w_eh, w_er, w_ed],
            [w_is, w_ie, w_ii, w_ih, w_ir, w_id],
            [w_hs, w_he, w_hi, w_hh, w_hr ,w_hd],
            [w_rs, w_re, w_ri, w_rh, w_rr, w_rd],
            [w_ds, w_de, w_di, w_dh, w_dr, w_dd]
             ])

        if int((np.sum(q))) != 0:
            print(q)
            print(np.sum(q, axis=0))
            raise "Invalid Transition Matrix"

        dx = np.dot(q, p)

        # ds = ny * n + delta() * r - 1 / n * beta(t) * s * i - my * s
        # de = 1 / n  * beta(t) * s * i - alpha() * e - my * e
        # di = alpha() * e - gamma() * i - zeta() * i - epsilon() * i - my * i
        # dih = zeta() * i - eta() * ih - epsilon_h(t, ih) * ih - my * ih
        # dr = gamma() * i + eta() * ih - delta() * r - my * r
        # dd = epsilon_h(ih) * ih + epsilon() * i

        dn = 0

        dc_s = 0
        dc_e = 0
        dc_i = 0
        dc_h = 0
        dc_r = 0
        dc_d = 0

        # dx = np.array([dn, ds, de, di, dih, dr, dd, dc_s, dc_e, dc_i, dc_ih, dc_r, dc_d])

        return np.concatenate((np.array([dn]), dx, np.array([dc_s, dc_e, dc_i, dc_h, dc_r, dc_d])))

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