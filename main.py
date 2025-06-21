import numpy as np

import visualization

# Move to separate file

# Approximates a function's values using the euler method: f(t + h) = f(t) + h * f'(t).
# f: The function that should be approximated
# x0: The initial value at which to start
# t: The final time point at which to end
# h: The delta value at which the function should be approximated
def euler_method(f, x0, t, h):
    x = x0
    result = [[0, x]]

    for k in range(1, int(1 + t / h)):
        t = k * h
        x = x + h * f(t, x)

        result.append([t, x])
    return result

# Move to separate file

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
def seir_model(alpha, beta, gamma, delta=0, zeta=0, eta=0, epsilon=0, ny=0, my=0):
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

def seasonality(parameter, amplitude=0, phi=0, period=365):
    def f(t):
        return parameter * (1 + amplitude * np.sin(2 * np.pi * (t + phi) / period))

    return f

def mortality(parameter, limit, magnitude=1):
    def f(x):
        return parameter if x <= limit else magnitude * parameter

    return f

def seir_simulation():
    healthcare_limit = 25000
    alpha = 1 / 5
    gamma = 1 / 10
    beta = seasonality(1.5 / 10, amplitude=0.1)
    delta = 1/365
    epsilon = mortality(1 / 100, limit=healthcare_limit, magnitude=2)
    zeta = 1 / 1000
    eta = 1 / 20
    ny = 0
    my = 0

    f = seir_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, zeta=zeta, eta=eta, ny=ny, my=my)

    n = 83200000
    s = n - 10000
    e = n - s
    ih = i = r = d = 0
    c_s = s
    c_e = e
    c_i = c_ih = c_r = c_d = 0

    x0 = np.array([n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d])

    t = 5 * 365
    h = 1

    results = euler_method(f, x0, t, h)

    return {
        "parameters": {
            "healthcare_limit": healthcare_limit,
            "alpha": alpha,
            "gamma": gamma,
            "beta": beta,
            "delta": delta,
            "epsilon": epsilon,
            "zeta": zeta,
            "eta": eta,
            "ny": ny,
            "my": my
        },
        "results": results
    }

visualization.diagram(seir_simulation())