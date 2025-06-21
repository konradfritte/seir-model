import numpy as np

import model

# TODO:: Create different simulations with different parameters
def seir_simulation():
    healthcare_limit = 25000
    alpha = 1 / 5
    gamma = 1 / 10
    beta = model.seasonal_contact_rate(1.5 / 10, amplitude=0.1)
    delta = 1/365
    epsilon = model.dynamic_mortality(1 / 100, limit=healthcare_limit, magnitude=2)
    zeta = 1 / 1000
    eta = 1 / 20
    ny = 0
    my = 0

    f = model.initialize_seir_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, zeta=zeta, eta=eta, ny=ny, my=my)

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