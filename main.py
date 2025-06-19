import numpy as np
import matplotlib.pyplot as plt

# Approximates a function's values using the euler method: f(t + h) = f(t) + h * f'(t).
# f: The function that should be approximated
# x0: The initial value at which to start
# t: The final time point at which to end
# h: The delta value at which the function should be approximated
def euler_method(f, x0, t, h):
    x = x0
    result = []

    for k in range(int(1 + t / h)):
        t = k * h
        x = x + h * f(t, x)

        result.append([t, x])
    return result

# Initializes the SEIR Model with its relevant parameters.
# Alpha: Determines the transition rate from exposed to infectious people. Its reciprocal 1/alpha is the incubation time in days.
# Beta: Determines the transmission rate from subsceptible to exposed people. It is usually connected to the reproduction number R0 and the recovery rate gamma.
# Gamma: Determines the recovery rate from infectious to recovered people. Its reciprocal 1/gamma is the recovery time in days.
def seir_model(alpha, beta, gamma, ny=0, my=0):
    def f(t, x):
        n, s, e, i, r = x

        dn = (ny - my) * n
        ds = ny * n - 1 / n * beta(t) * s * i - my * s
        de = 1 / n  * beta(t) * s * i - alpha * e - my * e
        di = alpha * e - gamma * i - my * i
        dr = gamma * i - my * r

        dx = np.array([dn, ds, de, di, dr])

        return dx

    return f

def beta_modulator(beta, amplitude=0, phi=0, period=365):
    def f(t):
        return beta * (1 + amplitude * np.sin(2 * np.pi * (t + phi) / period))

    return f

def seir_simulation():
    r0 = 2.4

    alpha = 1/5.5
    gamma = 1/3
    beta = beta_modulator(r0 * gamma)

    f = seir_model(alpha, beta, gamma)

    n = 83200000
    s = n - 40000
    e = n - s
    i = 0
    r = 0
    x0 = np.array([n, s, e, i, r])

    t = 365
    h = 1

    result = euler_method(f, x0, t, h)
    t, x = zip(*result)
    _, s, e, i, r = zip(*x)


    fig, ax = plt.subplots()


    plt.grid(linestyle=":", color="lightgrey")

    ax.plot(t, s, label="Subsceptible", color="blue")
    ax.plot(t, e, label="Exposed", color="orange")
    ax.plot(t, i, label="Infectious", color="red")
    ax.plot(t, r, label="Recovered", color="green")

    ax.legend()
    plt.show()

seir_simulation()