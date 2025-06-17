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
def seir_model(alpha, beta, gamma):
    def f(t, x):

        s, e, i, _ = x

        ds = -beta * s * i
        de = beta * s * i - alpha * e
        di = alpha * e - gamma * i
        dr = gamma * i

        dx = np.array([ds, de, di, dr])

        return dx

    return f

def seir_simulation():
    # Reproduction number R0: if R0 > 1, the pandemic starts, if R0 < 1, it ceases
    r0 = 2

    alpha = 1/5; gamma = 1/10; beta = r0 * gamma

    f = seir_model(alpha, beta, gamma)

    s = 0.999; e = 1 - s; i = 0; r = 0
    x0 = np.array([s, e, i, r])

    t = 365
    h = 1

    result = euler_method(f, x0, t, h)

    t, x = zip(*result)
    x = np.array(x)
    labels = ["s", "e", "i", "r"]

    plt.plot(t, x, label=labels)
    plt.legend()
    plt.show()

seir_simulation()