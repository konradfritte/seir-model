import numpy as np
import matplotlib.pyplot as plt

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

# Initializes the SEIR Model with its relevant parameters.
# Alpha: Determines the transition rate from exposed to infectious people. Its reciprocal 1/alpha is the incubation time in days.
# Beta: Determines the transition rate from susceptible to exposed people. It is usually connected to the reproduction number R0 and the recovery rate gamma.
# Gamma: Determines the recovery rate from infectious to recovered people. Its reciprocal 1/gamma is the recovery time in days.
# Delta : Determines the transition rate from recovered to susceptible people. Its reciprocal 1/delta is the immunity time in days.
# Ny: Determines the natural birth rate for the total population.
# My: Determines the natural death rate for the total population.
def seir_model(alpha, beta, gamma, epsilon=0, delta=0, ny=0, my=0):
    def f(t, x):
        n, s, e, i, r, _, _, _, _, _, _ = x

        dn = (ny - my) * n - epsilon * i

        ds = ny * n + delta * r - 1 / n * beta(t) * s * i - my * s
        de = 1 / n  * beta(t) * s * i - alpha * e - my * e
        di = alpha * e - gamma * i - epsilon * i - my * i
        dr = gamma * i - delta * r - my * r
        dd = epsilon * i

        dc_s = ny * n + delta * r
        dc_e = 1 / n  * beta(t) * s * i
        dc_i = alpha * e
        dc_r = gamma * i
        dc_d = epsilon * i

        dx = np.array([dn, ds, de, di, dr, dd, dc_s, dc_e, dc_i, dc_r, dc_d])

        return dx

    return f

def beta_modulator(beta, amplitude=0, phi=0, period=365):
    def f(t):
        return beta * (1 + amplitude * np.sin(2 * np.pi * (t + phi) / period))

    return f

def seir_simulation():
    r0 = 1.75

    alpha = 1 / 5
    gamma = 1 / 10
    beta = beta_modulator(beta=r0 * gamma, amplitude=0.1)
    delta = 1 / 365
    epsilon = 1 / 10000
    ny = 0
    my = 0

    f = seir_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, ny=ny, my=my)

    n = 83200000
    s = n - 10000
    e = n - s
    i = r = d = 0
    c_s = s
    c_e = e
    c_i = c_r = c_d = 0

    x0 = np.array([n, s, e, i, r, d, c_s, c_e, c_i, c_r, c_d])

    t = 10 * 365
    h = 1

    result = euler_method(f, x0, t, h)
    t, x = zip(*result)
    n, s, e, i, r, d, c_s, c_e, c_i, c_r, c_d = zip(*x)

    dc_s = np.insert(np.diff(c_s), 0, 0)
    dc_e = np.insert(np.diff(c_e), 0, 0)
    dc_i = np.insert(np.diff(c_i), 0, 0)
    dc_r = np.insert(np.diff(c_r), 0, 0)
    dc_d = np.insert(np.diff(c_d), 0, 0)

    plt.style.use("ggplot")

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)

    ax1.set_title("Development of SEIRD compartments")
    ax1.set_xlabel("Day")
    ax1.set_ylabel("Amount")

    ax1.plot(t, s, label="Susceptible", color="blue")
    ax1.plot(t, e, label="Exposed", color="orange")
    ax1.plot(t, i, label="Infectious", color="red")
    ax1.plot(t, r, label="Recovered", color="green")
    ax1.plot(t, d, label="Dead", color="grey")

    ax1.legend()

    ax2.set_title("New daily cases for SEIRD compartments")
    ax2.set_xlabel("Day")
    ax2.set_ylabel("Amount")

    ax2.plot(t, dc_s, label="Susceptible", color="blue")
    ax2.plot(t, dc_e, label="Exposed", color="orange")
    ax2.plot(t, dc_i, label="Infectious", color="red")
    ax2.plot(t, dc_r, label="Recovered", color="green")
    ax2.plot(t, dc_d, label="Dead", color="grey")

    ax2.legend()

    ax3.set_title("Cumulative cases for SEIRD compartments")
    ax3.set_xlabel("Day")
    ax3.set_ylabel("Amount")
    ax3.set_yscale("log")

    ax3.plot(t, c_s, label="Susceptible", color="blue")
    ax3.plot(t, c_e, label="Exposed", color="orange")
    ax3.plot(t, c_i, label="Infectious", color="red")
    ax3.plot(t, c_r, label="Recovered", color="green")
    ax3.plot(t, c_d, label="Dead", color="grey")

    ax3.legend()

    ax4.set_title("Transition rates for SEIRD compartments")
    ax4.set_xlabel("Day")
    ax4.set_ylabel("Rate")

    ax4.plot(t, np.full(len(t), alpha), label="Alpha", color="red")
    ax4.plot(t, beta(np.array(t)), label="Beta", color="orange")
    ax4.plot(t, np.full(len(t), gamma), label="Gamma", color="green")
    ax4.plot(t, np.full(len(t), delta), label="Delta", color="blue")
    ax4.plot(t, np.full(len(t), epsilon), label="Epsilon", color="grey")

    ax4.legend()

    ax5.set_title("Development of total population")
    ax5.set_xlabel("Day")
    ax5.set_ylabel("Amount")

    ax5.plot(t, n, label="Population", color="grey")

    ax5_twin = ax5.twinx()
    ax5_twin.set_ylabel("Rate")
    ax5_twin.plot(t, np.full(len(t), ny), linestyle=":", label="Ny", color="blue")
    ax5_twin.plot(t, np.full(len(t), my), linestyle=":", label="My", color="brown")

    ax5_twin.tick_params(axis="y", labelcolor="brown")
    ax5_twin.set_ylim(bottom=0, top=0.01)
    ax5_twin.grid(False)

    ax5.legend()
    ax5_twin.legend(loc="lower right")

    plt.show()


seir_simulation()