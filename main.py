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
        n, s, e, i, r, _, _, _ = x

        dn = (ny - my) * n

        ds = ny * n - 1 / n * beta(t) * s * i - my * s
        de = 1 / n  * beta(t) * s * i - alpha * e - my * e
        di = alpha * e - gamma * i - my * i
        dr = gamma * i - my * r

        d_c_e = 1 / n  * beta(t) * s * i
        d_c_i = alpha * e
        d_c_r = gamma * i

        dx = np.array([dn, ds, de, di, dr, d_c_e, d_c_i, d_c_r])

        return dx

    return f

def beta_modulator(beta, amplitude=0, phi=0, period=365):
    def f(t):
        return beta * (1 + amplitude * np.sin(2 * np.pi * (t + phi) / period))

    return f

def seir_simulation():
    r0 = 2.4

    alpha = 1/5
    gamma = 1/10
    beta = beta_modulator(r0 * gamma, 1)
    ny = 0.005
    my = 0.004

    f = seir_model(alpha, beta, gamma, ny, my)

    n = 83200000
    s = n - 40000
    e = n - s
    i = r = 0
    c_e = c_i = c_r = 0

    x0 = np.array([n, s, e, i, r, c_e, c_i, c_r])

    t = 1 * 365
    h = 1

    result = euler_method(f, x0, t, h)
    t, x = zip(*result)
    n, s, e, i, r, c_e, c_i, c_r = zip(*x)


    d_c_e = np.diff(c_e)
    d_c_i = np.diff(c_i)
    d_c_r = np.diff(c_r)

    plt.style.use("ggplot")

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)

    ax1.set_title("Development of SEIR compartments")
    ax1.set_xlabel("Day")
    ax1.set_ylabel("Amount")

    ax1.plot(t, s, label="Subsceptible", color="blue")
    ax1.plot(t, e, label="Exposed", color="orange")
    ax1.plot(t, i, label="Infectious", color="red")
    ax1.plot(t, r, label="Recovered", color="green")

    ax1.legend()

    ax2.set_title("New daily cases for EIR compartments")
    ax2.set_xlabel("Day")
    ax2.set_ylabel("Amount")

    ax2.plot(d_c_e, label="Exposed", color="orange")
    ax2.plot(d_c_i, label="Infectious", color="red")
    ax2.plot(d_c_r, label="Recovered", color="green")

    ax2.legend()

    ax3.set_title("Cumulative cases for EIR compartments")
    ax3.set_xlabel("Day")
    ax3.set_ylabel("Amount")

    ax3.plot(t, c_e, label="Exposed", color="orange")
    ax3.plot(t, c_i, label="Infectious", color="red")
    ax3.plot(t, c_r, label="Recovered", color="green")

    ax3.legend()

    ax4.set_title("Transition rates for EIR compartments")
    ax4.set_xlabel("Day")
    ax4.set_ylabel("Rate")

    ax4.axhline(alpha, label="Alpha", color="red")
    ax4.plot(t, beta(np.array(t)), label="Beta", color="orange")
    ax4.axhline(gamma, label="Gamma", color="green")

    ax4.legend()

    ax5.set_title("Development of total population")
    ax5.set_xlabel("Day")
    ax5.set_ylabel("Amount")

    ax5.plot(t, n, label="Population", color="grey")

    ax5_twin = ax5.twinx()
    ax5_twin.set_ylabel("Rate")
    ax5_twin.axhline(ny, linestyle=":", label="Ny", color="blue")
    ax5_twin.axhline(my, linestyle=":", label="My", color="brown")

    ax5_twin.tick_params(axis="y", labelcolor="brown")
    ax5_twin.set_ylim(bottom=0, top=0.01)
    ax5_twin.grid(False)

    ax5.legend()

    plt.show()

seir_simulation()