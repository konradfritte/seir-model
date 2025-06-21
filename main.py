import numpy as np
import matplotlib.pyplot as plt

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

def healthcare_limit(parameter, limit, magnitude=1):
    def f(x):
        return parameter if x <= limit else magnitude * parameter

    return f

def seir_simulation():
    alpha = 1 / 5
    gamma = 1 / 10
    beta = seasonality(1.5 / 10, amplitude=0.1)
    delta = 1/365
    epsilon = healthcare_limit(1 / 100, limit=25000, magnitude=2)
    zeta = 1 / 1000
    eta = 1 / 20
    ny = 0
    my = 0

    f = seir_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, zeta=zeta, eta=eta, ny=ny, my=my)

    n = 83200000
    s = n - 10000
    e = n - s
    ih = 0
    i = r = d = 0
    c_s = s
    c_e = e
    c_i = c_ih = c_r = c_d = 0

    x0 = np.array([n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d])

    t = 5 * 365
    h = 1

    result = euler_method(f, x0, t, h)
    t, x = zip(*result)
    n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d = zip(*x)

    dc_s = np.diff(c_s)
    dc_e = np.diff(c_e)
    dc_i = np.diff(c_i)
    dc_ih = np.diff(c_ih)
    dc_r = np.diff(c_r)
    dc_d = np.diff(c_d)

    # Move to separate file

    plt.style.use("ggplot")

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)

    ax1.set_title("Development of SEIRD compartments")
    ax1.set_xlabel("Day")
    ax1.set_ylabel("Amount")

    ax1.plot(t, s, label="Susceptible", color="blue")
    ax1.plot(t, e, label="Exposed", color="orange")
    ax1.plot(t, i, label="Infectious", color="red")
    ax1.plot(t, ih, label="Hospitalized", color="brown")
    ax1.plot(t, r, label="Recovered", color="green")
    ax1.plot(t, d, label="Dead", color="grey")

    ax2.set_title("New daily cases for SEIRD compartments")
    ax2.set_xlabel("Day")
    ax2.set_ylabel("Amount")

    ax2.plot(dc_s, label="Susceptible", color="blue")
    ax2.plot(dc_e, label="Exposed", color="orange")
    ax2.plot(dc_i, label="Infectious", color="red")
    ax2.plot(dc_ih, label="Hospitalized", color="brown")
    ax2.plot(dc_r, label="Recovered", color="green")
    ax2.plot(dc_d, label="Dead", color="grey")

    ax2.legend(loc="upper right")

    ax3.set_title("Cumulative cases for SEIRD compartments")
    ax3.set_xlabel("Day")
    ax3.set_ylabel("Amount")
    ax3.set_yscale("log")

    ax3.plot(t, c_s, label="Susceptible", color="blue")
    ax3.plot(t, c_e, label="Exposed", color="orange")
    ax3.plot(t, c_i, label="Infectious", color="red")
    ax3.plot(t, c_ih, label="Hospitalized", color="brown")
    ax3.plot(t, c_r, label="Recovered", color="green")
    ax3.plot(t, c_d, label="Dead", color="grey")

    ax3.legend(loc="upper right")

    # Basic Reproduction Number R0 with demographic parameters
    r0 = s[0]/n * alpha * beta(t[0]) / ((alpha + ny) * (gamma + zeta + ny))
    # Effective Reproduction Number Rq with demographic parameters
    rq = s / np.array(n) * alpha * beta(np.array(t)) / ((alpha + ny ) * (gamma + zeta + ny))

    ax4.set_title("Reproduction Coefficients")
    ax4.set_xlabel("Day")
    ax4.set_ylabel("Rate")

    ax4.plot(np.full(len(t), r0), label="Basic Reproduction Number", color="gray")
    ax4.plot(rq, label="Effective Reproduction Number", color="orange")
    ax4.plot(t, beta(np.array(t)), linestyle=":", label="Contact Rate", color="orange")
    ax4.axhline(1, linestyle="--", label="Balanced State", color="black")

    ax4_twin = ax4.twinx()
    ax4_twin.fill_between(t[1:], dc_e, linestyle=":", label="Daily Infections", color="grey", alpha=0.3)
    ax4_twin.grid(False)

    ax4_handles, ax4_labels = ax4.get_legend_handles_labels()
    ax4_twin_handles, ax4_twin_labels = ax4_twin.get_legend_handles_labels()

    ax4.legend(ax4_handles + ax4_twin_handles, ax4_labels + ax4_twin_labels, loc="upper right")

    ax5.set_title("Development of total population")
    ax5.set_xlabel("Day")
    ax5.set_ylabel("Amount")

    ax5.plot(t, n, label="Population", color="grey")

    ax5.legend(loc="upper right")

    ax6.set_title("Healthcare System Load")
    ax6.set_xlabel("Day")
    ax6.set_ylabel("Amount")

    ax6.fill_between(t, c_d, label="Cumulative death cases", color="grey", alpha=0.3)

    ax6_twin = ax6.twinx()
    ax6_twin.set_ylabel("Amount")
    ax6_twin.plot(t, ih, label="Hospitalized", color="brown")
    ax6_twin.axhline(25000, linestyle="--", label="Healthcare Limit", color="black")
    ax6_twin.grid(False)

    ax6_handles, ax6_labels = ax6.get_legend_handles_labels()
    ax6_twin_handles, ax6_twin_labels = ax6_twin.get_legend_handles_labels()

    ax6.legend(ax6_handles + ax6_twin_handles, ax6_labels + ax6_twin_labels, loc="upper right")

    print(np.argmax(rq))

    plt.show()


seir_simulation()