import matplotlib.pyplot as plt
import numpy as np

def diagram(data):
    healthcare_limit, alpha, gamma, beta, delta, epsilon, zeta, eta, ny, my = data["parameters"].values()

    t, x = zip(*data["results"])

    n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d = zip(*x)

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

    dc_s, dc_e, dc_i, dc_ih, dc_r, dc_d = np.diff([c_s, c_e, c_i, c_ih, c_r, c_d])

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
    ax6_twin.axhline(healthcare_limit, linestyle="--", label="Healthcare Limit", color="black")
    ax6_twin.grid(False)

    ax6_handles, ax6_labels = ax6.get_legend_handles_labels()
    ax6_twin_handles, ax6_twin_labels = ax6_twin.get_legend_handles_labels()

    ax6.legend(ax6_handles + ax6_twin_handles, ax6_labels + ax6_twin_labels, loc="upper right")

    plt.show()