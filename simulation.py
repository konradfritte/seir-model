import numpy as np

import model
import solver

# TODO:: Create different simulations with different parameters

def measles(duration=365):
    healthcare_limit = 100000

    alpha = model.alpha_moderator(1 / 10)
    beta = model.beta_moderator(4 / 5)
    gamma = model.gamma_moderator(1 / 5.5)
    delta = model.delta_moderator(0)
    epsilon = model.epsilon_moderator(1 / 1000000)
    epsilon_h = model.epsilon_moderator(5/1000, limit=healthcare_limit, magnitude=5)
    zeta = model.zeta_moderator(2/1000)
    eta = model.eta_moderator(1/10)

    ny=1/365/80
    my=1/365/80

    f = model.seihrsd_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, epsilon_h=epsilon_h, zeta=zeta, eta=eta, ny=ny, my=my)

    n = 83200000
    s = n - 10000
    e = n - s
    ih = i = r = d = 0
    c_s = s
    c_e = e
    c_i = c_ih = c_r = c_d = 0

    x0 = np.array([n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d])

    results = solver.euler_method(f, x0, duration)

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

def influenza(duration=365):
    healthcare_limit = 100000

    alpha = model.alpha_moderator(1 / 2)
    beta = model.beta_moderator(1.3 / 5, amplitude=0.5)
    gamma = model.gamma_moderator(1 / 5)
    delta = model.delta_moderator(1/365)
    epsilon = model.epsilon_moderator(1 / 100000)
    epsilon_h = model.epsilon_moderator(1/1000, limit=healthcare_limit, magnitude=5)
    zeta = model.zeta_moderator(1/1000)
    eta = model.eta_moderator(1/7)

    ny=1/365/80
    my=1/365/80

    f = model.seihrsd_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, epsilon_h=epsilon_h, zeta=zeta, eta=eta, ny=ny, my=my)

    n = 83200000
    s = n - 10000
    e = n - s
    ih = i = r = d = 0
    c_s = s
    c_e = e
    c_i = c_ih = c_r = c_d = 0

    x0 = np.array([n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d])

    results = solver.euler_method(f, x0, duration)

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

def covid(duration=365):
    healthcare_limit = 100000

    alpha = model.alpha_moderator(1 / 4)
    beta = model.beta_moderator(2 / 8.5, amplitude=0.1)
    gamma = model.gamma_moderator(1 / 8.5)
    delta = model.delta_moderator(1/365)
    epsilon = model.epsilon_moderator(1 / 100000)
    epsilon_h = model.epsilon_moderator(1/100, limit=healthcare_limit, magnitude=5)
    zeta = model.zeta_moderator(1/100)
    eta = model.eta_moderator(1/12)

    ny=1/365/80
    my=1/365/80

    f = model.seihrsd_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, epsilon_h=epsilon_h, zeta=zeta, eta=eta, ny=ny, my=my)

    n = 83200000
    s = n - 10000
    e = n - s
    ih = i = r = d = 0
    c_s = s
    c_e = e
    c_i = c_ih = c_r = c_d = 0

    x0 = np.array([n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d])

    results = solver.euler_method(f, x0, duration)

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

def seihrsd(duration=365):
    healthcare_limit = 100000

    alpha = model.alpha_moderator(1 / 5)
    beta = model.beta_moderator(2 / 10)
    gamma = model.gamma_moderator(1 / 10)
    delta = model.delta_moderator(1/365)
    epsilon = model.epsilon_moderator(0)
    epsilon_h = model.epsilon_moderator(1/100, limit=healthcare_limit, magnitude=5)
    zeta = model.zeta_moderator(1/1000)
    eta = model.eta_moderator(1/20)

    ny=0
    my=0

    f = model.seihrsd_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, epsilon_h=epsilon_h, zeta=zeta, eta=eta, ny=ny, my=my)

    n = 83200000
    s = n - 10000
    e = n - s
    ih = i = r = d = 0
    c_s = s
    c_e = e
    c_i = c_ih = c_r = c_d = 0

    x0 = np.array([n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d])

    results = solver.euler_method(f, x0, duration)

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

def seirsd_with_seasonality(duration=5*365):
    alpha = model.alpha_moderator(1 / 5)
    beta = model.beta_moderator(2 / 10, amplitude=0.2)
    gamma = model.gamma_moderator(1 / 10)
    delta = model.delta_moderator(1 / 365)
    epsilon = model.epsilon_moderator(1 / 1000000)
    epsilon_h = model.epsilon_moderator(0)
    zeta = model.zeta_moderator(0)
    eta = model.eta_moderator(0)

    ny = 0
    my = 0

    f = model.seihrsd_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, epsilon_h=epsilon_h, zeta=zeta, eta=eta, ny=ny, my=my)

    n = 83200000
    s = n - 10000
    e = n - s
    ih = i = r = d = 0
    c_s = s
    c_e = e
    c_i = c_ih = c_r = c_d = 0

    x0 = np.array([n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d])

    results = solver.euler_method(f, x0, duration)

    return {
        "parameters": {
            "healthcare_limit": 0,
            "alpha": alpha,
            "gamma": gamma,
            "beta": beta,
            "delta": delta,
            "epsilon": epsilon,
            "zeta": model.zeta_moderator(0),
            "eta": model.zeta_moderator(0),
            "ny": ny,
            "my": my
        },
        "results": results
    }

def seirsd(duration=5*365):
    alpha = model.alpha_moderator(1 / 5)
    beta = model.beta_moderator(2 / 10)
    gamma = model.gamma_moderator(1 / 10)
    delta = model.delta_moderator(1 / 365)
    epsilon = model.epsilon_moderator(1 / 1000000)
    epsilon_h = model.epsilon_moderator(0)
    zeta = model.zeta_moderator(0)
    eta = model.eta_moderator(0)

    ny = 0
    my = 0

    f = model.seihrsd_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, epsilon_h=epsilon_h, zeta=zeta, eta=eta, ny=ny, my=my)

    n = 83200000
    s = n - 10000
    e = n - s
    ih = i = r = d = 0
    c_s = s
    c_e = e
    c_i = c_ih = c_r = c_d = 0

    x0 = np.array([n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d])

    results = solver.euler_method(f, x0, duration)

    return {
        "parameters": {
            "healthcare_limit": 0,
            "alpha": alpha,
            "gamma": gamma,
            "beta": beta,
            "delta": delta,
            "epsilon": epsilon,
            "zeta": model.zeta_moderator(0),
            "eta": model.zeta_moderator(0),
            "ny": ny,
            "my": my
        },
        "results": results
    }

def seird(duration=365):
    alpha = model.alpha_moderator(1 / 5)
    beta = model.beta_moderator(2 / 10)
    gamma = model.gamma_moderator(1 / 10)
    delta = model.delta_moderator(0)
    epsilon = model.epsilon_moderator(1 / 1000000)
    epsilon_h = model.epsilon_moderator(0)
    zeta = model.zeta_moderator(0)
    eta = model.eta_moderator(0)

    ny = 0
    my = 0

    f = model.seihrsd_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, epsilon_h=epsilon_h, zeta=zeta, eta=eta, ny=ny, my=my)

    n = 83200000
    s = n - 10000
    e = n - s
    ih = i = r = d = 0
    c_s = s
    c_e = e
    c_i = c_ih = c_r = c_d = 0

    x0 = np.array([n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d])

    results = solver.euler_method(f, x0, duration)

    return {
        "parameters": {
            "healthcare_limit": 0,
            "alpha": alpha,
            "gamma": gamma,
            "beta": beta,
            "delta": delta,
            "epsilon": epsilon,
            "zeta": model.zeta_moderator(0),
            "eta": model.zeta_moderator(0),
            "ny": ny,
            "my": my
        },
        "results": results
    }

def seir_with_seasonality(duration=5*365):
    alpha = model.alpha_moderator(1 / 5)
    beta = model.beta_moderator(2 / 10, amplitude=0.2, phi=-365/4)
    gamma = model.gamma_moderator(1 / 10)
    delta = model.delta_moderator(0)
    epsilon = model.epsilon_moderator(0)
    epsilon_h = model.epsilon_moderator(0)
    zeta = model.zeta_moderator(0)
    eta = model.eta_moderator(0)

    ny = 0
    my = 0

    f = model.seihrsd_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, epsilon_h=epsilon_h, zeta=zeta, eta=eta, ny=ny, my=my)

    n = 83200000
    s = n - 10000
    e = n - s
    ih = i = r = d = 0
    c_s = s
    c_e = e
    c_i = c_ih = c_r = c_d = 0

    x0 = np.array([n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d])

    results = solver.euler_method(f, x0, duration)

    return {
        "parameters": {
            "healthcare_limit": 0,
            "alpha": alpha,
            "gamma": gamma,
            "beta": beta,
            "delta": delta,
            "epsilon": epsilon,
            "zeta": model.zeta_moderator(0),
            "eta": model.zeta_moderator(0),
            "ny": ny,
            "my": my
        },
        "results": results
    }

def seir_with_vitality_parameters(duration=365*40):
    alpha = model.alpha_moderator(1 / 5)
    beta = model.beta_moderator(2 / 10)
    gamma = model.gamma_moderator(1 / 10)
    delta = model.delta_moderator(0)
    epsilon = model.epsilon_moderator(0)
    epsilon_h = model.epsilon_moderator(0)
    zeta = model.zeta_moderator(0)
    eta = model.eta_moderator(0)

    ny = 1/365/40
    my = 1/365/80

    f = model.seihrsd_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, epsilon_h=epsilon_h, zeta=zeta, eta=eta, ny=ny, my=my)

    n = 83200000
    s = n - 10000
    e = n - s
    ih = i = r = d = 0
    c_s = s
    c_e = e
    c_i = c_ih = c_r = c_d = 0

    x0 = np.array([n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d])

    results = solver.euler_method(f, x0, duration)

    return {
        "parameters": {
            "healthcare_limit": 0,
            "alpha": alpha,
            "gamma": gamma,
            "beta": beta,
            "delta": delta,
            "epsilon": epsilon,
            "zeta": model.zeta_moderator(0),
            "eta": model.zeta_moderator(0),
            "ny": ny,
            "my": my
        },
        "results": results
    }

def seir(duration=365):
    alpha = model.alpha_moderator(1 / 5)
    beta = model.beta_moderator(2 / 10)
    gamma = model.gamma_moderator(1 / 10)
    delta = model.delta_moderator(0)
    epsilon = model.epsilon_moderator(0)
    epsilon_h = model.epsilon_moderator(0)
    zeta = model.zeta_moderator(0)
    eta = model.eta_moderator(0)

    ny = 0
    my = 0

    f = model.seihrsd_model(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, epsilon_h=epsilon_h, zeta=zeta, eta=eta, ny=ny, my=my)

    n = 83200000
    s = n - 10000
    e = n - s
    ih = i = r = d = 0
    c_s = s
    c_e = e
    c_i = c_ih = c_r = c_d = 0

    x0 = np.array([n, s, e, i, ih, r, d, c_s, c_e, c_i, c_ih, c_r, c_d])

    results = solver.euler_method(f, x0, duration)

    return {
        "parameters": {
            "healthcare_limit": 0,
            "alpha": alpha,
            "gamma": gamma,
            "beta": beta,
            "delta": delta,
            "epsilon": epsilon,
            "zeta": model.zeta_moderator(0),
            "eta": model.zeta_moderator(0),
            "ny": ny,
            "my": my
        },
        "results": results
    }