# Approximates a function's values using the euler method: f(t + h) = f(t) + h * f'(t).
# f: The function that should be approximated
# x0: The initial value at which to start
# t: The final time point at which to end
# h: The delta value at which the function should be approximated
def euler_method(f, x0, t, h=1):
    x = x0
    result = [[0, x]]

    for k in range(1, int(1 + t / h)):
        t = k * h
        x = x + h * f(t, x)

        result.append([t, x])
    return result