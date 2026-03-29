#!/usr/bin/env python3
"""interpolate - Lagrange, Newton, cubic spline, Chebyshev interpolation."""
import sys, argparse, math, json

def lagrange(xs, ys, x):
    n = len(xs); result = 0
    for i in range(n):
        term = ys[i]
        for j in range(n):
            if i != j: term *= (x - xs[j]) / (xs[i] - xs[j])
        result += term
    return result

def newton_divided_diff(xs, ys):
    n = len(xs); table = [list(ys)]
    for j in range(1, n):
        row = []
        for i in range(n - j):
            row.append((table[j-1][i+1] - table[j-1][i]) / (xs[i+j] - xs[i]))
        table.append(row)
    return [table[j][0] for j in range(n)]

def newton_eval(xs, coeffs, x):
    n = len(coeffs); result = coeffs[-1]
    for i in range(n-2, -1, -1):
        result = result * (x - xs[i]) + coeffs[i]
    return result

def cubic_spline(xs, ys):
    n = len(xs) - 1
    h = [xs[i+1] - xs[i] for i in range(n)]
    alpha = [0] + [3/h[i]*(ys[i+1]-ys[i]) - 3/h[i-1]*(ys[i]-ys[i-1]) for i in range(1, n)]
    l = [1] + [0]*n; mu = [0]*(n+1); z = [0]*(n+1)
    for i in range(1, n):
        l[i] = 2*(xs[i+1]-xs[i-1]) - h[i-1]*mu[i-1]
        mu[i] = h[i]/l[i]; z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i]
    l[n] = 1
    c = [0]*(n+1); b = [0]*n; d = [0]*n
    for j in range(n-1, -1, -1):
        c[j] = z[j] - mu[j]*c[j+1]
        b[j] = (ys[j+1]-ys[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3
        d[j] = (c[j+1]-c[j])/(3*h[j])
    def evaluate(x):
        for i in range(n):
            if xs[i] <= x <= xs[i+1]:
                dx = x - xs[i]
                return ys[i] + b[i]*dx + c[i]*dx**2 + d[i]*dx**3
        return None
    return evaluate

def chebyshev_nodes(a, b, n):
    return [(a+b)/2 + (b-a)/2 * math.cos((2*k+1)*math.pi/(2*n)) for k in range(n)]

def main():
    p = argparse.ArgumentParser(description="Interpolation methods")
    p.add_argument("--demo", action="store_true")
    args = p.parse_args()
    if args.demo:
        xs = [0, 1, 2, 3, 4]; ys = [1, 2.7, 7.4, 20.1, 54.6]
        print("=== Lagrange ===")
        for x in [0.5, 1.5, 2.5, 3.5]:
            print(f"  f({x}) = {lagrange(xs, ys, x):.4f}")
        print("\n=== Newton ===")
        coeffs = newton_divided_diff(xs, ys)
        for x in [0.5, 1.5, 2.5, 3.5]:
            print(f"  f({x}) = {newton_eval(xs, coeffs, x):.4f}")
        print("\n=== Cubic Spline ===")
        spline = cubic_spline(xs, ys)
        for x in [0.5, 1.5, 2.5, 3.5]:
            print(f"  f({x}) = {spline(x):.4f}")
        print("\n=== Chebyshev nodes (avoid Runge) ===")
        nodes = chebyshev_nodes(-5, 5, 11)
        runge = lambda x: 1/(1+x**2)
        cn_ys = [runge(x) for x in nodes]
        uniform = [i-5 for i in range(11)]
        un_ys = [runge(x) for x in uniform]
        test_x = 4.5
        print(f"  Chebyshev at {test_x}: {lagrange(nodes, cn_ys, test_x):.6f}")
        print(f"  Uniform at {test_x}: {lagrange(uniform, un_ys, test_x):.6f}")
        print(f"  Exact: {runge(test_x):.6f}")
    else: p.print_help()

if __name__ == "__main__":
    main()
