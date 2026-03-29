#!/usr/bin/env python3
"""interpolate - Lagrange, Newton, cubic spline, and linear interpolation."""
import sys, json

def linear_interp(points, x):
    points = sorted(points)
    for i in range(len(points)-1):
        x0, y0 = points[i]; x1, y1 = points[i+1]
        if x0 <= x <= x1:
            t = (x - x0)/(x1 - x0)
            return y0 + t*(y1 - y0)
    return None

def lagrange(points, x):
    n = len(points); result = 0
    for i in range(n):
        xi, yi = points[i]; basis = yi
        for j in range(n):
            if i != j:
                xj = points[j][0]
                basis *= (x - xj)/(xi - xj)
        result += basis
    return result

def newton_divided_diff(points):
    n = len(points)
    table = [[0]*n for _ in range(n)]
    for i in range(n): table[i][0] = points[i][1]
    for j in range(1, n):
        for i in range(n-j):
            table[i][j] = (table[i+1][j-1]-table[i][j-1])/(points[i+j][0]-points[i][0])
    return [table[0][j] for j in range(n)]

def newton_eval(points, coeffs, x):
    n = len(coeffs); result = coeffs[-1]
    for i in range(n-2, -1, -1):
        result = result*(x - points[i][0]) + coeffs[i]
    return result

def cubic_spline(points, x):
    n = len(points); points = sorted(points)
    h = [points[i+1][0]-points[i][0] for i in range(n-1)]
    alpha = [0]+[3*(points[i+1][1]-points[i][1])/h[i]-3*(points[i][1]-points[i-1][1])/h[i-1] for i in range(1,n-1)]
    l = [1]+[0]*(n-1); mu = [0]*n; z = [0]*n
    for i in range(1, n-1):
        l[i] = 2*(points[i+1][0]-points[i-1][0])-h[i-1]*mu[i-1]
        mu[i] = h[i]/l[i]; z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i]
    c = [0]*n; b = [0]*(n-1); d = [0]*(n-1)
    for j in range(n-2, -1, -1):
        c[j] = z[j]-mu[j]*c[j+1]
        b[j] = (points[j+1][1]-points[j][1])/h[j]-h[j]*(c[j+1]+2*c[j])/3
        d[j] = (c[j+1]-c[j])/(3*h[j])
    for i in range(n-1):
        if points[i][0] <= x <= points[i+1][0]:
            dx = x - points[i][0]
            return points[i][1]+b[i]*dx+c[i]*dx**2+d[i]*dx**3
    return None

def main():
    import math
    pts = [(0,0),(1,0.84),(2,0.91),(3,0.14),(4,-0.76),(5,-0.96)]
    print("Interpolation demo\n")
    x = 2.5
    print(f"At x={x}:")
    print(f"  Linear:   {linear_interp(pts, x):.6f}")
    print(f"  Lagrange: {lagrange(pts, x):.6f}")
    coeffs = newton_divided_diff(pts)
    print(f"  Newton:   {newton_eval(pts, coeffs, x):.6f}")
    print(f"  Spline:   {cubic_spline(pts, x):.6f}")
    print(f"  sin(2.5): {math.sin(x):.6f}")

if __name__ == "__main__":
    main()
