#!/usr/bin/env python3
"""Interpolation: linear, cubic spline, Lagrange."""
def linear(points,x):
    for i in range(len(points)-1):
        if points[i][0]<=x<=points[i+1][0]:
            t=(x-points[i][0])/(points[i+1][0]-points[i][0])
            return points[i][1]+t*(points[i+1][1]-points[i][1])
    return None
def lagrange(points,x):
    n=len(points);result=0
    for i in range(n):
        term=points[i][1]
        for j in range(n):
            if i!=j: term*=(x-points[j][0])/(points[i][0]-points[j][0])
        result+=term
    return result
def cubic_spline(points):
    n=len(points)-1;h=[points[i+1][0]-points[i][0] for i in range(n)]
    alpha=[0]+[3/h[i]*(points[i+1][1]-points[i][1])-3/h[i-1]*(points[i][1]-points[i-1][1]) for i in range(1,n)]
    l=[1]+[0]*n;mu=[0]*n;z=[0]*(n+1)
    for i in range(1,n):
        l[i]=2*(points[i+1][0]-points[i-1][0])-h[i-1]*mu[i-1]
        mu[i]=h[i]/l[i];z[i]=(alpha[i]-h[i-1]*z[i-1])/l[i]
    l.append(1)
    c=[0]*(n+1);b=[0]*n;d=[0]*n
    for j in range(n-1,-1,-1):
        c[j]=z[j]-mu[j]*c[j+1]
        b[j]=(points[j+1][1]-points[j][1])/h[j]-h[j]*(c[j+1]+2*c[j])/3
        d[j]=(c[j+1]-c[j])/(3*h[j])
    a=[points[i][1] for i in range(n)]
    def evaluate(x):
        for i in range(n):
            if points[i][0]<=x<=points[i+1][0]:
                dx=x-points[i][0]
                return a[i]+b[i]*dx+c[i]*dx**2+d[i]*dx**3
        return None
    return evaluate
if __name__=="__main__":
    import math
    pts=[(i,math.sin(i*0.5)) for i in range(8)]
    assert abs(linear(pts,2.5)-(-0.2474+1.2474)/2)<0.5
    l=lagrange(pts,2.5);print(f"Lagrange at 2.5: {l:.4f}")
    spline=cubic_spline(pts);s=spline(2.5);print(f"Spline at 2.5: {s:.4f}")
    print(f"Actual sin(1.25): {math.sin(1.25):.4f}")
    print("Interpolation OK")
