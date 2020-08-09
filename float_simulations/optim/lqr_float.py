from autolib import *

def draw_buoy(x):
    clear(ax)
    x=x.flatten()
    plot([-10,10],[0,0],'black',linewidth=1)    
    d=-x[0]
    P=array([[-ech,-1.8*ech],[ech,-1.8*ech],[ech,0],[-ech,0]])
    draw_polygon(ax,P,'blue')
    plot([   0,   L,  L,  L/2,   L/2,   L/2,  0,  0],
         [-L-d,-L-d, -d,   -d,   2-d,    -d, -d,-L-d],'black',linewidth=3)
    plot([-ech, ech],[zbar,zbar],'red',linewidth=1)
    b=-x[2]     
    P=array([[0,-L-d+L],[L,-L-d+L],[L,-L/2-L*b/2-d],[0,-L/2-L*b/2-d]])
    draw_polygon(ax,P,'white')
    
def lqr(A,B,Q,R):
    X = solve_continuous_are(A,B,Q,R)
    K = inv(R)@(B.T@X)
    return K

    
def f(X,u):
    return A@X+u*B
    
dt = 0.05
kmax=200

ech=5
p1,p2,p3=0.1,1,2
A=array([[0,1,0],[p1,-p2,p3],[0,0,0]])
B=array([[0],[0],[1]])
K=lqr(A,B,diag([1,1,0]),diag([1]))
zbar=-5
xbar=array([[zbar],[0],[-(p1/p3)*zbar]])
x = array([[-2],[0],[0]])
L=1
ax=init_figure(-ech,ech,-1.8*ech,0.2*ech)

U,Z=zeros(kmax),zeros(kmax)
for k in range(0,kmax):
    u = (-K @ (x-xbar))[0,0]
    U[k]=u
    Z[k]=x[0,0]
    x=x+dt*f(x,u)
    draw_buoy(x)
pause(1)

tmax=dt*kmax
ax=init_figure(0,tmax,-5,5)
T=arange(0,tmax,dt)
plot(T, U, 'red')
plot(T, Z, 'blue')
pause(10)

