#available at https://www.ensta-bretagne.fr/jaulin/autolib.py 
# used in AutoMOOC :  https://www.ensta-bretagne.fr/jaulin/automooc.html


import scipy.linalg
import matplotlib.pyplot as plt
from numpy import mean,pi,cos,sin,sqrt,tan,arctan,arctan2,tanh,arcsin,arccos,\
                    exp,dot,array,log,inf, eye, zeros, ones, inf,size,\
                    arange,reshape,vstack,hstack,diag,median,\
                    sign,sum,meshgrid,cross,linspace,append,round,trace,floor
from matplotlib.pyplot import *
from numpy.random import randn,rand
from numpy.linalg import inv, det, norm, eig,qr
from scipy.linalg import sqrtm,expm,logm,norm,block_diag,solve_continuous_are
from scipy.signal import place_poles

from mpl_toolkits.mplot3d import Axes3D
from math import factorial
from matplotlib.patches import Ellipse,Rectangle,Circle,Wedge,Polygon, Arc
from matplotlib.collections import PatchCollection

#from control.matlab import tf2ss,tf,bode # Not compatible with Pygame
# Unicode https://en.wikipedia.org/wiki/List_of_Unicode_characters
# αβδεθλΛμρτφψωΓ

   
    
def angle(x):
    x=x.flatten()
    return arctan2(x[1],x[0])

def tran2H(x,y):
    return array([[1,0,x],[0,1,y],[0,0,1]])


def rot2H(a):
    return array([[cos(a),-sin(a),0],[sin(a),cos(a),0],[0,0,1]])
    
    
def plot2D(M,col='black',w=1):
    plot(M[0, :], M[1, :], col, linewidth = w)         
        

def draw_segment(a,b,col='darkblue',w=1):
    plot2D(hstack((a,b)),col, w)
    #plot2D(a,'ro')
    #plot2D(b,'ro')      
  

def draw_disk(ax,c,r,col,alph=0.7,w=1):     #draw_disk(array([[1],[2]]),0.5,ax,"blue")
    e = Ellipse(xy=c, width=2*r, height=2*r, angle=0,linewidth = w)   
    ax.add_artist(e)
    e.set_clip_box(ax.bbox)
    e.set_alpha(alph)  # transparency
    e.set_facecolor(col)
        
def draw_box(ax,x1,x2,y1,y2,col): 
    c=array([[x1],[y1]])    
    rect = Rectangle(c, width=x2-x1, height=y2-y1, angle=0)
    rect.set_facecolor(array([0.4,0.3,0.6]))   
    ax.add_patch(rect)
    rect.set_clip_box(ax.bbox)
    rect.set_alpha(0.7)
    rect.set_facecolor(col)    

def draw_polygon(ax,P,col): 
    patches = []     
    patches.append(Polygon(P, True))    
    p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4, color=col)
    ax.add_collection(p)


def draw_arc(c,a,θ,col):
    s = arange(0,abs(θ),0.01)
    s = sign(θ) * s
    d = a-c
    r = norm(d)
    α = angle(d)
    w = c@ones((1,size(s)))+r*array([[cos(α),-sin(α)],[sin(α), cos(α)]])@array([cos(s),sin(s)])
    plot2D(w,col,3)  
    
    
def draw_arrow(x,y,θ,L,col):
    e=0.2
    M1=L*array([[0,1,1-e,1,1-e],[0,0,-e,0,e]])
    M=np.append(M1,[[1,1,1,1,1]],axis=0)
    R=array([[cos(θ),-sin(θ),x],[sin(θ),cos(θ),y],[0,0,1]])
    plot2D(R@M,col)    
    
def draw_sailboat(x,δs,δr,ψ,awind):
    x=x.flatten()
    θ=x[2]
    hull=array([[-1,5,7,7,5,-1,-1,-1],[-2,-2,-1,1,2,2,-2,-2],[1,1,1,1,1,1,1,1]])
    sail=array([[-7,0],[0,0],[1,1]])
    rudder=array([[-1,1],[0,0],[1,1]])
    R=array([[cos(θ),-sin(θ),x[0]],[sin(θ),cos(θ),x[1]],[0,0,1]])
    Rs=array([[cos(δs),-sin(δs),3],[sin(δs),cos(δs),0],[0,0,1]])
    Rr=array([[cos(δr),-sin(δr),-1],[sin(δr),cos(δr),0],[0,0,1]])
    draw_arrow(x[0]+5,x[1],ψ,5*awind,'red')
    plot2D(R@hull,'black');       
    plot2D(R@Rs@sail,'red',2);       
    plot2D(R@Rr@rudder,'red',2);

def draw_tank(x,col='darkblue',r=1,w=2):
    mx,my,θ=list(x[0:3,0])
    M = r*array([[1,-1,0,0,-1,-1,0,0,-1,1,0,0,3,3,0], [-2,-2,-2,-1,-1,1,1,2,2,2,2,1,0.5,-0.5,-1]])
    M=add1(M)
    plot2D(tran2H(mx,my)@rot2H(θ)@M,col,w)

def draw_wheel(x,y,θ,ρ,col='darkblue',w=1):    
    wheel = zeros((3, 0))
    for i in range(18):
        k = i*pi/8
        rayon = array([[ρ*cos(k)], [ρ*sin(k)], [1]])
        wheel = hstack((wheel, rayon, array([[0], [0], [1]]), rayon))
    wheel = array([[cos(θ),-sin(θ),x], [sin(θ),cos(θ),y], [0,0,1]]) @ wheel
    plot2D(wheel, col,w)    

    	
def draw_car(x,col='darkblue',L=1,w=2): # the car has a length L
    mx,my,θ,v,δ=list(x[0:5,0])
    M = add1(L*array([[-0.3, 1.3, 1.6,1.6,1.3,-0.3,-0.3,-0.3, 0, 0,-0.3, 0.3, 0,0,-0.3,0.3, 0, 0,1,1, 1],  
                      [-0.7,-0.7,-0.3,0.3,0.7, 0.7,-0.7,-0.7,-0.7,-1,-1,-1,-1,1, 1,1, 1, 0.7,0.7,1,-1]]))                
    R=tran2H(mx,my)@rot2H(θ)
    W = add1(L*array([[-0.3, 0.3], [0, 0]])) #Front Wheel                
    plot2D(R@M,col,w)          
    plot2D(R@tran2H(L,L)@rot2H(δ)@W,col,2)
    plot2D(R@tran2H(L,-L)@rot2H(δ)@W,col,2)
    
  

def tondarray(M):
    if type(M)==float:
        return array([[M]])
    elif type(M)==int:
        return array([[M]])        
    else:
        return M    


def mvnrnd(xbar,Γ,n): 
    X=randn(2,n)
    X = (xbar @ ones((1,n))) + sqrtm(Γ) @ X
    return(X)    



def mvnrnd2(x,G): 
    n=len(x)
    x1=x.reshape(n)
    y = np.random.multivariate_normal(x1,G).reshape(n,1)
    return(y)    

def mvnrnd1(G):
    G=tondarray(G)
    n=len(G)
    x=array([[0]] * n)
    return(mvnrnd2(x,G))  
    

def place(A,B,poles):
    return place_poles(A,B,poles).gain_matrix
    
    
def RegulKLH(A,B,C,E,pcom,pobs):
    K=place(A,B,pcom)
    L=place(A.T,C.T,pobs).T
    H=-inv(E@inv(A-B@K)@B)
    Ar=A-B@K-L@C
    Br=hstack((B@H,L))
    Cr=-K
    Dr=hstack((H,0*(B.T@C.T))) 
    return Ar,Br,Cr,Dr  


def lqr(A,B,Q,R):
    X = solve_continuous_are(A,B,Q,R)
    K = inv(R)@(B.T@X) 
    eigVals, eigVecs = eig(A-B*K) 
    return K

    
def RegulKLHlqr(A,B,C,E,p):
    na=A.shape[1]
    nb=B.shape[1]
    nc=C.shape[0]
    K=lqr(A, B, -p*eye(na), -p*eye(nb))
    L=lqr(A.T,C.T,-p*eye(na), -p*eye(nc))
    L=L.T
    H=-inv(E@inv(A-B@K)@B)
    Ar=A-B@K-L@C
    Br=hstack((B@H,L))
    Cr=-K
    Dr=hstack((H,0*(B.T@C.T))) 
    return Ar,Br,Cr,Dr  
    
    
def draw_tank(x,col='darkblue',r=1,w=2):
    mx,my,θ=list(x[0:3,0])
    M = r*array([[1,-1,0,0,-1,-1,0,0,-1,1,0,0,3,3,0], [-2,-2,-2,-1,-1,1,1,2,2,2,2,1,0.5,-0.5,-1]])
    M=add1(M)
    plot2D(tran2H(mx,my)@rot2H(θ)@M,col,w)  

def draw_invpend(ax,x,col='blue'): #inverted pendulum
    s,θ=x[0,0],x[1,0]
    draw_box(ax,s-0.7,s+0.7,-0.25,0,col)
    plot( [s,s-sin(θ)],[0,cos(θ)],'magenta', linewidth = 2)

def draw_segway(ax,x,col,w):
    s,θ,ds,dθ=list(x[0:4,0])
    ρ=1 #radius of the wheel
    M = add1(array([[0,0.6,1,0,-1,-0.6,0],
               [0,0.6,5,5.4,5,0.6,0]]))
    M = tran2H(-ρ*s,ρ)@rot2H(θ)@M
    ax.plot([-10,10],[0,0],col) 
    plot2D(M,col,w)
    draw_wheel(-ρ*s,ρ,s,ρ,'darkblue',w)
    pause(0.01)

  
  
def demo_draw():  
    ax=init_figure(-15,15,-15,15)
    
    c=array([[5],[0]])
    e = Ellipse(xy=c, width=13.0, height=2.0, angle=45)  
    ax.add_artist(e)
    e.set_clip_box(ax.bbox)
    e.set_alpha(0.9)
    e.set_facecolor(array([0.7,0.3,0.6]))   
    
    rect = Rectangle( (1,1), width=5, height=3)
    rect.set_facecolor(array([0.4,0.3,0.6]))   
    ax.add_patch(rect)    
        
    pause(0.2)    
    draw_tank(array([[-7],[5],[1]]))
    draw_tank(array([[-7],[5],[1]]),'red',0.2)

    
    draw_car(array([[1],[2],[3],[4],[0.5]]),'blue',3)   
    
    P=array([[5,-3],[9,-10],[7,-4],[7,-6]])
    draw_polygon(ax,P,'green')   
    draw_disk(ax,array([[-8],[-8]]),2,"blue")   
    draw_arc(array([[0],[5]]),array([[4],[6]]),2,'red')   
    show()  # only at the end. Otherwize, it closes the figure in a terminal mode

def loadcsv(file1):
    fichier = open(file1,'r')
    D = fichier.read().split("\n")
    fichier.close()
    for i in range(len(D)):
        D[i] = D[i].split(";")
    D = array([[float(elt) for elt in Ligne] for Ligne in D])
    return D


def init_figure(xmin,xmax,ymin,ymax): 
    fig = figure()
    ax = fig.add_subplot(111, aspect='equal')	
    ax.xmin=xmin
    ax.xmax=xmax
    ax.ymin=ymin
    ax.ymax=ymax
    clear(ax)
    return ax


def clear(ax):
    pause(0.001)
    cla()
    ax.set_xlim(ax.xmin,ax.xmax)
    ax.set_ylim(ax.ymin,ax.ymax)


def draw_field(ax,f,xmin,xmax,ymin,ymax,a,b=True):
    Mx    = arange(xmin,xmax,a)
    My    = arange(ymin,ymax,a)
    X1,X2 = meshgrid(Mx,My)
    VX,VY=f(X1,X2) 
    R=sqrt(VX**2+VY**2)
    if b: quiver(Mx,My,VX/R,VY/R)
    else: quiver(Mx,My,VX,VY)


############################################################    
############################################################    
############################################################    
#  ALL 3D is made with sketch in homogenous coordinates    
############################################################    
############################################################    
############################################################    



  
def clean3D(ax,x1=-10,x2=10,y1=-10,y2=10,z1=-10,z2=10):
    ax.clear()
    ax.set_xlim3d(x1,x2)
    ax.set_ylim3d(y1,y2)
    ax.set_zlim3d(z1,z2)

def draw_arrow3D(ax,x,y,z,wx,wy,wz,col):  # initial point : x ; final point x+w 
    ax.quiver(x,y,z,wx,wy,wz,color=col,lw=1,pivot='tail',length=norm([wx,wy,wz]))

def draw_axis3D(ax,x=0,y=0,z=0,R=eye(3),zoom=1):
    ax.scatter(x,y,z,color='magenta')
    R=zoom*R
    draw_arrow3D(ax,x,y,z,R[0,0],R[1,0],R[2,0],"red")
    draw_arrow3D(ax,x,y,z,R[0,1],R[1,1],R[2,1],"green")
    draw_arrow3D(ax,x,y,z,R[0,2],R[1,2],R[2,2],"blue")
    
    
  
def ToH(R): #transformation matrix to homogenous
    H=hstack((R,array([[0],[0],[0]])))
    V=vstack((H,array([0,0,0,1])))
    return V


def tran3H(x,y,z):
    return array([[1,0,0,x],[0,1,0,y],[0,0,1,z],[0,0,0,1]])


def rot3H(wx,wy,wz):
    return ToH(expm(adjoint(array([[wx],[wy],[wz]]))))



def eulerH(φ,θ,ψ):
    Ad_i = adjoint(array([1,0,0]))
    Ad_j = adjoint(array([0,1,0]))
    Ad_k = adjoint(array([0,0,1]))
    M=expm(ψ*Ad_k) @ expm(θ*Ad_j) @ expm(φ*Ad_i)
    return ToH(M)    

    
def adjoint(w):    
    w=w.flatten()
    return array([[0,-w[2],w[1]] , [w[2],0,-w[0]] , [-w[1],w[0],0]])
    

    
def draw3H(ax,M,col,shadow=False,mirror=1):                         #mirror=-1 in case z in directed downward
    ax.plot(mirror*M[0],M[1],mirror*M[2],color=col)
    if shadow: ax.plot(mirror*M[0],M[1],0*M[2],color='black')    
    

def add1(M):
    return vstack((M,ones(M.shape[1])))
    
def wheel3H(r):
    n = 12
    W0=[[0,0,0],[r,0,r],[0,0,0],[1,1,1]]
    W=W0
    R=rot3H(2*pi/n,0,0)
    for i in range(n+1):
        W0=R@W0
        W=hstack((W,W0))
    return W    

    
def circle3H(r):
    n = 10
    θ = linspace(0, 2*pi, n)
    x = r*cos(θ)+array(n*[0])
    y = r*sin(θ) + array(n*[0])
    z = zeros(n)
    return add1(array([x,y,z]))

    

def draw_segway3D(ax,x,y,θ,ψ,s1,s2):
    p=array([[x],[y],[0]])
    ρ=1 #radius of the wheel
    Rψ = rot3H(0,0,ψ)
    Mbody = array([[0,0,0,0,0],[-1,-0.1,0,0.1,1],[0,0,5,0,0],[1,1,1,1,1]])
    Rp=tran3H(x,y,ρ)
    draw3H(ax,Rp@rot3H(0,θ,ψ)@rot3H(0,θ,0)@Mbody,'red',False)
    W=wheel3H(ρ)
    W1=rot3H(0,0,pi/2)@rot3H(s1,0,0)@W
    W2=rot3H(0,0,pi/2)@rot3H(s2,0,0)@W
    T1=tran3H(0,1,0)
    T2=tran3H(0,-1,0)
    draw3H(ax,Rp@Rψ@T1@W1,'blue',True)
    draw3H(ax,Rp@Rψ@T2@W2,'blue',True)
    pause(0.01)


    
    

def demo_animation():    
    ax=init_figure(-15,15,-15,15)
    for t in arange(0,1,0.1) :
        clear(ax)
        draw_car(array([[t],[2],[3+t],[4],[5+t]]),'blue',3)    
#        if (t>50)&(k%2000==0):
#            fig.savefig('convoy'+str(k)+'.pdf', dpi=fig.dpi)


def demo_random():  
    N=1000
    xbar = array([[1],[2]])
    Γx = array([[3,1],[1,3]])
    X=randn(2,N)
    Y=rand(2,3)
    print("Y=",Y)
    X = (xbar @ ones((1,N))) + sqrtm(Γx) @ X
    xbar_ = mean(X,axis=1)
    Xtilde = X - xbar @ ones((1,N))
    Γx_ = (Xtilde @ Xtilde.T)/N
    ax=init_figure(-20,20,-20,20)
    ax.scatter(*X)    
    pause(0.1)
    plot()  


    
def sawtooth(x):
    return (x+pi)%(2*pi)-pi   # or equivalently   2*arctan(tan(x/2))

    
def clock_RK(f,x,u,dt):  # with Runge Kutta
    return x+dt*(0.25*f(x,u)+0.75*(f(x+dt*(2/3)*f(x,u),u))) 

def clock_Euler(f,x,u,dt):  # with Euler
    return x+dt*f(x,u) 






def clock_RK(f,x,dt):  # with Runge Kutta
    return x+dt*(0.25*f(x)+0.75*(f(x+dt*(2/3)*f(x)))) 

def clock_Euler(f,x,dt):  # with Euler
    return x+dt*f(x) 





    
    
if __name__ == "__main__":


    
    demo_random()
    demo_animation()
    demo_draw() 


    #     φ,θ,ψ=list(x[3:6])
    #    s,θ,ds,dθ =list(x[0:4,0])



#
#    M=array([[1,2],[5,6],[9,10]])
#    print(M)
#    x=array([[1], [2]])    
#    x2= M@x  #multiplication dans Python 3
#    
#    print (A)
##
#    G = array([[1, 0], [0, 1]])
#    x3=mvnrnd2(x,G)
#    print("x3=",x3)
#    
#    x4=mvnrnd1(G)
#    print(x4)
#    
#    
#    print(K)
#    
#    
