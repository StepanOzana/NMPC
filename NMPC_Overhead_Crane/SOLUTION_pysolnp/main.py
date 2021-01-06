#----------------------------------------------------------
#NMPC design for an overhead crane
#(C)Stepan Ozana (2020)
#email: stepan.ozana@vsb.cz
#https://www.youtube.com/user/StepanOzana
#www.stepan-ozana.com
#(C)Filip Krupa (2020)
#email: filip.krupa@vsb.cz

#----------------------------------------------------------
import pysolnp
import matplotlib.pyplot as plt
from matplotlib import rc
rc.usetex = True
import numpy as np

from timeit import default_timer as timer
import datetime

from math import sin as msin
from math import cos as mcos
from functools import reduce
from operator import add

xg=[]

#customized functions --------------------------------------------------

def rungeKutta(x_cur,u_next,u_cur):
    h=0.25
    xj=[0,0,0,0]
    
    up=(u_next+u_cur)*0.5

    k11 = x_cur[1]
    k12 = u_cur-x_cur[1]
    k12 = 10*k12
    k13 = x_cur[3]
    k14 = 0.0004*(-4905*msin(x_cur[2])-500*x_cur[3]+5000*(u_cur-x_cur[1])*mcos(x_cur[2]))

    x2p=x_cur[1]+k12*h*0.5
    x3p=x_cur[2]+k13*h*0.5
    x4p=x_cur[3]+k14*h*0.5
   
    k21 =x2p
    k22 =up-x2p
    k22=10*k22
    k23 =x4p
    k24=0.0004*(-4905*msin(x3p)-500*x4p+5000*(up-x2p)*mcos(x3p))
    
    x2q=x_cur[1]+k22*h*0.5
    x3q=x_cur[2]+k23*h*0.5
    x4q=x_cur[3]+k24*h*0.5
    
    k31 =x2q
    k32 =up-x2q
    k32=10*k32
    k33 =x4q
    k34=0.0004*(-4905*msin(x3q)-500*x4q+5000*(up-x2q)*mcos(x3q))
    
    x2e=x_cur[1]+k32*h
    x3e=x_cur[2]+k33*h
    x4e=x_cur[3]+k34*h
            
    k41 =x2e
    k42 =u_next-x2e
    k42=10*k42
    k43 =x4e
    k44=0.0004*(-4905*msin(x3e)-500*x4e+5000*(u_next-x2e)*mcos(x3e))


    xj[0]=x_cur[0]+h*(k11+2.0*k21+2.0*k31+k41)/6.0
    xj[1]=x_cur[1]+h*(k12+2.0*k22+2.0*k32+k42)/6.0
    xj[2]=x_cur[2]+h*(k13+2.0*k23+2.0*k33+k43)/6.0
    xj[3]=x_cur[3]+h*(k14+2.0*k24+2.0*k34+k44)/6.0
    return xj

#customized functions --------------------------------------------------

def custom_J(x,h):
    out=(x-h)*(x-h)
    return out

def custom_Sum(x,y):
    out=x+y
    return out   

def custom_Sub(x,y):    
    out=x-y
    return out

def custom_Mult(x,y):
    out=x*y
    return out

def custom_Reshape(h1,h2,h3,h4):
    out=[h1,h2,h3,h4]
    return out

def custom_SubL4(x,y):
    out=[x[0]-y[0],x[1]-y[1],x[2]-y[2],x[3]-y[3]]
    return out

def box_objective_function(x):
    cas=[0,1,2,3,4,5,6,7,8,9,10]
    j=0.0
    ul=x[0:N]
    x1l=x[N:2*N]
    x2l=x[2*N:3*N]
    x3l=x[3*N:4*N]
    x4l=x[4*N:5*N]
    Q=[1,0.1,50,0.1]
    R1=1
    R2=1
    j=Q[0]*sum(list(map(custom_J,x1l,x1)))+Q[1]*sum(list(map(custom_Mult,x2l,x2l)))+Q[2]*sum(list(map(custom_Mult,cas,map(custom_Mult,x3l,x3l))))+Q[3]*sum(list(map(custom_Mult,x4l,x4l)))+R1*(sum(list(map(custom_Mult,ul,ul))))+R2*(sum(list(map(custom_Mult,map(custom_Sub,ul[0:N-1],ul[1:N]),map(custom_Sub,ul[0:N-1],ul[1:N])))))
    return j


def box_equality_function(x):
    global xg
    N=11
    ul=x[0:N]
    x1l=x[N:2*N]
    x2l=x[2*N:3*N]
    x3l=x[3*N:4*N]
    x4l=x[4*N:5*N]
    c=[0.0]*45

    #RK4:
    rk_out=list(map(rungeKutta,map(custom_Reshape,x1l[0:N-1],x2l[0:N-1],x3l[0:N-1],x4l[0:N-1]),ul[1:N],ul[0:N-1]))
    c[0:40]=reduce(add,map(custom_SubL4,list(map(custom_Reshape,x1l[1:N],x2l[1:N],x3l[1:N],x4l[1:N])),rk_out))
    c[-5]=(x1l[-1]-x1[0])*(x1l[-1]-x1[0])
    c[-4]=x1l[0]-xg[0]
    c[-3]=x2l[0]-xg[1]
    c[-2]=x3l[0]-xg[2]
    c[-1]=x4l[0]-xg[3]
    return c


def solve_box():
    result = pysolnp.solve(
        obj_func=box_objective_function,
        par_start_value=z_init,
        par_lower_limit=z_lower,
        par_upper_limit=z_upper,
        eq_func=box_equality_function,
        eq_values=equality_values,
        max_major_iter=10,
        max_minor_iter=10,
        delta=1e-5,
        tolerance=1e-4)
        #debug=True
    return result
   
if __name__ == "__main__":
    N       =   11
    u_ub    =   1
    u_lb    =   -1

    x1_ub   =   4.0  
    x2_ub   =   10.0 
    x3_ub   =   0.5 
    x4_ub   =   10.0
    x1_lb   =   -4.0  
    x2_lb   =   -10.0
    x3_lb   =   -0.5
    x4_lb   =   -10.0

    x1=[]
    x2=[]
    x3=[]
    x4=[]
    
    z_init=[]
    z_upper=[]
    z_lower=[]
    equality_values=[]

    f=[0.0,0.0,0.0,0.0]

    uOK=[]
    xOK=[]

    u_act=0
    #reference
    for i in range(N):
        x1.append(2.0)  #position   
        x2.append(0.0)  #speed  
        x3.append(0.0)  #angle
        x4.append(0.0)  #angular speed 
    xg.append(0.0)
    xg.append(0.0)
    xg.append(0.0)
    xg.append(0.0)

    for i in range(N):
        z_init.append(0.0)
        z_lower.append(u_lb)
        z_upper.append(u_ub)

        equality_values.append(0.0)
        equality_values.append(0.0)
        equality_values.append(0.0)
        equality_values.append(0.0)

    equality_values.append(0.0)

    for i in range(N):
        z_init.append(0.0)
        z_lower.append(x1_lb)
        z_upper.append(x1_ub)

    for i in range(N):
        z_init.append(0.0)
        z_lower.append(x2_lb)
        z_upper.append(x2_ub)

    for i in range(N):
        z_init.append(0.0)
        z_lower.append(x3_lb)
        z_upper.append(x3_ub)

    for i in range(N):
        z_init.append(0.0)
        z_lower.append(x4_lb)
        z_upper.append(x4_ub)

    for smycka in range(101):
        print(smycka)

        if smycka==0:
            xg=[0.0,0.0,0.0,0.0]
        else:
            xg = rungeKutta(xg,u_act,u_act)
        
        result = solve_box()
        u_act=result.optimum[0]
        xOK.append(xg)
        uOK.append(result.optimum[0])
    #PLOTTING RESULTS----------------------
    cas=np.linspace(0.0, 100*0.25, num=101)
    plt.figure(1)
    plt.plot(cas,uOK)
    plt.title('Control variable (speed setpoint for the overhead crane)')
    plt.grid(True)
    plt.xlabel('time [s]')
    plt.ylabel('u [m/s]')
    
    plt.figure(2)
    plt.plot(cas,xOK)
    plt.title('State variables')
    plt.grid(True)
    plt.xlabel('time [s]')
    plt.ylabel('$x_1$,$x_2$,$x_3$,$x_4$')
    plt.legend(["$x_1 [m]$", "$x_2 [m/s]$","$x_3 [rad]$", "$x_4 [rad/s]$"], loc ="upper right")
    plt.show()





    
       
