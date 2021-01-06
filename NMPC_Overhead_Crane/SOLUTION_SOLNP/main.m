% NMPC design for an overhead crane
% (C)Stepan Ozana (2020)
% email: stepan.ozana@vsb.cz
% https://www.youtube.com/user/StepanOzana
% www.stepan-ozana.com
% (C)Filip Krupa (2020)
% email: filip.krupa@vsb.cz

clear all,close all,clc
cpk=101;
global x0 h x1ref x2ref x3ref x4ref N N1
ccount = 1;
N=10;
N1 = N+1;
x1_min = -4*ones(N1,1); x1_max = 4*ones(N1,1);
x2_min = -10*ones(N1,1); x2_max = 10*ones(N1,1);
x3_min = -0.5*ones(N1,1); x3_max = 0.5*ones(N1,1);
x4_min = -10*ones(N1,1); x4_max = 10*ones(N1,1);
u_min = -1*ones(N1,1);
u_max = 1*ones(N1,1);
zmin= [ u_min; x1_min; x2_min; x3_min; x4_min];
zmax= [ u_max; x1_max; x2_max; x3_max; x4_max];
zb_init=[zmin zmax];
xinit=[0;0;0;0]; %INITIAL STATE
h=0.25;
x1ref=2;
x2ref=0;
x3ref=0;
x4ref=0;
uOK=zeros(1,cpk);zOK=zeros(4,cpk);
tic
for m=1:cpk
m
            
    %I. MEASURE THE STATES
    if (m==1)
        x0=xinit;
    else %m>=2
        x0=rungeKutta(x0,u_current,u_current,h); %shoda se simulinkem, kdyz na vstup dame tvarovac (po dobu intervalu konstantni hodnota)
        %REAL:  x0 = MeasureStates();
    end;
    
    %II. SOLVE OPTIMIZATION PROBLEM
    %calling optimiation
    [z,oh,y,outh,outic,OK] = solnp(zb_init);
    
    %III. USE FIRST VALUE OF OPTIMAL SEQUENCE FOR THE FEEDBACK
    u_current=z(1);
    %REAL:  SendToRealPlant(u_current);
    
    %IV. STORE THE VALUES FOR POSTPROCESSING
    uOK(m)=u_current;
    zOK(1,m)=x0(1);
    zOK(2,m)=x0(2);
    zOK(3,m)=x0(3);
    zOK(4,m)=x0(4);
    
end;
toc
cas=(0:h:(cpk-1)*h);
figure,plot(cas,x1ref*ones(1,length(cas))),title('Reference position of the overhead crane'),xlabel('time [s]'),ylabel('r_x [m]'),grid on
figure,plot(cas,uOK),title('Control variable (speed setpoint for the overhead crane)'),xlabel('time [s]'),ylabel('u [m/s]'),grid on
figure,plot(cas,zOK'),title('State variables'),xlabel('time [s]'),ylabel('x_1,x_2,x_3,x_4'),legend('x_1 [m]','x_2 [m/s]','x_3 [rad]','x_4 [rad/s]'),grid on

