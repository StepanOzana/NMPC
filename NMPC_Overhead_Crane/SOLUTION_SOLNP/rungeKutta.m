function x_next = rungeKutta(x_current,u_next,u_current,h)
g=9.81;
m=100;
l=5;
tau=0.1;
I=m*l^2;
b=500;


up=(u_next+u_current)/2;

        k11 = x_current(2);
        k12 = 1/tau*u_current-1/tau*x_current(2);
        k13 = x_current(4);
        k14 = 1/I*(-m*g*l*sin(x_current(3))-b*x_current(4)+m*l*(1/tau*u_current-1/tau*x_current(2))*cos(x_current(3)));
		
        x1p=x_current(1)+k11*h/2.0;
		x2p=x_current(2)+k12*h/2.0;
		x3p=x_current(3)+k13*h/2.0;
        x4p=x_current(4)+k14*h/2.0;

        
        
        k21 =x2p;
        k22 =1/tau*up-1/tau*x2p;
        k23 =x4p;
        k24=1/I*(-m*g*l*sin(x3p)-b*x4p+m*l*(1/tau*up-1/tau*x2p)*cos(x3p));
        
        		
		x1q=x_current(1)+k21*h/2.0;
		x2q=x_current(2)+k22*h/2.0;
		x3q=x_current(3)+k23*h/2.0;
        x4q=x_current(4)+k24*h/2.0;
		
	    k31 =x2q;
        k32 =1/tau*up-1/tau*x2q;
        k33 =x4q;
        k34=1/I*(-m*g*l*sin(x3q)-b*x4q+m*l*(1/tau*up-1/tau*x2q)*cos(x3q));
		
	    x1e=x_current(1)+k31*h;
		x2e=x_current(2)+k32*h;
		x3e=x_current(3)+k33*h;
        x4e=x_current(4)+k34*h;
				
	    k41 =x2e;
        k42 =1/tau*u_next-1/tau*x2e;
        k43 =x4e;
        k44=1/I*(-m*g*l*sin(x3e)-b*x4e+m*l*(1/tau*u_next-1/tau*x2e)*cos(x3e));
	
		
		x1 = x_current(1)+h*(k11+2.0*k21+2.0*k31+k41)/6.0;
		x2 = x_current(2)+h*(k12+2.0*k22+2.0*k32+k42)/6.0;
		x3 = x_current(3)+h*(k13+2.0*k23+2.0*k33+k43)/6.0;
        x4 = x_current(4)+h*(k14+2.0*k24+2.0*k34+k44)/6.0;

       x_next=[x1;x2;x3;x4];
       
end