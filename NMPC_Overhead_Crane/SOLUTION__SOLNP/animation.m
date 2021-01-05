
g=9.81;
m=100;
l=5;
tau=0.1;
I=m*l^2;
b=500;

close all
O=[0 0];
% P=[0 0];
% axis(gca,'equal')
axis([0 3 -6 0])
grid on
% pause(10)
for i=1:length(cas)

O=[zOK(1,i) 0];
P=[O(1)-l*sin(zOK(3,i)) O(2)-l*cos(zOK(3,i))];
O_circ=viscircles(O,0.01);
pend=line([O(1) P(1)],[O(2) P(2)]);
ball=viscircles(P,0.1);
pause(0.1)
 
    
    if (i<length(cas))
        delete(pend)
        delete(ball)
        delete(O_circ)
    end
    
end;