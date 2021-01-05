function [cost] = costF(z,par)
global x0 h x1ref  x2ref x3ref x4ref N N1
% N=30;
% t0=0; % Initial time
% tf=30; % Final time
% N1 = N+1; % Number of grid points
%h = (tf-t0)/N; % Integration step
u = z(1:N1); % Control input
% States
x1=z(N1+1:2*N1);
x2=z(2*N1+1:3*N1);
x3=z(3*N1+1:4*N1);
x4=z(4*N1+1:5*N1);
%x=[x1 x2 x3 x4];
Q=1*[1;0.1;50;0.1];
R1=1;
R2=1;

cas=(0:length(x3)-1)';

%fobj = Q(1)*sum((x1-x1ref).^2)+0.01*sum(x3.^2)+R*sum(abs(diff(u)));
%fobj = Q(1)*sum((x1-x1ref).^2)+0.01*sum(x3.^2)+R*sum(abs(diff(u)));
%fobj = 1*sum((x1-x1ref).^2)+100*sum(cas.*x3.^2)+1*sum(abs(diff(u)));

%fobj = 5*(abs(max(x1)-x1ref))^2+Q(1)*sum((x1-x1ref).^2)+Q(2)*sum(x2.^2)+Q(3)*sum(x3.^2)+Q(4)*sum(x4.^2)+R*sqrt(sum(u.^2))+M*sqrt(sum((diff(u)).^2));
%fobj = 25*(abs(max(x1)-x1ref))^2+Q(1)*sum((x1-x1ref).^2)+Q(2)*sum(x2.^2)+Q(3)*sum(cas.*x3.^2)+Q(4)*sum(x4.^2)+R*sqrt(sum(u.^2))+M*sqrt(sum((diff(u)).^2));
%fobj = 5*(x1(end)-x1ref)^2+Q(1)*sum((x1-x1ref).^2)+Q(2)*sum(x2.^2)+Q(3)*sum(cas.*x3.^2)+Q(4)*sum(x4.^2)+R*sqrt(sum(u.^2))+M*sqrt(sum((diff(u)).^2));
%fobj = 100*(x1(end)-x1ref)^2+Q(1)*sum((x1-x1ref).^2)+Q(2)*sum(x2.^2)+Q(3)*sum(cas.*x3.^2)+Q(4)*sum(x4.^2)+R*sqrt(sum(u.^2))+M*sqrt(sum((diff(u)).^2));
%fobj = Q(1)*sum((x1-x1ref).^2)+Q(2)*sum((x2-x2ref).^2)+Q(3)*sum((x3-x3ref).^2)+Q(4)*sum((x4-x4ref).^2)+R*sqrt(sum(u.^2))+M*sqrt(sum((diff(u)).^2));
%fobj=0.01*sum(u.^2);

%fobj = Q(1)*sum((x1-x1ref).^2)+Q(2)*sum(x2.^2)+Q(3)*sum(cas.*x3.^2)+Q(4)*sum(x4.^2)+R*sqrt(sum(u.^2))+M*sqrt(sum((diff(u)).^2));

fobj = Q(1)*sum((x1-x1ref).^2)+Q(2)*sum(x2.^2)+Q(3)*sum(cas.*x3.^2)+Q(4)*sum(x4.^2)+R1*sum(u.^2)+R2*sum((diff(u)).^2);

%fobj =(x1(end)-x1ref)^2+x2(end)^2+x3(end)^2+x4(end)^2+0.01*sum(u.^2);   %N=15, E(5) disabled, tolerance a delta 1e-3


%l=5;
%O=[z(1,:) 0];
%P=[O(1)-l*sin(z(3,:)) O(2)-l*cos(z(3,:))];
%fobj = Q(1)*sum((x1-x1ref).^2)+(P(1)-x1ref)^2;

% Differential defect constraints...
Z = [];
E = zeros(4,1);
E(1) = x1(1)- x0(1);
E(2) = x2(1)- x0(2);
E(3) = x3(1)- x0(3);
E(4) = x4(1)- x0(4);
E(5) = (x1(end)-x1ref)^2;


Z = E;

for k=1:N1-1
xcurrent = [x1(k), x2(k), x3(k), x4(k)]';
xnext = [x1(k+1),x2(k+1),x3(k+1),x4(k+1)]';
unext=u(k+1);
ucurrent=u(k);

 condyn = xnext-rungeKutta(xcurrent,unext,ucurrent,h);


Z = [Z; condyn];
end;
cost=[fobj; Z];






