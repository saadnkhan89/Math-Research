clc
clear all

%% Transport Equation.

M=12;
I=-M:M;

c=20;   % For second plot use c=20

N=250;

dt=[.1 .1/2 .1/4 .1/8 .1/16 .1/32 .1/64];
T=1./dt;

x=linspace(-.5,.5,N)';
x=x(2:end-1);


for p=1:size(dt,2)
D2_f=[(1+c*dt(p)*1i*pi*I).*exp(1i*pi*x*I);exp(1i*pi*-0.5*I)];
U0=sin(x);

for z=1:T(p)
F=[U0;sin(-0.5-c*dt(p)*z)];
C=D2_f\F;
U1=real(exp(1i*pi*x*I)*C);
% U_num(:,z)=U1;
% E(z)=norm(U1-sin(x-c*dt*z),2);
U0=U1;
end


E(p)=norm(U1-sin(x-c),2);
end

p=polyfit(log2(dt),log2(E),1);

figure(1)
plot(log2(dt),log2(E),'o')
hold on
plot(log2(dt),p(1)*log2(dt)+p(2))
title('c=20','Interpreter','latex')
xlabel('$\log_2(\delta_t)$','interpreter','latex')
ylabel('$\log_2(E)$','interpreter','latex')
legend('Error','Line of best fit','interpreter','latex','location','northwest')
text(log2(dt(3)),log2(dt(3))*p(1)+p(2),'\leftarrow slope 0.8097')
set(gca,'fontsize',18)
grid on