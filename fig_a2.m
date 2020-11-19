clc 
clear all

%% Error Oscillatory behaviour Figure(4)
M=12;
I=-M:M;

c=4;

N=250;

dt=[.1 .01 .001 .0001];
T=2./dt;

x=linspace(-.5,.5,N)';
x=x(2:end-1);
M=cell(size(dt,2),1);

for p=1:size(dt,2)
D2_f=[(1+c*dt(p)*1i*pi*I).*exp(1i*pi*x*I);exp(1i*pi*-0.5*I)];
U0=sin(x);
E=zeros(T(p),1);


for z=1:T(p)
F=[U0;sin(-0.5-c*dt(p)*z)];
C=D2_f\F;
U1=real(exp(1i*pi*x*I)*C);
% U_num(:,z)=U1;
% E(z)=norm(U1-sin(x-c*dt*z),2);
E(z)=norm(U1-sin(x-c*dt(p)*z),2);

U0=U1;
end

M{p,1}=E;
end



figure(1)
hold on
plot(dt(end-3)*(1:T(end-3))',log10(M{end-3,1}),'kx')
plot(dt(end-2)*(1:T(end-2))',log10(M{end-2,1}),'bx')
plot(dt(end-1)*(1:T(end-1))',log10(M{end-1,1}),'gx')
plot(dt(end)*(1:T(end))',log10(M{end,1}),'cx')
hold off
title('Error as a function of time','Interpreter','latex')
xlabel('$Time (T)$','interpreter','latex')
ylabel('$\log_{10}(Error)$','interpreter','latex')
legend('\delta_t=0.1','\delta_t=0.01','\delta_t=0.001','\delta_t=0.0001',...
   'location','southeast')
set(gca,'fontsize',18)
grid on