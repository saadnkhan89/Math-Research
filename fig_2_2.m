clc
clear all


eps=10^(-5);
N=[5 10 15 20 25 35];
M=[2*N+1;10*N;25*N];
Cel=cell(3,1);

for l=1:size(M,1)
E=zeros(size(N,2),1);
for k=1:size(N,2)

    I=-N(k):N(k);
    x=linspace(-.5,.5,M(l,k))'; 
    x=x(2:end-1);
    
    t=x.*exp(x);

    A=[exp(1i*pi*x*I)*(-pi^2*diag(-N(k):N(k)).^2-2*eye(2*N(k)+1));exp(1i*pi*-.5*I)*(1i*pi*diag(-N(k):N(k)));exp(1i*pi*.5*I)];

    [U,S,V]=svd(A);
    S=diag(S);
    S=S(S>eps);
    S=1./S;
    S=[S;zeros((2*N(k)+1)-size(S,1),1)];
    S=[diag(S) zeros(size(S,1),M(l,k)-(2*N(k)+1))];
    c=V*S*U'*[2*exp(x)-x.*exp(x);0.5*exp(-0.5);.5*exp(.5)];

    Sol=real(exp(1i*pi*x*I)*c);
    E(k)=norm(t-Sol,2);
end
Cel{l,1}=E;
end


pt=linspace(-0.5,.5,20);
figure(1)
subplot(1,2,1)
plot(x,Sol,'--',pt,pt.*exp(pt),'x')
title('Numerical Solution','interpreter','latex')
xlabel('x','interpreter','latex')
ylabel('$u(x)$','interpreter','latex')
legend('Numerical Solution','True Solution','interpreter','latex')
grid on
subplot(1,2,2)
hold on
plot(2*N'+1,log10(Cel{1,1}),'x-')
plot(2*N'+1,log10(Cel{2,1}),'o-')
plot(2*N'+1,log10(Cel{3,1}),'*-')
hold off
title('Error Convergence for Different Sampling','interpreter','latex')
xlabel('N','interpreter','latex')
ylabel('$log_{10}$(Error)','interpreter','latex')
legend('M=N','M=10N','M=25N','interpreter','latex')
grid on
