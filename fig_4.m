clc
clear all

eps=10^(-4);

N=[3 4 5 10 15 20 24];
M_int=[26 34 42 60 100 120 136];
M_bd=[20 30 40 100 130 160 180];

Cat=zeros(2,size(N,2));
Cel=cell(1,size(N,2));

for k=1:size(N,2)
[I_1,I_2]=meshgrid(-N(k):N(k));
I=[I_1(:) I_2(:)];

[x,y]=meshgrid(linspace(-0.5,0.5,M_int(k)));

A=x.*(y>-.5-x & y<0 & (-0.5<=x)&(x<=0))+x.*(y>-.5+x & y<0 & (x<=0.5)&(x>0))+...
    x.*(y<(.25-x.^2).^.5 & x<.5 & x>-.5 & 0<=y);
B=y.*(y>-.5-x & y<0 & (-0.5<=x)&(x<=0))+y.*(y>-.5+x & y<0 & (x<=0.5)&(x>0))+...
    y.*(y<(.25-x.^2).^.5 & x<.5 & x>-.5& 0<=y);

Int=[nonzeros(A) nonzeros(B)];

x_bd=linspace(-.5,0,3*M_bd(k))';
theta=linspace(0,pi,5*M_bd(k))';
bd=[.5*cos(theta) .5*sin(theta);x_bd -x_bd-.5;-x_bd(1:end-1) -x_bd(1:end-1)-.5];

f=-5*ones(size(Int,1),1);

A=[(-pi^2*(I(:,1).^2+I(:,2).^2)').*exp(pi*1i*(Int(:,1)*I(:,1)'+Int(:,2)*I(:,2)'));...
     exp(pi*1i*( bd(:,1)*I(:,1)'+bd(:,2)*I(:,2)' ) )];

 M=size(A,1);
    
    [U,S,V]=svd(A);
    S=diag(S);
    S=S(S>eps);
    S=1./S;
    S=[S;zeros((2*N(k)+1)^2-size(S,1),1)];
    S=[diag(S) zeros(size(S,1),M-(2*N(k)+1)^2)];
    c=V*S*U'*[f;zeros(size(bd,1),1)];

    Cel{1,k}=c;
    Cat(1,k)=size(Int,1);
    Cat(2,k)=size(bd,1);
    Cat(3,k)=Cat(1,k)+Cat(2,k);
end

%%
pt=120;
EC=cell(2,size(N,2));
[x,y]=meshgrid(linspace(-0.5,0.5,pt));
for k=1:size(N,2)
    
[I_1,I_2]=meshgrid(-N(k):N(k));
I=[I_1(:) I_2(:)];

Plt=real(exp(pi*1i*(x(:)*I(:,1)'+y(:)*I(:,2)'))*Cel{1,k}); 
Plt=reshape(Plt,[pt pt]);
Plt=Plt.*(y>-.5-x & y<0 & (-0.5<=x)&(x<=0))+Plt.*(y>-.5+x & y<0 & (x<=0.5)&(x>0))+...
    Plt.*(y<(.25-x.^2).^.5 & x<.5 & x>-.5 & 0<=y);
EC{1,k}=Plt;
EC{2,k}=nonzeros(Plt(:));
end

E=zeros(1,size(N,2)-1);

for k=1:size(N,2)-1
    E(k)=1/sqrt(size(EC{2,end},1))*norm(EC{2,k}-EC{2,end},2);
end

%%
Cat(4,:)=(2*N+1).^2;
EC{1,end}(EC{1,end}==0)=inf;

figure(1)
surf(x,y,EC{1,end})
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$u(x,y)$','interpreter','latex')
title('Numerical Solution','interpreter','latex')
grid on
figure(2)
hold on
plot(Cat(4,1:end-1),log10(E(1,:)),'x--')
hold off
xlabel('N','interpreter','latex')
ylabel('$log_{10}$(Error)','interpreter','latex')
title('Error Convergence for Fixed Sampling rate of 5','interpreter','latex')
legend('M\geq 5N','interpreter','latex')
grid on

