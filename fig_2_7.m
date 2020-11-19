clc
clear all


N=[3 4 5 10 15 20 24];
M_int=[26 34 42 60 100 120 136];
M_bd=[20 30 40 100 130 160 180];

epsl=10^(-4);

Cat=zeros(2,size(N,2));
Cel=cell(1,size(N,2));

for k=1:size(N,2)
[I_1,I_2]=meshgrid(-N(k):N(k));
I=[I_1(:) I_2(:)];

[x,y]=meshgrid(linspace(-0.5,0.5,M_int(k)));


A=x.*(x<.25 & x>-.25 & y<.25 & y>-.25)+...
    x.*(x<=-.25 & x>-.5 & (x+.25).^2+y.^2<1/16 & -.25<y & y<.25)+...
    x.*(x.^2+(y-.25).^2<1/16 & (x<0.25)&(x>-0.25) & 0.25<y & y<.5)+...
    x.*(x.^2+(y+.25).^2<1/16 & (x<0.25)&(x>-0.25) & -.5<y & y<-.25)+...
    x.*((x-.25).^2+y.^2<1/16 & x>=.25 &-.25<y & y<.25);

B=y.*(x<.25 & x>-.25 & y<.25 & y>-.25)+...
    y.*(x<=-.25 & x>-.5 & (x+.25).^2+y.^2<1/16 & -.25<y & y<.25)+...
    y.*(x.^2+(y-.25).^2<1/16 & (x<0.25)&(x>-0.25) & 0.25<y & y<.5)+...
    y.*(x.^2+(y+.25).^2<1/16 & (x<0.25)&(x>-0.25) & -.5<y & y<-.25)+...
    y.*((x-.25).^2+y.^2<1/16 & x>=.25 & -.25<y & y<.25);

Int=[nonzeros(A) nonzeros(B)];

x_bd=linspace(-.5,0,2*M_bd(k))';
theta=linspace(0,pi,2*M_bd(k))';
theta1=linspace(0,pi/2,M_bd(k))';
theta2=linspace(pi/2,pi,M_bd(k))';
bd=[.25*cos(theta) .25+.25*sin(theta);.25*cos(theta2)-.25 .25*sin(theta2);...
    .25*cos(theta1)+.25 .25*sin(theta1);.25*cos(theta) -.25-.25*sin(theta);...
    .25*cos(theta2(1:end-1))-.25 -.25*sin(theta2(1:end-1));...
    .25*cos(theta1(2:end))+.25 -.25*sin(theta1(2:end))];

f=Int(:,1).*Int(:,2);

A=[(-1+10^(-3)*pi^2*(I(:,1).^2+I(:,2).^2)').*exp(pi*1i*(Int(:,1)*I(:,1)'+Int(:,2)*I(:,2)'));...
     exp(pi*1i*( bd(:,1)*I(:,1)'+bd(:,2)*I(:,2)' ) )];

 M=size(A,1);
    
    [U,S,V]=svd(A);
    S=diag(S);
    S=S(S>epsl);
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
Plt=Plt.*(x<.25 & x>-.25 & y<.25 & y>-.25)+...
    Plt.*(x<=-.25 & x>-.5 & (x+.25).^2+y.^2<1/16 & -.25<y & y<.25)+...
    Plt.*(x.^2+(y-.25).^2<1/16 & (x<0.25)&(x>-0.25) & 0.25<y & y<.5)+...
    Plt.*(x.^2+(y+.25).^2<1/16 & (x<0.25)&(x>-0.25) & -.5<y & y<-.25)+...
    Plt.*((x-.25).^2+y.^2<1/16 & x>=.25 &-.25<y & y<.25);
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
subplot(1,2,1)
surf(x,y,EC{1,end})
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$u(x,y)$','interpreter','latex')
title('Numerical Solution','interpreter','latex')
grid on
subplot(1,2,2)
hold on
plot(Cat(4,1:end-1),log10(E(1,:)),'x--')
hold off
xlabel('N','interpreter','latex')
ylabel('$log_{10}$(Error)','interpreter','latex')
title('Error Convergence for Fixed Sampling rate of 5','interpreter','latex')
legend('M\geq 5N','interpreter','latex')
grid on
