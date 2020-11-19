eps=10^(-4);
N=[3 10 15];

N_int=[23 31 41 51 60;23 31 41 51 60;31 41 51 61 70];
N_bd=[100 175 375 575 2000;100 175 375 575 2000;275 575 775 1075 1500];


pt=zeros(size(N_int,1),size(N_int,2));
E=zeros(size(N_int,1),size(N_int,2));

for j=1:size(N_int,1)
    [I_1,I_2]=meshgrid(-N(j):N(j));
    I=[I_1(:) I_2(:)];

for k=1:size(N_int,2)
[x,y]=meshgrid(linspace(-.5,.5,N_int(j,k)));
Int=zeros(size(x(x.^2+y.^2<.25),1),2);
Int(:,1)=x(x.^2+y.^2<.25);
Int(:,2)=y(x.^2+y.^2<.25);

bd=zeros(N_bd(j,k),2);
bd(:,1)=.5*cos(linspace(0,2*pi,N_bd(j,k))');
bd(:,2)=.5*sin(linspace(0,2*pi,N_bd(j,k))');

f=4*Int(:,1)-2*Int(:,2)-12;


A=[(-pi^2*(2*I(:,1).^2+4*I(:,2).^2)'+1i*pi*(-2*I(:,1)+I(:,2))').*exp(pi*1i*(Int(:,1)*I(:,1)'+Int(:,2)*I(:,2)'));...
     exp(pi*1i*( bd(:,1)*I(:,1)'+bd(:,2)*I(:,2)' ) )];

 M=size(A,1);
%  
    [U,S,V]=svd(A);
    S=diag(S);
    S=S(S>eps);
    S=1./S;
    S=[S;zeros((2*N(j)+1)^2-size(S,1),1)];
    S=[diag(S) zeros(size(S,1),M-(2*N(j)+1)^2)];
    c=V*S*U'*[f;zeros(N_bd(j,k),1)];

    Sol=real(exp(pi*1i*(Int(:,1)*I(:,1)'+Int(:,2)*I(:,2)'))*c);
    
    t=.25-(Int(:,1).^2+Int(:,2).^2);

    E(j,k)=1/sqrt(size(t,1))*norm(t-Sol,2);
    pt(j,k)=M;
end
end

%%

Plt=real(exp(pi*1i*(x(:)*I(:,1)'+y(:)*I(:,2)'))*c);
Plt=reshape(Plt,[70 70]);
Plt(x.^2+y.^2>=.25)=0;

figure(1)
subplot(1,2,1)
surf(x,y,Plt)
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$u(x,y)$','interpreter','latex')
title('Numerical Solution','interpreter','latex')
grid on
subplot(1,2,2)
hold on
plot(pt(1,:),log10(E(1,:)),'x--')
plot(pt(2,:),log10(E(2,:)),'x--')
plot(pt(3,:),log10(E(3,:)),'x--')
hold off
xlabel('M','interpreter','latex')
ylabel('$log_{10}$(Error)','interpreter','latex')
title('Error Convergence for Increasing Sampling rate','interpreter','latex')
legend('N=49','N=441','N=961','interpreter','latex')
grid on