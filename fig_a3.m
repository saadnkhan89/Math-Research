%% Number of Modes BDF2
M=5;
I=(-M:M);

%% Spatial discretization
N=21;
x=linspace(-0.5,0.5,N)';
x=x(2:end-1);

%% Time step

k=.001;
T=200;

%% Forcing and derivative operator
% f_0=1-4*pi^2*sin(2*pi*x);
D2_f=[(1+(2*k/3)*pi^2*I.^2).*exp(1i*pi*x*I);exp(1i*pi*[-0.5;.5]*I)];

%% Initialization

U=zeros(N,T);
U0=[4*x(1:5)+2;-4*x(6:10);20*x(11:15);-20*x(16:end)+10];
% U0=[1;0;1;0;1;0;1];

C=exp(1i*pi*[-.5;x;.5]*I)\[0;U0;0];
% f_0=real((-pi^2*I.^2).*exp(1i*pi*[-.5;x;.5]*I)*C);
f_0=zeros(N,1);
f_1=U0+(k)*f_0(2:end-1);

C1=[(1+(k/2)*pi^2*I.^2).*exp(1i*pi*x*I);exp(1i*pi*[-0.5;.5]*I)]\[f_1;0;0];
U1=real(exp(1i*pi*x*I)*C);

U(2:end-1,1)=U0;
U(2:end-1,2)=U1;

E=zeros(T,1);
E(1)=norm(U0,2);
E(2)=norm(U1,2);

for z=3:T
    C2=D2_f\[4*U1/3-U0/3;0;0];
    U2=real(exp(1i*pi*x*I)*C2);
    U(2:end-1,z)=U2;
    E(z)=norm(U2,2);
    U0=U1;
    U1=U2;
end


%% First Plot
figure(1)
plot((1:200)'*k,E,'x-')
xlabel('Time','Interpreter','latex')
ylabel('Energy','Interpreter','latex')
title('Energy of Solution','Interpreter','latex')
set(gca,'fontsize',18)
grid on


% figure(1)
% plot([-0.5;x;0.5],U(:,1))
% xlabel('x','Interpreter','latex')
% ylabel('$u(t,x)$','Interpreter','latex')
% title('Initial Condition','Interpreter','latex')
% grid on


%% Animation

N_p=50;
x_p=linspace(-0.45,.45,N_p-2);


U_interp=zeros(N_p,T);
for z=2:T
    U_interp(2:end-1,z)=interp1(x,U(2:end-1,z),x_p,'cubic');
end

[t_p1,x_p1]=meshgrid(0:k:k*(T-1),linspace(-0.5,0.5,N));
U_p1=zeros(N,T);
U_p1=1./U_p1;
U_p1(:,1)=U(:,1);
[t_p,x_p]=meshgrid(0:k:k*(T-1),linspace(-0.5,0.5,N_p));
U_p=zeros(N_p,T);
U_p=1./U_p;


for z=1:size(U,2)
    if z==1
        figure(2)
        surf(t_p1,x_p1,U_p1)
    end
    U_p(:,z)=U_interp(:,z);
    figure(2)
    hold on
    surf(t_p,x_p,U_p)
    xlim([0,.2])
    set(gca,'fontsize',18)
    grid on
    hold off
    
end

%% 2nd plot
figure(3)
surf(t_p1,x_p1,U_p1)
hold on
surf(t_p,x_p,U_p)
    xlabel('t','Interpreter','latex')
        ylabel('$x$','Interpreter','latex')
        zlabel('$u(t,x)$','Interpreter','latex')
        title('Evolution of solution','Interpreter','latex')
    xlim([0,.2])
    set(gca,'fontsize',18)
