clc
clear all
close all

%%% Eigenvalue convergence for 1D Unit Interval

tic

% This vector contains the nodal discretization value that is multiplied to
% the number of modes.
ND=[1 3 10];

% Vector containing the absolute value of the largest fourier mode. The
% number of modes at each iteration is 2M+1.
M=4:4:52;

% This is a cell that collects the numerical error. The first column
% correspond to the first and second eigenvalue error. The third and fourth
% column correspond to the fifth and sixth eigenvalue error.
Cel=cell(1,4);

% Vector containing the true solution.
Tru=pi^2*((1:1:6).^2)';

% First loop for different nodal discretizations. Each iteration increases
% the number of nodal discretization
for z=1:size(ND,2)
N=ND(z)*(2*M+1);

% Main loop increasing the number of modes
for k=1:size(M,2)    
x_p=linspace(0,1,N(k));    

% vector of Modes
y=-M(k):1:M(k);

% Vector of Interior points.
x_int=linspace(-0.5,0.5,N(k));
x_int=x_int(2:end-1);

% Constructing the matrix A and B
A_int=pi^2*y.^2.*exp(1i*pi*x_int'*y);
A_bd=zeros(2,size(y,2));
A_bd(1,:)=exp(1i*y*pi*-0.5);
A_bd(2,:)=exp(1i*y*pi*0.5);
A=[A_int ;A_bd];
B=[exp(1i*pi*x_int'*y) ;zeros(2,size(y,2))];

% Computing the eigenvalues
E=sort(abs(1./eigs(A\B)));

% Computing the error
Cel{1,1}(z,k)=log10(norm(E(1)-Tru(1),2));
Cel{1,2}(z,k)=log10(norm(E(2)-Tru(2),2));
Cel{1,3}(z,k)=log10(norm(E(5)-Tru(5),2));
Cel{1,4}(z,k)=log10(norm(E(6)-Tru(6),2));

end
end

%% Plotting

figure(1)
plot(2*M+1,Cel{1,1}(1,:),'c--x')
hold on
plot(2*M+1,Cel{1,1}(2,:),'r--o')
hold on
plot(2*M+1,Cel{1,1}(3,:),'g--*')
hold off
title('$1^{st}$ Eigenvalue Convergence on Unit Interval','interpreter','latex')
xlabel('Number of Modes','interpreter','latex')
ylabel('$\log_{10}(Error)$','interpreter','latex')
legend('M=1\times (2N+1)','M=3\times (2N+1)','M=10\times (2N+1)','interpreter','latex')

figure(2)
plot(2*M+1,Cel{1,1}(1,:),'c--x')
hold on
plot(2*M+1,Cel{1,1}(2,:),'r--o')
hold on
plot(2*M+1,Cel{1,1}(3,:),'g--*')
hold off
title('$2^{nd}$ Eigenvalue Convergence on Unit Interval','interpreter','latex')
xlabel('Number of Modes','interpreter','latex')
ylabel('$\log_{10}(Error)$','interpreter','latex')
legend('M=1\times (2N+1)','M=3\times (2N+1)','M=10\times (2N+1)','interpreter','latex')

figure(3)
plot(2*M+1,Cel{1,1}(1,:),'c--x')
hold on
plot(2*M+1,Cel{1,1}(2,:),'r--o')
hold on
plot(2*M+1,Cel{1,1}(3,:),'g--*')
hold off
title('$5^{th}$ Eigenvalue Convergence on Unit Interval','interpreter','latex')
xlabel('Number of Modes','interpreter','latex')
ylabel('$\log_{10}(Error)$','interpreter','latex')
legend('M=1\times (2N+1)','M=3\times (2N+1)','M=10\times (2N+1)','interpreter','latex')

figure(4)
plot(2*M+1,Cel{1,4}(1,:),'c--x')
hold on
plot(2*M+1,Cel{1,4}(2,:),'r--o')
hold on
plot(2*M+1,Cel{1,4}(3,:),'g--*')
hold off
title('$6^{th}$ Eigenvalue Convergence on Unit Interval','interpreter','latex')
xlabel('Number of Modes','interpreter','latex')
ylabel('$\log_{10}(Error)$','interpreter','latex')
legend('M=1\times (2N+1)','M=3\times (2N+1)','M=10\times (2N+1)','interpreter','latex')

toc