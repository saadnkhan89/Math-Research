clear all
clc
close all

%%% Eigenvalue Convergence for Circle centered at zero.

tic

% Vector containing the largest mode at each iteration. The total number of
% modes at each iteration will be 2M+1
M=4:2:16;

% Nodal discretization at each iteration is 3 times the number of
% modes for interior points and 25 times half the number of boundary
% points.
N_int=3*(2*M+1);
N_bd=25*(2*M+1);

% This cell collects the eigenvalues. The first row corresponds to the
% numerically computed first eigenvalue as the number of modes increase.
% The second row corresponds to the numerically computed second eigenvalue
% as the number of modes increase.
Cel=cell(2,size(M,2));

% Matrix that collects the error. The first row is the error for the first
% eigenvalue. The second is the error for the second eigenvalue.
Err=zeros(2,size(M,2)-1);

% Main loop increasing the number of modes
for k=1:size(M,2)

% Creating the index for the modes in 2D      
[i_1,i_2]=meshgrid(-M(k):1:M(k));
I=[i_1(:) i_2(:)];

% Creating Interior Points. The first column is the x-value wheras the
% second column is the corresponding y value.
[x,y]=meshgrid(linspace(-0.5,0.5,N_int(k)));
Int=[x(:),y(:)];
Int=Int(x(:).^2+y(:).^2<.5^2,:);

% Creating the matrix of boundary points.
bd=[linspace(-0.5,0.5,N_bd(k))' sqrt(0.5^2-linspace(-0.5,0.5,N_bd(k)).^2)';...
    linspace(-0.5,0.5,N_bd(k))' -sqrt(0.5^2-linspace(-0.5,0.5,N_bd(k)).^2)'];

% Constructing the matrix A and B
A_int=(pi^2*(I(:,1).^2+I(:,2).^2))'.*exp(1i*pi*Int*I');
A_bd=exp(1i*pi*bd*I');
A=[A_int ;A_bd];
B=[exp(1i*pi*Int*I') ;zeros(size(bd,1),size(I,1))];

% Collecting the smallest two eigenvalues eigenvalue
E=eigs(A\B);
E=sort(abs(1./E));
Cel{1,k}=E(1);
Cel{2,k}=E(2);

end

% Computing the Error

for j=1:2
     for k=1:size(M,2)-1
       Err(j,k)=norm(Cel{j,k}-Cel{j,7},2);
     end
end

%% Plotting
figure(1)
plot((2*M(1:end-1)+1).^2,log10(Err(1,:)),'r--x');
hold on
plot((2*M(1:end-1)+1).^2,log10(Err(2,:)),'b--*');
title('Convergence of Eigenvalues on Unit circle $\Omega$','interpreter','latex')
xlabel('Number of Modes','interpreter','latex')
ylabel('$\log_{10}(Error)$','interpreter','latex')
legend('\lambda_1','\lambda_2','interpreter','latex')

toc