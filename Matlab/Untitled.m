clear; clc; close all;
%
% Consider n coupled second order dynamics 
% d^2 xi/dt^2 = -gamma dxi/dt + SUM_j A(i,j) * (x(j)-x(i)) + Pi;
%
% To use ode45, we convert the n second order equations 2n first order equations 


% Dimension of the system 
n =10;

% Create a (random) matrix 
[coupling ,P]=getNetwork('network3');


n=size(coupling,1);

% Create the Laplacian
L = diag(sum(coupling,2)) - coupling;

gamma=0.9

% Perform diagonalization
[V,Lambda] = eig(L);
Time = 0:.0018:20;

% Run the diagonalized dynamics
%
%  To diagonalize, pre-multiply by [ inv(V)   O_n   ]
%                                  [  O_n    inv(V) ]
%
Q = V \ P;
diagonalized_dynamics = @(t,x) [zeros(n), eye(n); -Lambda, -gamma*eye(n)] * x + [zeros(n,1); Q];

% Run ode45 for diagonalized dynamics
x0 = zeros(2*n,1);

[T, X] = ode45(diagonalized_dynamics,Time, x0);


couplingFail= coupling;
couplingFail(3,1)=0;
couplingFail(1,3)=0;
        
% Create the Laplacian new
LF = diag(sum(couplingFail,2)) - couplingFail;
[VF,LambdaF] = eig(LF);
        
QF= VF \ P;
        
%simulation will start with these inital conditions. Line
%has failed.
Init =X(size(X,1),:)

     
diagonalized_dynamics = @(t,x) [zeros(n), eye(n); -LambdaF, -gamma*eye(n)] * x + [zeros(n,1); QF];
                
[T2,X2] = ode45(diagonalized_dynamics, Time, Init);
                
%append data from before failing and after failing.
X3 = [X;X2];
T3 = [T;T2+T(size(T,1),1)];
                
k=1                
%save data in cell.
val{1,k}=coupling;
val{2,k}=V;
val{3,k}=Lambda;
val{4,k}=Q;
val{5,k}=Init;
        
val{6,k}=couplingFail;
val{7,k}=VF;
val{8,k}=LambdaF;
val{9,k}=QF;
val{10,k} =X2(size(X2,1),:);
val{11,k}=X3;
val{12,k}=T3;
val{13,k}=strcat('Line',num2str(j),', ',num2str(i));


val{14,k}=[X(:,1:n)*V(2,:)'-X(:,1:n)*V(1,:)';X2(:,1:n)*VF(2,:)'-X2(:,1:n)*VF(1,:)'];



figure('Name','Measured Data');
plot(val{12,1},abs(val{14,1}));


figure('Name','Measured Data');
plot(val{12,1},val{11,1});