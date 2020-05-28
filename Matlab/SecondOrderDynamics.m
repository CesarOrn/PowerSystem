clear; clc; close all;
%
% Consider n coupled second order dynamics 
% d^2 xi/dt^2 = -gamma dxi/dt + SUM_j A(i,j) * (x(j)-x(i)) + Pi;
%
% To use ode45, we convert the n second order equations 2n first order equations 


% Dimension of the system 
n = 6;

KOUPLING=140;

% Create a (random) matrix 

A =  zeros(2,2);

A(1,2)=KOUPLING;
A(2,1)=KOUPLING;

A(1,3)=KOUPLING;
A(3,1)=KOUPLING;

A(1,5)=KOUPLING;
A(5,1)=KOUPLING;

A(4,3)=KOUPLING;
A(3,4)=KOUPLING;

A(4,5)=KOUPLING;
A(5,4)=KOUPLING;

A(6,5)=KOUPLING;
A(5,6)=KOUPLING;

% Create the Laplacian
L = diag(sum(A,2)) - A;

[V,Lambda] = eig(L);

% Damping
gamma = 0.9;

% Forcing vector 
P = zeros(n,1);
P(1)= 40;
P(2)= -20;
P(3)= -25;
P(4)= 40;
P(5)= -25;
P(6)= -10;
% Perform diagonalization
[V,Lambda] = eig(L);

lookat=[zeros(n), eye(n); -L, -gamma*eye(n)];

% Define the linear, second order, dynamics
%
% [ dx/dt ] = [  O     I    ] [ x ]
% [ dv/dt ] = [ -L -gamma*I ] [ v ]
%
linear_dynamics = @(t,x) [zeros(n), eye(n); -L, -gamma*eye(n)] * x + [zeros(n,1); P];

% Random initial conditions
x0 = zeros(1,2*n);

% Run ode45 
[T,X] = ode45(linear_dynamics, [0,30], x0);

% Make some plots
figure(1)
title('Coupled Dynamics');
subplot(2,1,1)
plot(T,X(:,1:n))
xlabel('Time')
ylabel('Positions')

subplot(2,1,2)
plot(T,X(:,n+1:2*n))
xlabel('Time')
ylabel('Velocities')

% Run the diagonalized dynamics
%
%  To diagonalize, pre-multiply by [ inv(V)   O_n   ]
%                                  [  O_n    inv(V) ]
%
Q = V \ P;
diagonalized_dynamics = @(t,x) [zeros(n), eye(n); -Lambda, -gamma*eye(n)] * x + [zeros(n,1); eye(n)*Q];

% Run ode45 for diagonalized dynamics
y0 = zeros(1,2*n);
%y0 = [V, zeros(n); zeros(n), V] \ x0
[T, Y] = ode45(diagonalized_dynamics, T, y0);

figure(2)
title('Diagonalized Dynamics');
subplot(2,1,1);
plot(T,Y(:,1:n));
xlabel('Time');
ylabel('Diagonalized Positions');

subplot(2,1,2);
plot(T,Y(:,1+n:2*n))
xlabel('Time');
ylabel('Diagonalized Velocities');

%
% Demonstrate that we can map Y back to X, xmapped = V*Y
%

xmapped = transpose([V, zeros(n); zeros(n), V] * transpose(Y));

% Define the error to be |xi - (Vy)i]
err = abs(X - xmapped);

figure(3)
title('Numerical Error Introduced by Diagonalization');
subplot(2,1,1)
plot(T,log10(err(:,1:n)));
xlabel('Time');
ylabel('Position Error');

subplot(2,1,2);
plot(T,log10(err(:,1+n:2*n)));
xlabel('Time')
ylabel('Velocity Error');