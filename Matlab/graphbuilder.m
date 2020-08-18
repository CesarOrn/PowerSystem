

clc
clear
n =3;
gamma=0.9

% get network and power
[A, P] = symmetricNetwork('starquo')
%A = A-diag(diag(A));
%A= A+transpose(A);

% Create the Laplacian

L = diag(sum(A,2)) - A
coupling=A;
dataSize=size(coupling);

[V,Lambda] = eig(L);
[line, ~]=findLines(coupling,V);
popo=size(line');
numLines=popo(2);
Time = 0:.0018:20;
Q = V \ P;
linear_dynamics = @(t,x) [zeros(n), eye(n); -L, -gamma*eye(n)] * x + [zeros(n,1); P];

% Random initial conditions
x0 = zeros(2*n,1);

% Run ode45 
[T,X] = ode45(linear_dynamics, Time ,x0);

    % looks a upper triangule of the matrix. If A(i,j) is not equal to 0. We have
    % a connection to j,i.
    q=1;
    for i = 1:dataSize
        upper=dataSize-i; 
        for j = dataSize-upper:dataSize
            if(coupling(i,j)~=0)
                X1=X(:,j)-X(:,i);
                X2(:,q)=X1;
                q=q+1;
            end
        end       
    end
    
figure(1)
plot(T,X2(:,1:n))
xlabel('Time')
ylabel('Flow');

figure (2)
plot(T,X(:,1:n))
xlabel('Time')
ylabel('Displacement')
