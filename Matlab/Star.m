
clc;clear;
n=10
gamma=0.9

A(1,2)=1;
        A(2,1)=1;
        
        A(1,3)=1;
        A(3,1)=1;
        
        A(1,4)=1;
        A(4,1)=1;
        
        A(2,5)=1;
        A(5,2)=1;
        
        A(5,6)=1;
        A(6,5)=1;
        A(6,2)=1;
        A(2,6)=1;
        
        A(6,7)=1;
        A(7,6)=1;
        A(7,3)=1;
        A(3,7)=1;
        A(8,3)=1;
        A(3,8)=1;
        A(8,7)=1;
        A(7,8)=1;
        A(8,9)=1;
        A(9,8)=1;
        
        
        A(4,9)=1;
        A(9,4)=1;
        A(9,10)=1;
        A(10,9)=1;
        
        A(10,4)=1;
        A(4,10)=1;
        A(5,10)=1;
        A(10,5)=1;
        
    
        
        E(1,1)=1;
E(6,3)=1;

E(10,3)=1;
E(9,3)=1;
E(8,3)=1;
E(7,3)=1;

E(5,3)=1;
E(4,2)=1;
E(3,2)=1;
E(2,2)=1;

 L = diag(sum(A,2)) - A;
   
[V,Lambda] = eig(L)


Lambda1=zeros([n-1,1]);
for j=1:n
    Lambda1(j,1)=Lambda(j,j)
end
Lambda1=Lambda1(2:n)
omega=sqrt(Lambda1)


zeta=gamma./(2.*omega)
OS=exp((-zeta.*pi)./sqrt(1-(zeta).^2))*100

EE=(E'*E)^-1*E'

Q=EE*A*E

P(1,1)=3;
        P(2,1)=1;
        
        P(3,1)=1;
        P(4,1)=1;
        P(5,1)=-1;
        P(6,1)=-1;
        
        P(7,1)=-1;
        P(8,1)=-1;
        P(9,1)=-1;
        P(10,1)=-1;
        
        QP= EE*P
        %Set time steps;
Time = 0:.0018:10;

linear_dynamics = @(t,x) [zeros(n), eye(n); -L, -gamma*eye(n)] * x + [zeros(n,1); P];

% Random initial conditions
x0 = zeros(2*n,1);

% Run ode45 
[T,X] = ode45(linear_dynamics, Time ,x0);

figure (7)
plot(T,X(:,1:10))
title('Displacement over time for all nodes')
xlabel('time')
ylabel('displacement')


X1=X(:,1);
X2=(X(:,5)+X(:,6)+X(:,7)+X(:,8)+X(:,9)+X(:,10))./6;
X3=(X(:,2)+X(:,3)+X(:,4))./3;

figure (8)
plot(T,[X1,X2,X3])
title('Displacement over time for all clusters')
xlabel('time')
ylabel('displacement')
