function [val]= LDwFailure(coupling,gamma,power,n,Time,numLines)
% Create the Laplacian
L = diag(sum(coupling,2)) - coupling;

% Define the linear, second order, dynamics
%
% [ dx/dt ] = [  O     I    ] [ x ]
% [ dv/dt ] = [ -L -gamma*I ] [ v ]
%

linear_dynamics = @(t,x) [zeros(n), eye(n); -L, -gamma*eye(n)] * x + [zeros(n,1); power];

% Random initial conditions
x0 = zeros(2*n,1);

% Run ode45 
[T,X] = ode45(linear_dynamics, Time ,x0);
val = {}


    k=0;
    dataSize=size(coupling,1);
    % looks a upper triangule of the matrix. If A(i,j) is not equal to 0. We have
    % a connection to j,i.
    for i = 1:dataSize
        upper=dataSize-i;
        for j = dataSize-upper:dataSize
            if(coupling(i,j)~=0)
            couplingFail= coupling;
            couplingFail(i,j)=0;
            couplingFail(j,i)=0;
        
            % Create the Laplacian new
            LF = diag(sum(couplingFail,2)) - couplingFail;
            Init =X(size(X,1),:);
            linear_dynamics = @(t,x) [zeros(n), eye(n); -LF, -gamma*eye(n)] * x + [zeros(n,1); power];
            [T2,X2] = ode45(linear_dynamics, Time, Init);
            X3 = [X;X2];
            T3 = [T;T2+T(size(T,1),1)];
            
            k=k+1;
        
            val{1,k}=coupling;
            val{2,k}=L;
            val{3,k}=Init;
            val{4,k}=couplingFail;
            val{5,k}=LF;
            val{6,k} =X2(size(X2,1),:);
            val{7,k}=X3;
            val{8,k}=T3;
            
            if(k>=numLines)
               return; 
            end
            

            end
        end
    end


end