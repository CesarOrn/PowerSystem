%coupling: Coupling matrix.
%gamma:
%power:
%Time set time to calculate ODE.
%numLines: number of lined failures to look at.

%val: cell containing data
%k are the diffrent simulation

%val{1,k}=coupling
%val{2,k}=Laplacian
%val{3,k}=Last Values if simulaiton before line failure

%val{4,k}=coupling after failing Line
%val{5,k}=Laplacian After Failure
%val{6,k}=Last values after line failure

%val{7,k}=All data of velocities and position
%val{8,k}=Time for whole simulation

%val{9,k}=Line Failed.


function [val]= LDwFailure(coupling,gamma,power,Time,numLines)

% number of elements
n=size(coupling,1);

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
val = {};


    k=0;
    dataSize=size(coupling,1);
    % looks a upper triangule of the matrix. If A(i,j) is not equal to 0. We have
    % a connection to j,i.
    for i = 1:dataSize
        upper=dataSize-i;
        for j = dataSize-upper:dataSize
            if(coupling(i,j)~=0)
            % fail a line once we have founda line to fail.
            couplingFail= coupling;
            couplingFail(i,j)=0;
            couplingFail(j,i)=0;
        %val{13,k}=Line Failed.
            % Create the Laplacian new
            LF = diag(sum(couplingFail,2)) - couplingFail;
            
            % start simulaiton  with initial condition before line cut.
            Init =X(size(X,1),:);
            linear_dynamics = @(t,x) [zeros(n), eye(n); -LF, -gamma*eye(n)] * x + [zeros(n,1); power];
            [T2,X2] = ode45(linear_dynamics, Time, Init);
            
            %append new simulation data.
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
            
            val{9,k}=strcat('Line',num2str(j),', ',num2str(i));
            
            if(k>=numLines)
               return; 
            end
            

            end
        end
    end


end