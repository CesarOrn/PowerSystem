%coupling: Coupling matrix.
%gamma:
%power:
%Time set time to calculate ODE.
%numLines: number of lines of line failures.

%val: cell containing data
%k are the diffrent simulation

%val{1,k}=coupling;
%val{2,k}=EigenVector;
%val{3,k}=Lambdas;
%val{4,k}=Q;
%val{5,k}=Last values first part simulation;

%val{6,k}=coupling after failure;
%val{7,k}=eigenVector after failure;
%val{8,k}=Lambda  after failure;
%val{9,k}=Q After failure;
%val{10,k} =Last values second part simulation;
        
%val{11,k}=Values of all simulation;
%val{12,k}=Time Values.

function [val]=DDwFailure(coupling,gamma,power,Time,numLines)

% number of elements
n=size(coupling,1);

% Create the Laplacian
L = diag(sum(coupling,2)) - coupling;

% Perform diagonalization
[V,Lambda] = eig(L);

% Run the diagonalized dynamics
%
%  To diagonalize, pre-multiply by [ inv(V)   O_n   ]
%                                  [  O_n    inv(V) ]
%
Q = V \ power;
diagonalized_dynamics = @(t,x) [zeros(n), eye(n); -Lambda, -gamma*eye(n)] * x + [zeros(n,1); Q];

% Run ode45 for diagonalized dynamics
x0 = zeros(2*n,1);
[T, X] = ode45(diagonalized_dynamics,Time, x0);

val = {}

    k=0;
    dataSize=size(coupling,1);
    % looks a upper triangule of the matrix. If A(i,j) is not equal to 0. We have
    % a connection to j,i.
    for i = 1:dataSize
        upper=dataSize-i;
        for j = dataSize-upper:dataSize
            %we find a line that is connected. If it is connected. Fail the
            %line and simulate the new ODE with inital conditions.
            if(coupling(i,j)~=0)
                k=k+1;
                %fail a line by setting A(i,j)=A(j,i)=0
                couplingFail= coupling;
                couplingFail(i,j)=0;
                couplingFail(j,i)=0;
        
                 % Create the Laplacian new
                 LF = diag(sum(couplingFail,2)) - couplingFail;
                [VF,LambdaF] = eig(LF);
        
                QF= VF \ power;
        
                %simulation will start with these inital conditions. Line
                %has failed.
                Init =X(size(X,1),:);
                
                diagonalized_dynamics = @(t,x) [zeros(n), eye(n); -LambdaF, -gamma*eye(n)] * x + [zeros(n,1); QF];
                
                [T2,X2] = ode45(diagonalized_dynamics, Time, Init);
                
                %append data from before failing and after failing.
                X3 = [X;X2];
                T3 = [T;T2+T(size(T,1),1)];
                
                
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
                if(k>=numLines)
                    return; 
                end

            end
        end
    end    
        
    
    end