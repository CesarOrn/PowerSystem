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


%
%eta to function
for i=1:size(V,2)
    Lambda=val{3,1}(i,i);
    omega=sqrt(val{3,1}(i,i));
    zeta=0.9/(2*omega);
    
    Qs=val{4,1}(i,1);
    
    const=Qs/Lambda;
    
   leftFunc{1,i}=@YF;
   leftFunc{2,i}=omega;
   leftFunc{3,i}=zeta;
   leftFunc{4,i}=const;
end

for i=1:size(V,2)
    Lambda=val{8,1}(i,i);
    omega=sqrt(val{8,1}(i,i));
    zeta=0.9/(2*omega);
    
    Qs=val{9,1}(i,1);
    
    const=Qs/Lambda;
    
   rightFunc{1,i}=@YF;
   rightFunc{2,i}=omega;
   rightFunc{3,i}=zeta;
   rightFunc{4,i}=const;
   
end
%

figure('Name','Measured Data');
plot(val{12,1},abs(val{14,1}));


figure('Name','Measured Data');
plot(val{12,1},val{11,1}(:,1:n));
for j =2:size(V)
for i = 1:size(Time,2)
    data(i,j)=YF(Time(i),leftFunc{2,j},leftFunc{3,j},leftFunc{4,j});
end
end
figure('Name','eta from func');
plot(Time,data);

for i = 1:size(Time,2)
    fData(i)=compostion(Time(i),leftFunc,val{2,1},[2 1]);
end

figure('Name','Plot flow');
plot(Time,fData);

GoldenSectionSearch(@compostion,1,2,1e-5,leftFunc,val{2,1},[2 1])

function Y = YF(t,omega,zeta,const) 
    expo=exp(-1*zeta*omega*t)/sqrt(1-zeta^2);
    trig=sin(omega*sqrt(1-zeta^2)*t-acos(zeta));
    Y=const*(1+(expo*trig));
end



function Y = compostion(t,funLeft,eig,line) 
Y=0;
for i=2:size(eig)
    v=eig(line(2),i);
    omega=funLeft{2,i};
    zeta=funLeft{3,i};
    const=funLeft{4,i};
    expo=exp(-1*zeta*omega*t)/sqrt(1-zeta^2);
    trig=sin(omega*sqrt(1-zeta^2)*t-acos(zeta));
    Y=Y+v*const*(1+(expo*trig));
end

for i=2:size(eig)
    v=eig(line(1),i);
    omega=funLeft{2,i};
    zeta=funLeft{3,i};
    const=funLeft{4,i};
    expo=exp(-1*zeta*omega*t)/sqrt(1-zeta^2);
    trig=sin(omega*sqrt(1-zeta^2)*t-acos(zeta));
    Y=Y+v*-1*const*(1+(expo*trig));
end
line(1);
eig(line(1),:);

line(2);
eig(line(2),:);

end
