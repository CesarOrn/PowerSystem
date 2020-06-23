clear; clc; close all;
%
% Consider n coupled second order dynamics 
% d^2 xi/dt^2 = -gamma dxi/dt + SUM_j A(i,j) * (x(j)-x(i)) + Pi;
%
% To use ode45, we convert the n second order equations 2n first order equations 


% Dimension of the system 
n =6;

% Create a (random) matrix 
A =  getNetwork(network1);
%A = A-diag(diag(A));
%A= A+transpose(A);

% Create the Laplacian
L = diag(sum(A,2)) - A;

% Damping
gamma = 0.9;

% Forcing vector 
P = rand(n,1);
P = P -mean(P);

%Set time steps;
Time = 0:.0018:10;



% Perform diagonalization

val=LDwFailure(A,gamma,P,Time,1);

% Make some plots
figure(1)
T=val{8,1};
X=val{7,1};
subplot(2,1,1)
plot(T,X(:,1:n))
xlabel('Time')
ylabel('Positions')
title('Coupled Dynamics');

subplot(2,1,2)
plot(T,X(:,n+1:2*n))
xlabel('Time')
ylabel('Velocities')

% Run the diagonalized dynamics
%
%  To diagonalize, pre-multiply by [ inv(V)   O_n   ]
%                                  [  O_n    inv(V) ]
%

valDD=DDwFailure(A,gamma,P,Time);
pop=size(valDD);
A=valDD{1,1}
V=valDD{2,1};
Lambda=valDD{3,1};
Q=valDD{4,1};

Y=valDD{9,1};
T=valDD{12,1};
Q1=Q(2:n); % Removing the zero in the first positons so that we don't get infinity in the later calculations
Lambda1=zeros([n-1,1]);
for j=1:n-1
    Lambda1(j,1)=Lambda(j+1,j+1);
end
[Line1, Coeff1]=findLines(A,V); %to find flow coefficients and labels
Line=Line1' % Transpose to plug into error approximation function
Coeff=Coeff1'
omega=sqrt(Lambda1);

zeta=gamma./(2.*omega);

SS1= Q1./Lambda1;
SSbeforefal=[0;SS1];  % Steady state
pt1=pi./(omega.*sqrt(1-zeta.^2));
ptbeforefal=[0;pt1]; % Peak time
peak1=Q1./Lambda1.*(1+exp(-pi.*zeta./sqrt(1-zeta.^2)));
peakbeforefal=[0;peak1]; % Peak
peakbeforefal=peakbeforefal
writematrix(peakbeforefal,'peakbeforefal.xlsx')
SSbeforefal=SSbeforefal
writematrix(SSbeforefal,'SSbeforefal.xlsx')
ptbeforefal=ptbeforefal
writematrix(ptbeforefal,'ptbeforefal.xlsx')

for i =1:pop(2)
Ind=valDD{13,i}
A=valDD{6,i}
V=valDD{7,i};
Lambda=valDD{8,i};
Q=valDD{9,i};

Y=valDD{10,i};
T=valDD{12,i};
Q1=Q(2:n); % Removing the zero in the first positons so that we don't get infinity in the later calculations
Lambda1=zeros([n-1,1]);
for j=1:n-1
    Lambda1(j,1)=Lambda(j+1,j+1);
end
[Line1, Coeff1]=findLines(A,V); %to find flow coefficients and labels
Line=Line1' ;% Transpose to plug into error approximation function
Coeff=Coeff1'
Coeff2{1,i}=Coeff
omega=sqrt(Lambda1);
omegareal=[0;omega];
zeta=gamma./(2.*omega);
zetareal=[0;zeta];
SS1= Q1./Lambda1;
SS(:,i)=[0;SS1];  % Steady state
pt1=pi./(omega.*sqrt(1-zeta.^2));
pt(:,i)=[0;pt1]; % Peak time
peak1=Q1./Lambda1.*(1+exp(-pi.*zeta./sqrt(1-zeta.^2)));
peak(:,i)=[0;peak1]; % Peak
dim=size(Coeff)
n=dim(1,1);
m=dim(1,2);


end

peak=peak
writematrix(peak,'peak.xlsx')
SS=SS
writematrix(SS,'SS.xlsx')
pt=pt
writematrix(pt,'pt.xlsx')



    q=1
    r=1
for i =1:length(Coeff2)
    Coeffmat(q:q+n-1,r:r+m-1)=Coeff2{1,i};
    q=n+q+2;
    r=m+r+2;
end

%Multiply coefficient to steady State value

for k = 1:m
        Linsm=zeros([n,m]);
        newCoeff=cell2mat(Coeff2(1,k))
for j=1:m
    Linsm(:,j)=newCoeff(:,j).*SS(:,k)
end
for l = 1:m
    SSsum(1,l)=sum(Linsm(:,l))
end
    SSsum2{1,k}=SSsum
end

   
    
for i =1:length(SSsum2)
    SSmat(i,1:m)=SSsum2{1,i};
end

writematrix(Coeffmat,'Coeff.xlsx')
writematrix(SSmat, 'SSsum.xlsx')







