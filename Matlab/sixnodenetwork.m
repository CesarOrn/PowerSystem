clear; clc; close all;
%
% Consider n coupled second order dynamics 
% d^2 xi/dt^2 = -gamma dxi/dt + SUM_j A(i,j) * (x(j)-x(i)) + Pi;
%
% To use ode45, we convert the n second order equations 2n first order equations 


% Dimension of the system 
n =6;

% get network and power
[A, P] =  getNetwork('network1');
%A = A-diag(diag(A));
%A= A+transpose(A);

% Create the Laplacian
L = diag(sum(A,2)) - A;

% Damping
gamma = 0.9;


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

%BEFORE FALIURE
valDD=DDwFailure(A,gamma,P,Time); %to get coupling matrix, Eigen value and vectors and Q of non failed grid
A=valDD{1,1} %Coupling matrix
V=valDD{2,1};%Eigen vectos
Lambda=valDD{3,1};%Eigen values
Q=valDD{4,1};%Q values


%Finding omega, zeta
Q1=Q(2:n); % Removing the zero in the first positons so that we don't get infinity in the later calculations

Lambda1=zeros([n-1,1]);
for j=1:n-1
    Lambda1(j,1)=Lambda(j+1,j+1);
end
omega=sqrt(Lambda1);

zeta=gamma./(2.*omega);

% Coefficient and line indexes
[Line1, Coeff1]=findLines(A,V); %to find flow coefficients and labels
Line=Line1' % Transpose to plug into error approximation function
Coeff=Coeff1'
writematrix(Coeff,'coeffbeforefail.xlsx')
dim=size(Coeff);
n=dim(1,1);
m=dim(1,2);




%find steady state error, peak, peak time
SS1= Q1./Lambda1;
SSbeforefal=[0;SS1]  % Steady state
pt1=pi./(omega.*sqrt(1-zeta.^2));
ptbeforefal=[0;pt1] % Peak time
peak1=Q1./Lambda1.*(1+exp(-pi.*zeta./sqrt(1-zeta.^2)));
peakbeforefal=[0;peak1] % Peak


%write in excel files
writematrix(peakbeforefal,'peakbeforefal.xlsx')
writematrix(SSbeforefal,'SSbeforefal.xlsx')
writematrix(ptbeforefal,'ptbeforefal.xlsx')



%Find linear sum of steady state to get thetaj-thetai
Linsm=zeros([n,m]);
for j=1:m
    Linsm(:,j)=Coeff(:,j).*SSbeforefal;
end

% Sum each of the columns to attain summation value

SSsum=zeros([1,m]);
for i = 1:m
    SSsum(1,i)=sum(Linsm(:,i));             
end

 writematrix(SSsum,'SSsumbefofal.xlsx')   

 
 %AFTER FALIURE
for i =1:m
Ind=valDD{13,i} %Index of lines
A=valDD{6,i} % Coupling matrix
V=valDD{7,i}; %Eigen vectors
Lambda=valDD{8,i}; %Eigen values
Q=valDD{9,i}; %new Q values
Y=valDD{11,i};
T=valDD{12,i};
Q1=Q(2:n); % Removing the zero in the first positons so that we don't get infinity in the later calculations


figure(i+1)
subplot(2,1,1);
plot(T,Y(:,1:n));
xlabel('Time');
ylabel('Diagonalized Positions');
title('Diagonalized Dynamics');
subplot(2,1,2);
plot(T,Y(:,1+n:2*n))
xlabel('Time');
ylabel('Diagonalized Velocities');


Lambda1=zeros([n-1,1]);
for j=1:n-1
    Lambda1(j,1)=Lambda(j+1,j+1);
end

%Find coefficents and line indexes using find Lines
[Line1, Coeff1]=findLines(A,V); %to find flow coefficients and labels
Line=Line1' ;% Transpose to plug into error approximation function
Coeff=Coeff1'

Coeff2{1,i}=Coeff % Store all coefficients in a cell array

%Find omega, zeta, peak, steady state and teak time
omega=sqrt(Lambda1);
zeta=gamma./(2.*omega);
SS1= Q1./Lambda1;
SS(:,i)=[0;SS1];  % Steady state
pt1=pi./(omega.*sqrt(1-zeta.^2));
pt(:,i)=[0;pt1]; % Peak time
peak1=Q1./Lambda1.*(1+exp(-pi.*zeta./sqrt(1-zeta.^2)));
peak(:,i)=[0;peak1]; % Peak
dim=size(Coeff);
n=dim(1,1);
m=dim(1,2);
end

%display and write values in excel files
peak=peak
writematrix(peak,'peak.xlsx')
SS=SS
writematrix(SS,'SS.xlsx')
pt=pt
writematrix(pt,'pt.xlsx')



%Writing coefficient cell array in exel file
    q=1;
    r=1;
for i =1:length(Coeff2)
    Coeffmat(q:q+n-1,r:r+m-1)=Coeff2{1,i};
    q=n+q+2;
    r=m+r+2;
end

%Finding linear sum of etas to find thetaj-thetai after faliures

%Multiply coefficient to steady State value

for k = 1:m+1
        Linsm=zeros([n,m]);
        newCoeff=cell2mat(Coeff2(1,k));
    for j=1:m
        Linsm(:,j)=newCoeff(:,j).*SS(:,k);
    end

        SSsum=zeros([1 m]);
%Sum the columns together
    for l = 1:m
        SSsum(1,l)=sum(Linsm(:,l))
    end

%store value in array
    SSsum2{1,k}=SSsum;
end

   
    %turn array to matrix to write file
for i =1:length(SSsum2);
    SSmat(i,1:m)=SSsum2{1,i};
end

writematrix(Coeffmat,'Coeff.xlsx')
writematrix(SSmat, 'SSsum.xlsx')






