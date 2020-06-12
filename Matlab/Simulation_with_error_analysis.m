clear; clc; close all;
%
% Consider n coupled second order dynamics 
% d^2 xi/dt^2 = -gamma dxi/dt + SUM_j A(i,j) * (x(j)-x(i)) + Pi;
%
% To use ode45, we convert the n second order equations 2n first order equations 


% Dimension of the system 
n =1000;

% Create a (random) matrix 
A =  full(BAgraph(n));
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

val=LDwFailure(A,gamma,P,n,Time,1);

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
valDD=DDwFailure(A,gamma,P,n,Time,1);
valDD{1,1}

A=valDD{1,1}
V=valDD{2,1}
Lambda=valDD{3,1}
Q=valDD{4,1}

Y=valDD{10,1};
T=valDD{11,1};

[Line1, Coeff1]=findLines(A,V); %to find flow coefficients and labels
Line=Line1' % Transpose to plug into error approximation function
Coeff=Coeff1'
Q1=Q(2:n); % Removing the zero in the first positons so that we don't get infinity in the later calculations
Lambda1=zeros([n-1,1]);
for i=1:n-1
    Lambda1(i,1)=Lambda(i+1,i+1);
end

omega=sqrt(Lambda1);

zeta=1./(2.*omega);
SS1= Q1./Lambda1;
SS=[0;SS1]  % Steady state
pt1=pi./(omega.*sqrt(1-zeta.^2));
pt=[0;pt1] % Peak time
peak1=Q1./Lambda1.*(1+exp(-pi.*zeta./sqrt(1-zeta.^2)));
peak=[0;peak1] % Peak


 
[SSerror,SSerror2,l]=Errorapproxfunc(size(Coeff),Coeff,SS) % error approx
for i= 1:l
    Maxerror(:,i)=max(SSerror2(:,i)); %matrix of max errors of each flow
end

    


figure(2)
subplot(2,1,1);
plot(T,Y(:,1:n));
xlabel('Time');
ylabel('Diagonalized Positions');
title('Diagonalized Dynamics');
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
subplot(2,1,1)
plot(T,log10(err(:,1:n)));
xlabel('Time');
ylabel('Position Error');
title('Numerical Error Introduced by Diagonalization');
subplot(2,1,2);
plot(T,log10(err(:,1+n:2*n)));
xlabel('Time')
ylabel('Velocity Error');





(sum(sum(A)))/2;
k=0;

for i = 1:n
    upper=n-i;
   for j = n-upper:n
       if(A(i,j)~=0)
           k=k+1;
           etaTheta(k,:)=V(j,:)-V(i,:);
           etaSort(k,:)= sort(abs(etaTheta(k,:)),2,'descend');
           k;
           [etaSort(k,1)/etaSort(k,2)];
           ratio{k,1}=[etaSort(k,1)/etaSort(k,2)];
           strcat('link:',num2str(j),num2str(i));
           ratio{k,2}=strcat('link: ',num2str(j),' , ',num2str(i));
       end        
   end
end
k=0;
size(etaTheta,1);
for i=1:size(etaTheta,1)
    for j=1:n
        k=k+1;
        data(k,:)= [i,etaTheta(i,j)];
    end
end



figure(4);
scatter(data(:,1),data(:,2));
title('Plot Theta_j-Theta_i');
xlabel('Line');
ylabel('Eta coefficient (Theta_j-Theat_i)');    

figure(5);
plot(SSerror2,1:l)
title('Error spread over number of dominant etas used')
xlabel('error percentage');
ylabel('number of dominant etas used');


figure(6);
bar(1:l,Maxerror)
title('Max error VS number of dominant etas used');
xlabel('number of dominant etas used')
ylabel('error percentage')




