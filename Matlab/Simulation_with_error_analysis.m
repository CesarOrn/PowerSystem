clear; clc; close all;
%
% Consider n coupled second order dynamics 
% d^2 xi/dt^2 = -gamma dxi/dt + SUM_j A(i,j) * (x(j)-x(i)) + Pi;
%
% To use ode45, we convert the n second order equations 2n first order equations 


% Dimension of the system 
n =10;

% Create a (random) matrix 
[LK,Power]=getNetwork('network1')
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
Time = 0:.0018:20;

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


valDD=DDwFailure(A,gamma,P,Time,100);

Coff={};
fprintf('error');
AB=valDD{1,1};
VB=valDD{2,1};
[BFLine1, BFCoeff1]=findLines(AB,VB);
Coff{1,1}=BFLine1';
Coeff=BFCoeff1';
Coff{2,1}=BFCoeff1';
       
Lambda=valDD{3,1};
Q=valDD{4,1};
    
trans=transientResponse(Lambda,Q);
BFOmega=trans{2,1};
Coff{3,1}=trans{4,1};
Coff{4,1}=trans{5,1};
SS=trans{4,1};
Coff{5,1}=trans{6,1};

BFtrans=cell2mat({trans{4,1},trans{5,1},trans{6,1}});
    
[SSerror,SSerror2,l]=Errorapproxfunc(size(Coeff),Coeff,SS) ;% error approx

BFMaxerror=max(SSerror2);
Coff{6,1}= BFMaxerror;
    
for i=1:size(valDD,2)
    A=valDD{6,i};
    VF=valDD{7,i};
    [Line1, Coeff1]=findLines(A,VF);
    Coff{1,i+1}=Line1';
    Coeff=Coeff1';
    Coff{2,i+1}=Coeff1';
    
    
    Lambda=valDD{8,i};
    Q=valDD{9,i};
    
    trans=transientResponse(Lambda,Q);
    Omega=trans{2,1};
    Coff{3,i+1}=trans{4,1};
    Coff{4,i+1}=trans{5,1};
    SS=trans{4,1};
    Coff{5,i+1}=trans{6,1};
    trans=cell2mat({trans{4,1},trans{5,1},trans{6,1}});
    
    [SSerror,SSerror2,l]=Errorapproxfunc(size(Coeff),Coeff,SS) ;% error approx

    Maxerror=max(SSerror2);
    Coff{6,i+1}= Maxerror;
    name= strcat('myData',num2str(i),'.xlsx');
    exportTable(VB,BFtrans ,BFLine1,BFCoeff1,BFMaxerror,BFOmega,VF,trans ,Line1,Coeff1,Maxerror,Omega,valDD{13,i},name);
end



Y=valDD{11,1};
T=valDD{12,1};

 

fprintf('plotting');


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

V=VB;

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
l
plot(SSerror2,1:l)
title('Error spread over number of dominant etas used')
xlabel('error percentage');
ylabel('number of dominant etas used');


figure(6);
l
bar(1:l,Maxerror)
title('Max error VS number of dominant etas used');
xlabel('number of dominant etas used')
ylabel('error percentage')




