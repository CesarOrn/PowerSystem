clear; clc; close all;
%
% Consider n coupled second order dynamics 
% d^2 xi/dt^2 = -gamma dxi/dt + SUM_j A(i,j) * (x(j)-x(i)) + Pi;
%
% To use ode45, we convert the n second order equations 2n first order equations 



% Create a (random) matrix 
[coupling,P]=getNetwork('network1');

% Dimension of the system 
n =size(coupling,1);


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

P(1)=2;

        
% Create the Laplacian new
LF = diag(sum(coupling,2)) - coupling;
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
        
val{6,k}=coupling;
val{7,k}=VF;
val{8,k}=LambdaF;
val{9,k}=QF;
val{10,k} =X2(size(X2,1),:);
val{11,k}=X3;
val{12,k}=T3;
val{13,k}=strcat('Line',num2str(2),', ',num2str(1));


val{14,k}=[X(:,1:n)*V(2,:)'-X(:,1:n)*V(1,:)';X2(:,1:n)*VF(2,:)'-X2(:,1:n)*VF(1,:)'];

Y=val{11,1};
T=val{12,1};




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
    
    Qs=val{9,1}(i,1)-val{4,1}(i,1);
    
    const=Qs/Lambda;
    
   rightFunc{1,i}=@YF;
   rightFunc{2,i}=omega;
   rightFunc{3,i}=zeta;
   rightFunc{4,i}=const;
   
end

figure('Name','Measured Data');
plot(val{12,1},val{11,1}(:,1:n));
for j =2:size(V)
    for i = 1:size(Time,2)
        Left(i,1)=YF(Time(i),leftFunc{2,j},leftFunc{3,j},leftFunc{4,j});
        Right(i,1)=YF(Time(i),rightFunc{2,j},rightFunc{3,j},rightFunc{4,j});
        
    end

    data(:,j)=[Left;Left(size(Left,1))+Right];
end

for j =2:size(V)
    for i = 1:3
        [peakT,peak]=MaxTrig(i,leftFunc{2,j},leftFunc{3,j},leftFunc{4,j});
        LeftM(i,1)=peakT;
        LeftM(i,2)=peak;
        
        [peakT,peak]=MaxTrig(i,rightFunc{2,j},rightFunc{3,j},rightFunc{4,j});
        RightM(i,1)=peakT;
        RightM(i,2)=peak;
    end
    
    syms x
    S1=(LeftM(1,2))/(LeftM(1,1));
    S2=(LeftM(2,2)-LeftM(1,2))/(LeftM(2,1)-LeftM(1,1));
    S3=(LeftM(3,2)-LeftM(2,2))/(LeftM(3,1)-LeftM(2,1));
    LeftP=piecewise(x<LeftM(1,1),S1*x,LeftM(1,1)<x<LeftM(2,1),S2*(x-LeftM(1,1))+LeftM(1,2),LeftM(2,1)<x<LeftM(3,1),S3*(x-LeftM(2,1))+LeftM(2,2),leftFunc{4,j});
    leftFunc{4,j}
    
    S1=(RightM(1,2))/(RightM(1,1));
    S2=(RightM(2,2)-RightM(1,2))/(RightM(2,1)-RightM(1,1));
    S3=(RightM(3,2)-RightM(2,2))/(RightM(3,1)-RightM(2,1));
    RightP=piecewise(x<RightM(1,1),S1*x,RightM(1,1)<x<RightM(2,1),S2*(x-RightM(1,1))+RightM(1,2),RightM(2,1)<x<RightM(3,1),S3*(x-RightM(2,1))+RightM(2,2),rightFunc{4,j});
    rightFunc{4,j}
    rightFunc{4,j}+leftFunc{4,j}
    LTD=subs(LeftP, x, Time)';
    RTD=subs(RightP, x, Time)';
    trigDataS(:,j)=[LTD;LTD(size(LTD,1))+RTD];
end


trigData=double(trigDataS);

figure('Name','eta from func');
plot([Time Time+Time(size(Time,2))],[data,trigData]);

figure('Name','eta from tri');
plot([Time Time+Time(size(Time,2))],trigData);


fprintf('plotting');


figure('Name','Eta w Trig');

subplot(3,3,1)
plot(T,[Y(:,1),trigData(:,1)]);
xlabel('Time');
ylabel('Eta');
title('EtaWTri');
subplot(3,3,2)
plot(T,[Y(:,2),trigData(:,2)]);
xlabel('Time');
ylabel('Eta');
title('EtaWTri');

subplot(3,3,3)
plot(T,[Y(:,3),trigData(:,3)]);
xlabel('Time');
ylabel('Eta');
title('EtaWTri');

subplot(3,3,4)
plot(T,[Y(:,4),trigData(:,4)]);
xlabel('Time');
ylabel('Eta');
title('EtaWTri');

subplot(3,3,5)
plot(T,[Y(:,5),trigData(:,5)]);
xlabel('Time');
ylabel('Eta');
title('EtaWTri');


subplot(3,3,6)
plot(T,[Y(:,6),trigData(:,6)]);
xlabel('Time');
ylabel('Eta');
title('EtaWTri');



subplot(3,3,7)
plot(T,[Y(:,1:n),trigData]);
xlabel('Time');
ylabel('PowerData');
title('Power');



[~,~,Pairs]=findLines(coupling,V);

size(Pairs,1)
size(Time,2)
for i =1:size(Pairs,1)
    LD(:,1)=sum(X(:,1:n).*val{2,1}(Pairs(i,2),:),2)-sum(X(:,1:n).*val{2,1}(Pairs(i,1),:),2);
        
        
    LDA(:,1)=sum(X2(:,1:n).*val{7,1}(Pairs(i,2),:),2)-sum(X2(:,1:n).*val{7,1}(Pairs(i,1),:),2);
        
    lineData(:,i)=[LD;LDA];
end

SizesF=size(trigData,1)/2;
SizesF2=size(trigData,1);
[sum(trigData(1:SizesF,:).*val{2,1}(:,Pairs(1,2))',2),sum(trigData(1:SizesF,:).*val{2,1}(:,Pairs(1,1))',2)]

for i =1:size(Pairs,1)
    
    SizesF=size(trigData,1)/2;
    SizesF2=size(trigData,1);
    
    LD2(:,1)=sum(trigData(1:SizesF,:).*val{2,1}(Pairs(i,2),:),2)-sum(trigData(1:SizesF,:).*val{2,1}(Pairs(i,1),:),2);
        
    
    LDA2(:,1)=sum(trigData(SizesF+1:SizesF2,:).*val{7,1}(Pairs(i,2),:),2)-sum(trigData(SizesF+1:SizesF2,:).*val{7,1}(Pairs(i,1),:),2);
    
    lineData2(:,i)=[LD2;LDA2];
end

figure('Name','Powers');
plot(T,lineData);
xlabel('Time');
ylabel('PowerData');
title('Power');

figure('Name','Powers wTrig');

subplot(3,3,1)
plot(T,[lineData(:,1),lineData2(:,1)]);
xlabel('Time');
ylabel('PowerData');
title('Power');

subplot(3,3,2)
plot(T,[lineData(:,2),lineData2(:,2)]);
xlabel('Time');
ylabel('PowerData');
title('Power');

subplot(3,3,3)
plot(T,[lineData(:,3),lineData2(:,3)]);
xlabel('Time');
ylabel('PowerData');
title('Power');

subplot(3,3,4)
plot(T,[lineData(:,4),lineData2(:,4)]);
xlabel('Time');
ylabel('PowerData');
title('Power');

subplot(3,3,5)
plot(T,[lineData(:,5),lineData2(:,5)]);
xlabel('Time');
ylabel('PowerData');
title('Power');


subplot(3,3,6)
plot(T,[lineData(:,6),lineData2(:,6)]);
xlabel('Time');
ylabel('PowerData');
title('Power');

subplot(3,3,7)
plot(T,[lineData(:,7),lineData2(:,7)]);
xlabel('Time');
ylabel('PowerData');
title('Power');

subplot(3,3,8)
plot(T,[lineData,lineData2]);
xlabel('Time');
ylabel('PowerData');
title('Power');


figure('Name','Powers Trig');
plot(T,lineData2);
xlabel('Time');
ylabel('PowerData');
title('Power');

err = abs(lineData - lineData2)./lineData;

flowER=(log10(err(:,1:n)*100));

figure('Name','Power Error');
plot(T,flowER);
xlabel('Time');
ylabel('Error');
title('Flow Error');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [maxt,Y] = MaxTrig(n,omega,zeta,const) 
    maxt= (pi*n)/(omega*sqrt(1-zeta^2));
    expo=exp(-1*zeta*omega*maxt)/sqrt(1-zeta^2);
    trig=sin(omega*sqrt(1-zeta^2)*maxt-acos(zeta));
    Y=const*(1+(expo*trig));
end

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
