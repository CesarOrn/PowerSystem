clear; clc; close all;
%
% Consider n coupled second order dynamics 
% d^2 xi/dt^2 = -gamma dxi/dt + SUM_j A(i,j) * (x(j)-x(i)) + Pi;
%
% To use ode45, we convert the n second order equations 2n first order equations 


% Dimension of the system 
%n =10;

% Create a (random) matrix
[A,P]=symmetricNetwork('star');
A
%A = A-diag(diag(A));
%A= A+transpose(A);

n=size(A,1)

% Create the Laplacian
L = diag(sum(A,2)) - A;

% Damping
gamma = 0.9;

% Forcing vector 
%P = rand(n,1);
%P = P -mean(P);

%Set time steps;
Time = 0:.0018:20;

% Perform diagonalization

initPower=P;
numPert=10;
powerlist=zeros(n,10);
for i=1:numPert
    size(powerlist);
    powerlist(:,i)=P;
    powerlist(randi(n),i)= randi(3);
end

initCon=zeros(2*n,1);
initlist=rand(numPert,2*n);

[pert,initsize,psize]=perturbation(A,Time,initPower,powerlist,initCon,initlist);
q_i=pert{3,1}./sum(pert{2,1},2);

%[linePert,dataPert,SSP]= findLines(A,pert{1,1},q_i);

for i=1:size(A,1)
   etaSt{i,:}=strcat('eta',num2str(i,'%i'));
end

%val{1,k}=EigenVector;
%val{2,k}=Lambdas;
%val{3,k}=Q;
%val{4,k}=LastValues;
      
%val{5,k}=Condition
%val{6,k}=Power

%val{7,k}=Values;
%val{8,k}=Time Values.

%val{7,k}=Values;
%val{8,k}=Time Values.

OL=sum(pert{2,1},2)
OQ=sum(pert{3,1},2)

IC=pert{5,1}
IP=pert{6,1}

writecell(etaSt,'Data.xlsx','Sheet',1,'Range','A2');
%writematrix(cell2mat(SSP),'Data.xlsx','Sheet',1,'Range','D1');

writematrix("Lambda",'Data.xlsx','Sheet',1,'Range',' B1');
writematrix(sum(pert{2,1},2),'Data.xlsx','Sheet',1,'Range',' B2');

writematrix("Q",'Data.xlsx','Sheet',1,'Range',' D1');
writematrix(pert{3,1},'Data.xlsx','Sheet',1,'Range',' D2');

writematrix("InitCondition",'Data.xlsx','Sheet',1,'Range',' F1');
writematrix(pert{5,1},'Data.xlsx','Sheet',1,'Range',' F2');

writematrix("InitPower",'Data.xlsx','Sheet',1,'Range',' H1');
writematrix(pert{6,1},'Data.xlsx','Sheet',1,'Range',' H2');

writematrix("SS",'Data.xlsx','Sheet',1,'Range',' J1');
writematrix(pert{11,1},'Data.xlsx','Sheet',1,'Range',' J2');

writecell(etaSt,'Data.xlsx','Sheet',1,'Range','L2');
%writematrix(cell2mat(SSP),'Data.xlsx','Sheet',1,'Range','D1');

writematrix("Lambda",'Data.xlsx','Sheet',1,'Range',' M1');
writematrix(sum(pert{2,1},2),'Data.xlsx','Sheet',1,'Range',' M2');

writematrix("Q",'Data.xlsx','Sheet',1,'Range',' O1');
writematrix(pert{3,1},'Data.xlsx','Sheet',1,'Range',' O2');

writematrix("InitCondition",'Data.xlsx','Sheet',1,'Range',' Q1');
writematrix(pert{5,1},'Data.xlsx','Sheet',1,'Range',' Q2');

writematrix("InitPower",'Data.xlsx','Sheet',1,'Range',' S1');
writematrix(pert{6,1},'Data.xlsx','Sheet',1,'Range',' S2');

writematrix("SS",'Data.xlsx','Sheet',1,'Range',' U1');
writematrix(pert{12,1},'Data.xlsx','Sheet',1,'Range',' U2');


ylim([-2 5])
plot([pert{8,1};pert{10,1}],[pert{7,1}(:,1:n); pert{9,1}(:,1:n)]);
title('No Perturbation')
xlabel('Time(S)') 
ylabel('Eta') 
saveas(gcf,strcat('Data','.png'));




initsize
for i=1:initsize
    index=1+i;
    baseName=strcat('Data',num2str(i,'%i'),'I');
    name= strcat(baseName,'.xlsx');

    writecell(etaSt,name,'Sheet',1,'Range','A2');
%writematrix(cell2mat(SSP),'Data.xlsx','Sheet',1,'Range','D1');

writematrix("Lambda",name,'Sheet',1,'Range',' B1');
writematrix(OL,name,'Sheet',1,'Range',' B2');

writematrix("Q",name,'Sheet',1,'Range',' D1');
writematrix(OQ,name,'Sheet',1,'Range',' D2');

writematrix("InitCondition",name,'Sheet',1,'Range',' F1');
writematrix(IC,name,'Sheet',1,'Range',' F2');

writematrix("InitPower",name,'Sheet',1,'Range',' H1');
writematrix(IP,name,'Sheet',1,'Range',' H2');

writematrix("SS",name,'Sheet',1,'Range',' J1');
writematrix(pert{11,index},name,'Sheet',1,'Range',' J2');

writecell(etaSt,name,'Sheet',1,'Range','L2');
%writematrix(cell2mat(SSP),'Data.xlsx','Sheet',1,'Range','D1');

writematrix("Lambda",name,'Sheet',1,'Range',' M1');
writematrix(sum(pert{2,index},2),name,'Sheet',1,'Range',' M2');

writematrix("Q",name,'Sheet',1,'Range',' O1');
writematrix(pert{3,index},name,'Sheet',1,'Range',' O2');

writematrix("Condition",name,'Sheet',1,'Range',' Q1');
writematrix(pert{5,index}',name,'Sheet',1,'Range',' Q2');

writematrix("PowerPert",name,'Sheet',1,'Range',' S1');
writematrix(pert{6,index}',name,'Sheet',1,'Range',' S2');

writematrix("SS",name,'Sheet',1,'Range',' U1');
writematrix(pert{12,index},name,'Sheet',1,'Range',' U2');
    
    ylim([-2 5])
    plot([pert{8,index}; pert{10,index}],[pert{7,index}(:,1:n); pert{9,index}(:,1:n)]);
    title('Inital Condition Perturbation')
    xlabel('Time(S)') 
    ylabel('Eta') 
    saveas(gcf,strcat(baseName,'.png'));

    
end
psize
for i=1:psize
    index=1+initsize+i;
    baseName=strcat('Data',num2str(i,'%i'),'P');
    name= strcat(baseName,'.xlsx');

    
    writecell(etaSt,name,'Sheet',1,'Range','A2');
%writematrix(cell2mat(SSP),'Data.xlsx','Sheet',1,'Range','D1');

writematrix("Lambda",name,'Sheet',1,'Range',' B1');
writematrix(OL,name,'Sheet',1,'Range',' B2');

writematrix("Q",name,'Sheet',1,'Range',' D1');
writematrix(OQ,name,'Sheet',1,'Range',' D2');

writematrix("InitCondition",name,'Sheet',1,'Range',' F1');
writematrix(IC,name,'Sheet',1,'Range',' F2');

writematrix("InitPower",name,'Sheet',1,'Range',' H1');
writematrix(IP,name,'Sheet',1,'Range',' H2');

writematrix("SS",name,'Sheet',1,'Range',' J1');
writematrix(pert{11,index},name,'Sheet',1,'Range',' J2');

writecell(etaSt,name,'Sheet',1,'Range','L2');
%writematrix(cell2mat(SSP),'Data.xlsx','Sheet',1,'Range','D1');

writematrix("Lambda",name,'Sheet',1,'Range',' M1');
writematrix(sum(pert{2,index},2),name,'Sheet',1,'Range',' M2');

writematrix("Q",name,'Sheet',1,'Range',' O1');
writematrix(pert{3,index},name,'Sheet',1,'Range',' O2');

writematrix("Condition",name,'Sheet',1,'Range',' Q1');
writematrix(pert{5,index}',name,'Sheet',1,'Range',' Q2');

writematrix("PowerPert",name,'Sheet',1,'Range',' S1');
writematrix(pert{6,index},name,'Sheet',1,'Range',' S2');

writematrix("SS",name,'Sheet',1,'Range',' U1');
writematrix(pert{12,index},name,'Sheet',1,'Range',' U2');
    
    plot([pert{8,index}; pert{10,index}],[pert{7,index}(:,1:n); pert{9,index}(:,1:n)]);
    ylim([-2 5])
    title('Power Perturbation')
    xlabel('Time(S)') 
    ylabel('Eta') 
    saveas(gcf,strcat(baseName,'.png'));
    
end

T=[pert{8,1}; pert{10,1}];
I=[pert{7,2}(:,1:n); pert{9,2}(:,1:n)]
    
PP=[pert{7,3+psize}(:,1:n); pert{9,3+psize}(:,1:n)];

G = table(T,I);

writetable(G,'powerPert.txt');
type powerPert.txt;

K = table(T,PP);

writetable(K,'InitPert.txt');
type InitPert.txt;

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

E

EE=(E'*E)^-1*E'

Q=EE*A*E

P
%Quotion Power. Should it be EE *P *E ???
QP= EE*P


val=DDwFailure(Q,gamma,QP,Time,1);
initPowerq=QP;
numPert=10;
for i=1:numPert
    size(powerlist);
    powerlist(:,i)
    powerlistq(:,i)=EE*powerlist(:,i)
end

size(initlist,1)
initConq=[EE*initCon(1:n); EE*initCon(n+1:2*n)]
for i=1:size(initlist,1)
   initlistq(i,:)=[EE*initlist(i,1:n)'; EE*initlist(i,n+1:2*n)'] ;
end

[pert,initsizeq,psizeq]=perturbation(Q,Time,initPowerq,powerlistq,initConq,initlistq);
q_i=pert{3,1}./sum(pert{2,1},2);

%[linePert,dataPert,SSP]= findLines(Q,pert{1,1},q_i);

for i=1:size(Q,1)
   etaStQ{i,:}=strcat('eta',num2str(i,'%i'))
end
nq=size(Q,1);


OL=sum(pert{2,1},2)
OQ=pert{3,1}

IC=pert{5,1}
IP=pert{6,1}




initsizeq
for i=1:initsizeq
    index=1+i;
    %[linePert,dataPert,SSP]= findLines(A,pert{1,index},pert{3,index}./sum(pert{2,index},2));
    baseName=strcat('Data',num2str(i,'%i'),'IQ');
    name= strcat(baseName,'.xlsx');

    writecell(etaSt,name,'Sheet',1,'Range','A2');
%writematrix(cell2mat(SSP),'Data.xlsx','Sheet',1,'Range','D1');

writematrix("Lambda",name,'Sheet',1,'Range',' B1');
writematrix(OL,name,'Sheet',1,'Range',' B2');

writematrix("Q",name,'Sheet',1,'Range',' D1');
writematrix(OQ,name,'Sheet',1,'Range',' D2');

writematrix("InitCondition",name,'Sheet',1,'Range',' F1');
writematrix(IC,name,'Sheet',1,'Range',' F2');

writematrix("InitPower",name,'Sheet',1,'Range',' H1');
writematrix(IP,name,'Sheet',1,'Range',' H2');

writematrix("SS",name,'Sheet',1,'Range',' J1');
writematrix(pert{11,index},name,'Sheet',1,'Range',' J2');

writecell(etaSt,name,'Sheet',1,'Range','L2');
%writematrix(cell2mat(SSP),'Data.xlsx','Sheet',1,'Range','D1');

writematrix("Lambda",name,'Sheet',1,'Range',' M1');
writematrix(sum(pert{2,index},2),name,'Sheet',1,'Range',' M2');

writematrix("Q",name,'Sheet',1,'Range',' O1');
writematrix(pert{3,index},name,'Sheet',1,'Range',' O2');

writematrix("Condition",name,'Sheet',1,'Range',' Q1');
writematrix(pert{5,index}',name,'Sheet',1,'Range',' Q2');

writematrix("PowerPert",name,'Sheet',1,'Range',' S1');
writematrix(pert{6,index}',name,'Sheet',1,'Range',' S2');

writematrix("SS",name,'Sheet',1,'Range',' U1');
writematrix(pert{12,index},name,'Sheet',1,'Range',' U2');
    
    ylim([-2 5])
    plot([pert{8,index}; pert{10,index}],[pert{7,index}(:,1:nq); pert{9,index}(:,1:nq)]);
        title('Inital Condition Perturbation')
    xlabel('Time(S)') 
    ylabel('Eta') 
    saveas(gcf,strcat(baseName,'.png'));

    
end
psizeq
for i=1:psizeq
    index=1+initsize+i;
    %[linePert,dataPert,SSP]= findLines(A,pert{1,index},pert{3,index}./sum(pert{2,index},2));
    baseName=strcat('Data',num2str(i,'%i'),'PQ');
    name= strcat(baseName,'.xlsx');

    writecell(etaSt,name,'Sheet',1,'Range','A2');
%writematrix(cell2mat(SSP),'Data.xlsx','Sheet',1,'Range','D1');

writematrix("Lambda",name,'Sheet',1,'Range',' B1');
writematrix(OL,name,'Sheet',1,'Range',' B2');

writematrix("Q",name,'Sheet',1,'Range',' D1');
writematrix(OQ,name,'Sheet',1,'Range',' D2');

writematrix("InitCondition",name,'Sheet',1,'Range',' F1');
writematrix(IC,name,'Sheet',1,'Range',' F2');

writematrix("InitPower",name,'Sheet',1,'Range',' H1');
writematrix(IP,name,'Sheet',1,'Range',' H2');

writematrix("SS",name,'Sheet',1,'Range',' J1');
writematrix(pert{11,index},name,'Sheet',1,'Range',' J2');

writecell(etaSt,name,'Sheet',1,'Range','L2');
%writematrix(cell2mat(SSP),'Data.xlsx','Sheet',1,'Range','D1');

writematrix("Lambda",name,'Sheet',1,'Range',' M1');
writematrix(sum(pert{2,index},2),name,'Sheet',1,'Range',' M2');

writematrix("Q",name,'Sheet',1,'Range',' O1');
writematrix(pert{3,index},name,'Sheet',1,'Range',' O2');

writematrix("Condition",name,'Sheet',1,'Range',' Q1');
writematrix(pert{5,index}',name,'Sheet',1,'Range',' Q2');

writematrix("PowerPert",name,'Sheet',1,'Range',' S1');
writematrix(pert{6,index},name,'Sheet',1,'Range',' S2');

writematrix("SS",name,'Sheet',1,'Range',' U1');
writematrix(pert{12,index},name,'Sheet',1,'Range',' U2');
    
    plot([pert{8,index}; pert{10,index}],[pert{7,index}(:,1:nq); pert{9,index}(:,1:nq)]);
    ylim([-2 5])
        title('Power Perturbation')
    xlabel('Time(S)') 
    ylabel('Eta') 
    saveas(gcf,strcat(baseName,'.png'));
    
end

[pert{8,index}; pert{10,index}],[pert{7,index}(:,1:nq); pert{9,index}(:,1:nq)]
[pert{8,index}; pert{10,index}],[pert{7,index}(:,1:nq); pert{9,index}(:,1:nq)]

T=[pert{8,1}; pert{10,1}];
I=[pert{7,2}(:,1:nq); pert{9,2}(:,1:nq)]
    
PP=[pert{7,2+psize}(:,1:nq); pert{9,2+psize}(:,1:nq)];



G = table(T,I);

writetable(G,'powerQPert.txt');
type powerQPert.txt;

K = table(T,PP);

writetable(K,'InitQPert.txt');
type InitQPert.txt;
