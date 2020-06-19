%two componets

path=pwd;
path=split(path,"/");
pathSize=size(path,1);
path(pathSize,1)={'Chart'};
path(pathSize+1,1)={'Setup4_components'};

pathComp2=path;
pathSize=size(pathComp2,1);
pathComp2(pathSize+1,1)={'Components2'};
pathComp2(pathSize+2,1)={'Data.xlsm'};
pathComp2=string(join(pathComp2,"/"));

pathComp1=path;
pathSize=size(pathComp1,1);
pathComp1(pathSize+1,1)={'Components1'};
pathComp1(pathSize+2,1)={'Data.xlsm'};
pathComp1=string(join(pathComp1,"/"));



[A_2,P_2]=compNetwork('twoComp')
L_2 = diag(sum(A_2,2)) - A_2;
[eig_2,Lambda_2]=eig(L_2)
Q_2=eig_2\P_2
Q_2=round(Q_2,3);
eig_2=round(eig_2,3);
Lambda_2=round(Lambda_2,3);

writematrix(eig_2,pathComp2,'Sheet',1,'Range','A1');
writematrix("Eigen Vector",pathComp2,'Sheet',1,'Range','A8');

writematrix(Lambda_2,pathComp2,'Sheet',1,'Range','I1');
writematrix("Eigen Values",pathComp2,'Sheet',1,'Range','I8');

writematrix(P_2,pathComp2,'Sheet',1,'Range','R1');
writematrix("Power",pathComp2,'Sheet',1,'Range','R8');

writematrix(Q_2,pathComp2,'Sheet',1,'Range','T1');
writematrix("Q",pathComp2,'Sheet',1,'Range','T8');


[A_1,P_1]=compNetwork('oneComp');
L_1 = diag(sum(A_1,2)) - A_1;
[eig_1,Lambda_1]=eig(L_1)
Q_1=eig_1\P_1
Q_1=round(Q_1,3);
eig_1=round(eig_1,3);
Lambda_1=round(Lambda_1,3);

writematrix(eig_1,pathComp1,'Sheet',1,'Range','A1');
writematrix("Eigen Vector",pathComp1,'Sheet',1,'Range','A8');

writematrix(Lambda_1,pathComp1,'Sheet',1,'Range','I1');
writematrix("Eigen Values",pathComp1,'Sheet',1,'Range','I8');

writematrix(P_1,pathComp1,'Sheet',1,'Range','R1');
writematrix("Power",pathComp1,'Sheet',1,'Range','R8');

writematrix(Q_1,pathComp1,'Sheet',1,'Range','T1');
writematrix("Q",pathComp1,'Sheet',1,'Range','T8');