clear;
[A,P]=symmetricNetwork('fullNetwork')

%indecator matrix
E(1,1)=1;
E(2,1)=1;

E(3,2)=1;
E(4,2)=1;
E(5,2)=1;
E(6,2)=1;

E(7,3)=1;
E(8,3)=1;
E(9,3)=1;
E(10,3)=1;
E(11,3)=1;
E(12,3)=1;

%E

EE=(E'*E)^-1*E'

Q=EE*A*E

P
%Quotion Power. Should it be EE *P *E ???
QP= EE*P

[BA,BP]=symmetricNetwork('smallNetwork')

sum(P)