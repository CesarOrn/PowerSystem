clear;
[A,P]=toyNetwork('toyfull')

[V,e]=eig(A)

%indecator matrix
E(1,1)=1;
E(5,1)=1;

E(2,2)=1;
E(4,2)=1;

E(3,3)=1;
E(6,4)=1;
E(7,5)=1;
E(10,7)=1;

E(9,6)=1;
E(8,6)=1;
E

E'*V(:,1)


%E
E'*A*E


EE=(E'*E)^-1*E'

Q=EE*A*E
[V,e]=eig(Q)

P
%Quotion Power. Should it be EE *P *E ???
QP= EE*P

[BA,BP]=symmetricNetwork('smallNetwork')

sum(P)