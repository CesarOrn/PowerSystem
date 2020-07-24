function[val,initsize,psize]=perturbation(coupling,Time,initPower,powerlist,initCon,initlist)


%val{1,k}=EigenVector;
%val{2,k}=Lambdas;
%val{3,k}=Q;
%val{4,k}=LastValues;
      
%val{5,k}=Condition
%val{6,k}=Power

%val{8,k}=Values;
%val{9,k}=Time Values.

k=1;
gamma=0.9;

% number of elements
n=size(coupling,1);

% Create the Laplacian
L = diag(sum(coupling,2)) - coupling;


% Perform diagonalization
[V,Lambda] = eig(L);
%V(V<1e-6)=0;
%Lambda(Lambda<1e-6)=0;


% Run the diagonalized dynamics
%
%  To diagonalize, pre-multiply by [ inv(V)   O_n   ]
%                                  [  O_n    inv(V) ]
%
Q = V \ initPower;
Q(Q<1e-6)=0;
diagonalized_dynamics = @(t,x) [zeros(n), eye(n); -Lambda, -gamma*eye(n)] * x + [zeros(n,1); Q];

% Run ode45 for diagonalized dynamics
x0 = initCon;
[T, X] = ode45(diagonalized_dynamics,Time, x0);

val{1,k}=V;
val{2,k}=Lambda;
val{3,k}=Q;



val{4,k}=X(size(X,1),:);

      
val{5,k}=x0;
val{6,k}=initPower;

val{7,k}=X;
val{8,k}=T;


diagonalized_dynamics = @(t,x) [zeros(n), eye(n); -Lambda, -gamma*eye(n)] * x + [zeros(n,1); Q];

x0=X(size(X,1),:);

[T2, X2] = ode45(diagonalized_dynamics,Time, x0);

val{9,k}=X2;
val{10,k}=[T2+T(size(T,1),1)];
val{11,k}=sum(V.*(Q./sum(Lambda,2)),2);
val{12,k}=sum(V.*(Q./sum(Lambda,2)),2);


% inital condition changes
initsize=size(initlist,1)
for i=1:initsize
    k=k+1;
    val{1,k}=V;
    val{2,k}=Lambda;
    val{3,k}=Q;
    val{4,k}=X(size(X,1),:);

      
    val{5,k}=initlist(i,:);
    val{6,k}=initPower';
    
    val{7,k}=X;
    val{8,k}=T;
    
    
    x0=X(size(X,1),:)+initlist(i,:);
    [T2, X2] = ode45(diagonalized_dynamics,Time, x0);

    val{9,k}=X2;
    val{10,k}=[T2+T(size(T,1),1)];
    val{11,k}=sum(V.*(Q./sum(Lambda,2)),2);
    val{12,k}=sum(V.*(Q./sum(Lambda,2)),2);
    
end


% power changes
psize=size(powerlist,2)
for i=1:psize
    k=k+1;
    
    power= powerlist(:,i);
    QN = V \ power;
    QN(QN<1e-6)=0;
    
    diagonalized_dynamics = @(t,x) [zeros(n), eye(n); -Lambda, -gamma*eye(n)] * x + [zeros(n,1); QN];
    
    x0=X(size(X,1),:);
    
    
    val{1,k}=V;
    val{2,k}=Lambda;
    val{3,k}=QN;
    val{4,k}=X(size(X,1),:);

      
    val{5,k}=x0;
    val{6,k}=powerlist(:,i);

    val{7,k}=X;
    val{8,k}=T;
    
    
    [T2, X2] = ode45(diagonalized_dynamics,Time, x0);

    val{9,k}=X2;
    val{10,k}=[T2+T(size(T,1),1)];
    
val{11,k}=sum(V.*(Q./sum(Lambda,2)),2);
val{12,k}=sum(V.*(QN./sum(Lambda,2)),2);
end


end