A=[2 3;0 1];
S=[12 14;13 15];
BLA=S*A;
SI= inv(S);
lamda = BLA\S;
[V,D]= eig(lamda)