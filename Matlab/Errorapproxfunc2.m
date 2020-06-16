function [SSerror,SSerror2,l]=Errorapproxfunc2(dim,Coeff,SS)
%To determine the error in taking dominant etas VS summation of all etas
%SSerror is the final steady state error when the while loop stops
%SSerror2 is the matrix of all steady state errors arranged in cloumns i.e
%1 first column has SS error when l=1, 2nd column when l=2 and so on.
%l is the index on which loop stops i.e the number of dominant etas taken
%for error to be less than 5%

n=dim(1,1);
m=dim(1,2);


                %           STEADY STATE
%Multiply coefficient to steady State value


Linsm=zeros([n,m]);
for j=1:m
    Linsm(:,j)=Coeff(:,j).*SS;
end

% Sum each of the columns to attain summation value

SSsum=zeros([1,m]);
for i = 1:m
    SSsum(1,i)=sum(Linsm(:,i));
end
    
%For dominant eta
SSerror=10; %arbitrary value above 5 to start the loop
l=1; %index for number of dominant 
while max(SSerror)>5   
   Sdom=zeros([1,m]);
   
       
    for k = 1:m 
        Summ(:,k)=SS.*Coeff(:,k); %multiply steady stae value and coefficient
    end
    
    for i =1:m
        [~,idx]=sort(abs(Summ(:,i)), 'Descend');
        Summ2(:,i)=Summ(idx,i);
    end
    
    for p=1:m
        Sdom(:,p)=sum(Summ2(1:l,p)); %sum all the multiples
    end
    
    SSerror= abs(((SSsum-Sdom)./SSsum).*100);%Steady state error
    SSerror2(:,l)=SSerror;
    l=l+1;  
end
l=l-1; %end of loop has 1 added so we have to adjust


