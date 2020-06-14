function [SSerror,SSerror2,l]=Errorapproxfunc(dim,Coeff,SS)
%To determine the error in taking a dominant eta VS summation of all etas


n=dim(1,1);
m=dim(1,2);


                %           STEADY STATE
%Multiply coefficient to steady State value
Linsm=Coeff.*SS;


% Sum each of the columns to attain summation value
SSsum =sum(Linsm) ;

%For dominant eta, take just the coefficient of eta multiplied by the
%steady state
SSerror=10;
l=0;
while max(SSerror)>5   
    l=l+1;
   Sdom=zeros([1,m]);
   
   [~,idx]=sort(abs(SS), 'Descend');    
        SS2=SS(idx);
        Coeff2=Coeff(idx,:);

        
    for k = 1:m 
        Summ(1:l,k)=SS2(1:l).*Coeff2(1:l,k);
    end
    
    for p=1:m
        Sdom(:,p)=sum(Summ(:,p));
    end
    
    SSerror= abs(((SSsum-Sdom)./SSsum).*100);%Steady state error
    SSerror2(:,l)=SSerror;
      
end

end


