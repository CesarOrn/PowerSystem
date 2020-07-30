function [val]=GoldenSectionSearch(f,a,b,tol,fun,eig,lines)

gr=(sqrt(5)+1)/2;



c=b-((b-a)/gr);
d=a+((b-a)/gr);
total=0;
while(abs(c-d)>tol && total<1000)
    total=total+1;
    [c d]
    %[f(c,fun,eig,lines)     f(d,fun,eig,lines)]
    %pause(1)
    if(-1*f(c,fun,eig,lines)<-1*f(d,fun,eig,lines))
       b=d; 
    else
       a=c;
    end
    
    c=b-((b-a)/gr);
    d=a+((b-1)/gr);

end

val= (b+a)/2;
end