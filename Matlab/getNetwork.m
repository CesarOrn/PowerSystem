function [A,P]=getNetwork(name)
%name: Selecte network to get.
%option:
    %network1
    %network2
    %default: network1

%A: return coupling matrix.
KOUPLING=140;
P=zeros(6,1);
        
P(1,1)=  40;
P(2,1)= -20;
P(3,1)= -25;
P(4,1)=  40;
P(5,1)= -25;
P(6,1)= -10;

    if(name=='network1')
        A=zeros(6);
        A(1,2)=KOUPLING;
        A(2,1)=KOUPLING;
        
        A(3,2)=KOUPLING;
        A(2,3)=KOUPLING;

        A(1,3)=KOUPLING;
        A(3,1)=KOUPLING;

        A(1,5)=KOUPLING;
        A(5,1)=KOUPLING;

        A(4,3)=KOUPLING;
        A(3,4)=KOUPLING;

        A(4,5)=KOUPLING;
        A(5,4)=KOUPLING;

        A(6,5)=KOUPLING;
        A(5,6)=KOUPLING;
       
        return ;
    end
    
        if(name=='network2')
        A=zeros(6);
        A(1,2)=KOUPLING;
        A(2,1)=KOUPLING;
        
        A(3,2)=KOUPLING;
        A(2,3)=KOUPLING;

        A(1,3)=KOUPLING;
        A(3,1)=KOUPLING;

        A(1,5)=KOUPLING;
        A(5,1)=KOUPLING;

        A(4,3)=KOUPLING;
        A(3,4)=KOUPLING;

        A(4,5)=KOUPLING;
        A(5,4)=KOUPLING;

        A(6,5)=KOUPLING;
        A(5,6)=KOUPLING;
        
        A(6,1)=KOUPLING;
        A(1,6)=KOUPLING;
        return
        end
    
        A=zeros(6);
        A(1,2)=KOUPLING;
        A(2,1)=KOUPLING;
        
        A(3,2)=KOUPLING;
        A(2,3)=KOUPLING;

        A(1,3)=KOUPLING;
        A(3,1)=KOUPLING;

        A(1,5)=KOUPLING;
        A(5,1)=KOUPLING;

        A(4,3)=KOUPLING;
        A(3,4)=KOUPLING;

        A(4,5)=KOUPLING;
        A(5,4)=KOUPLING;

        A(6,5)=KOUPLING;
        A(5,6)=KOUPLING;

end
