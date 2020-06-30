function [A,P]=compNetwork(name)

%A: return coupling matrix
%P: returns power vector.

%Name: options
%woCompReduceR: Right side reduce to one node
%twoCompReduceM: bottle neck middle combined
%twoComp: break after two components
%oneComp,Default: Bottle neck with middle inteact

switch name
    
    case 'twoCompReduceR'
        P(1,1)=3;
        P(2,1)=-1;
        P(3,1)=-1;
        P(4,1)=1;
        P(5,1)=-2;
        
        A(1,2)=1;
        A(2,1)=1;
        
        A(1,3)=1;
        A(3,1)=1;
        
        A(2,3)=1;
        A(3,2)=1;
        
        A(2,4)=1;
        A(4,2)=1;
        
        A(3,4)=1;
        A(4,3)=1;
        
        A(2,5)=1;
        A(5,2)=1;     
        
    case 'twoCompReduceM'
        P(1,1)=3;
        P(2,1)=-2;
        P(3,1)=-1;
        P(4,1)=1;
        P(5,1)=-2;
        P(6,1)=1;
        
        A(1,2)=1;
        A(2,1)=1;
        
        A(1,3)=1;
        A(3,1)=1;
        
        A(2,3)=1;
        A(3,2)=1;
        
        A(2,4)=1;
        A(4,2)=1;
        
        A(3,4)=1;
        A(4,3)=1;
   
        A(2,5)=1;
        A(5,2)=1;
        
        A(6,2)=1;
        A(2,6)=1;
        
        A(6,5)=1;
        A(5,6)=1;
   


    
    case 'twoComp'
        
        P(1,1)=3;
        P(2,1)=-1;
        P(3,1)=-1;
        P(4,1)=1;
        P(5,1)=-1;
        P(6,1)=-2;
        P(7,1)=1;

        
        A(1,2)=1;
        A(2,1)=1;
        
        A(1,3)=1;
        A(3,1)=1;
        
        A(2,3)=1;
        A(3,2)=1;
        
        A(2,4)=1;
        A(4,2)=1;
        
        A(3,4)=1;
        A(4,3)=1;
        
        A(2,5)=0;
        A(5,2)=0;
        
        A(6,5)=1;
        A(5,6)=1;
        
        A(7,5)=1;
        A(5,7)=1;
        
        A(7,6)=1twoComp;
        A(6,7)=1;
        
    otherwise
        P(1,1)=3;
        P(2,1)=-1;
        P(3,1)=-1;
        P(4,1)=1;
        P(5,1)=-1;
        P(6,1)=-2;
        P(7,1)=1;

        
        A(1,2)=1;
        A(2,1)=1;
        
        A(1,3)=1;
        A(3,1)=1;
        
        A(2,3)=1;
        A(3,2)=1;
        
        A(2,4)=1;
        A(4,2)=1;
        
        A(3,4)=1;
        A(4,3)=1;
        
        A(2,5)=1;
        A(5,2)=1;
        
        A(6,5)=1;
        A(5,6)=1;
        
        A(7,5)=1;
        A(5,7)=1;
        
        A(7,6)=1;
        A(6,7)=1;

end