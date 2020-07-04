function [A,P]=symmetricNetwork(option)

%A: Coupling matrix
%P: Power vector

%option:
%fullNetwork: no quotion matrix

%smallNetwork: quotiont matrix
%P = EE*P(Full network)


switch option
    case 'fullNetwork'
        
        A(1,3)=1;
        A(3,1)=1;
        
        A(1,5)=1;
        A(5,1)=1;
        
        A(4,2)=1;
        A(2,4)=1;
        
        A(2,6)=1;
        A(6,2)=1;
        
        A(4,3)=1;
        A(3,4)=1;
        
        A(3,5)=1;
        A(5,3)=1;
        
        A(7,3)=1;
        A(3,7)=1;
        
        A(3,8)=1;
        A(8,3)=1;
        
        A(9,3)=1;
        A(3,9)=1;
        
        A(4,6)=1;
        A(6,4)=1;
        
        A(4,7)=1;
        A(7,4)=1;
        
        A(8,4)=1;
        A(4,8)=1;
        
        A(4,9)=1;
        A(9,4)=1;
        
        A(6,5)=1;
        A(5,6)=1;
        
        A(10,5)=1;
        A(5,10)=1;
        
        A(5,11)=1;
        A(11,5)=1;
        
        A(5,12)=1;
        A(12,5)=1;
        
        A(6,10)=1;
        A(10,6)=1;
        
        A(6,11)=1;
        A(11,6)=1;
        
        A(6,12)=1;
        A(12,6)=1;
        
        A(7,8)=1;
        A(8,7)=1;
        
        A(7,9)=1;
        A(9,7)=1;
        
        A(8,9)=1;
        A(9,8)=1;
        
        A(10,11)=1;
        A(11,10)=1;
        
        A(10,12)=1;
        A(12,10)=1;
        
        A(11,12)=1;
        A(12,11)=1;
        
        %power
        P(1,1)=1;
        P(2,1)=1;
        
        P(3,1)=-2;
        P(4,1)=-2;
        P(5,1)=-2;
        P(6,1)=-2;
        
        P(7,1)=1;
        P(8,1)=1;
        P(9,1)=1;
        P(10,1)=1;
        P(11,1)=1;
        P(12,1)=1;
        
    case 'smallNetwork'
        
        A(1,2)=2;
        
        A(2,1)=1;
        A(2,2)=2;
        A(2,3)=3;
        
        A(3,2)=2;
        A(3,3)=2;
        
        P(1,1)=1;
        P(2,1)=-2;
        P(3,1)=1;
        
end



end
