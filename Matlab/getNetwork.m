function [A,P]=getNetwork(name)
%name: Selecte network to get.
%option:
    %network1
    %network2
    %network3 Bottle necknetwork.
    %default: network1

%A: return coupling matrix.
KOUPLING=1;
    % 7 Connection network.
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

        P=zeros(6,1);

        P(1,1)=  4;
        P(2,1)= -1;
        P(3,1)= -2;
        P(4,1)=  2;
        P(5,1)= -2;
        P(6,1)= -1;


        return ;
    end
        % 8 Connection network.
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

        P=zeros(6,1);

        P(1,1)=  40;
        P(2,1)= -20;
        P(3,1)= -25;
        P(4,1)=  40;
        P(5,1)= -25;
        P(6,1)= -10;

        return;

        end
        % Bottle neck network
        if(name=='network3')

        A=zeros(7);
        A(1,2)=1;
        A(2,1)=1;

        A(3,2)=1;
        A(2,3)=1;

        A(1,3)=1;
        A(3,1)=1;

        A(2,4)=1;
        A(4,2)=1;

        A(4,3)=1;
        A(3,4)=1;

        A(2,5)=1;
        A(5,2)=1;

        A(6,5)=1;
        A(5,6)=1;

        A(7,5)=1;
        A(5,7)=1;

        A(7,6)=1;
        A(6,7)=1;

        P=zeros(7,1);

        P(1,1)=  3.0;
        P(2,1)= -1.0;
        P(3,1)= -1.0;
        P(4,1)=  1.0;
        P(5,1)= -1.0;
        P(6,1)= -2.0;
        P(7,1)=  1.0;


        return;
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

        P=zeros(6,1);

        P(1,1)=  40;
        P(2,1)= -20;
        P(3,1)= -25;
        P(4,1)=  40;
        P(5,1)= -25;
        P(6,1)= -10;


end
