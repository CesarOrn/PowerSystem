%ode solver for power grid

K=140;

[t,y] = ode45(@ThetaOnly,[0,20],[0 0 0 0 0 0 0 0 0 0 0 0]);
power= zeros(size(y,1),6);
size(y,1)
size(y,2)
for i=1:size(y,1)
    a=[K*sin(y(i,1)-y(i,3)) K*sin(y(i,1)-y(i,5)) K*sin(y(i,7)-y(i,5)) K*sin(y(i,7)-y(i,9)) K*sin(y(i,1)-y(i,9)) K*sin(y(i,9)-y(i,11))];
    power(i,:)=a;
end

figure();
plot(t,y(:,[1,3,5,7,9,11]));
legend({'Generator1','Consumer1','Consumer2','Generator2','Consumer3','Consumer4'},'Location','southwest');
xlabel('Time');
ylabel('Theta') 

figure();
plot(t,power);
legend({'Line1','Line2','Line3','Line4','Line5','Line6'},'Location','southwest')
xlabel('Time');
ylabel('Power') 

function dy = Theta(t,y)
    D=.9;
    K=140;
    P=40;
    
    connected=1;
    
    %if (t>30)
      %  connected=0;
    %end
        
   
    % Generator 1
    dy(1) = y(2);
    dy(2) = P-D*y(2)+K*sin(y(3)-y(1))+K*sin(y(5)-y(1))+connected*K*sin(y(9)-y(1));
    
    % Consumer 1
    dy(3)= y(4);
    dy(4) = -20-D*y(4)+K*sin(y(1)-y(3));
    
    % Consumer 2
    dy(5)= y(6);
    dy(6) = -25-D*y(6)+K*sin(y(1)-y(5))+K*sin(y(7)-y(5));
    
    % Generator 2
    dy(7) = y(8);
    dy(8) = P-D*y(8)+K*sin(y(5)-y(7))+K*sin(y(9)-y(7));
    
    % Consumer 3
    dy(9) = y(10);
    dy(10) = -25-D*y(10)+connected*K*sin(y(1)-y(9))+K*sin(y(7)-y(9))+K*sin(y(11)-y(9));
    
    % Consumer 4
    dy(11) = y(12);
    dy(12) = -10-D*y(12)+K*sin(y(9)-y(11));
    dy = dy(:);
    
end

function dy = ThetaOnly(t,y)
    D=0.9;
    K=140;
    P=40;
    
    connected=1;
    
    %if (t>30)
     %   connected=0;
    %end
        
   
    % Generator 1
    dy(1) = y(2);
    dy(2) = P-D*y(2)+K*(y(3)-y(1))+K*(y(5)-y(1))+connected*K*(y(9)-y(1));
    
    % Consumer 1
    dy(3)= y(4);
    dy(4) = -20-D*y(4)+K*(y(1)-y(3));
    
    % Consumer 2
    dy(5)= y(6);
    dy(6) = -25-D*y(6)+K*(y(1)-y(5))+K*(y(7)-y(5));
    
    % Generator 2
    dy(7) = y(8);
    dy(8) = P-D*y(8)+K*(y(5)-y(7))+K*(y(9)-y(7));
    
    % Consumer 3
    dy(9) = y(10);
    dy(10) = -25-D*y(10)+connected*K*(y(1)-y(9))+K*sin(y(7)-y(9))+K*(y(11)-y(9));
    
    % Consumer 4
    dy(11) = y(12);
    dy(12) = -10-D*y(12)+K*(y(9)-y(11));
    dy = dy(:);
    
end



