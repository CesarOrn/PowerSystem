
coupling = zeros(6,6);
power = zeros(6,1);
KOUPLING=140;

coupling(1,2)=KOUPLING;
coupling(2,1)=KOUPLING;

coupling(1,3)=KOUPLING;
coupling(3,1)=KOUPLING;

coupling(1,5)=KOUPLING;
coupling(5,1)=KOUPLING;

coupling(4,3)=KOUPLING;
coupling(3,4)=KOUPLING;

coupling(4,5)=KOUPLING;
coupling(5,4)=KOUPLING;

coupling(6,5)=KOUPLING;
coupling(5,6)=KOUPLING;



power(1)= 40;
power(2)= -20;
power(3)= -25;
power(4)= 40;
power(5)= -25;
power(6)= -10;

test=laplacian(coupling);

[V,D]=eig(test);
Q= power\V;
damping=0.9;

[t,y] = ode45(@(t,y) f(t,y,damping,D,Q),[0 20],[0 0 0 0 0 0 0 0 0 0 0 0]);

save('test')
size(y,1)
theta=zeros(size(y,1),size(y,2)/2);
for i=1:size(y,1)
    theta(i,1)=dot(y(i,[1,3,5,7,9,11]),V(1,:));
    theta(i,2)=dot(y(i,[1,3,5,7,9,11]),V(2,:));
    theta(i,3)=dot(y(i,[1,3,5,7,9,11]),V(3,:));
    theta(i,4)=dot(y(i,[1,3,5,7,9,11]),V(4,:));
    theta(i,5)=dot(y(i,[1,3,5,7,9,11]),V(5,:));
    theta(i,6)=dot(y(i,[1,3,5,7,9,11]),V(6,:));
end
figure(4)
plot(t,y(:,1),t,y(:,3),t,y(:,5),t,y(:,7),t,y(:,9),t,y(:,11),'LineWidth',1.25)
%plot(t,theta(:,1),t,theta(:,2),t,theta(:,3),t,theta(:,4),t,theta(:,5),t,theta(:,6),'LineWidth',1.25)

%plot(t,KOUPLING*sin(theta(:,1)-theta(:,2)),t,KOUPLING*sin(theta(:,1)-theta(:,3)), t,KOUPLING*sin(theta(:,1)-theta(:,5)), t,KOUPLING*sin(theta(:,4)-theta(:,3)),t,KOUPLING*sin(theta(:,4)-theta(:,5)),t,KOUPLING*sin(theta(:,5)-theta(:,6)),'LineWidth',1.25);

%plot(t,KOUPLING*sin(theta(:,6)-theta(:,5)));


function dydt = f(t,y,damp,eig,Q)


dydt(1) = y(2);
dydt(2) = -1*damp*y(2)-eig(1,1)*y(1)+Q(1);

dydt(3) =y(4);
dydt(4) = -1*damp*y(4)-eig(2,2)*y(3)+Q(2);

dydt(5) = y(6);
dydt(6) = -1*damp*y(6)-eig(3,3)*y(5)+Q(3);

dydt(7) = y(8);
dydt(8) = -1*damp*y(8)-eig(4,4)*y(7)+Q(4);

dydt(9) = y(10);
dydt(10) = -1*damp*y(10)-eig(5,5)*y(9)+Q(5);

dydt(11) = y(12);
dydt(12) = -1*damp*y(12)-eig(6,6)*y(11)+Q(6);

dydt=dydt(:);
end

function lap=laplacian(x)
    lap = zeros(size(x,1),size(x,2));

    for i=1:size(x,1)
        for j=1:size(x,2)
            sum = -1*x(i,j);
            if i==j
              
                for k=1:size(x,1)
                sum = x(i,k)+sum;
                end
                
            end
            lap(i,j)=sum;
        end
    end
end


