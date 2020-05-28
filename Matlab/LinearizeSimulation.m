clear; %clc; close all;
coupling = zeros(6,6);
power = zeros(6,1);
KOUPLING=140;

%coupling matrix
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


%power of each node.
power(1,1)= 40;
power(2,1)= -20;
power(3,1)= -25;
power(4,1)= 40;
power(5,1)= -25;
power(6,1)= -10;

L=laplacian(coupling);

%calcualte eigenvalues and eigenvector Matrix
[V,eig]=eig(L);
Q= V\power;
damping=0.9;

%summulate the ddot(eta)= gamma *dot(eta) + eigenvalues *eta +Q

[t,y] = ode45(@(t,y) f(t,y,damping,eig,Q),[0 20],[0.01 0.00 0.10 0 0.01 0 0.01 0 0 0.01 0 0]);

save('test')
size(y,1)

%calculates each of the theta by summing the etas
theta=zeros(size(y,1),size(y,2)/2);
for i=1:size(y,1)
    theta(i,1)=dot(y(i,[1,3,5,7,9,11]),V(1,:));
    theta(i,2)=dot(y(i,[1,3,5,7,9,11]),V(2,:));
    theta(i,3)=dot(y(i,[1,3,5,7,9,11]),V(3,:));
    theta(i,4)=dot(y(i,[1,3,5,7,9,11]),V(4,:));
    theta(i,5)=dot(y(i,[1,3,5,7,9,11]),V(5,:));
    theta(i,6)=dot(y(i,[1,3,5,7,9,11]),V(6,:));
end

color=get(gca,'colororder');

% plot each eta individualy and save as png.
for i=0:5

    fig =figure(4);
    name='eta';
    %plot(t,y(:,1),t,y(:,3),t,y(:,5),t,y(:,7),t,y(:,9),t,y(:,11),'LineWidth',1.25);
    plot(t,y(:,2*i+1),'LineWidth',1.25,'Color',color(i+1,:));
    axis([0 20 -0.6 0.4])
    xlabel('Time(s)') ;
    ylabel('Eta') ;
    title(strcat('Eta ',num2str(i+1),' Plot'));
    name =strcat(name,num2str(i+1),'.png');
    saveas(fig,name);
end

% plot all eta  and save as png.
fig =figure(4);
name='etas';
plot(t,y(:,1),t,y(:,3),t,y(:,5),t,y(:,7),t,y(:,9),t,y(:,11),'LineWidth',1.25);
xlabel('Time(s)') ;
ylabel('Eta') ;
title(strcat('Eta Plot'));
name =strcat(name,'.png');
saveas(fig,name);

% plot each theta individualy
for i=1:6
    %save figure to output to file later;
    fig=figure(5);

    name='theta';

    %plot the theta and configure plot
    plot(t,theta(:,i),'LineWidth',1.25,'Color',color(i,:));
    axis([0 20 -0.4 0.4])
    xlabel('Time(s)') ;
    ylabel('Theta') ;
    title(strcat('Theta ',num2str(i),' Plot'));

    %name of output file.
    name =strcat(name,num2str(i),'.png');
    saveas(fig,name);
end

% plot all theta and save as png.
name='theta';
fig=figure(5);
plot(t,theta(:,1),t,theta(:,2),t,theta(:,3),t,theta(:,4),t,theta(:,5),t,theta(:,6),'LineWidth',1.25);
xlabel('Time(s)') ;
ylabel('Theta') ;
title(strcat('Theta Plot'));
name =strcat(name,'.png');
saveas(fig,name);


%Define an edge in software with names.
edge =["theta 2 - theta 1" "theta 3 - theta 1" "theta 5 - theta 1" "theta 3 - theta 4" "theta 5 - theta 4" "theta 5 - theta 6"];

lines(1,:)= theta(:,1)-theta(:,2);
lines(2,:)= theta(:,1)-theta(:,3);
lines(3,:)= theta(:,1)-theta(:,5);
lines(4,:)= theta(:,4)-theta(:,3);
lines(5,:)= theta(:,4)-theta(:,5);
lines(6,:)= theta(:,5)-theta(:,6);

% plot each Theta_j-Theta_j individualy and save as png.
for i= 1:6
fig=figure(6);

name='Theta';

plot(t,sin(lines(i,:)),'LineWidth',1.25,'Color',color(i,:))
  xlabel('Time(s)') ;
  ylabel('Theta') ;
  title(edge(i));
  axis([0 20 -0.2 0.6])
  name =strcat(edge(i),'.png');
  saveas(fig,name);
end

% plot all Theta_j-Theta_j  and save as png.
fig=figure(6);
plot(t,sin(theta(:,1)-theta(:,2)),t,sin(theta(:,1)-theta(:,3)), t,sin(theta(:,1)-theta(:,5)), t,sin(theta(:,4)-theta(:,3)),t,sin(theta(:,4)-theta(:,5)),t,sin(theta(:,5)-theta(:,6)),'LineWidth',1.25);
title('Theta_j-Theta_i');
name =strcat('Theta','.png');
saveas(fig,name);



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
