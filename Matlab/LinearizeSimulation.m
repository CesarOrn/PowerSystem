clear; %clc; close all;
n=6
coupling = zeros(n,n);
power = zeros(n,1);
KOUPLING=140;

%coupling = full(BAgraph(n));

%coupling matrix
coupling(1,2)=KOUPLING;
coupling(2,1)=KOUPLING;

coupling(1,3)=0;
coupling(3,1)=0;

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

%power = 25*rand(n,1);
%power = power -mean(power)

L=laplacian(coupling);

%calcualte eigenvalues and eigenvector Matrix
[V,eig]=eig(L);

Q= V\power;
damping=0.9;

%summulate the ddot(eta)= gamma *dot(eta) + eigenvalues *eta +Q



% Run ode45 for diagonalized dynamics
x0 = zeros(2*n,1);
diagonalized_dynamics = @(t,x) [zeros(n), eye(n); -eig, -damping*eye(n)] * x + [zeros(n,1); Q];

[t, y] = ode45(diagonalized_dynamics, [0 20], x0);

save('test')
size(y,1)

%calculates each of the theta by summing the etas
theta=zeros(size(y,1),size(y,2)/2);
for i=1:size(y,1)
    theta(i,1)=dot(y(i,1:n),V(1,:));
    theta(i,2)=dot(y(i,1:n),V(2,:));
    theta(i,3)=dot(y(i,1:n),V(3,:));
    theta(i,4)=dot(y(i,1:n),V(4,:));
    theta(i,5)=dot(y(i,1:n),V(5,:));
    theta(i,6)=dot(y(i,1:n),V(6,:));
end

color=get(gca,'colororder');

% plot each eta individualy and save as png.
saveFigs('Eta','Eta',t,y,[0 20 -0.35 0.1],4,n)


% plot each theta individualy
saveFigs('Theta','Theta',t,theta,[0 20 -0.2 0.25],5,n)



%Define an edge in software with names.
edge =["theta 2 - theta 1" "theta 3 - theta 1" "theta 5 - theta 1" "theta 3 - theta 4" "theta 5 - theta 4" "theta 5 - theta 6"];
lines(:,1)= (theta(:,1)-theta(:,2));
lines(:,2)= (theta(:,1)-theta(:,3));
lines(:,3)= (theta(:,1)-theta(:,5));
lines(:,4)= (theta(:,4)-theta(:,3));
lines(:,5)= (theta(:,4)-theta(:,5));
lines(:,6)= (theta(:,5)-theta(:,6));
size(lines,2)

% plot each Theta_j-Theta_j individualy and save as png.
saveFigs('Theta',edge,t,lines,[0 20 -0.15 0.3],6,size(lines,2));



function saveFigs(Ylabel,Name,t,data,AX,numfig,numofplots)
type ='.fig'
% list of all colors
colors=get(gca,'colororder');

    num=numofplots;
    for i=1:num
        %save figure to output to file later;
        if (isstring(Name))
            name=Name(i);
            titleName=strcat(name,' Plot');
            name =strcat(name,type);
        else
            name=Name;
            titleName= strcat(name,num2str(i),' Plot');
            name =strcat(name,num2str(i),type);
        end
        fig=figure(numfig);     

        %plot the theta and configure plot
        i
        plot(t,data(:,i),'LineWidth',1.25,'Color',colors(mod(i,size(colors,1)),:));
        axis(AX)
        xlabel('Time(s)') ;
        ylabel(Ylabel) ;
        title(titleName);

        %name of output file.
        
        saveas(fig,name);
    end
    
    % plot all the lines together.
    fig=figure(numfig);
    plot(t,data(:,1:numofplots),'LineWidth',1.25)
    xlabel('Time(s)') ;
    ylabel(Ylabel) ;
    
    if isstring(Name)
        name ='theta_j-theta_i';
        title(name);
        name =strcat(name,type)
    else
         
         name =strcat(Name,type);
         title(Name);
    end
    saveas(fig,name);

end

function lap=laplacian(A)
  lap =diag(sum(A,2)) - A;
end
