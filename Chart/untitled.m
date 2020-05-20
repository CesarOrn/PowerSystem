fid = fopen('/home/cesar/Code/C++/PowerSystem/Chart/Setup1/k140Line2/Linear/LinePower.txt','rt');
LinearData = textscan(fid,'%f%f%f%f%f%f%f%f%f','HeaderLines',2,'CollectOutput',1);
LinearData = cell2mat(LinearData);
fclose(fid);


fid = fopen('/home/cesar/Code/C++/PowerSystem/Chart/Setup1/k140Line2/nonLinear/LinePower.txt','rt');
nonLinear = textscan(fid,'%f%f%f%f%f%f%f%f%f','HeaderLines',2,'CollectOutput',1);
nonLinear = cell2mat(nonLinear);
fclose(fid);

line=8;

figure('Name','Non-Linear');
plot(nonLinear(:,1),nonLinear(:,2),nonLinear(:,1),nonLinear(:,3),nonLinear(:,1),nonLinear(:,4),nonLinear(:,1),nonLinear(:,5),nonLinear(:,1),nonLinear(:,6),nonLinear(:,1),nonLinear(:,7),nonLinear(:,1),nonLinear(:,8),nonLinear(:,1),nonLinear(:,9),'LineWidth',1.1);
xlabel('Time(s)') 
ylabel('Power') 

figure('Name','Linear');
plot(LinearData(:,1),LinearData(:,2),LinearData(:,1),LinearData(:,3),LinearData(:,1),LinearData(:,4),LinearData(:,1),LinearData(:,5),LinearData(:,1),LinearData(:,6),LinearData(:,1),LinearData(:,7),LinearData(:,1),LinearData(:,8),LinearData(:,1),LinearData(:,9),'LineWidth',1.1);

figure('Name','Both');
plot(LinearData(:,1),LinearData(:,2),nonLinear(:,1),nonLinear(:,2));
xlabel('Time(s)') 
ylabel('Power') 