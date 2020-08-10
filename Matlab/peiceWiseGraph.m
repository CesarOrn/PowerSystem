function [data]=peiceWiseGraph(PeakTime,Peaks,Time,SS)
Slopes(1,1)=Peaks(1,1)/PeakTime(1,1);

for i=2:size(Peaks)
    Slopes(i,1)=(Peaks(i,1)-Peaks(i-1,1))/(PeakTime(i,1)-PeakTime(i-1,1));
end
Slopes
for i=1:size(Time)
if(Time(i,1)<PeakTime(1,1))
    data(i,1)=Slopes(1,1)*Time(i,1);
elseif(Time(i,1)<PeakTime(size(PeakTime,1),1))
    for j=2:size(Peaks,1)
        if(PeakTime(j-1,1)<Time(i,1)&&Time(i,1)<PeakTime(j,1))
           data(i,1)= Slopes(j,1)*(Time(i,1)-PeakTime(j-1,1))+Peaks(j-1,1);
           break;
        end
        
    end
    
else
    data(i,1)=SS;
end
end