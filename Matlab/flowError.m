function [Times,errors]=flowError(flow,flowAprox,Time)
flowSizeHalf= size(flow,1)/2
flowSize= size(flow,1)

%real

for i=1:size(flow,2)
    [RMF,RIF]=max(abs(flow(1:flowSizeHalf,i)));
    [RMS,RIS]=max(abs(flow(flowSizeHalf+1:flowSize,i)));
    
    IndexesR(i,:)=[RIF RIS+flowSizeHalf];
end

%aprox

for i=1:size(flowAprox,2)
    [AMF,AIF]=max(abs(flowAprox(1:flowSizeHalf,i)));
    [AMS,AIS]=max(abs(flowAprox(flowSizeHalf+1:flowSize,i)));
    
    IndexesA(i,:)=[AIF AIS+flowSizeHalf];
end


for i=1:size(IndexesR,1)
    RealTimeF=Time(IndexesR(i,1))
    RealFlowF=flow(IndexesR(i,1),i)

    
    RealTimeS=Time(IndexesR(i,2))
    RealFlowS=flow(IndexesR(i,2),i)

    
    
    AproxTimeF=Time(IndexesA(i,1))
    AproxFlowF=flow(IndexesA(i,1),i)

    
   
    AproxTimeS=Time(IndexesA(i,2))
    AproxFlowS=flow(IndexesA(i,2),i)
    fprintf('new\n')
    
    TimeErrorF=((RealTimeF-AproxTimeF)/RealTimeF)*100;
    TimeErrorS=((RealTimeS-AproxTimeS)/RealTimeS)*100;
    
    FlowErrorF=((RealFlowF-AproxFlowF)/RealFlowF)*100;
    FlowErrorS=((RealFlowS-AproxFlowS)/RealFlowS)*100;
    

    errors(i,:)=[TimeErrorF, FlowErrorF,TimeErrorS,FlowErrorS];
end
Times=[IndexesR,IndexesA]


end