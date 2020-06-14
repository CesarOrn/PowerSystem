function [var]=transientResponse(Lambda,Q)
    lambda= Lambda*ones(size(Lambda,1),1);

    omega =sqrt(lambda);
    zeta=1./(2.*omega);
    
    SS= Q./lambda;
    SS(1,1)=0;  % Steady state
    
    pt=pi./(omega.*sqrt(1-zeta.^2));
    pt(1,1)=0; % Peak time
    
    peak=Q./lambda.*(1+exp(-pi.*zeta./sqrt(1-zeta.^2)));
    peak(1,1)=0; % Peak

    var{1,1}=lambda;
    
    var{2,1}=omega;
    var{3,1}=zeta;
    
    var{4,1}=SS;
    
    var{5,1}=pt;
    
    var{6,1}=peak;

end