function []=exportTable(eigBF,TransBF,LineNameBF,LineCoffBF,ErrorBF,OmegaBF,eigAF,TransAF,LineNameAF,LineCoffAF,ErrorAF,OmegaAF,LineFailed,name)

for i =1:size(eigBF)
   lineEta(i,:)= strcat("eta:",num2str(i));
   space(i,:)= "";
end

%eigBF
char(65)
char(90)
%
place=0;
title=[char(66),'1']
index=[char(65+place),'2']


writematrix("EigenVector",name,'Sheet',1,'Range',title);
writematrix(lineEta,name,'Sheet',1,'Range',index);



writematrix(eigBF,name,'Sheet',1,'Range','B2');

place =size(eigBF,1)
index=[char(65+place),'2']

writematrix(space,name,'Sheet',1,'Range','K2');

writematrix(lineEta,name,'Sheet',1,'Range','L2');
writematrix(TransBF,name,'Sheet',1,'Range','M2');
writematrix("SS",name,'Sheet',1,'Range','M1');
writematrix("PeakTime",name,'Sheet',1,'Range','N1');
writematrix("Peak",name,'Sheet',1,'Range','O1');
writematrix(OmegaBF,name,'Sheet',1,'Range','P2');
writematrix("Freq",name,'Sheet',1,'Range','p1');


writematrix(space,name,'Sheet',1,'Range','Q2');
writematrix(string(LineNameBF),name,'Sheet',1,'Range','R2');
writematrix(LineCoffBF,name,'Sheet',1,'Range','S2');



writematrix("EigenVector",name,'Sheet',1,'Range','AD1');
bla=convertCharsToStrings(LineFailed);
writematrix(bla,name,'Sheet',1,'Range','AE1');

writematrix(lineEta,name,'Sheet',1,'Range','AD2');
writematrix(eigAF,name,'Sheet',1,'Range','AE2');


writematrix(space,name,'Sheet',1,'Range','AO2');


writematrix(lineEta,name,'Sheet',1,'Range','AP2');
writematrix(TransAF,name,'Sheet',1,'Range','AQ2');
writematrix("SS",name,'Sheet',1,'Range','AQ1');
writematrix("PeakTime",name,'Sheet',1,'Range','AR1');
writematrix("Peak",name,'Sheet',1,'Range','AS1');
writematrix(OmegaAF,name,'Sheet',1,'Range','AT2');
writematrix("Freq",name,'Sheet',1,'Range','AT1');

writematrix(space,name,'Sheet',1,'Range','AU2');
writematrix(string(LineNameAF),name,'Sheet',1,'Range','AV2');
writematrix(LineCoffAF,name,'Sheet',1,'Range','AW2');

writematrix(space,name,'Sheet',1,'Range','BG2');

writematrix("Before Failure Error",name,'Sheet',1,'Range','BH1');
writematrix(ErrorBF,name,'Sheet',1,'Range','BH2');
writematrix("After Failure Error",name,'Sheet',1,'Range','BH3');
writematrix(ErrorAF,name,'Sheet',1,'Range','BH4');
end