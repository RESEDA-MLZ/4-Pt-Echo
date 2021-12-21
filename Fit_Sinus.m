function [FittedCurve Contrast Phase Average] = Fit_Sinus(XData, YData)

Emnid = zeros(3,1);
[Emnid(1) Emnid2Ind] = max(YData);

Emnid(2) = XData(Emnid2Ind)/max(XData)*2*pi;
    
Emnid(3) = sum(YData)/max(XData);


model = @(Emnid)sum((Emnid(1)*sin(XData/max(XData)*2*pi-Emnid(2))+Emnid(3)-YData).^2);
    
[p,fval] = fminsearch(model,Emnid);
steps = (max(XData)-min(XData)) ./ (length(XData) - 1 );
    
FineXData = min(XData):steps/10:max(XData);
FittedCurve = [FineXData ; p(1)*sin(FineXData/max(XData)*2*pi - p(2))+p(3)];
    
Contrast = p(1)/p(3);
Phase = p(2);
Average = p(3);


        
    
