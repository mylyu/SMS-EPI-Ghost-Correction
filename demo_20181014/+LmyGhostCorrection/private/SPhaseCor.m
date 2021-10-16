% correct the raw data in x-Ky space by input parameter (constant and linear)
function CorxKy = SPhaseCor(xKy, PhasPara)
ConPhase = PhasPara(1);
LinPhase = PhasPara(2);
Phase = ConPhase+(1:size(xKy,2))*LinPhase;
Phasemap = repmat([Phase;-Phase],[size(xKy,1)/2 1 size(xKy,3)]);
CorxKy = xKy.*exp(1i.*Phasemap);
end