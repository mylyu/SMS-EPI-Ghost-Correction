function CorxKy = IPhaseCor(xKy, PhasPara)


ConPhase = PhasPara(1);
LinPhase = PhasPara(2);

Phase = ConPhase+(1:size(xKy,2))*LinPhase;
Phasemap = repmat([Phase;-Phase],[size(xKy,1)/4 1]);
CorxKy = zeros(size(xKy));
CorxKy(1:2:end,:) = xKy(1:2:end,:).*exp(1i.*Phasemap);
CorxKy(2:2:end,:) = xKy(2:2:end,:).*exp(1i.*Phasemap);

end