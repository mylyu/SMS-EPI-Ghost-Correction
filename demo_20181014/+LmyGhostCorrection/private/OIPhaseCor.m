%evaluated by overall entropy
%20130325 by victor
%20130414 add parameters of Nshot
function CorxKy = OIPhaseCor(xKy, PhasPara, Nshot)



ConPhase = PhasPara(1);
LinPhase = PhasPara(2);

Phase = ConPhase+(1:size(xKy,2))*LinPhase;
Phasemap = repmat([Phase;-Phase],[size(xKy,1)/2/Nshot 1 size(xKy,3)]);
CorxKy = complex(zeros(size(xKy)));

% CorxKy(1:2:end,:,:) = xKy(1:2:end,:,:).*exp(1i.*Phasemap);
% CorxKy(2:2:end,:,:) = xKy(2:2:end,:,:).*exp(1i.*Phasemap);


for ii = 1:Nshot
    CorxKy(ii:Nshot:end,:,:) = xKy(ii:Nshot:end,:,:).*exp(1i.*Phasemap);
end