function [image_coilcombined_usingCSM, unmixing_map_coilCombined, unmixing_map_coilWise] = ...
    GRAPPA_interpolate_imageSpaceCoilCombined(undersampled_kspace_kykxc, header, GRAPPA_weights, coil_sensitivity_map, unmixing_map_coilCombined)
import GRAPPA.*
[Npe,Nfe,Ncoil]=size(undersampled_kspace_kykxc);
header.Npe=Npe;
header.Nfe=Nfe;
if nargin<5 || isempty(unmixing_map_coilCombined)
    disp('recalculating coil combined GRAPPA unmixing map')
    unmixing_map_coilWise = GRAPPA.getGrappaImageSpaceCoilCoeff(header,GRAPPA_weights);
    [unmixing_map_coilCombined] = GRAPPA.getGrappaImageSpaceFinalCoeff(unmixing_map_coilWise,coil_sensitivity_map);
end
I_aliased=ifft2c(undersampled_kspace_kykxc(:,:,:));
image_coilcombined_usingCSM=sum(I_aliased.*unmixing_map_coilCombined, 3);

% I_coils = zeros(Npe, Nfe, header.consider_coils);
% for ii = 1:header.consider_coils
%     I_coils(:,:,ii)=sum(I_aliased(:,:,COIL_index(ii,:)).*squeeze(w(:,:,ii,:)),3);
% end
% I_final=sos(I_coils);

return