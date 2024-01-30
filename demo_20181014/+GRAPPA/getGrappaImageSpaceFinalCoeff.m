function [unmixing_coeff] = getGrappaImageSpaceFinalCoeff(GRAPPA_unmixing_map,coil_sensitivity_map_xyc)
%%% If you have coil sensitivity maps, you can use this function to
%%% directly obtained the coil-combined image using GRAPPA. However, I note
%%% this type of reconstruction is less robust than SOS type of GRAPPA,
%%% because coil sensitivity map sometimes cannot capture very rapid phase
%%% change, resulting ringing artifact in poorly shimmed areas.
[Npe, Nfe, Ncoil,~]=size(GRAPPA_unmixing_map);
unmixing_coeff = complex(zeros(Npe, Nfe, Ncoil));
COIL_index=zeros(Ncoil,Ncoil);
for index = 1 : Ncoil
        COIL_index(index, :) = 1: Ncoil;
end
for ii = 1:Ncoil
    for jj = 1:Ncoil
        unmixing_coeff(:,:,COIL_index(jj,ii)) = unmixing_coeff(:,:,COIL_index(jj,ii)) + GRAPPA_unmixing_map(:,:,jj,ii).*conj(coil_sensitivity_map_xyc(:,:,jj));
    end
end
unmixing_coeff=unmixing_coeff./repmat(sum(abs(coil_sensitivity_map_xyc).^2,3),[1 1 size(coil_sensitivity_map_xyc,3)]);
% I_aliased=ifft2c(subsampled_kspace_data(:,:,:));
% I_final=sum(I_aliased.*ucomp, 3);