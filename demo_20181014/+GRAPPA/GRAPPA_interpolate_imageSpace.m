function [image_coilcombined_sos, kspace_coils_kykxc, unmixing_map_coilWise] = ...
    GRAPPA_interpolate_imageSpace(undersampled_kspace_kykxc, header, GRAPPA_weights, unmixing_map_coilWise)
%%% header must contain the following fields
%%% header.blocks;
%%% header.column;
%%% header.Npe;
%%% header.Nfe;
import GRAPPA.*
subsampled_kspace_data_org=undersampled_kspace_kykxc;
[Npe,Nfe,Ncoil]=size(undersampled_kspace_kykxc);
R=header.subsampling_factor;
header.Npe=Npe;
header.Nfe=Nfe;
if nargin<4 || isempty(unmixing_map_coilWise)
    disp('recalculate GRAPPA unmixing map')
    unmixing_map_coilWise = GRAPPA.getGrappaImageSpaceCoilCoeff(header,GRAPPA_weights);
end

% remove ACS lines if any before deconvolution, we will fill
% them back later.
firstAcquireLine=find(sum(abs(undersampled_kspace_kykxc(:,:)),2),1);
temp=undersampled_kspace_kykxc(firstAcquireLine:R:end,:,:);
undersampled_kspace_kykxc(:)=0;
undersampled_kspace_kykxc(firstAcquireLine:R:end,:,:)=temp;
% deconvolution
I_aliased=ifft2c(undersampled_kspace_kykxc(:,:,:));
I_coils = zeros(Npe, Nfe, Ncoil);
for ii = 1:Ncoil
    I_coils(:,:,ii)=sum(I_aliased(:,:,:).*squeeze(unmixing_map_coilWise(:,:,ii,:)),3);
end
image_coilcombined_sos=sos(I_coils);
kspace_coils_kykxc=fft2c(I_coils);
% refilling ACS lines
acquired_positions=subsampled_kspace_data_org~=0;
kspace_coils_kykxc(acquired_positions)=subsampled_kspace_data_org(acquired_positions);
