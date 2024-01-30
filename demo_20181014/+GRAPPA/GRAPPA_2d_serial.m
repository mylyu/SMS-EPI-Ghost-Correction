function [kspace_coils_kykxc] = GRAPPA_2d_serial(undersampled_kspace_kykxc,ACS_kykxc,header, regularization_factor)
%%% Note that serial 2D GRAPPA is suboptimal because synthesized low SNR
%%% data is used in step 2.
import GRAPPA.*
if nargin<3 || isempty(regularization_factor)
    regularization_factor=0.001;
end
disp(['regularization_factor ' num2str(regularization_factor)])

if ~isfield(header,'blocks')|| isempty(header.blocks)
    header.blocks=2;
end
if ~isfield(header,'column')|| isempty(header.column)
    header.column=3;
end

Ry=header.subsampling_factor(1);
Rx=header.subsampling_factor(2);

% sampling pattern
sampling_mask=sum(undersampled_kspace_kykxc,3)~=0;
sampling_mask_ACS=sum(crop(undersampled_kspace_kykxc,size(ACS_kykxc)),3)~=0;
% sampling_mask_ky=squeeze(sampling_mask(:,1));
sampling_mask_kx=squeeze(sampling_mask(1,:));
sampling_mask_ACS_kx=squeeze(sampling_mask_ACS(1,:));
% recon ky
header_1=header;
header_1.subsampling_factor=Ry;
header_1.blocks=header.blocks;
header_1.column=header.column;
ACS_kykxc_1=ACS_kykxc(:,sampling_mask_ACS_kx,:);
undersampled_kspace_kykxc_1=undersampled_kspace_kykxc(:,sampling_mask_kx,:);
[GRAPPA_weights_1] = GRAPPA_calibrate_weights(ACS_kykxc_1, header_1, regularization_factor);
[~, kspace_temp_1, ~] = ...
    GRAPPA_interpolate_imageSpace(undersampled_kspace_kykxc_1, header_1, GRAPPA_weights_1);
undersampled_kspace_kykxc(:,sampling_mask_kx,:)= kspace_temp_1;

% recon kx
header_2=header;
header_2.subsampling_factor=Rx;
header_2.blocks=header.column-1;% swap
header_2.column=header.blocks+1;
ACS_kykxc_2=permute(ACS_kykxc,[2 1 3]);
undersampled_kspace_kykxc_2=permute(undersampled_kspace_kykxc,[2 1 3]);
[GRAPPA_weights_2] = GRAPPA_calibrate_weights(ACS_kykxc_2, header_2, regularization_factor);
[~, kspace_temp_2, ~] = ...
    GRAPPA_interpolate_imageSpace(undersampled_kspace_kykxc_2, header_2, GRAPPA_weights_2);
kspace_coils_kykxc=permute(kspace_temp_2,[2 1 3]);

