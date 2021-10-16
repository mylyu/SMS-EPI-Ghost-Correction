clear;
close all;
clc;
import LmyParallelImaging.*
%% Data setting
import LmyUtility.*
Data_set_option='invivo';
MB_factor=4;
%% CSM setting
% CSM_option.source='segEPI LPC';
CSM_option.source='segEPI SAKE';
%% PEM setting
% PEM_option.model='1dLpc';
PEM_option.model='2dFreeBART';
%% Fixed settings
LOAD_ROOT='./data';
nFE=128;
nPE=128;
nCoil=32;
nSeg=2;% even or odd echo
R_factor=2*nSeg;
load_path=fullfile(LOAD_ROOT);
%% load data and info according to MB_factor
[MB_kxkyzsc,correspoding_single_slices,CAIPI_pattern]...
    =MBGC_load_dataAndParmeters_20161102(MB_factor,load_path);
MB_kxkyzsc=fft2c(LmyUtility.windows_filter_2d(ifft2c(MB_kxkyzsc),[nFE nPE],@(x)tukeywin(x,0.2)));
nSlice=size(MB_kxkyzsc,3)*MB_factor;
cmrr_kxkyzc_pos=squeeze(MB_kxkyzsc(:,:,:,1,:));
cmrr_kxkyzc_neg=squeeze(MB_kxkyzsc(:,:,:,2,:));
% aliased_images_xyzsc=ifft2c_MN(cmrr_kxkyzsc,nFE*2,nPE);
%% load CSM and trim it
disp('loading CSM')
CSM_xyzc=MBGC_load_CSM_xyzc(CSM_option,load_path);
disp('done')

load(fullfile(load_path,['CSM_130y_xyzc_epi_sakeCor_bart']),'CSMweight_xyz_epi_sakeCor')%
CSMweight_used=crop(CSMweight_xyz_epi_sakeCor,nFE,nPE,size(CSMweight_xyz_epi_sakeCor,3));
CSM_mask_threshold=0.2;
mask=CSMweight_used>CSM_mask_threshold;

% caution
% mask=mask*0+1;% no mask with this line
%
CSM_xyzc(isnan(CSM_xyzc))=0;
CSM_xyzc(~repmat(mask,[1 1 1 nCoil]))=0;
%% Get phase error maps
disp('loading PEM')
phase_error_to_use=MBGC_load_phase_error_maps(PEM_option,load_path, [nFE nPE]);
% phase_error_to_use=conj(phase_error_to_use);
disp('done')
%% reshape beforehand for each SMS package for easier parallel computation
CSM_forEachSmsPack_xyzcp=zeros(nFE,nPE,MB_factor,nCoil,nSlice/MB_factor);
PEM_forEachSmsPack_xyzp=zeros(nFE,nPE,MB_factor,nSlice/MB_factor);
for iSmsPack=1:nSlice/MB_factor
    this_correspoding_single_slices=correspoding_single_slices{iSmsPack};
    CSM_forEachSmsPack_xyzcp(:,:,:,:,iSmsPack)=CSM_xyzc(:,:,this_correspoding_single_slices,:);
    PEM_forEachSmsPack_xyzp(:,:,:,iSmsPack)=phase_error_to_use(:,:,this_correspoding_single_slices);
end
%% start PEC-SENSE
Shift_distance_fromNegToPos=-nSeg;%
recon_slices_PEC=zeros(nFE,nPE,nSlice);
recon_PEC_xyzp=zeros(nFE,nPE,MB_factor,nSlice/MB_factor);
for iSmsPack=1:nSlice/MB_factor
    disp(iSmsPack)
    % adjust positive echoes
    MB_kxkyc_pos=squeeze(cmrr_kxkyzc_pos(:,:,iSmsPack,:));
    pseudoPi_kxkyc_pos=MBGC_readoutCascadeSmsData(MB_kxkyc_pos, MB_factor);
    CSM_xyzc_thisPack_pos=CSM_forEachSmsPack_xyzcp(:,:,:,:,iSmsPack);
    adjustedCSM_xyc_pos=MBGC_readoutCascadeAdjustCSM(CSM_xyzc_thisPack_pos, CAIPI_pattern);
    % adjust negative echoes
    MB_kxkyc_neg=circshift(squeeze(cmrr_kxkyzc_neg(:,:,iSmsPack,:)),[0 Shift_distance_fromNegToPos]);% shift to match sampling location
    pseudoPi_kxkyc_neg=MBGC_readoutCascadeSmsData(MB_kxkyc_neg, MB_factor);
    kspaceShifted_CSM=LmyUtility.addPhaseRamp(CSM_xyzc_thisPack_pos,[0 -nSeg]);% shift CSM as well
    for iSlice=1:MB_factor
        kspaceShifted_CSM(:,:,iSlice,:)=kspaceShifted_CSM(:,:,iSlice,:).*exp(-1i*pi*(iSlice-1));
    end
    %   kspaceShifted_CSM=CSM_xyzc_thisPack_pos; % for debug
    CSM_xyzc_thisPack_neg=bsxfun(@times,kspaceShifted_CSM,PEM_forEachSmsPack_xyzp(:,:,:,iSmsPack));
    %     CSM_xyzc_thisPack_neg= kspaceShifted_CSM; % for debug
    adjustedCSM_xyc_neg=MBGC_readoutCascadeAdjustCSM(CSM_xyzc_thisPack_neg, CAIPI_pattern);
    % combine two parts of data and do as conventional PI recon
    pseudoPi_kxkyc_full=cat(3, pseudoPi_kxkyc_pos, pseudoPi_kxkyc_neg);
    adjustedCSM_xyc_full=cat(3, adjustedCSM_xyc_pos, adjustedCSM_xyc_neg);
    reg_factor=0.0;
%     reg_factor=0.0001;
        recon_xy_temp=LmyParallelImaging.ismrm_cartesian_SENSE_general2(pseudoPi_kxkyc_full,adjustedCSM_xyc_full,[MB_factor R_factor],[],[],reg_factor);
%             recon_xy_temp=bart('pics -i 100 -R Q:0.01 -R L:7:7:0.0015',reshape(pseudoPi_kxkyc_full,[nFE*MB_factor,nPE,1,nCoil*2]),reshape(adjustedCSM_xyc_full,[nFE*MB_factor,nPE,1,nCoil*2]));
%             recon_xy_temp=bart('pics -i 100 -l1 -r 0.01',reshape(pseudoPi_kxkyc_full,[nFE*MB_factor,nPE,1,nCoil*2]),reshape(adjustedCSM_xyc_full,[nFE*MB_factor,nPE,1,nCoil*2]));
    recon_PEC_xyzp(:,:,:,iSmsPack)=MBGC_readoutCascadeRestoreReconSlices(recon_xy_temp, CAIPI_pattern);
end
%% Reorder slices
for iSmsPack=1:nSlice/MB_factor
    this_correspoding_single_slices=correspoding_single_slices{iSmsPack};
    recon_slices_PEC(:,:,this_correspoding_single_slices)=recon_PEC_xyzp(:,:,:,iSmsPack);
end

%% Plot CSM absolute values
figure(411);
MY_montage(abs(CSM_xyzc(:,:,:,[21])),'size',[4 size(recon_slices_PEC,3)/4],'displayrange',[]),title('CSM of channel 21')
%% for m/s
slice2show=[1:20];
im2show=recon_slices_PEC(:,:,slice2show);
disp_func=@(x)rot90(x,1);
disp_range_max=prctile(abs(im2show(:)),99);
figure;MY_montage(disp_func(abs(im2show)),'displayrange',[0 disp_range_max],'size',[3 8]);
figure;MY_montage(disp_func(abs(im2show)*5),'displayrange',[0 disp_range_max],'size',[3 8]);
im2show=phase_error_to_use(:,:,slice2show);
figure;MY_montage(disp_func(angle(im2show)),'displayrange',[-0.3 0.3],'size',[3 8]);colormap(gca,jet);

