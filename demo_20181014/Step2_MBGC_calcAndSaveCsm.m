clear
load_path='./data';
% 'SbData_48kx48ky_sake_ksize3_Thres4ExceptS1T3_Iter50' is vc-sake
% corrected single-band data with optimized thresholds.
load(fullfile(load_path,'SbData_48kx48ky_sake_ksize3_Thres4ExceptS1T3_Iter50'),'epi_kxkyzc_2shot_sakeCor')
% load(fullfile(load_path,'SbData_48kx48ky_sake'),'epi_kxkyzc_2shot_sakeCor')%%
import LmyUtility.*
iNpe = 130; % has to use this value for this particular dataset to ensure FOV match
epi_xyzc_2shot_sakeCor=LmyUtility.ifft2c_MN(epi_kxkyzc_2shot_sakeCor, 128, iNpe);
[nFE,nPE,nSlice,nCoil]=size(epi_xyzc_2shot_sakeCor);
CSM_xyzc_epi_sakeCor=zeros(nFE,nPE,nSlice,nCoil);
CSMweight_xyz_epi_sakeCor=zeros(nFE,nPE,nSlice);
calib_size=[48,48];
%% using BART
parfor iSlice=1:nSlice
    disp(iSlice)
    [CSM_xyzc_epi_sakeCor(:,:,iSlice,:), CSMweight_xyz_epi_sakeCor(:,:,iSlice)]    ...
        = bart(['ecalib -c0 -m 1 -r ' num2str(calib_size(1))],fft2c(epi_xyzc_2shot_sakeCor(:,:,iSlice,:)));
end
CSM_xyzc_epi_sakeCor=single(CSM_xyzc_epi_sakeCor);
CSMweight_xyz_epi_sakeCor=single(CSMweight_xyz_epi_sakeCor);
% warning('results not saved')
save(fullfile(load_path,['CSM_' num2str(nPE) 'y_xyzc_epi_sakeCor_bart']),'CSM_xyzc_epi_sakeCor','CSMweight_xyz_epi_sakeCor','calib_size','-v7.3')%
disp('done')