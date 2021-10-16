%% defined parameters
clear
import LmyFilterPackage.*
import LmyParallelImaging.*
import LmyUtility.*
load_folder='./data';
nRep=2;
R=2*nRep;
nPE_extended=130; %k-space size need to pad to to match MB data resolution
%% load EPI data
load(fullfile(load_folder,'SbData_128kx'),'epi_kxkyzc_2shot','epi_2shot_sampling_mask_kys','epi_2shot_fovY','epi_2shot_rangeKy')
% epi_kxkyzc_2shot=double(epi_kxkyzc_2shot(:,2:end-1,:,:));
epi_kxkyzc_2shot=double(epi_kxkyzc_2shot(:,1:end-2,:,:));
[nFE,nPE,nSlice,nCoil]=size(epi_kxkyzc_2shot);
assert(nFE==128&&nPE==128)
epi_kxkyzcr_pos=zeros(nFE,nPE,nSlice,nCoil,nRep);
epi_kxkyzcr_neg=zeros(nFE,nPE,nSlice,nCoil,nRep);
for iRep=1:nRep
    epi_kxkyzcr_pos(:,iRep:2*nRep:end,:,:,iRep)=epi_kxkyzc_2shot(:,iRep:2*nRep:end,:,:);
    epi_kxkyzcr_neg(:,iRep+nRep:2*nRep:end,:,:,iRep)=epi_kxkyzc_2shot(:,iRep+nRep:2*nRep:end,:,:);
end
%% load CSM
load(fullfile(load_folder,['CSM_130y_xyzc_epi_sakeCor_bart']),'CSM_xyzc_epi_sakeCor','CSMweight_xyz_epi_sakeCor')%
CSM_used=crop(CSM_xyzc_epi_sakeCor,128,128,size(CSM_xyzc_epi_sakeCor,3),size(CSM_xyzc_epi_sakeCor,4));
CSMweight_used=CSMweight_xyz_epi_sakeCor;
%% trim CSM
CSM_mask_threshold=0;
CSM_mask=CSMweight_used<CSM_mask_threshold;
CSM_used(repmat(CSM_mask,[1 1 1 nCoil]))=0;
%% SENSE recon on even/odd only echoes
recon_xyzr_pos=zeros(nFE,nPE,nSlice,nRep);
recon_xyzr_neg=recon_xyzr_pos;
for iRep=1:nRep
    for iSlice=1:nSlice
        disp(iSlice)
        reg=0.001;
        CsmThisSlice=squeeze(CSM_used(:,:,iSlice,:));
        KspTemp=squeeze(epi_kxkyzcr_pos(:,:,iSlice,:,iRep));
        [recon_xyzr_pos(:,:,iSlice,iRep)] = ismrm_cartesian_SENSE_1D(KspTemp,CsmThisSlice,R,[],[],reg);
        KspTemp=squeeze(epi_kxkyzcr_neg(:,:,iSlice,:,iRep));
        [recon_xyzr_neg(:,:,iSlice,iRep)] = ismrm_cartesian_SENSE_1D(KspTemp,CsmThisSlice,R,[],[],reg);
        figure(12);
        MY_montage(abs(cat(3,recon_xyzr_pos(:,:,iSlice,iRep),recon_xyzr_neg(:,:,iSlice,iRep))));title(['slice' num2str(iSlice)])
    end
end
%% Resize recon images to match MB images by padding zeros in kspace
recon_xyzr_pos=LmyUtility.ifft2c_MN(fft2c(recon_xyzr_pos),nFE,nPE_extended);
recon_xyzr_neg=LmyUtility.ifft2c_MN(fft2c(recon_xyzr_neg),nFE,nPE_extended);
figure(15);MY_montage(abs(recon_xyzr_pos(:,:,1:end/4,:)),'displayrange',[0 50]);title('pos, 2 shots')
figure(16);MY_montage(abs(recon_xyzr_neg(:,:,1:end/4,:)),'displayrange',[0 50]);title('neg, 2 shots')
figure(115);MY_montage(angle(recon_xyzr_pos(:,:,1:end/4,:)),'displayrange',[]);title('pos')
figure(116);MY_montage(angle(recon_xyzr_neg(:,:,1:end/4,:)),'displayrange',[]);title('neg')
%% do phase alignment if with multiple input
recon_xyz_pos_RepSum=squeeze(recon_xyzr_pos(:,:,:,1));
recon_xyz_neg_RepSum=squeeze(recon_xyzr_neg(:,:,:,1));
for iRep=2:nRep
    recon_xyz_pos_RepSum=pos_neg_add(recon_xyz_pos_RepSum,recon_xyzr_pos(:,:,:,iRep))/2;
    recon_xyz_neg_RepSum=pos_neg_add(recon_xyz_neg_RepSum,recon_xyzr_neg(:,:,:,iRep))/2;
end
figure(25);MY_montage(abs(recon_xyz_pos_RepSum(:,:,1:end/4)),'displayrange',[]);title('pos sum')
figure(26);MY_montage(abs(recon_xyz_neg_RepSum(:,:,1:end/4)),'displayrange',[]);title('neg sum')
%% save result if requested

save_filename='SbData_PosNeg_Recon';
save(fullfile(load_folder,save_filename),'recon_xyz_pos_RepSum','recon_xyz_neg_RepSum','recon_xyzr_pos','recon_xyzr_neg','-v7.3');