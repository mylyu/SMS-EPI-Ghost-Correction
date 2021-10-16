function [epi_kxkyzc_sakeCor,epi_kxkyzcs_sakeCor_nSeg]=twoDimNonLinCorr_vcSake(epi_kxkyzc_lpcCor, nSeg, nCalib, ksize, wnthresh_list, sakeIter, plotIter)
%%% This function performs virtual coil SAKE reconstruction on EPI data to
%%% remove Nyquist ghosts and inter-shot phase variation. Ideally, the data
%%% should be firstly corrected by 1d lpc; Note: pos_neg_add from Hoge's
%%% toolbox is needed.

% ncalib = 48;
% ksize = [3,3]; % ESPIRiT kernel-window-size
% sakeIter = 100;
import LmyGhostCorrection.*
if nargin<7
    plotIter=[]; % no plot
end

org_size=size(epi_kxkyzc_lpcCor);

if ndims(epi_kxkyzc_lpcCor)==3 % assume single-slice
epi_kxkyzc_lpcCor=reshape(epi_kxkyzc_lpcCor,size(epi_kxkyzc_lpcCor,1),size(epi_kxkyzc_lpcCor,2),1,size(epi_kxkyzc_lpcCor,3));
end
[nFE,nPE,nSlice,nCoil]=size(epi_kxkyzc_lpcCor);

if numel(wnthresh_list)==1
    wnthresh_list=repmat(wnthresh_list(1),nSlice);
end

epi_kxkyzc_sakeCor=zeros([nCalib,nCalib,nSlice,nCoil]);
epi_kxkyzcs_sakeCor_nSeg=zeros([nCalib,nCalib,nSlice,nCoil,nSeg]);
disp('Performing SAKE recovery of calibration');
for iSlice=1:nSlice
    disp(['VC-SAKE on slice ' num2str(iSlice)])
    DATA=squeeze(epi_kxkyzc_lpcCor(:,:,iSlice,:));
    % convert shot to VCC
    DATA_org=DATA;
    for iSeg=2:nSeg
        DATA=cat(3,DATA,circshift(DATA_org,-iSeg+1,2));
    end
    mask = zeros(size(DATA,1),size(DATA,2));
    mask(:,1:nSeg:end) = 1;
    DATAc = DATA;
    calibc = crop(DATAc,[nCalib,nCalib,size(DATA,3)]);
    
    %% Perform SAKE reconstruction to recover the calibration area
    tic; 
    calib_sake = LmyGhostCorrection.SAKEwithInitialValue(calibc, ksize, wnthresh_list(iSlice),sakeIter,plotIter,...
        repmat(crop(mask,[nCalib,nCalib]),[1 1 size(DATA,3)]));
    toc;
    
    %% combine shots
    calib_sake=reshape(calib_sake,[nCalib,nCalib,nCoil,nSeg]);
    for iSeg=2:nSeg
        calib_sake(:,:,:,iSeg)=circshift(calib_sake(:,:,:,iSeg),[0 iSeg-1]);
    end
    epi_kxkyzcs_sakeCor_nSeg(:,:,iSlice,:,:)=calib_sake;
    temp_img=ifft2c(calib_sake);
    % every time combine the 1st half with 2nd half.
    %     for iRound=1:log2(nSeg)
    %         for iSeg=1:nSeg/(2.^iRound)
    %             temp_img(:,:,:,iSeg)=pos_neg_add(double(temp_img(:,:,:,iSeg)),double(temp_img(:,:,:,iSeg+nSeg/(2.^iRound))))/2;
    %         end
    %     end
    for iSeg=2:nSeg
        temp_img(:,:,:,1)=pos_neg_add(double(temp_img(:,:,:,1))*(nSeg-1),double(temp_img(:,:,:,iSeg)))/nSeg;
    end
    
    epi_kxkyzc_sakeCor(:,:,iSlice,:)=fft2c(temp_img(:,:,:,1));
end

epi_kxkyzc_sakeCor=reshape(epi_kxkyzc_sakeCor,[nCalib,nCalib,org_size(3:end)]);