clear
import LmyGhostCorrection.*
import LmyUtility.*
load_folder=fullfile('.','data');
load(fullfile(load_folder,'SbData_128kx128ky_lpc'),'epi_kxkyzc_2shot_lpcCor')

nSlice=size(epi_kxkyzc_2shot_lpcCor,3);
nCoil=size(epi_kxkyzc_2shot_lpcCor,4);
%%
ncalib = 48;
% threshold_list=[4*ones(1,10),4.5*ones(1,15),5*ones(1,30),4.5*ones(1,15),4*ones(1,10)]; % good for phantom
% threshold_list=[linspace(4,5,40),linspace(5,4,40)];
threshold_list=[linspace(3,4,5),linspace(4,4,15)];
ksize = [3,3]; % ESPIRiT kernel-window-size
sakeIter = 50;
% wnthresh = 4.5; % 3 or 4 good for brain

epi_kxkyzc_2shot_sakeCor=zeros(ncalib,ncalib,nSlice,nCoil);
epi_kxkyzc_2shot_sakeCor_fullCoils=...
    cat(4,epi_kxkyzc_2shot_sakeCor,epi_kxkyzc_2shot_sakeCor,epi_kxkyzc_2shot_sakeCor,epi_kxkyzc_2shot_sakeCor);
for iSlice=1:20
    % for iSlice=5
    disp(iSlice)
    DATA=squeeze(epi_kxkyzc_2shot_lpcCor(:,:,iSlice,:));
    nSeg=4;
    % convert shot to VCC
    DATA_org=DATA;
    for iSeg=2:nSeg
        DATA=cat(3,DATA,circshift(DATA_org,-iSeg+1,2));
    end
    if exist('threshold_list','var')
        wnthresh = threshold_list(iSlice); % Window-normalized number of singular values to threshold
    end
    [sx,sy,Nc] = size(DATA);
    mask = zeros(size(DATA,1),size(DATA,2));
    mask(:,1:nSeg:end) = 1;
    DATA2recon=DATA.* repmat(mask,[1,1,size(DATA,3)]);
    DATAc = DATA;
    calibc = crop(DATAc,[ncalib,ncalib,size(DATA,3)]);
    
    %% Perform SAKE reconstruction to recover the calibration area
    im = ifft2c(DATAc);
    disp('Performing SAKE recovery of calibration');
    tic; calib_sake = SAKEwithInitialValue(calibc, [ksize], wnthresh,sakeIter, 0,repmat(crop(mask,[ncalib,ncalib]),[1 1 size(DATA,3)]));toc
    calib_sake(:,:,end/4+1:end)=circshift(calib_sake(:,:,end/4+1:end),[0 1]);
    calib_sake(:,:,end/2+1:end)=circshift(calib_sake(:,:,end/2+1:end),[0 1]);
    calib_sake(:,:,end*3/4+1:end)=circshift(calib_sake(:,:,end*3/4+1:end),[0 1]);
    epi_kxkyzc_2shot_sakeCor_fullCoils(:,:,iSlice,:)= calib_sake;
    a=pos_neg_add(ifft2c(calib_sake(:,:,1:end/4)),ifft2c(calib_sake(:,:,end/4+1:end/2)));
    b=pos_neg_add(ifft2c(calib_sake(:,:,end/2+1:end*3/4)),ifft2c(calib_sake(:,:,end*3/4+1:end)));
    epi_kxkyzc_2shot_sakeCor(:,:,iSlice,:)=fft2c(pos_neg_add(a,b))/4;
end
disp('Done')
figure();LmyUtility.MY_montage((sos(ifft2c(crop(epi_kxkyzc_2shot_lpcCor,[48,48,20,32])))),'displayrange',[0 200]);title('LPC')
figure();LmyUtility.MY_montage((sos(ifft2c(epi_kxkyzc_2shot_sakeCor))),'displayrange',[0 200]);title('VC SAKE')
figure();LmyUtility.MY_montage((abs(ifft2c(crop(epi_kxkyzc_2shot_lpcCor(:,:,:,21),[48,48,20])))),'displayrange',[0 20]);title('LPC')
figure();LmyUtility.MY_montage((abs(ifft2c(epi_kxkyzc_2shot_sakeCor(:,:,:,21)))),'displayrange',[0 20]);title('VC SAKE')
figure();LmyUtility.MY_montage((abs(ifft2c(crop(epi_kxkyzc_2shot_lpcCor(:,:,:,11),[48,48,20])))),'displayrange',[0 20]);title('LPC')
figure();LmyUtility.MY_montage((abs(ifft2c(epi_kxkyzc_2shot_sakeCor(:,:,:,11)))),'displayrange',[0 20]);title('VC SAKE')
figure();LmyUtility.MY_montage((sos(ifft2c(crop(epi_kxkyzc_2shot_lpcCor,[48,48,20,32])))),'displayrange',[0 25]);title('LPC')
figure();LmyUtility.MY_montage((sos(ifft2c(epi_kxkyzc_2shot_sakeCor))),'displayrange',[0 25]);title('VC SAKE')
figure();LmyUtility.MY_montage(db(sos(ifft2c(crop(epi_kxkyzc_2shot_lpcCor,[48,48,20,32])))));title('LPC')
figure();LmyUtility.MY_montage(db(sos(ifft2c(epi_kxkyzc_2shot_sakeCor))));title('VC SAKE')
%%
save_name = 'SbData_48kx48ky_sake.mat';
warning('results not saved')
% save(fullfile(load_folder,save_name),'epi_kxkyzc_2shot_sakeCor')
