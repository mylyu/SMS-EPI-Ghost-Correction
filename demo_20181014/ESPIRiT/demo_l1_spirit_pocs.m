if exist('spirit.m')==0
	addpath('utils')
end

%load phantom.mat
load brain_8ch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Reconstruction Parameters %%%%%%%%%%%%%%%%%%%
%						     %	 			
kSize = [5,5];  % SPIRiT kernel size
nIter = 20; % number of iteration; phantom requires twice as much as the brain.
mask_type = 'random4'; % options are: 'unif','random4','random3'
CalibTyk = 0.01;  % Tykhonov regularization in the calibration
wavWeight = 0.0015;  % Wavelet soft-thresholding regularization in the reconstruction (SPIRiT only)
addNoiseSTD = 0.0; % add noise. Use only for phantom. Brain has enough noise!
skipGRAPPA = 1; % skip GRAPPA recon.







DATA = DATA+randn(size(DATA))*addNoiseSTD + i*randn(size(DATA))*addNoiseSTD;
im = ifft2c(DATA);

switch mask_type
    case 'unif'
           mask = mask_unif;
           disp('Using uniform undersampling.')
           disp('Change mask_type in the code for other sampling patterns.');
           disp(' ');
           disp(' ');
           
    case 'random3'
            mask = mask_randm_x3;
            if skipGRAPPA==0
                disp('Warning:GRAPPA recon is slow for random sampling.')
                disp('change to skipGRAPPA=1 in the code to skipp.')
                disp(' ');
                disp(' ');
                
                gSkip=1;
            end
            
            
    case 'random4'
            mask = mask_randm_x4;
            if skipGRAPPA==0
                disp('Warning:GRAPPA recon is slow for random sampling.')
                disp('change to skipGRAPPA=1 in the code to skipp.')
                disp(' ')
                disp(' ')
                
                gSkip=1;
            end
    otherwise
        mask = mask_unif'
end

[CalibSize, dcomp] = getCalibSize(mask);  % get size of calibration area from mask
pe = size(DATA,2); fe = size(DATA,1); coils = size(DATA,3); % get sizes
DATA = DATA.*repmat(mask,[1,1,coils]); % multiply with sampling matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scale the data such that the zero-filled density compensated      %%%%%%%%%
% k-space norm is 1. This is useful in order to use similar         %%%%%%%%%
% regularization penalty values for different problems.             %%%%%%%%%

DATAcomp = DATA.*repmat(dcomp,[1,1,coils]);
scale_fctr = norm(DATAcomp(:))/sqrt(coils)/20;
DATA = DATA/scale_fctr;
DATAcomp = DATAcomp/scale_fctr;

im_dc = ifft2c(DATAcomp);
im = im/scale_fctr;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GRAPPA                                   %%%%%%%%%%%
if skipGRAPPA
    disp('skipping grappa, replacing with zero filling');
    res_grappa = DATA;
else
    disp('performing traditional GRAPPA reconstruction');
    kCalib = crop(DATA,[CalibSize,coils]);
    res_grappa = GRAPPA(DATA,kCalib,kSize,CalibTyk);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%		 Perform Calibration                       %%%%%%%%%%
disp('performing calibration for SPIRiT')
kCalib = crop(DATA,[CalibSize,coils]);
kernel = zeros([kSize,coils,coils]);

[AtA] = corrMatrix(kCalib,kSize);
for n=1:coils
	kernel(:,:,:,n) = calibrate(AtA,kSize,coils,n,CalibTyk);
end
GOP = SPIRiT(kernel, 'fft',[fe,pe]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                  Reconstruction                        %%%%%%%%%


disp('performing pocs reconstruction')
tic;
[res_pocs] = pocsSPIRiT(DATA,GOP,nIter,DATA,wavWeight,0);
toc

im_pocsspirit = ifft2c(res_pocs);
im_grappa = ifft2c(res_grappa);

im_pocsspirit_err = im_pocsspirit - im;
im_grappa_err = im_grappa - im;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%		Display                                  %%%%%%%%%%%
im_pocsspirit_sqr = sos(im_pocsspirit);
im_grappa_sqr = sos(im_grappa);
im_dc_sqr = sos(im_dc);
im_sqr = sos(im);
im_pocsspirit_err_sqr = sos(im_pocsspirit_err);
im_grappa_err_sqr = sos(im_grappa_err);


figure, imshow(cat(1,cat(2,im_sqr,im_dc_sqr),cat(2,im_pocsspirit_sqr,im_grappa_sqr)),[],'InitialMagnification',150);
title ('top-left:Full		top-right:zero-fill w/dc	bottom-left:SPIRiT 	 	bottom-right:GRAPPA');
figure, imshow(cat(2,im_pocsspirit_err_sqr,im_grappa_err_sqr),[],'InitialMagnification',150);
title ('Difference images: SPIR-iT               GRAPPA');



