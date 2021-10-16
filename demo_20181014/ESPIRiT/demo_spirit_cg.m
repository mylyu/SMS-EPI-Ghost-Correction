if exist('spirit.m')==0
	addpath('utils')
end

%load phantom.mat
load brain_8ch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Reconstruction Parameters %%%%%%%%%%%%%%%%%%%
%						     %	 			
kSize = [5,5];  % SPIRiT kernel size
nIterCG = 30; % number of iteration; phantom requires twice as much as the brain.
mask_type = 'unif'; % options are: 'unif','random4','random3'
CalibTyk = 0.01;  % Tykhonov regularization in the calibration
ReconTyk = 1e-5;  % Tykhovon regularization in the reconstruction (SPIRiT only)
addNoiseSTD = 0.0; % add noise. Use only for phantom. Brain has enough noise!
skipGRAPPA = 0; % skip GRAPPA recon.







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


%kernel = zeros([kSize,coils,coils]);

%[AtA] = corrMatrix(kCalib,kSize);
%for n=1:coils
%	kernel(:,:,:,n) = calibrate(AtA,kSize,coils,n,CalibTyk);
%end

kernel = calibSPIRiT(kCalib, kSize, coils, CalibTyk);
GOP = SPIRiT(kernel, 'fft',[fe,pe]);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                  Reconstruction                        %%%%%%%%%


disp('performing CG reconstruction')
tic;
[res_cg, RESVEC] = cgSPIRiT(DATA,GOP,nIterCG,ReconTyk, DATA);
toc

im_cgspirit = ifft2c(res_cg);
im_grappa = ifft2c(res_grappa);

im_cgspirit_err = im_cgspirit - im;
im_grappa_err = im_grappa - im;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%		Display                                  %%%%%%%%%%%
im_cgspirit_sqr = sos(im_cgspirit);
im_grappa_sqr = sos(im_grappa);
im_dc_sqr = sos(im_dc);
im_sqr = sos(im);
im_cgspirit_err_sqr = sos(im_cgspirit_err);
im_grappa_err_sqr = sos(im_grappa_err);


figure, imshow(cat(1,cat(2,im_sqr,im_dc_sqr),cat(2,im_cgspirit_sqr,im_grappa_sqr)),[],'InitialMagnification',150);
title ('top-left:Full		top-right:zero-fill w/dc	bottom-left:SPIRiT 	 	bottom-right:GRAPPA');
figure, imshow(cat(2,im_cgspirit_err_sqr,im_grappa_err_sqr),[],'InitialMagnification',150);
title ('Difference images: SPIR-iT CG              GRAPPA');



