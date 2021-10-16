%% Speedup of Geometric Decomposition Coil Compression
% This is a demo of the Geometric-decomposition Coil Compression (GCC) based on  
% Zhang et. al MRM 2013;69(2):571-82. The demo shows how to use the compression
% and in addition it shows the use of alignment of the compression
% matrices. The example is on 2D data, but can also be applied on 3D as
% well. GCC removes correlations between channels along a fully sampled
% readout dimension. This allows for significant reduction in data, so
% parallel imaging and other multi-channel algorithms are much more
% computationally efficient.  This demo demonstrate effective compression
% from 32->5 virtual channels. 
%
% This demo aims to show speedup of computation for parallel imaging



%% 
% Set parameters:

load brain_32ch.mat
im = ifft2c(DATA);
DATA = DATA/max(abs(im(:)));

% number of compressed virtual coils
ncc = 5;
dim = 2;
ncalib = 20; % use 24 calibration lines to compute compression
[sx,sy,Nc] = size(DATA);
slwin = 5; % sliding window length for GCC

mask = zpad(ones(ncalib,sy,Nc),[sx,sy,Nc]);
mask(1:2:end,:,:) = 1;
dens = -zpad(ones(ncalib,sy,Nc),[sx,sy,Nc]) + 2; 
DATAc = DATA.*mask;
% crop calibration data
calib = crop(DATAc,[ncalib,sy,Nc]);

%%

disp('Performing GCC')
tic
gccmtx = calcGCCMtx(calib,dim,slwin);
gccmtx_aligned = alignCCMtx(gccmtx(:,1:ncc,:));
DATAcc = CC(DATAc,gccmtx_aligned,dim);
calibcc = CC(calib,gccmtx_aligned,dim);
gcct = toc;

disp('perfomring GRAPPA on GCC coils');
tic; resGRAPPAcc = GRAPPA(DATAcc,calibcc,[5,5],0.01);grappacct = toc;
disp('performing GRAPPA on physical coils');
tic; resGRAPPA = GRAPPA(DATAc,calib,[5,5],0.01);grappat = toc;

errGcc = unCC(resGRAPPAcc,dim,gccmtx_aligned)-DATA;
errG = resGRAPPA - DATA;

%%
% GCC has many order of magnitude improvement in recon time. Note thst the
% error image does not have any aliasing patterns in it and is coming from the compression
% error, rather than the parallel imaging recon

figure, imshow(cat(2, sos(ifft2c(resGRAPPAcc)), sos(ifft2c(resGRAPPA)), sos(ifft2c(DATAc.*dens))),[]);
title('Left: GRAPPA with GCC, middle GRAPPA, right: Zero-filling');

figure, imshow(cat(2, sos(ifft2c(errGcc)), sos(ifft2c(errG))),[]);
title('Reconstruction error in %. Left: GRAPPA with GCC, right GRAPPA');
colorbar

disp(sprintf('Speedup is: %f times!', grappat/(grappacct + gcct)));
