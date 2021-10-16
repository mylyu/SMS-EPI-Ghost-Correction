if exist('spirit.m')==0
	addpath('utils')
end

if exist('nufft.m')==0
    addpath ('nufft_files');
end


% Initial Parameters

kSize = [6,6];                  % SPIRiT Kernel size
CalibSize = [30,30];            % size of the calibration region
N = [260,360];                  % size of the target image
nIterCG = 15;                   % number of reconstruction iterations
CalibTyk = 0.02;                % Tykhonov regularization for calibration 0.01-0.05 recommended
accel = 3;                      % Acceleration factor

eigThresh_k = 0.02; % threshold of eigenvectors in k-space
eigThresh_im = 0.95; % threshold of eigenvectors in image space

% load data
load spiral.mat

% perform SVD coil compression
disp('perform coil compression at  5\% tollerance.')
D = reshape(data,size(data,1)*size(data,2),size(data,3));
[U,S,V] = svd(D,'econ');
nCoil = max(find(diag(S)/S(1)>0.05));
data = reshape(D*V(:,1:nCoil),size(data,1),size(data,2),nCoil);
disp(sprintf('Using %d virtual channels',nCoil'));
%nCoil = 8;

% Calibration
% In this example, we calibrate from a fully sampled acquisition
% reconstructed with gridding and density compensation.
% Then undersample it to see how the reconstruction performs.

disp('Starting Calibration');
disp('Generating NUFFT object for calibration')

% Gridding and density compensation reconstruction
GFFT1 = NUFFT(k,w, [0,0] , N);
im = GFFT1'*(data.*repmat(sqrt(w),[1,1,nCoil]));
kData = fft2c(im);


kCalib = crop(kData,[CalibSize,nCoil]); % crop center k-space for calibration

[kernel,s] = dat2Kernel(kCalib,kSize);
idx = max(find(S >= S(1)*eigThresh_k));
[M,W] = kernelEig(kernel(:,:,:,1:idx),N);
maps = M(:,:,:,end);
weights = double(W(:,:,end) >  eigThresh_im);
ESP = ESPIRiT(maps,weights);

disp('Done Calibrating');


% Undersample the data and prepare new NUFFT operator for it
idx = (1:accel:size(k,2));
k_u = k(:,idx);
w_u = w(:,idx);  % use w_u = w(:,idx)*0+1; if you don't want to use density weighting
                 % this may improve noise, but will converge slower. use
                 % larger lambda.
kData_u = data(:,idx,:);

disp('generating nufft object for recon')
GFFT_u = NUFFT(k_u,w_u, [0,0], N);

im_dc = GFFT_u'*(kData_u.*repmat(sqrt(w_u),[1,1,nCoil]))*accel;
res = im_dc;
disp('initiating reconstruction')
tic

XOP = Wavelet('Daubechies',4,6);
[res] = cgL1ESPIRiT(kData_u.*repmat(sqrt(w_u),[1,1,nCoil]), ESP'*im_dc, GFFT_u, ESP, nIterCG,XOP,1,0.5,20);
toc
disp('done!');

figure(100), imshow(cat(2,sos(im),sos(im_dc),abs(res)),[]), 
