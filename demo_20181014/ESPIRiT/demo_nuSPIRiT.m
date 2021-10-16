if exist('spirit.m')==0
	addpath('utils')
end

if exist('nufft.m')==0
    addpath ('nufft_files');
end


% Initial Parameters

kSize = [7,7];                  % SPIRiT Kernel size
CalibSize = [30,30];            % size of the calibration region
N = [260,360];                  % size of the target image
nIterCG = 15;                   % number of reconstruction iterations
CalibTyk = 0.02;                % Tykhonov regularization for calibration 0.01-0.05 recommended
lambda = 1;                     % Ratio between data and calibration consistency. 1 recommended when density compensation is used.
accel = 3;                      % Acceleration factor

% load data
load spiral.mat

% perform SVD coil compression
disp('perform coil compression at  5\% tollerance.')
D = reshape(data,size(data,1)*size(data,2),size(data,3));
[U,S,V] = svd(D,'econ');
nCoil = max(find(diag(S)/S(1)>0.05));
data = reshape(D*V(:,1:nCoil),size(data,1),size(data,2),nCoil);
disp(sprintf('Using %d virtual channels',nCoil'));


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


kernel = zeros([kSize,nCoil,nCoil]);
kCalib = crop(kData,[CalibSize,nCoil]); % crop center k-space for calibration

%prepare calibration matrix for calibration ( Phil Beatty calls this
%correlation values matrix. See thesis by Phil Beatty.)
%[AtA,] = corrMatrix(kCalib,kSize);
%for n=1:nCoil
%    disp(sprintf('Calibrating coil %d',n));
%	kernel(:,:,:,n) = calibrate(AtA,kSize,nCoil,n,CalibTyk);
%end

kernel = calibSPIRiT(kCalib, kSize, nCoil, CalibTyk);
GOP = SPIRiT(kernel, 'image',N); % init the SPIRiT Operator
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
res = cgNUSPIRiT(kData_u.*repmat(sqrt(w_u),[1,1,nCoil]),res,GFFT_u,GOP,nIterCG,lambda);
toc
disp('done!');


figure(100), imshow(cat(2,sos(im),sos(im_dc),sos(res)),[]), 
