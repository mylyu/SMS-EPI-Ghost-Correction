function [res,LIST, KEY] = GRAPPA_lustig(kData,kCalib,kSize,lambda,selected_coils2recon,LIST, KEY)
% res = GRAPPA(kData,kCalib,kSize,lambda [, disp)
%
% This is a GRAPPA reconstruction algorithm that supports
% arbitrary Cartesian sampling. However, the implementation
% is highly inefficient in Matlab because it uses for loops.
% This implementation is very similar to the GE ARC implementation.
%
% The reconstruction looks at a neighborhood of a point and
% does a calibration according to the neighborhood to synthesize
% the missing point. This is a k-space varying interpolation.
% A sampling configuration is stored in a list, and retrieved
% when needed to accelerate the reconstruction (a bit)
%
% Inputs:
%       kData     -   [Size x, Size y, num coils] 2D multi-coil k-space data to reconstruct from.
%                   Make sure that the missing entries have exact zeros in them.
%       kCalib    -   calibration data (fully sampled k-space)
%       kSize     - size of the 2D GRAPPA kernel [kx, ky]
%       lambda    - Tykhonov regularization for the kernel calibration.
%       dispp      - Figure number to display images as they are
%                   reconstructed
% Outputs:
%       res       - k-space data where missing entries have been filled in.
%
% Example:
%   [x,y] = meshgrid(linspace(0,1,128));
%   % Generate fake Sensitivity maps
%   sMaps = cat(3,x.^2,1-x.^2,y.^2,1-y.^2);
%   % generate 4 coil phantom
%   imgs = repmat(phantom(128),[1,1,4]).*sMaps;
%   DATA = fft2c(imgs);
%   % crop 20x20 window from the center of k-space for calibration
%   kCalib = crop(DATA,[20,20,4]);
%
%   %calibrate a kernel
%   kSize = [5,5];
%   coils = 4;
%
%   % undersample by a factor of 2
%   DATA(1:2:end,2:2:end,:) = 0;
%   DATA(2:2:end,1:2:end,:) = 0;
%
%   %reconstruct:
%   [res] = GRAPPA(DATA,kCalib, kSize, 0.01);
%   figure, imshow(cat(2,sos(imgs), 2*sos(ifft2c(DATA)), sos(ifft2c(res))),[]);
%   title('full,  zero-fill,   result')
%
%
%
% (c) Michael Lustig 2008

    flag_do_calib=0;
if nargin < 7||isempty(KEY)||isempty(LIST)
    KEY=[];LIST=[];
    flag_do_calib=1;
    disp('No input of GRAPPA weights, calibration will be done')
end
if nargin < 5||isempty(selected_coils2recon)
    %     dispp = 0;
    selected_coils2recon=1:size(kData,3);
end
pe = size(kData,2); fe = size(kData,1); coils = size(kData,3); % get sizes
res = kData*0;

if flag_do_calib==1
    %[AtA] = corrMatrix(kCalib,kSize); % build coil correlation matrix
    AtA = dat2AtA(kCalib, kSize); % build coil calibrating matrix
    for n=1:numel(selected_coils2recon)
%         disp(sprintf('calibrate and recon coil %d',selected_coils2recon(n)));
        [res(:,:,n), LIST(:,:,n), KEY(:,:,n)] = ARC(kData, AtA,kSize, selected_coils2recon(n),lambda, [], []); % reconstruct single coil image
    end
else
    for n=1:numel(selected_coils2recon)
%         disp(sprintf('recon coil %d',selected_coils2recon(n)));
        [res(:,:,n)] = ARC(kData, [],kSize, selected_coils2recon(n),lambda, squeeze(LIST(:,:,n)), squeeze(KEY(:,:,n))); % reconstruct single coil image
    end
end


function [res, LIST, KEY] = ARC(kData, AtA, kSize, c, lambda, LIST, KEY)
[sx,sy,nCoil] = size(kData);
kData = zpad(kData,[sx+kSize(1)-1, sy+kSize(2)-1,nCoil]);

dummyK = zeros(kSize(1),kSize(2),nCoil); dummyK((end+1)/2,(end+1)/2,c) = 1;
idxy = find(dummyK);

res = zeros(sx,sy);

MaxListLen = 100;
if isempty(KEY)||isempty(LIST)
    LIST = zeros(kSize(1)*kSize(2)*nCoil,MaxListLen);
    KEY =  zeros(kSize(1)*kSize(2)*nCoil,MaxListLen);
end
count = 0;

%H = waitbar(0);
for y = 1:sy
    for x=1:sx
        %	waitbar((x + (y-1)*sx)/sx/sy,H);
        tmp = kData(x:x+kSize(1)-1,y:y+kSize(2)-1,:);
        pat = abs(tmp)>0;
        if pat(idxy) | sum(pat)==0
            res(x,y) = tmp(idxy);
        else
            key = pat(:);
            idx = 0;
            for nn=1:size(KEY,2);
                if sum(key==KEY(:,nn))==length(key)
                    idx = nn;
                    break;
                end
            end
            if idx == 0
                if isempty(KEY)||isempty(LIST)||isempty(AtA)
                    warning('No GRAPPA weights are provided for one of the patterns, pls check input')
                end
                count = count + 1;
                kernel = calibrate(AtA,kSize,nCoil,c,lambda,pat);
                KEY(:,mod(count,MaxListLen)+1) = key(:);
                LIST(:,mod(count,MaxListLen)+1) = kernel(:);
%                 disp('add another key');%size(KEY,2)
            else
                kernel = LIST(:,idx);
            end
            res(x,y) = sum(kernel(:).*tmp(:));
        end
        
    end
end

%close(H);



