function res = SAKE(DATA, kSize, wnRank,nIter, show)
%
% res = SAKE(DATA, kSize, wnRank,nIter)
%
% function performs an autocalibration reconstrution without calibration
% lines based on the SAKE method (P. Shin et. al, "Calibrationless Parallel 
% Imaging Reconstruction Based on Structured Low-Rank Matrix Completion"
% 2013, submitted to MRM. 
%
% It is recommended that this method be used to complete a calibration
% region and then use ESPIRiT to generate ESPIRiT maps. 
%
% INPUTS:
%       DATA [nx X ny X nc]  - 2D multi-coil k-space data with missing
%                               entries, preferably random! missing entries
%                               should be EXACTLY = 0!
%       kSize [ 1 x 2]       - sliding window size
%       wnRank [ 1 x 1]       - the window-normalized rank to enforce 
%                               (# of coils = full rank, 2.0 would be typical 
%                               for 8 coils with window size [5 5]
%       nIter                - number of iTerations.
%       show                 - show intermediate reconstruction show=0 will
%                               skip. Show = 100 will plot in figure 100
%
%
% Outputs:
%       res -                - 2D multi coil k-space data where the missing
%                               data was filled. 
%
% (c) Michael Lustig 2012
%


[sx,sy,nc] = size(DATA);
mask = abs(DATA)>0;

res = DATA;
for n=1:nIter
    
    % reorder data to get Hankel structure. 
    tmp = im2row(res,kSize); [tsx,tsy,tsz] = size(tmp);
    A = reshape(tmp,tsx,tsy*tsz);
    
    % SVD thresholding
    [U,S,V] = svd(A,'econ');
    keep = 1:floor(wnRank*prod(kSize));
    A = U(:,keep)*S(keep,keep)*V(:,keep)';
    
    % Enforce Hankel structure
           A = reshape(A,tsx,tsy,tsz);
    tmp = row2im(A,[sx,sy,nc],kSize);
    
    % enforce data consistency
    res = tmp.*(1-mask) + DATA;
    
    if show
        if max(size(res)) < 100
            figure(show), imshow(sos(ifft2c(res)),[],'InitialMagnification',200);
            drawnow;
        else
            figure(show), imshow(sos(ifft2c(res)),[])
        drawnow;
        n
        end
    end
end



