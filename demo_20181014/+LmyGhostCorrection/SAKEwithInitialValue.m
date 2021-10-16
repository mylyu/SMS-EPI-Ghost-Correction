function res = SAKEwithInitialValue(DATA, kSize, wnRank,nIter, show, mask)
import LmyUtility.*
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
% mask = abs(DATA)>0;
LmyUtility.textprogressbar(['totally ' num2str(nIter) ' iterations '])
res = DATA;
for n=1:nIter
    LmyUtility.textprogressbar(n/nIter*100)
% %     % shift back, temporary use
% %     res=cat(3,res(:,:,1:end/2),circshift(res(:,:,end/2+1:end),[0 1]));
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
    % enforce magnitude of virtual coils
    img_1=ifft2c(tmp(:,:,1:end/2));
    img_2=ifft2c(tmp(:,:,end/2+1:end));
    mean_abs=(abs(img_1)+abs(img_2))/2;
    tmp=fft2c(cat(3,mean_abs.*exp(1i*angle(img_1)),mean_abs.*exp(1i*angle(img_2))));
% %     % enforce both magnitude and pem consistency
% %     img_1=ifft2c(tmp(:,:,1:end/2));
% %     img_2=ifft2c(tmp(:,:,end/2+1:end));
% %     mean_abs=(abs(img_1)+abs(img_2))/2;
% %     img_1_chComb=sos(img_1);
% %     img_2_chComb=comb_ch(img_2,img_1./repmat(img_1_chComb,[1 1 nc/2]));
% %     phase_diff_chComb=repmat(img_2_chComb.*conj(img_1_chComb),[1 1 nc/2]);
% %     tmp=fft2c(cat(3,mean_abs.*exp(1i*angle(img_1)),mean_abs.*exp(1i*angle(img_1)).*exp(1i*angle(phase_diff_chComb))));
% % %     % shift to virtual coil pattern, temporary use
% % %     tmp=cat(3,tmp(:,:,1:end/2),circshift(tmp(:,:,end/2+1:end),[0 -1]));
    % enforce data consistency
    res = tmp.*(1-mask) + DATA.*mask;
   
    
    if show
        disp_max=prctile(abs(res(:)),99);
        if size(res,1) < 300
            figure(show), imshow(reshape((abs(ifft2c(res(:,:,1:ceil(end/4/4):ceil(end/4))))),size(res,1),[]),[0 disp_max/3],'InitialMagnification',40000/size(res,1));
%             figure(show), imshow(reshape((abs(ifft2c(res(:,:,[9 13 17 22])))),size(res,1),[]),[0 disp_max/3],'InitialMagnification',30000/size(res,1));

    drawnow;
        else
            figure(show), imshow((abs(ifft2c(res(:,:,1)))),[0 disp_max/3])
        drawnow;
        n
        end
    end
end
LmyUtility.textprogressbar('finished')



