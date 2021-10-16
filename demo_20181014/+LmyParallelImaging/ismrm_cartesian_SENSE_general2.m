function [img,gmap,snr,snr_pseudo,gmap_pseudo,noise_psf_pseudo] = ismrm_cartesian_SENSE_general2(inp,csm,acc_factor,replicas,samp_mat,reg)
%
%   [img,gmap,snr_pseudo,gmap_pseudo,noise_psf_pseudo] = ismrm_cartesian_SENSE(inp,csm,acc_factor,replicas)
%
%   Cartesian SENSE reconstruction.
%
%   It is recommended to input data with pre-whitened noise scale to sd=1.
%
%
%   INPUT:
%     - inp         [kx,ky,coils]        : Undersampled input k-space data
%     - acc+factor  [acc_x,acc_y]        : 2D or 1D Parallel Imaging Acceleration factor
%     - csm         [x,y,coil]           : Optional coil sensitivity map (for coil combination)
%     - replicas    scalar (dafault 100) : Number of replicas to run if SNR
%                                          is requested
%     - samp_mat    [kx, ky]             : Sampling matrix 1 or 0.
%
%   OUTPUT:
%     - img              [x,y]                : Output image
%     - gmap             [x,y]                : g-map (calulated from unmixing coefficients)
%     - snr              [x,y]                : An image in SNR units. From direct recon.
%     - snr_pseudo       [x,y]                : An image in SNR units. Using psudo replica method.
%     - gmap_pseudo      [x,y]                : A g-map (Pseudo replica)
%     - noise_psf_pseudo [x,y]                : Point spread function of the noise
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%
import LmyParallelImaging.*
if nargin<4,
    replicas = 100;
end

if nargin<5||isempty(samp_mat),
   samp_mat = (sum(abs(inp),3) > 0); 
end
% added by mengye
if nargin<6,
   reg = []; 
end

[unmix_sense, gmap]   = ismrm_calculate_sense_unmixing_general2(acc_factor, csm,[],reg);


img_alias = sqrt(prod(acc_factor))*ismrm_transform_kspace_to_image(inp,[1,2]);
img = sum(img_alias .* unmix_sense,3);
snr = img ./ sqrt(sum(abs(unmix_sense).^2,3));


if (nargout > 3),
    image_formation_func = @(x) sum(ismrm_transform_kspace_to_image(sqrt(prod(acc_factor))*x .* repmat(samp_mat,[1 1 size(csm,3)]),[1,2]).*unmix_sense,3);
    [snr_pseudo,gmap_pseudo,noise_psf_pseudo] = ismrm_pseudo_replica(inp, image_formation_func,replicas);
    gmap_pseudo = gmap_pseudo .* sqrt(sum(abs(csm).^2,3));
end

