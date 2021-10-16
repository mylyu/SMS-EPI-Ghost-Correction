function [unmix, gmap] = ismrm_calculate_sense_unmixing_general2(acc_factor, csm, noise_matrix, regularization_factor)
%
% function [unmix, gmap] = ismrm_calculate_sense_unmixing(acc_factor, csm, noise_matrix, regularization_factor)
%
% Calculates the unmixing coefficients for a 2D image
%
% INPUT:
%       acc_factor  2D             : Acceleration factor, e.g. [2 3]
%       csm         [x, y, coil]       : Coil sensitivity map
%       noise_matrix [coil, coil]      : noise covariance matrix
%       regularization_factor scaler   : adds Tychonov regularization.
%                                        0 = no regularization
%                                        0.001 = default
%                                        set higher for more aggressive
%                                        regularization.
%
% OUTPUT:
%       unmix       [x, y, coil] : Image unmixing coefficients for a single x location
%       gmap        [x, y]       : Noise enhancement map
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
%
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip Beatty (philip.beatty@sri.utoronto.ca)
%
import LmyParallelImaging.*
assert(nargin > 1, 'At least 2 arguments needed');
assert(ndims(csm)==3, 'coil sensitivity map must have 3 dimensions');

if nargin < 3,
    noise_matrix = [];
end

if nargin < 4,
    regularization_factor = [];
end

if isempty(noise_matrix),
    noise_matrix = eye(size(csm,3));
end
if isempty(regularization_factor),
    regularization_factor = 0.00;
end

noise_matrix_inv = pinv(noise_matrix);

if length(acc_factor)==1
    unmix=zeros(size(csm));
    for x=1:size(csm,1),
        unmix(x,:,:) = ismrm_calculate_sense_unmixing_1d(acc_factor, squeeze(csm(x,:,:)), noise_matrix_inv, regularization_factor);
    end
elseif length(acc_factor)==2
    unmix = ismrm_calculate_sense_unmixing_2d_actual(acc_factor, csm, noise_matrix_inv, regularization_factor);
else
    error('acceleration factor cannot have more than 3 elements')
end

if (nargout > 1),
    gmap = sqrt(sum(abs(unmix).^2,3)).*sqrt(sum(abs(csm).^2,3));
end

function [unmix2d] = ismrm_calculate_sense_unmixing_2d_actual(acc_factor, csm2d, noise_matrix_inv, regularization_factor)

[nx, ny, nc] = size(csm2d);
if mod(nx, acc_factor(1)) ~= 0,
    error('nx must be a multiple of acc_factor(1)');
end
if mod(ny, acc_factor(2)) ~= 0,
    error('ny must be a multiple of acc_factor(2)');
end

acc_factor_x = acc_factor(1);
acc_factor_y = acc_factor(2);
unmix2d = zeros(nx,ny, nc);

nx_blocks = nx/acc_factor_x;
ny_blocks = ny/acc_factor_y;
for ix_blocks=1:nx_blocks
    for iy_blocks = 1:ny_blocks
        A = reshape(csm2d(ix_blocks:nx_blocks:nx,iy_blocks:ny_blocks:ny,:),acc_factor_x*acc_factor_y,[]).';
        if max(abs(A(:))) > 0.01
            %        unmix1d(index:n_blocks:ny, :) = pinv(A);
            pixels_unaliased=zeros(size(A,2),size(A,1));
            nonzeros_index=(sum(abs(A),1)>0.01);
            
            A_reduced = A(:,nonzeros_index);
            AHA = A_reduced' * A_reduced;
            reduced_eye = diag(abs(diag(AHA))>0.01);
            n_alias = sum(reduced_eye(:));
            scaled_reg_factor = regularization_factor * trace(AHA)/n_alias;
            
            pixels_unaliased(nonzeros_index,:)=reshape((AHA + reduced_eye .* scaled_reg_factor)...
                \ A_reduced',...
                size(A_reduced,2),nc);
            unmix2d(ix_blocks:nx_blocks:nx,iy_blocks:ny_blocks:ny, :) = reshape(pixels_unaliased,acc_factor_x, acc_factor_y, nc);
        end
    end
end

function [unmix1d] = ismrm_calculate_sense_unmixing_1d(acc_factor, csm1d, noise_matrix_inv, regularization_factor)

[ny, nc] = size(csm1d);

if mod(ny, acc_factor) ~= 0,
    error('ny must be a multiple of acc_factor');
end

unmix1d = zeros(ny, nc);

n_blocks = ny/acc_factor;
for index = 1:n_blocks
    A = csm1d(index:n_blocks:ny,:).';
    if max(abs(A(:))) > 0,  
%        unmix1d(index:n_blocks:ny, :) = pinv(A);
        AHA = A'*noise_matrix_inv * A;
        reduced_eye = diag(abs(diag(AHA))>0);
        n_alias = sum(reduced_eye(:));
        scaled_reg_factor = regularization_factor * trace(AHA)/n_alias;
        
        unmix1d(index:n_blocks:ny, :) = pinv(AHA + reduced_eye .* scaled_reg_factor) * A' * noise_matrix_inv;
    end
end

