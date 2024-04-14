function [I_recon_SOS, kspace_coils] = GRAPPA_interpolate_kSpace_2d(data_undersampled_KyKxC, header, GRAPPA_weights)
% Reconstruct k-space by input GRAPPA weights in 2D
% Assumes GRAPPA_weights are calibrated for both kx and ky

Nblock = header.blocks;
Ncolumn = header.column;
Rpe = header.subsampling_factor(1);
Rfe = header.subsampling_factor(2);
[Npe, Nfe, Ncoil] = size(data_undersampled_KyKxC);

% Define padding margins
margin_top_pe = Rpe * (Nblock / 2 + 1);
margin_bottom_pe = margin_top_pe;
margin_left_fe = Rfe * (Ncolumn / 2 + 1);
margin_right_fe = margin_left_fe;

% Determine which PE and FE lines are acquired
% isLineAcquiredPE = any(any(data_undersampled_KyKxC ~= 0, 3), 2);
% isLineAcquiredFE = any(any(data_undersampled_KyKxC ~= 0, 3), 1);
firstAcquirePoint_ky=find(sum(abs(data_undersampled_KyKxC(:,:,1)),2),1);
firstAcquirePoint_kx=find(sum(abs(data_undersampled_KyKxC(:,:,1)),1),1);
allAcquiredLinesPE = firstAcquirePoint_ky:Rpe:Npe;
allAcquiredLinesFE = firstAcquirePoint_kx:Rfe:Nfe;
% allMissingLinesPE = find(~isLineAcquiredPE);
% allMissingLinesFE = find(~isLineAcquiredFE);

% Pad data to handle borders effectively
paddedData_KyKxC = padarray(data_undersampled_KyKxC, [margin_top_pe, margin_left_fe], 'both');

%% Interpolation in both PE and FE dimensions
for iCoil = 1:Ncoil
    for iPattern = 1:(Rpe * Rfe - 1) % For each interpolation pattern excluding the fully sampled
        [x, y] = ind2sub([Rpe, Rfe], iPattern + 1);
        iPattern_pe = x - 1; % Pattern shifts in PE
        iPattern_fe = y - 1; % Pattern shifts in FE

        % Intersect the target ranges with the missing lines
%         targetPE_range = intersect(allAcquiredLinesPE + iPattern_pe, allMissingLinesPE);
%         targetFE_range = intersect(allAcquiredLinesFE + iPattern_fe, allMissingLinesFE);
        targetPE_range = allAcquiredLinesPE + iPattern_pe;
        targetFE_range = allAcquiredLinesFE + iPattern_fe;
        
        NtargetPE = length(targetPE_range);
        NtargetFE = length(targetFE_range);
        SourceLines_thisPattern = complex(zeros(Nblock, Ncolumn, NtargetPE, NtargetFE, Ncoil));

        % Extract source blocks for interpolation
        for iBlock = 1:Nblock
            for iColumn = 1:Ncolumn
                iBlock_offset = -iPattern_pe - Rpe * (Nblock / 2 - 1) + (iBlock - 1) * Rpe;
                iColumn_offset = -iPattern_fe - Rfe * (Ncolumn / 2 - 1) + (iColumn - 1) * Rfe;
                SourceLines_thisPattern(iBlock, iColumn, :, :, :) = paddedData_KyKxC(...
                    targetPE_range + margin_top_pe + iBlock_offset, ...
                    targetFE_range + margin_left_fe + iColumn_offset, :);
            end
        end

        % Reshape source pattern matrix for GRAPPA weight application
        SourceMatrix_thisPattern = reshape(permute(SourceLines_thisPattern, [3, 4, 1, 2, 5]), ...
                                           NtargetPE * NtargetFE, Nblock * Ncolumn * Ncoil);

        % Apply GRAPPA weights to interpolate missing data
        interpolated_kSpace = SourceMatrix_thisPattern * squeeze(GRAPPA_weights(iPattern, iCoil, :));

        % Place interpolated data back into the padded k-space array
        paddedData_KyKxC(targetPE_range + margin_top_pe, targetFE_range + margin_left_fe, iCoil) = ...
            reshape(interpolated_kSpace, [NtargetPE, NtargetFE]);
    end
end

% Crop the padded k-space to its original size
kspace_coils = paddedData_KyKxC(1 + margin_top_pe:end - margin_bottom_pe, 1 + margin_left_fe:end - margin_right_fe, :);
% data consistency enforcement for the central ACS lines if any
kspace_coils(data_undersampled_KyKxC~=0)=data_undersampled_KyKxC(data_undersampled_KyKxC~=0);
% Compute the reconstructed image using Sum of Squares (SOS)
I_recon_SOS = sqrt(sum(abs(ifft2c(kspace_coils)).^2, 3));

end
