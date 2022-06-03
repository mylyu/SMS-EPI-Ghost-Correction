function [I_recon_SOS, kspace_coils] = GRAPPA_interpolate_kSpace(data_undersampled_KyKxC, header,GRAPPA_weights)
% Reconstruct each interleave by input grappa weights
% savedcoeff: Grappa weights
% subsampled_kspace_data: data from each interleave
% ReconKspace: Kspace data reconstructed by grappa
% I_recon_SOS: reconstructed image
% Orginal Version by PULSA toolbox and Victor Bin Xie
% Clearer implementation by Lyu Mengye
Nblock=header.blocks;
Ncolumn=header.column;
R=header.subsampling_factor;
[Npe,Nfe,Ncoil]=size(data_undersampled_KyKxC);
%%
margin_top_pe=R*(Nblock/2+1);
margin_bottom_pe=R*(Nblock/2+1);
margin_left_fe=(Ncolumn-1)/2+1;
margin_right_fe=(Ncolumn-1)/2+1;
isThisLineAcquired=any(any(data_undersampled_KyKxC~=0,3),2);
allAcquiredLines=find(isThisLineAcquired);
allMissingLines=find(~isThisLineAcquired);
% firstMissingLine=find(~isThisLineAcquired,1);% find the first missing line;
paddedData_KyKxC=padarray(data_undersampled_KyKxC,[margin_top_pe margin_left_fe],'both');
% to deal with situations where the first several lines are not acquired
allAcquiredLines_extended=[allAcquiredLines(1)-R;allAcquiredLines];
%% interpolation
for iCoil = 1:Ncoil
    for iPattern = 1:R - 1 %i.e. how far this pattern is away from last sample in PE
        % handle the lines around the margin and skip ACS lines and lines out of range
    	  targetPE_range = intersect(allAcquiredLines_extended+iPattern,allMissingLines);
          targetFE_range = 1:Nfe;
          NtargetPE=length(targetPE_range);
          SourceLines_thisPattern = complex(zeros(Nblock,Ncolumn,NtargetPE,Nfe,Ncoil));
        for iBlock=1:Nblock
            for iColumn=1:Ncolumn
                iBlock_offset=-iPattern-R*(Nblock/2-1)+(iBlock-1)*R;
                iColumn_offset=-(Ncolumn-1)/2+iColumn-1;  
                SourceLines_thisPattern(iBlock,iColumn,:,:,:)=paddedData_KyKxC...
                    (targetPE_range+margin_top_pe+iBlock_offset, targetFE_range+margin_left_fe+iColumn_offset, :);
            end
        end
        SourceMatrix_thisPattern=reshape(permute(SourceLines_thisPattern,[3 4 1 2 5]),NtargetPE*Nfe,Nblock*Ncolumn*Ncoil);
        % use \ in matlab is fast
        paddedData_KyKxC(targetPE_range+margin_top_pe,targetFE_range+margin_left_fe,iCoil) ...
            = reshape(SourceMatrix_thisPattern*squeeze(GRAPPA_weights(iPattern, iCoil, :)),NtargetPE,Nfe);
    end
end
kspace_coils=paddedData_KyKxC(1+margin_top_pe:end-margin_bottom_pe,1+margin_left_fe:end-margin_right_fe,:); % crop
%%
I_recon_SOS=sqrt(sum(abs(ifft2c(kspace_coils)).^2,3));