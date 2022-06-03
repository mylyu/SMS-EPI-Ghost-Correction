function [GRAPPA_weights] = GRAPPA_calibrate_weights(ACS_kykxc, header, regularization_factor)
% Calculate GRAPPA weights (both kx and ky used)
% with L2 norm relguarization
% by Lyu Mengye and Victor Bin Xie
% based on PULSA toolbox
%Nblocks   Number of reference blocks for GRAPPA reconstruction
%Ncolumn   Number of column blocks (when equals to 1, 1D GRAPPA)
import GRAPPA.*
if nargin<3 || isempty(regularization_factor)
    regularization_factor=0.001;
end
disp(['regularization_factor ' num2str(regularization_factor)])
Nblock=header.blocks;
Ncolumn=header.column;
R=header.subsampling_factor;
% Npe=header.Npe;
% Nfe=header.Nfe;
Npe=size(ACS_kykxc,1);
Nfe=size(ACS_kykxc,2);
Ncoil=size(ACS_kykxc,3);

margin_top_pe=R*(Nblock/2+1);
margin_bottom_pe=R*(Nblock/2+1);
margin_left_fe=(Ncolumn-1)/2+1;
margin_right_fe=(Ncolumn-1)/2+1;
targetPE_range=margin_top_pe:Npe-margin_bottom_pe+1;
targetFE_range=margin_left_fe:Nfe-margin_right_fe+1;
GRAPPA_weights=complex(zeros(R-1,Ncoil,Nblock*Ncolumn*Ncoil));
NtargetPE=length(targetPE_range);
NtargetFE=length(targetFE_range);
% Coefficient calculation using ACS lines
for iCoil = 1:Ncoil
    %     disp(['coil = ' num2str(z_index)]);
    for iPattern = 1:R - 1 %i.e. how far this pattern is away from last sample in PE
        TargetLines = ACS_kykxc(targetPE_range, targetFE_range, iCoil);
        TargetLines =  TargetLines(:);
        SourceLines_thisPattern=complex(zeros(Nblock,Ncolumn,NtargetPE,NtargetFE,Ncoil));
        for iBlock=1:Nblock
            for iColumn=1:Ncolumn
                iBlock_offset=-iPattern-R*(Nblock/2-1)+(iBlock-1)*R;
                iColumn_offset=-(Ncolumn-1)/2+iColumn-1;
                SourceLines_thisPattern(iBlock,iColumn,:,:,:)=ACS_kykxc...
                    (targetPE_range+iBlock_offset, targetFE_range+iColumn_offset, :);
            end
        end
        SourceMatrix_thisPattern=reshape(permute(SourceLines_thisPattern,[3 4 1 2 5]),NtargetPE*NtargetFE,Nblock*Ncolumn*Ncoil);
        % use pinv is easy to understand, but slow
        % coefficient = pinv(SourceMatrix_thisPattern) * TargetLines;
        % use \ in matlab is fast
        %         coefficient = SourceMatrix_thisPattern\TargetLines;
        % L2 norm regularize
        A=SourceMatrix_thisPattern;
        AHA = A'* A;
        reduced_eye = diag(abs(diag(AHA))>0);
        n_alias = sum(reduced_eye(:));
        scaled_reg_factor = regularization_factor * trace(AHA)/n_alias;
%         coefficient = pinv(AHA + reduced_eye .* scaled_reg_factor) * A'*TargetLines;
        coefficient = (AHA + reduced_eye .* scaled_reg_factor) \ (A'*TargetLines);
        GRAPPA_weights(iPattern, iCoil, :) = coefficient;
    end
end