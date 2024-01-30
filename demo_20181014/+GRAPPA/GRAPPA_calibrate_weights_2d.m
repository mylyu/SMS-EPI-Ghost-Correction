function [GRAPPA_weights] = GRAPPA_calibrate_weights_2d(ACS_kykxc, header, regularization_factor)
% Calculate 2D GRAPPA weights (both kx and ky used)
% with L2 norm relguarization
% by Lyu Mengye and Victor Bin Xie
% based on PULSA toolbox

import GRAPPA.*
if nargin<3 || isempty(regularization_factor)
    regularization_factor=0.001;
end
disp(['regularization_factor ' num2str(regularization_factor)])
nBlocks=ceil(header.blocks/2)*2;
nColumns=ceil(header.column/2)*2;
Rpe=header.subsampling_factor(1);
Rfe=header.subsampling_factor(2);
% Npe=header.Npe;
% Nfe=header.Nfe;
Npe=size(ACS_kykxc,1);
Nfe=size(ACS_kykxc,2);
Ncoil=size(ACS_kykxc,3);

margin_top_pe=Rpe*(nBlocks/2+1);
margin_bottom_pe=Rpe*(nBlocks/2+1);
margin_left_fe=Rfe*(nColumns/2+1);
margin_right_fe=Rfe*(nColumns/2+1);
targetPE_range=margin_top_pe:Npe-margin_bottom_pe+1;
targetFE_range=margin_left_fe:Nfe-margin_right_fe+1;
% GRAPPA_weights=complex(zeros(R-1,Ncoil,Kpe*Kfe*Ncoil));
GRAPPA_weights=complex(zeros(Rpe*Rfe-1,Ncoil,nBlocks*nColumns*Ncoil,class(ACS_kykxc)));
NtargetPE=length(targetPE_range);
NtargetFE=length(targetFE_range);
% Coefficient calculation using ACS lines
for iCoil = 1:Ncoil
    %     disp(['coil = ' num2str(z_index)]);
    for iType = 1:Rpe*Rfe-1 %i.e. index away from last sample in uppper left direction
        [x, y] = ind2sub([Rpe,Rfe],iType+1);
        iPattern_pe=x-1;iPattern_fe=y-1;
        TargetLines = ACS_kykxc(targetPE_range, targetFE_range, iCoil);
        TargetLines =  TargetLines(:);
        SourceLines_thisPattern=complex(zeros(nBlocks,nColumns,NtargetPE,NtargetFE,Ncoil,class(ACS_kykxc)));
        for iBlock=1:nBlocks
            for iColumn=1:nColumns
                iBlock_offset=-iPattern_pe-Rpe*(nBlocks/2-1)+(iBlock-1)*Rpe;
%                 iColumn_offset=-(Kfe-1)/2+iColumn-1;
                iColumn_offset=-iPattern_fe-Rfe*(nColumns/2-1)+(iColumn-1)*Rfe;
                SourceLines_thisPattern(iBlock,iColumn,:,:,:)=ACS_kykxc...
                    (targetPE_range+iBlock_offset, targetFE_range+iColumn_offset, :);
            end
        end
        SourceMatrix_thisPattern=reshape(permute(SourceLines_thisPattern,[3 4 1 2 5]),NtargetPE*NtargetFE,nBlocks*nColumns*Ncoil);
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
        GRAPPA_weights(iType, iCoil, :) = coefficient;
    end
end