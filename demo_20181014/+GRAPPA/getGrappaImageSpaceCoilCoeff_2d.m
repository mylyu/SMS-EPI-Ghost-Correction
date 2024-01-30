function GrappaUnmixingMap = getGrappaImageSpaceCoilCoeff_2d(header,GRAPPA_weights)
%%% convert GRAPPA weights into image space unmixing maps
Nblock=header.blocks;
Ncolumn=header.column;
Npe=header.Npe;
Nfe=header.Nfe;
% Rpe=size(GRAPPA_weights,1)+1;
Rpe=header.subsampling_factor(1);
Rfe=header.subsampling_factor(2);
Ncoil=size(GRAPPA_weights,2);

% new_weights_full=complex(zeros(Npe, Nfe, Ncoil,  Ncoil, Rpe*Rfe-1));
new_weights_full_sumPattern=complex(zeros(Npe, Nfe, Ncoil,  Ncoil,'single'));
center_ky=floor(Npe/2)+1;
center_kx=floor(Nfe/2)+1;
new_weights = reshape(GRAPPA_weights, [Rpe*Rfe-1, Ncoil,  Nblock, Ncolumn, Ncoil]);
new_weights = permute(new_weights, [1 3 4 2 5]);
% direct flip the kernel by index
ky2use_closest2Lastsampled=center_ky+1+Rpe*(Nblock/2-1):-Rpe:center_ky+1-Rpe*Nblock/2;
kx2use_closest2Lastsampled=center_kx+1+Rfe*(Ncolumn/2-1):-Rfe:center_kx+1-Rfe*Ncolumn/2;
% kx2use=center_kx+(Ncolumn-1)/2:-1:center_kx-(Ncolumn-1)/2;
%% totally R-1 types of patterns of different PE lines
for iTypes=1:Rpe*Rfe-1
    [x, y] = ind2sub([Rpe,Rfe],iTypes+1);
    iTypes_pe=x-1;iTypes_fe=y-1;
    shift_relative2firstType_pe=iTypes_pe-1;
    shift_relative2firstType_fe=iTypes_fe-1;
    ky2use=ky2use_closest2Lastsampled+shift_relative2firstType_pe;
    kx2use=kx2use_closest2Lastsampled+shift_relative2firstType_fe;
    new_weights_full_sumPattern(ky2use,kx2use,:,:)=new_weights_full_sumPattern(ky2use,kx2use,:,:)+...
        shiftdim(new_weights(iTypes,:,:,:,:));
end
% self to self kernel should be one in the center
for iCoil=1:Ncoil
%     new_weights_full_sumPattern(center_ky,center_kx,iCoil, iCoil)=1/(Rpe*Rfe-1);
        new_weights_full_sumPattern(center_ky,center_kx,iCoil, iCoil)=1;
end

GrappaUnmixingMap = ifft2c(new_weights_full_sumPattern)*sqrt(Npe*Nfe);
% deconvolution
% I_aliased=ifft2c(subsampled_kspace_data(:,:,:));
% I_coils = complex(zeros(Npe, Nfe, Ncoils));
% for ii = 1:Ncoils
%     I_coils(:,:,ii)=sum(I_aliased(:,:,:).*squeeze(w(:,:,ii,:)),3);
% end
% I_final=sos(I_coils);
return