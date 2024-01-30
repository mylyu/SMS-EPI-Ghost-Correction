function GrappaUnmixingMap = getGrappaImageSpaceCoilCoeff(header,GRAPPA_weights)
%%% convert GRAPPA weights into image space unmixing maps
Nblock=header.blocks;
Ncolumn=header.column;
Npe=header.Npe;
Nfe=header.Nfe;
R=size(GRAPPA_weights,1)+1;
Ncoil=size(GRAPPA_weights,2);

new_weights_full=complex(zeros(Npe,Nfe, Ncoil,  Ncoil,R-1));
center_ky=floor(Npe/2)+1;
center_kx=floor(Nfe/2)+1;
new_weights = reshape(GRAPPA_weights, [R-1, Ncoil,  Nblock, Ncolumn, Ncoil]);
new_weights = permute(new_weights, [1 3 4 2 5]);
% direct flip the kernel by index
ky2use_closest2Lastsampled=center_ky+1+R*(Nblock/2-1):-R:center_ky+1-R*Nblock/2;
kx2use=center_kx+(Ncolumn-1)/2:-1:center_kx-(Ncolumn-1)/2;
%% totally R-1 types of patterns of different PE lines
for iTypes=1:R-1
    shift_relative2firstType=iTypes-1;
    ky2use=ky2use_closest2Lastsampled+shift_relative2firstType;
    new_weights_full(ky2use,kx2use,:,:,iTypes)=new_weights(iTypes,:,:,:,:);
end
% self to self kernel should be one in the center
for iCoil=1:Ncoil
    new_weights_full(center_ky,center_kx,iCoil, iCoil,:)=1/(R-1);
end

GrappaUnmixingMap = ifft2c(sum(new_weights_full,5))*sqrt(Npe*Nfe);
% deconvolution
% I_aliased=ifft2c(subsampled_kspace_data(:,:,:));
% I_coils = complex(zeros(Npe, Nfe, Ncoils));
% for ii = 1:Ncoils
%     I_coils(:,:,ii)=sum(I_aliased(:,:,:).*squeeze(w(:,:,ii,:)),3);
% end
% I_final=sos(I_coils);
return