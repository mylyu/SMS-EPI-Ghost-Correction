function restoredSlice_xyz=MBGC_readoutCascadeRestoreReconSlices(readoutCascadeRecon_xy,CAIPI_pattern)
% restore slice to original locations out of immediate readout-cascaded SMS
% recon.
nPE=size(readoutCascadeRecon_xy,2);
MB_factor=length(CAIPI_pattern);
nFE=size(readoutCascadeRecon_xy,1)/MB_factor;
restoredSlice_xyz=zeros(nFE,nPE,MB_factor);
readoutCascadeRecon_xy=circshift(readoutCascadeRecon_xy,+0.5*nFE*(mod(MB_factor,2)-1));
readoutCascadeRecon_xy=reshape(readoutCascadeRecon_xy,nFE,MB_factor,nPE);
for iSlice=1:MB_factor
    restoredSlice_xyz(:,:,iSlice)=circshift(readoutCascadeRecon_xy(:,iSlice,:),[0 0 -1*CAIPI_pattern(iSlice)*nPE]);
end
end