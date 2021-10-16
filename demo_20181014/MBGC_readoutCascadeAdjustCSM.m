function adjustedCSM_xyzc=MBGC_readoutCascadeAdjustCSM(CSM_xyzc, CAIPI_pattern)
%%% This function convert the CSM into a inplane parallel
%%% imaging problem by cascade in readout direction
[nFE,nPE,MB_factor,nCoil]=size(CSM_xyzc);
adjustedCSM_xyzc=zeros(nFE,MB_factor,nPE,nCoil);
for iSlice=1:MB_factor
    adjustedCSM_xyzc(:,iSlice,:,:)=circshift(CSM_xyzc(:,:,iSlice,:),[0 CAIPI_pattern(iSlice)*nPE]);
end
adjustedCSM_xyzc=reshape(adjustedCSM_xyzc,nFE*MB_factor,nPE,nCoil);
adjustedCSM_xyzc=circshift(adjustedCSM_xyzc,-nFE*0.5*(mod(MB_factor,2)-1));% shift if even MB factor
end