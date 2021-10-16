function [epi_kxkyzc_lpcCor]=oneDimLinearCorr_parameter(epi_kxkyzc_raw,nShot,phasepara_zp)
%%% This function corrects Nyquist ghosts with 1d linear phase model from 
%%% input parameters "phasepara_zp", whose dimesions should be
%%% slice-paramters.
import LmyGhostCorrection.*
org_size=size(epi_kxkyzc_raw);
if ndims(epi_kxkyzc_raw)<4
    epi_kxkyzc_raw=reshape(epi_kxkyzc_raw,size(epi_kxkyzc_raw,1),size(epi_kxkyzc_raw,2),1,size(epi_kxkyzc_raw,3));
end
epi_kxkyzc_lpcCor=zeros(size(epi_kxkyzc_raw));
nSlice=size(epi_kxkyzc_lpcCor,3);
for iSlice = 1:nSlice
    Kyxc=fftshift(fft(permute(squeeze(epi_kxkyzc_raw(:,:,iSlice,:)),[2 1 3]),[],2),2);
    CorxKy_en = OIPhaseCor(Kyxc,phasepara_zp(iSlice,:),nShot);
    epi_kxkyzc_lpcCor(:,:,iSlice,:) = permute(ifft(fftshift(CorxKy_en,2),[],2),[2 1 3]);
end
epi_kxkyzc_lpcCor=reshape(epi_kxkyzc_lpcCor,org_size);
