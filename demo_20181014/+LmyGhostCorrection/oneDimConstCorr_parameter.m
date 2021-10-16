function [epi_kxkyzc_cpcCor]=oneDimConstCorr_parameter(epi_kxkyzc_raw,nShot,phasepara_const_zp)
%%% This function corrects inter-shot constant phase variation from input
%%% paramters. The paramters should have slices as its 1st dimenions.
import LmyGhostCorrection.*
org_size=size(epi_kxkyzc_raw);
if ndims(epi_kxkyzc_raw)==3 % assume single-slice
epi_kxkyzc_raw=reshape(epi_kxkyzc_raw,size(epi_kxkyzc_raw,1),size(epi_kxkyzc_raw,2),1,size(epi_kxkyzc_raw,3));
end
epi_kxkyzc_cpcCor=zeros(size(epi_kxkyzc_raw));
nSlice=size(epi_kxkyzc_cpcCor,3);
for iSlice = 1:nSlice
    Kyxc=fftshift(fft(permute(squeeze(epi_kxkyzc_raw(:,:,iSlice,:)),[2 1 3]),[],2),2);
    for iShot=2:Nshot
        Kyxc(iShot:nShot:end,:,:) = Kyxc(iShot:nShot:end,:,:)*exp(1i*phasepara_const_zp(iShot-1));
    end
    epi_kxkyzc_cpcCor(:,:,iSlice,:) = permute(ifft(fftshift(Kyxc,2),[],2),[2 1 3]);
end
epi_kxkyzc_cpcCor=reshape(epi_kxkyzc_cpcCor,org_size);
