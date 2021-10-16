function [epi_kxkyzc_cpcCor,phasepara_const]=oneDimConstCorr_entropy(epi_kxkyzc_raw,nShot)
%%% This function corrects inter-shot constant phase variation by entropy
%%% minimization. It packs the old functions by Victor Xie to be more
%%% easily usable.
import LmyGhostCorrection.*
org_size=size(epi_kxkyzc_raw);
if ndims(epi_kxkyzc_raw)==3 % assume single-slice
epi_kxkyzc_raw=reshape(epi_kxkyzc_raw,size(epi_kxkyzc_raw,1),size(epi_kxkyzc_raw,2),1,size(epi_kxkyzc_raw,3));
end
epi_kxkyzc_cpcCor=zeros(size(epi_kxkyzc_raw));
phasepara_const=zeros(size(epi_kxkyzc_cpcCor,3),nShot);
nSlice=size(epi_kxkyzc_cpcCor,3);
middleSliceIndex=ceil(nSlice/2);
for iSlice = [middleSliceIndex:nSlice,middleSliceIndex-1:-1:1]
    disp(iSlice)
    Kyxc=fftshift(fft(permute(squeeze(epi_kxkyzc_raw(:,:,iSlice,:)),[2 1 3]),[],2),2);
    if iSlice==middleSliceIndex
        [CorxKy_en, phasepara_const(iSlice,:)] = OCEntropyBasedCor_forCompile(Kyxc,nShot,[]);
    elseif iSlice>middleSliceIndex
        [CorxKy_en, phasepara_const(iSlice,:)] = OCEntropyBasedCor_forCompile(Kyxc,nShot,phasepara_const(iSlice-1,:));
    else
        [CorxKy_en, phasepara_const(iSlice,:)] = OCEntropyBasedCor_forCompile(Kyxc,nShot,phasepara_const(iSlice+1,:));
    end
    epi_kxkyzc_cpcCor(:,:,iSlice,:) = permute(ifft(fftshift(CorxKy_en,2),[],2),[2 1 3]);
end
epi_kxkyzc_cpcCor=reshape(epi_kxkyzc_cpcCor,org_size);