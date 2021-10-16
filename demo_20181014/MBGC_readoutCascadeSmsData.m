function [pseudoPi_data_kxkyc]=MBGC_readoutCascadeSmsData(MB_kxkyc, MB_factor)
%%% This function convert a MB factor data set into a inplane parallel
%%% imaging problem by cascade in readout direction
[nFE,nPE,nCoil]=size(MB_kxkyc);
pseudoPi_data_kxkyc=zeros(nFE*MB_factor,nPE,nCoil);
pseudoPi_data_kxkyc(1:MB_factor:end,:,:)=MB_kxkyc;
end