%%% inline script to load data and sort them, only for data
%%% from_Markus_20161102
% load_path='./';
function [cmrr_kxkyzsc,correspoding_single_slices,CAIPI_shift_of_slices]...
    =MBGC_load_dataAndParmeters_20161102(MB_factor,load_path,iFrameToLoad)
% load_path='G:\Mengye\MB_RAW_DATA\from_Markus_20161102\processed\oblique';
if nargin<3 || isempty (iFrameToLoad)
    iFrameToLoad=1;
end
switch MB_factor
    case 2
        if iFrameToLoad>1
            error('MB2 data has only 1 frame.')
        end
        load(fullfile(load_path,'Mb2Data_128kx'),'cmrr_kxkyzc_2x2','cmrr_2x2_sampling_mask_kys')
        cmrr_kxkyzsc=repmat(permute(double(cmrr_kxkyzc_2x2),[1 2 3 5 4]),[1 1 1 size(cmrr_2x2_sampling_mask_kys,2) 1]);
        nSlice=size(cmrr_kxkyzsc,3)*MB_factor;
        for iPolarity=1:size(cmrr_kxkyzsc,4)
            cmrr_kxkyzsc(:,cmrr_2x2_sampling_mask_kys(:,iPolarity)==0,:,iPolarity,:)=0;
        end
        metaZorder=[2:2:nSlice/2,1:2:nSlice/2];
        slice_order_temp=[metaZorder; nSlice/MB_factor+metaZorder ];
        correspoding_single_slices=num2cell(slice_order_temp,1);
        CAIPI_shift_of_slices=[0 0.25];
    case 4
        load(fullfile(load_path,'Mb4Data_128kx'),'cmrr_kxkyzc_2x4','cmrr_2x4_sampling_mask_kys')
        cmrr_kxkyzsc=repmat(permute(double(cmrr_kxkyzc_2x4),[1 2 3 5 4]),[1 1 1 size(cmrr_2x4_sampling_mask_kys,2) 1]);
        nSlice=size(cmrr_kxkyzsc,3)*MB_factor;
        for iPolarity=1:size(cmrr_kxkyzsc,4)
            cmrr_kxkyzsc(:,cmrr_2x4_sampling_mask_kys(:,iPolarity)==0,:,iPolarity,:)=0;
        end
        metaZorder=[2:2:nSlice/4,1:2:nSlice/4];
        slice_order_temp=[metaZorder; nSlice/MB_factor+metaZorder;  nSlice/MB_factor*2+metaZorder;  nSlice/MB_factor*3+metaZorder; ];
        correspoding_single_slices=num2cell(slice_order_temp,1);
        CAIPI_shift_of_slices=[0 0.25 0.5 0.75];
end
import LmyUtility.*
cmrr_kxkyzsc=LmyUtility.do2d(@(x)zpad(x,[128 128]),cmrr_kxkyzsc);
