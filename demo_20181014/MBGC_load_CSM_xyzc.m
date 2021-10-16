%%%inline script for getting CSM_xyzc, requiring CSM_source
function CSM_xyzc=MBGC_load_CSM_xyzc(CSM_option,load_path)
% % % CSM_option.source='segEPI LPC';
% % % CSM_option.source='segEPI SAKE';
% % % CSM_option.filter='ESPIRiT';
switch CSM_option.source
    case 'segEPI LPC'
        disp(['CSM source from LPC corrected EPI'])
        raw_data_filepath=fullfile(load_path,'SbData_128kx128ky_lpc.mat');
        raw_data_var_name='epi_kxkyzc_2shot_lpcCor';
        ESPIRiT_filepath=fullfile(load_path,'CSM_130y_xyzc_epi_lpcCor_bart.mat');
        ESPIRiT_var_name='CSM_xyzc_epi_lpcCor';
    case 'segEPI SAKE'
        disp(['CSM source from VC SAKE corrected EPI'])
        %         raw_data_filepath=fullfile(load_path,'SbData_128kx128ky_sake');
        %         raw_data_var_name='epi_kxkyzc_2shot_vcSAKE_48x48';
        %         raw_data_filepath=fullfile(load_path,'SbData_48kx48ky_sake_ksize3Thres4Iter50');
        if strfind(load_path,'oblique')
            raw_data_filepath=fullfile(load_path,'SbData_48kx48ky_sake_ksize3_ThreSlice5Every20ToEnd4.1 4.6 4.9 4.4 Iter50.mat');
            disp('oblique data')
        else
            raw_data_filepath=fullfile(load_path,'SbData_48kx48ky_sake_ksize3_Thres4ExceptS1T3_Iter50');
            disp('human data')
        end
        raw_data_var_name='epi_kxkyzc_2shot_sakeCor';
        ESPIRiT_filepath=fullfile(load_path,'CSM_130y_xyzc_epi_sakeCor_bart.mat');
        ESPIRiT_var_name='CSM_xyzc_epi_sakeCor';
    otherwise
        error('unknown CSM option')
end

disp('CSM calculated using ESPIRiT')
s=load(ESPIRiT_filepath,ESPIRiT_var_name);
CSM_xyzc=s.(ESPIRiT_var_name);% dynamic refering to the field
CSM_xyzc=crop(CSM_xyzc,[128 128 size(CSM_xyzc,3) size(CSM_xyzc,4)]);
