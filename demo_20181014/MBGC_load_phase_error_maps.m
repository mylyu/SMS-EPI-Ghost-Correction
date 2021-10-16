function phase_error_to_use=MBGC_load_phase_error_maps(PEM_option,load_path,PEM_size)
nFE=PEM_size(1);nPE=PEM_size(2);
switch PEM_option.model
    case '1dLpc'
        %%  Pruessmann's group; generate 2D phase error map that is equilavent to linear phase error
        load(fullfile(load_path,['SbData_' num2str(nFE) 'kx' num2str(nPE) 'ky_lpc']),'phasepara_evenOddMajar_zp')
        %comupte moving average of phase error from adjacent slices, as better intial estimation
        phasepara_allslices_ordered= phasepara_evenOddMajar_zp;
        nSlice=size(phasepara_allslices_ordered,1);
        phase_error_to_use=zeros(nFE,nPE,nSlice);
        for iSlice=1:nSlice
            ConPhase = phasepara_allslices_ordered(iSlice,1);
            LinPhase = phasepara_allslices_ordered(iSlice,2);
            PhaseAlongX = ConPhase+(1:nFE)*LinPhase;
            %             if strfind(load_path,'oblique')
            %                 Phasemap = repmat(2*PhaseAlongX(:),[1 nPE]);if iSlice==1,disp('temp Siemens Phantom type of 1dc'),end
            %             else
            Phasemap = flip(repmat(2*PhaseAlongX(:),[1 nPE]),1);if iSlice==1,disp('Bruker/Siemens brain type of 1d LPC'),end
            %             end
            phase_error_to_use(:,:,iSlice)= exp(1i*Phasemap);
        end
        return % done with this function
    case '2dFreeBART'
        load(fullfile(load_path,'SbData_PosNeg_Recon'),'recon_xyz_neg_RepSum','recon_xyz_pos_RepSum');
        vc=cat(4,recon_xyz_neg_RepSum,recon_xyz_pos_RepSum);
        vc_csm=zeros(size(vc));
        if strcmp(PEM_option.model,'2dFreeBART')
            disp('BART ESPIRIT filtered PEM')
            %              vc= MBS_addnoise(vc,15);
            %              warning('adding noise to pem')
            parfor iSlice=1:size(vc,3)
                %                     [temp] = squeeze(bart('ecalib -c0 -r 80 -k 9',fft2c(vc(:,:,iSlice,:)) ));% for oblique phantom
                [temp] = squeeze(bart('ecalib -c0 -r 48',fft2c(vc(:,:,iSlice,:)) ));
                vc_csm(:,:,iSlice,:)=temp(:,:,:,1);
            end
        end
        PEM_xyz_NegFromPos_filtered=vc_csm(:,:,:,1).*conj(vc_csm(:,:,:,2));
        PEM_xyz_NegFromPos_filtered=PEM_xyz_NegFromPos_filtered./abs(PEM_xyz_NegFromPos_filtered);
        phase_error_to_use=crop(PEM_xyz_NegFromPos_filtered,[PEM_size size(PEM_xyz_NegFromPos_filtered,3)]);
        return
end
