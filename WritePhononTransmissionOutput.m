function WritePhononTransmissionOutput(DataFilesDir,wvec,LeftPhonon,RightPhonon,PhononData)

CurrDir = pwd;

nwmax = numel(wvec);

filename_t_L = 'Output_Left_Transmission.dat';
filename_r_L = 'Output_Left_Reflection.dat';
filename_t_R = 'Output_Right_Transmission.dat';
filename_tran = 'Output_Transmission.dat';

cd(DataFilesDir);
fid_t_L = fopen(filename_t_L,'w');
fid_r_L = fopen(filename_r_L,'w');
fid_t_R = fopen(filename_t_R,'w');
fid_tran = fopen(filename_tran,'w');
cd(CurrDir);

for nw = 1:1:nwmax
    w = wvec(nw);

    qx_L_temp = [LeftPhonon(nw).VecQx_minus_adv]; % wave vector (x-direction) of left lead phonons  
    qy_L_temp = [LeftPhonon(nw).VecQy_minus_adv]; % wave vector (y-direction) of left lead phonons  
    qz_L_temp = [LeftPhonon(nw).VecQz_minus_adv]; % wave vector (z-direction) of left lead phonons  
    b_L_temp =  [LeftPhonon(nw).VecB_minus_adv];  % band index of left lead phonons  
    velx_L_temp = [LeftPhonon(nw).VecVelx_minus_adv]; % group velocity (x-direction) of left lead phonons  
    vely_L_temp = [LeftPhonon(nw).VecVely_minus_adv]; % group velocity (y-direction) of left lead phonons  
    velz_L_temp = [LeftPhonon(nw).VecVelz_minus_adv]; % group velocity (z-direction) of left lead phonons  

    qx_R_temp = [RightPhonon(nw).VecQx_plus];     % wave vector (x-direction) of right lead phonons
    qy_R_temp = [RightPhonon(nw).VecQy_plus];     % wave vector (y-direction) of right lead phonons
    qz_R_temp = [RightPhonon(nw).VecQz_plus];     % wave vector (z-direction) of right lead phonons
    b_R_temp =  [RightPhonon(nw).VecB_plus];      % band index of right lead phonons  
    velx_R_temp = [RightPhonon(nw).VecVelx_plus];     % group velocity (x-direction) of right lead phonons
    vely_R_temp = [RightPhonon(nw).VecVely_plus];     % group velocity (y-direction) of right lead phonons
    velz_R_temp = [RightPhonon(nw).VecVelz_plus];     % group velocity (z-direction) of right lead phonons

    q_L_temp = [LeftPhonon(nw).VecQ_minus_adv]; % wave vector (long.) of left lead phonons (b/w 0 and 2 pi)  
    q_R_temp = [RightPhonon(nw).VecQ_plus];     % wave vector (long.) of right lead phonons (b/w 0 and 2 pi)
    t_L_temp = [PhononData(nw).tt_L_in];           % transmission coeffs. of left lead phonons (b/w 0 and 1)
    t_R_temp = [PhononData(nw).tt_L_out];           % transmission coeffs. of right lead phonons (b/w 0 and 1)
    q_tran1_L_temp = [LeftPhonon(nw).VecQ_tran1];  % wave vector (transv. 1) of left lead phonons (b/w 0 and 2 pi)
    q_tran2_L_temp = [LeftPhonon(nw).VecQ_tran2];  % wave vector (transv. 2) of left lead phonons (b/w 0 and 2 pi)
    q_tran1_R_temp = [RightPhonon(nw).VecQ_tran1]; % wave vector (transv. 1) of right lead phonons (b/w 0 and 2 pi)
    q_tran2_R_temp = [RightPhonon(nw).VecQ_tran2]; % wave vector (transv. 2) of right lead phonons (b/w 0 and 2 pi)
    v_L_temp = [LeftPhonon(nw).VecV_minus_adv]; % group velocity of left lead phonons
    v_R_temp = [RightPhonon(nw).VecV_plus];     % group velocity of right lead phonons
    p_L_temp = [LeftPhonon(nw).VecT_minus_adv]; % maximum transmission of propagating left lead phonons
    p_R_temp = [RightPhonon(nw).VecT_plus];     % maximum transmission of propagating right lead phonons

    nrange_L_temp = gt(abs(b_L_temp),0);
    nrange_R_temp = gt(abs(b_R_temp),0);

    qx_L_temp = qx_L_temp(nrange_L_temp);
    qy_L_temp = qy_L_temp(nrange_L_temp);
    qz_L_temp = qz_L_temp(nrange_L_temp);
    b_L_temp  = b_L_temp(nrange_L_temp);
    velx_L_temp = velx_L_temp(nrange_L_temp);
    vely_L_temp = vely_L_temp(nrange_L_temp);
    velz_L_temp = velz_L_temp(nrange_L_temp);

    qx_R_temp = qx_R_temp(nrange_R_temp);
    qy_R_temp = qy_R_temp(nrange_R_temp);
    qz_R_temp = qz_R_temp(nrange_R_temp);
    b_R_temp  = b_R_temp(nrange_R_temp);
    velx_R_temp = velx_R_temp(nrange_R_temp);
    vely_R_temp = vely_R_temp(nrange_R_temp);
    velz_R_temp = velz_R_temp(nrange_R_temp);

    t_L_temp = t_L_temp(nrange_L_temp);
    t_R_temp = t_R_temp(nrange_R_temp);
    q_tran1_L_temp = q_tran1_L_temp(nrange_L_temp);
    q_tran2_L_temp = q_tran2_L_temp(nrange_L_temp);
    q_tran1_R_temp = q_tran1_R_temp(nrange_R_temp);
    q_tran2_R_temp = q_tran2_R_temp(nrange_R_temp);   
    v_L_temp = v_L_temp(nrange_L_temp);
    v_R_temp = v_R_temp(nrange_R_temp);   
    p_L_temp = p_L_temp(nrange_L_temp);
    p_R_temp = p_R_temp(nrange_R_temp);   

    q_L_temp = q_L_temp(nrange_L_temp);
    q_R_temp = q_R_temp(nrange_R_temp);


    if gt(numel(b_R_temp),0)
        for nq = 1:numel(b_R_temp)
            % fprintf(fid_t_R,'%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n ', ...
            %         w, q_R_temp(nq), q_tran1_R_temp(nq), q_tran2_R_temp(nq), ... 
            %         v_R_temp(nq), p_R_temp(nq), t_R_temp(nq));            
            fprintf(fid_t_R,'%14.6e %4d %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \n', ...
                    w, b_R_temp(nq), qx_R_temp(nq), qy_R_temp(nq), qz_R_temp(nq), v_R_temp(nq), ... 
                    p_R_temp(nq), t_R_temp(nq), velx_R_temp(nq), vely_R_temp(nq), velz_R_temp(nq));       
        end
    end

    if gt(numel(b_L_temp),0)
        for nq = 1:numel(b_L_temp)
            % fprintf(fid_t_L,'%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n ', ...
            %         w, q_L_temp(nq), q_tran1_L_temp(nq), q_tran2_L_temp(nq), ... 
            %         v_L_temp(nq), p_L_temp(nq), t_L_temp(nq));            
            fprintf(fid_t_L,'%14.6e %4d %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \n', ...
                    w, b_L_temp(nq), qx_L_temp(nq), qy_L_temp(nq), qz_L_temp(nq), v_L_temp(nq), ... 
                    p_L_temp(nq), t_L_temp(nq), velx_L_temp(nq), vely_L_temp(nq), velz_L_temp(nq));            
        end
    end
    
    qxr_L_temp = [LeftPhonon(nw).VecQx_minus]; % wave vector (x-direction) of reflected left lead phonons  
    qyr_L_temp = [LeftPhonon(nw).VecQy_minus]; % wave vector (y-direction) of reflected left lead phonons  
    qzr_L_temp = [LeftPhonon(nw).VecQz_minus]; % wave vector (z-direction) of reflected left lead phonons  
    br_L_temp  = [LeftPhonon(nw).VecB_minus];  % band index of left lead phonons  
    velxr_L_temp = [LeftPhonon(nw).VecVelx_minus]; % group (x-direction) of reflected left lead phonons  
    velyr_L_temp = [LeftPhonon(nw).VecVely_minus]; % group (y-direction) of reflected left lead phonons  
    velzr_L_temp = [LeftPhonon(nw).VecVelz_minus]; % group (z-direction) of reflected left lead phonons  

    qr_L_temp = [LeftPhonon(nw).VecQ_minus]; % wave vector (long.) of reflected left lead phonons (b/w 0 and 2 pi)  
    tr_L_temp = [PhononData(nw).rr_L_out];        % reflection coeffs. of left lead phonons (b/w 0 and 1)
    qr_tran1_L_temp = [LeftPhonon(nw).VecQ_tran1];  % wave vector (transv. 1) of left lead phonons (b/w 0 and 2 pi)
    qr_tran2_L_temp = [LeftPhonon(nw).VecQ_tran2];  % wave vector (transv. 2) of left lead phonons (b/w 0 and 2 pi)
    vr_L_temp = [LeftPhonon(nw).VecV_minus]; % group velocity of reflected left lead phonons
    pr_L_temp = [LeftPhonon(nw).VecT_minus]; % maximum transmission of propagating reflected left lead phonons

    nrrange_L_temp = gt(abs(br_L_temp),0);

    qxr_L_temp = qxr_L_temp(nrrange_L_temp);
    qyr_L_temp = qyr_L_temp(nrrange_L_temp);
    qzr_L_temp = qzr_L_temp(nrrange_L_temp);
    br_L_temp  = br_L_temp(nrrange_L_temp);
    velxr_L_temp = velxr_L_temp(nrrange_L_temp);
    velyr_L_temp = velyr_L_temp(nrrange_L_temp);
    velzr_L_temp = velzr_L_temp(nrrange_L_temp);

    tr_L_temp = tr_L_temp(nrrange_L_temp);
    qr_tran1_L_temp = qr_tran1_L_temp(nrrange_L_temp);
    qr_tran2_L_temp = qr_tran2_L_temp(nrrange_L_temp);
    vr_L_temp = vr_L_temp(nrrange_L_temp);
    pr_L_temp = pr_L_temp(nrrange_L_temp);

    qr_L_temp = qr_L_temp(nrrange_L_temp);

    if gt(numel(br_L_temp),0)
        for nq = 1:numel(br_L_temp)
            % fprintf(fid_r_L,'%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n ', ...
            %         w, qr_L_temp(nq), qr_tran1_L_temp(nq), qr_tran2_L_temp(nq), ... 
            %         vr_L_temp(nq), pr_L_temp(nq), tr_L_temp(nq));            
            fprintf(fid_r_L,'%14.6e %4d %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \n', ...
                    w, br_L_temp(nq), qxr_L_temp(nq), qyr_L_temp(nq), qzr_L_temp(nq), vr_L_temp(nq), ... 
                    pr_L_temp(nq), tr_L_temp(nq), velxr_L_temp(nq), velyr_L_temp(nq), velzr_L_temp(nq));            
        end
    end

    Xi_L = real(LeftPhonon(nw).Xi_negf);
    Xi_R = real(RightPhonon(nw).Xi_negf);
    Xi_C_negf = real(PhononData(nw).Xi_negf);
    Xi_C_mode = real(PhononData(nw).Xi_mode);
    
    fprintf(fid_tran,'%16.8e %16.8e %16.8e %16.8e %16.8e \n', w, Xi_L, Xi_R, Xi_C_negf, Xi_C_mode);                
end

fprintf(1,'\t  <%s> \n', filename_t_L);
fprintf(1,'\t  <%s> \n', filename_r_L);
fprintf(1,'\t  <%s> \n', filename_t_R);
fprintf(1,'\t  <%s> \n', filename_tran);

fclose(fid_t_L);
fclose(fid_r_L);
fclose(fid_t_R);
fclose(fid_tran);







