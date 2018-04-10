function SavePhononTransmissionData(DataFilesDir,wvec,LeftPhonon,RightPhonon,PhononData)

nwmax = length(wvec);

filename_t_L = 'Output_Left_Transmission.mat';
filename_t_R = 'Output_Right_Transmission.mat';
filename_tran = 'Output_Transmission.mat';

data_t_L = [];
data_t_R = [];
data_tran = [];

for nw = 1:1:nwmax
    w = wvec(nw);
    q_L_temp = -[LeftPhonon(nw).QR_adv];
    q_R_temp = [RightPhonon(nw).QR];
    t_L_temp = [PhononData(nw).tt_L];
    t_R_temp = [PhononData(nw).tt_R];
    n_L_temp = [LeftPhonon(nw).N_sub];
    n_R_temp = [RightPhonon(nw).N_sub];
    v_L_temp = [LeftPhonon(nw).VecVR_adv];
    v_R_temp = [RightPhonon(nw).VecVR];
    p_L_temp = [LeftPhonon(nw).TR_adv];
    p_R_temp = [RightPhonon(nw).TR];

    t_L_temp = t_L_temp(gt(abs(q_L_temp),0));
    t_R_temp = t_R_temp(gt(abs(q_R_temp),0));
    n_L_temp = n_L_temp(gt(abs(q_L_temp),0));
    n_R_temp = n_R_temp(gt(abs(q_R_temp),0));   
    v_L_temp = v_L_temp(gt(abs(q_L_temp),0));
    v_R_temp = v_R_temp(gt(abs(q_R_temp),0));   
    p_L_temp = p_L_temp(gt(abs(q_L_temp),0));
    p_R_temp = p_R_temp(gt(abs(q_R_temp),0));   

    q_L_temp = q_L_temp(gt(abs(q_L_temp),0));
    q_R_temp = q_R_temp(gt(abs(q_R_temp),0));

    if gt(length(q_L_temp),0)
        for nq = 1:length(q_L_temp)
            data_t_L = [data_t_L; w q_L_temp(nq) t_L_temp(nq) n_L_temp(nq) v_L_temp(nq) p_L_temp(nq)];
        %     fprintf(fid_t_L,'%16.8e %16.8e %16.8e %5d %16.8e \n ', w, q_L_temp(nq), t_L_temp(nq), n_L_temp(nq), v_L_temp(nq));                          
        end        
    end

    if gt(length(q_R_temp),0)
        for nq = 1:length(q_R_temp)
            data_t_R = [data_t_R; w q_R_temp(nq) t_R_temp(nq) n_R_temp(nq) v_R_temp(nq) p_R_temp(nq)];
        %     fprintf(fid_t_R,'%16.8e %16.8e %16.8e %5d %16.8e \n ', w, q_R_temp(nq), t_R_temp(nq), n_R_temp(nq), v_R_temp(nq));            
        end        
    end
    
    Xi_L = real(LeftPhonon(nw).Xi_negf);
    Xi_R = real(RightPhonon(nw).Xi_negf);
    Xi_C = real(PhononData(nw).Xi_negf);
    
    % fprintf(fid_tran,'%16.8e %16.8e %16.8e %16.8e \n ', w, Xi_L, Xi_R, Xi_C);                
    data_tran = [data_tran; w Xi_L Xi_R Xi_C];
end

cd(DataFilesDir);
save(filename_t_L,'data_t_L');
save(filename_t_R,'data_t_R');
save(filename_tran,'data_tran');
cd('..');

fprintf(1,'\t  <%s> \n', filename_t_L);
fprintf(1,'\t  <%s> \n', filename_t_R);
fprintf(1,'\t  <%s> \n', filename_tran);








