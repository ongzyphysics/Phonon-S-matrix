function WritePhononLDOSOutput(DataFilesDir,wvec,PhononData,LDOS_output_flag)
% Plots the local phonon density of states in the scattering region

cd(DataFilesDir);

nwmax = length(wvec);

for nw = 1:1:nwmax
    LDOS(nw).w   = wvec(nw);
    LDOS(nw).rho = PhononData(nw).rho;
end

rho = [LDOS(:).rho]'; % density of states

if eq(LDOS_output_flag,1) 
    % output in matlab format
    filename_LDOS = 'Output_LDOS.mat';
    save(filename_LDOS,'wvec','rho');
    fprintf(1,'\t  <%s> \n', filename_LDOS);
else
    % output in text format 
    filename_LDOS = 'Output_LDOS.txt';
    fid_LDOS = fopen(filename_LDOS,'w');

    fprintf(fid_LDOS,'%11.4e', wvec);
    fprintf(fid_LDOS,'\n');
    for ns = 1:size(rho,2)
        fprintf(fid_LDOS,'%11.4e', rho(:,ns));
        fprintf(fid_LDOS,'\n');
    end

    fclose(fid_LDOS);
    fprintf(1,'\t  <%s> \n', filename_LDOS);
end

cd('..');

