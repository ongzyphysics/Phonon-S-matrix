function WritePhononDispersionOutput(DataFilesDir,Left,Right)

filename_phon_L = 'Output_Left_Dispersion.dat';
filename_phon_R = 'Output_Right_Dispersion.dat';

cd(DataFilesDir);
fid_ph_L = fopen(filename_phon_L,'w');
fid_ph_R = fopen(filename_phon_R,'w');
cd('..');

% ----------------------------------------------------------------
% Left phonon dispersion
% ----------------------------------------------------------------

[kphon, wphon, vphon] = GetPhononDispersion(Left);
kphon = kphon';
wphon = wphon';
nkmax = numel(kphon);

for nk = 1:1:nkmax
    fprintf(fid_ph_L,'%14.6e ', [kphon(nk,1) wphon(nk,:)]);    
    fprintf(fid_ph_L,'\n');
end

% ----------------------------------------------------------------
% Right phonon dispersion
% ----------------------------------------------------------------

[kphon, wphon, vphon] = GetPhononDispersion(Right);
kphon = kphon';
wphon = wphon';
nkmax = numel(kphon);

for nk = 1:1:nkmax
    fprintf(fid_ph_R,'%14.6e ', [kphon(nk,1) wphon(nk,:)]);    
    fprintf(fid_ph_R,'\n');
end

% ----------------------------------------------------------------

fprintf(1,'\t  <%s> \n', filename_phon_L);
fprintf(1,'\t  <%s> \n', filename_phon_R);
fclose(fid_ph_L);
fclose(fid_ph_R);
