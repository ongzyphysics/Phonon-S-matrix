function wvec_scatter = ReadScatteringFrequency(DataFilesDir)

CurrDir = pwd;

cd(DataFilesDir);

filename = 'Scattering_Output_Frequency.agf';
fid = fopen(filename,'r');

if gt(fid,0)
    data = importdata(filename,' ',0); % read file containing range of frequency values  
    wvec_scatter = data(:,1);
else
    warning(sprintf('<!!!> %s not found.',filename));
    wvec_scatter = 0;
end
fclose(fid);



cd(CurrDir);
