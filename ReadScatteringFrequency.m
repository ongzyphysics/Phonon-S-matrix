function wvec_scatter = ReadScatteringFrequency(DataFilesDir)
% NOTE: Reads the frequency points for outputting channel scattering data
% NOTE: This is optional. If file is not found, then we assume no channel scattering data output is wanted

  CurrDir = pwd;

  cd(DataFilesDir);

  filename = 'Scattering_Output_Frequency.agf'; % data file name 
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
end
