function [wvec, wwvec] = ReadTransmissionFrequency(DataFilesDir)
% NOTE: Read phonon frequency points for transmission/reflection calculation 
  CurrDir = pwd;

  % -----------------------------------------------------------------
  % Read data file
  % -----------------------------------------------------------------
  cd(DataFilesDir);
  filename = 'Input_Frequency.agf';
  fid = fopen(filename,'r');

  if gt(fid,0)
      % data = load('Input_Frequency.txt'); % read file containing range of frequency values
      data = importdata(filename,' ',0); % read file containing range of frequency values  
  else
      error(sprintf('<!!!> %s not found.',filename));
  end
  fclose(fid);

  cd(CurrDir);

  % -----------------------------------------------------------------
  % Set input frequency and frequency-square variables with small imaginary part
  % -----------------------------------------------------------------
  wvec = data(:,1);  % this is the range of frequencies
  delta = data(:,2); % this is the small imaginary values needed for the calculations
  wwvec = wvec.*wvec.*(1.0+delta*1i); % this is the frequency square + small imaginary part

  clear data;

end
