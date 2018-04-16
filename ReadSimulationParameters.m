function [LeftLeadParam, RightLeadParam, CenterParam] = ReadSimulationParameters(DataFilesDir)
% NOTE:
% NOTE:

  angstrom = 1.0E-10;

  CurrDir = pwd;
  cd(DataFilesDir);

  % -------------------------------------------------------------------
  % Read parameters for left lead (KL, KC, KR & M matrices and parameters)
  % -------------------------------------------------------------------

  % === Left KL matrix ===
  filename = 'Left_KL.agf';
  fid = fopen(filename,'r');
  if gt(fid,0)
      fprintf(1,'\t  <%s> ',filename); % name of text file containing KL matrix
      fclose(fid);

      TempMat = importdata(filename,' ',0);
      if eq(size(TempMat,2),2*size(TempMat,1)) % i.e. if there are imaginary components
          LeftLeadParam.MatKL = ...
              TempMat(1:size(TempMat,1),1:size(TempMat,1)) ...
              + 1i*TempMat(1:size(TempMat,1),(1:size(TempMat,1))+size(TempMat,1));
      else % i.e. if there are no imaginary components
          LeftLeadParam.MatKL = TempMat;
      end    
  else
      error(sprintf('<!!!> %s not found.',filename));
  end

  % === Left KC matrix ===
  filename = 'Left_KC.agf';
  fid = fopen(filename,'r');
  if gt(fid,0)
      fprintf(1,'   <%s> ',filename); % name of text file containing KC matrix
      fclose(fid);

      TempMat = importdata(filename,' ',0);
      if eq(size(TempMat,2),2*size(TempMat,1)) % i.e. if there are imaginary components
          LeftLeadParam.MatKC = ...
              TempMat(1:size(TempMat,1),1:size(TempMat,1)) ...
              + 1i*TempMat(1:size(TempMat,1),(1:size(TempMat,1))+size(TempMat,1));
      else % i.e. if there are no imaginary components
          LeftLeadParam.MatKC = TempMat;
      end    
  else
      error(sprintf('<!!!> %s not found.',filename));
  end

  % === Left KR matrix ===
  filename = 'Left_KR.agf';
  fid = fopen(filename,'r');
  if gt(fid,0)
      fprintf(1,'   <%s> ',filename); % name of text file containing KR matrix
      fclose(fid);

      TempMat = importdata(filename,' ',0);
      if eq(size(TempMat,2),2*size(TempMat,1)) % i.e. if there are imaginary components
          LeftLeadParam.MatKR = ...
              TempMat(1:size(TempMat,1),1:size(TempMat,1)) ...
              + 1i*TempMat(1:size(TempMat,1),(1:size(TempMat,1))+size(TempMat,1));
      else % i.e. if there are no imaginary components
          LeftLeadParam.MatKR = TempMat;
      end    
  else
      error(sprintf('<!!!> %s not found.',filename));
  end

  % === Left M matrix ===
  filename = 'Left_M.agf';
  fid = fopen(filename,'r'); % name of text file containing M matrix
  if gt(fid,0)
      fprintf(1,'   <%s> ',filename);
      fclose(fid);

      TempMat = importdata(filename,' ',0);
      if eq(size(TempMat,2),2*size(TempMat,1)) % i.e. if there are imaginary components
          LeftLeadParam.MatM = ...
              TempMat(1:size(TempMat,1),1:size(TempMat,1)) ...
              + 1i*TempMat(1:size(TempMat,1),(1:size(TempMat,1))+size(TempMat,1));
      else % i.e. if there are no imaginary components
          LeftLeadParam.MatM = TempMat;
      end    
  else
      error(sprintf('<!!!> %s not found.',filename));
  end

  % === Left lead parameters ===
  filename = 'Left_Parameters.agf';
  fid = fopen(filename,'r');
  if gt(fid,0)
      fprintf(1,'   <%s> \n',filename);
    
      while 1
          InputLine = fgetl(fid); % read text line

          if ~ischar(InputLine) % end of file
              break;
          end

          CheckBeginLine = [InputLine repmat(' ',1,30)]; % padded input text line

          % === Read lattice (non-primitive) vectors ===
          if strcmp(CheckBeginLine(1:numel('BEGIN_LATTICE')),'BEGIN_LATTICE')
              InputLine = fgetl(fid); % read text line
              InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
              LeftLeadParam.rvec_long = double([InputText{1} InputText{2} InputText{3}])*angstrom;
              LeftLeadParam.a_long = norm(LeftLeadParam.rvec_long);

              InputLine = fgetl(fid); % read text line
              InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
              LeftLeadParam.rvec_tran1 = double([InputText{1} InputText{2} InputText{3}])*angstrom;

              InputLine = fgetl(fid); % read text line
              InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
              LeftLeadParam.rvec_tran2 = double([InputText{1} InputText{2} InputText{3}])*angstrom;
          end

          % === Read tranverse span ===
          if strcmp(CheckBeginLine(1:numel('BEGIN_TRANSVERSE_SPAN')),'BEGIN_TRANSVERSE_SPAN')
              InputLine = fgetl(fid); % read text line
              InputText = textscan(InputLine,'%d %d %d ','CommentStyle','#');
              LeftLeadParam.n_tran1 = InputText{2};
              LeftLeadParam.n_tran2 = InputText{3};
          end

          % === Read tranverse cell information ===
          if strcmp(CheckBeginLine(1:numel('BEGIN_CELL_INDEX')),'BEGIN_CELL_INDEX')
              while 1
                  InputLine = fgetl(fid); % read text line
                  CheckEndLine = [InputLine repmat(' ',1,30)];
  
                  if strcmp(CheckEndLine(1:numel('END_CELL_INDEX')),'END_CELL_INDEX')
                      break;
                  else 
                      InputText = textscan(InputLine,'%d %d %f %f %f ','CommentStyle','#');
                      ncell = InputText{1};  % transverse cell index
                      natoms = InputText{2}; % number of atoms in transverse cell
                      LeftLeadParam.TransCell(ncell).nrange = (1:(3*natoms)) + (ncell-1)*3*natoms;
                      frac_long = InputText{3};
                      frac_tran1 = InputText{4};
                      frac_tran2 = InputText{5};
                      LeftLeadParam.TransCell(ncell).rvec = frac_long*LeftLeadParam.rvec_long ...
                        + frac_tran1*LeftLeadParam.rvec_tran1 + frac_tran2*LeftLeadParam.rvec_tran2;
                  end
              end
          end 
      end
      fclose(fid);
  else
      error(sprintf('<!!!> %s not found.',filename));
  end


  % -------------------------------------------------------------------
  % Read parameters for right lead (KL, KC, KR & M matrices and parameters)
  % -------------------------------------------------------------------

  % === Right KL matrix ===
  filename = 'Right_KL.agf';
  fid = fopen(filename,'r');
  if gt(fid,0)
      fprintf(1,'\t  <%s> ',filename); % name of text file containing KL matrix
      fclose(fid);

      TempMat = importdata(filename,' ',0);
      if eq(size(TempMat,2),2*size(TempMat,1)) % i.e. if there are imaginary components
          RightLeadParam.MatKL = ...
              TempMat(1:size(TempMat,1),1:size(TempMat,1)) ...
              + 1i*TempMat(1:size(TempMat,1),(1:size(TempMat,1))+size(TempMat,1));
      else % i.e. if there are no imaginary components
          RightLeadParam.MatKL = TempMat;
      end    
  else
      error(sprintf('<!!!> %s not found.',filename));
  end

  % === Right KC matrix ===
  filename = 'Right_KC.agf';
  fid = fopen(filename,'r');
  if gt(fid,0)
      fprintf(1,'   <%s>',filename); % name of text file containing KC matrix
      fclose(fid);

      TempMat = importdata(filename,' ',0);
      if eq(size(TempMat,2),2*size(TempMat,1)) % i.e. if there are imaginary components
          RightLeadParam.MatKC = ...
              TempMat(1:size(TempMat,1),1:size(TempMat,1)) ...
              + 1i*TempMat(1:size(TempMat,1),(1:size(TempMat,1))+size(TempMat,1));
      else % i.e. if there are no imaginary components
          RightLeadParam.MatKC = TempMat;
      end    
  else
      error(sprintf('<!!!> %s not found.',filename));
  end

  % === Right KR matrix ===
  filename = 'Right_KR.agf';
  fid = fopen(filename,'r');
  if gt(fid,0)
      fprintf(1,'   <%s>',filename); % name of text file containing KR matrix
      fclose(fid);

      TempMat = importdata(filename,' ',0);
      if eq(size(TempMat,2),2*size(TempMat,1)) % i.e. if there are imaginary components
          RightLeadParam.MatKR = ...
              TempMat(1:size(TempMat,1),1:size(TempMat,1)) ...
              + 1i*TempMat(1:size(TempMat,1),(1:size(TempMat,1))+size(TempMat,1));
      else % i.e. if there are no imaginary components
          RightLeadParam.MatKR = TempMat;
      end    
  else
      error(sprintf('<!!!> %s not found.',filename));
  end

  % === Right M matrix ===
  filename = 'Right_M.agf';
  fid = fopen(filename,'r');
  if gt(fid,0)
      fprintf(1,'   <%s>',filename); % name of text file containing M matrix
      fclose(fid);

      TempMat = importdata(filename,' ',0);
      if eq(size(TempMat,2),2*size(TempMat,1)) % i.e. if there are imaginary components
          RightLeadParam.MatM = ...
              TempMat(1:size(TempMat,1),1:size(TempMat,1)) ...
              + 1i*TempMat(1:size(TempMat,1),(1:size(TempMat,1))+size(TempMat,1));
      else % i.e. if there are no imaginary components
          RightLeadParam.MatM = TempMat;
      end    
  else
      error(sprintf('<!!!> %s not found.',filename));
  end

  % === Right lead parameters ===
  filename = 'Right_Parameters.agf';
  fid = fopen(filename,'r');
  if gt(fid,0)
      fprintf(1,'   <%s> \n',filename);
    
      while 1
          InputLine = fgetl(fid); % read text line
 
          if ~ischar(InputLine) % end of file
              break;
          end

          CheckBeginLine = [InputLine repmat(' ',1,30)]; % padded input text line

          % === Read lattice (non-primitive) vectors ===
          if strcmp(CheckBeginLine(1:numel('BEGIN_LATTICE')),'BEGIN_LATTICE')
              InputLine = fgetl(fid); % read text line
              InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
              RightLeadParam.rvec_long = double([InputText{1} InputText{2} InputText{3}])*angstrom;
              RightLeadParam.a_long = norm(RightLeadParam.rvec_long);

              InputLine = fgetl(fid); % read text line
              InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
              RightLeadParam.rvec_tran1 = double([InputText{1} InputText{2} InputText{3}])*angstrom;

              InputLine = fgetl(fid); % read text line
              InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
              RightLeadParam.rvec_tran2 = double([InputText{1} InputText{2} InputText{3}])*angstrom;
          end

          % === Read tranverse span ===
          if strcmp(CheckBeginLine(1:numel('BEGIN_TRANSVERSE_SPAN')),'BEGIN_TRANSVERSE_SPAN')
              InputLine = fgetl(fid); % read text line
              InputText = textscan(InputLine,'%d %d %d ','CommentStyle','#');
              RightLeadParam.n_tran1 = InputText{2};
              RightLeadParam.n_tran2 = InputText{3};
          end


          % === Read tranverse cell information ===
          if strcmp(CheckBeginLine(1:numel('BEGIN_CELL_INDEX')),'BEGIN_CELL_INDEX')
              while 1
                  InputLine = fgetl(fid); % read text line
                  CheckEndLine = [InputLine repmat(' ',1,30)];

                  if strcmp(CheckEndLine(1:numel('END_CELL_INDEX')),'END_CELL_INDEX')
                      break;
                  else 
                      InputText = textscan(InputLine,'%d %d %f %f %f ','CommentStyle','#');
                      ncell = InputText{1};  % transverse cell index
                      natoms = InputText{2}; % number of atoms in transverse cell
                      RightLeadParam.TransCell(ncell).nrange = (1:(3*natoms)) + (ncell-1)*3*natoms;
                      frac_long = InputText{3};
                      frac_tran1 = InputText{4};
                      frac_tran2 = InputText{5};
                      RightLeadParam.TransCell(ncell).rvec = frac_long*RightLeadParam.rvec_long ...
                        + frac_tran1*RightLeadParam.rvec_tran1 + frac_tran2*RightLeadParam.rvec_tran2;
                  end
              end
          end 
      end
      fclose(fid);
  else
      error(sprintf('<!!!> %s not found.',filename));
  end

  % -------------------------------------------------------------------
  % Read parameters for channel/scattering region 
  % (KL, KC, KR and M submatrices)  
  %
  %     To facilitate the use of the recursive Green's function (RCF) 
  %     algorithm, we partition, if possible, the scattering region 
  %     into successive layers with their own KL, KC, KR and M 
  %     submatrices. This partitioning  converts the giant matrices 
  %     associated with the scattering region into block-tridiagonal 
  %     matrices such that the diagonal and off-diagonal submatrices 
  %     correspond to the aforementioned KL, KC, KR and M matrices. The 
  %     partitioning also reduce the costs of the matrix inversion 
  %     operations.
  %
  % -------------------------------------------------------------------

  nlayer = 1;
  filename = sprintf('Center_M%d.agf',nlayer);
  fid = fopen(filename);

  if gt(fid,0)   
      while( gt(fid,0) )
          % === M submatrix for layer 'nlayer' ===
          M_filename = sprintf('Center_M%d.agf',nlayer);

          TempMat = importdata(M_filename,' ',0);
          if eq(size(TempMat,2),2*size(TempMat,1))
              Lyr(nlayer).MatM = ...
                  TempMat(1:size(TempMat,1),1:size(TempMat,1)) ...
                  + 1i*TempMat(1:size(TempMat,1),(1:size(TempMat,1))+size(TempMat,1));
          else
              Lyr(nlayer).MatM = TempMat;
          end    

          % fprintf(1,'\t  <%s>', KL_filename);
          % fprintf(1,'   <%s>', KC_filename);
          % fprintf(1,'   <%s>', KR_filename);
          % fprintf(1,'\t  <%s> \n', M_filename);            
        
          fclose(fid);        
          nlayer = nlayer + 1;
        	
          filename = sprintf('Center_M%d.agf',nlayer);
          fid = fopen(filename);
      end
  end

  % -------------------------------------------------------------------

  nlayer = 1;
  filename = sprintf('Center_KC%d.agf',nlayer);
  fid = fopen(filename);

  if gt(fid,0)   
      while( gt(fid,0) )
          if eq(nlayer,1)
              nrows_L = size(LeftLeadParam.MatM,1);
          else
              nrows_L = size(Lyr(nlayer-1).MatM,1);
          end

          nrows_C = size(Lyr(nlayer).MatM,1);

          if eq(nlayer,numel(Lyr))
              nrows_R = size(RightLeadParam.MatM,1);
          else
              nrows_R = size(Lyr(nlayer+1).MatM,1);
          end

          % === KL submatrix for layer 'nlayer' ===
          KL_filename = sprintf('Center_KL%d.agf',nlayer);

          TempMat = importdata(KL_filename,' ',0);
          if eq(size(TempMat,2),2*nrows_L)
              Lyr(nlayer).MatKL = ...
                  TempMat(1:nrows_C,1:nrows_L) ...
                  + 1i*TempMat(1:nrows_C,(1:nrows_L)+nrows_L);
          else
              Lyr(nlayer).MatKL = TempMat;
          end    
          %{
          if eq(size(TempMat,2),2*size(TempMat,1))
              Lyr(nlayer).MatKL = ...
                  TempMat(1:size(TempMat,1),1:size(TempMat,1)) ...
                  + 1i*TempMat(1:size(TempMat,1),(1:size(TempMat,1))+size(TempMat,1));
          else
              Lyr(nlayer).MatKL = TempMat;
          end    
          %}
          % === KC submatrix for layer 'nlayer' ===
          KC_filename = sprintf('Center_KC%d.agf',nlayer);
 
          TempMat = importdata(KC_filename,' ',0);
          if eq(size(TempMat,2),2*nrows_C)
              Lyr(nlayer).MatKC = ...
                  TempMat(1:nrows_C,1:nrows_C) ...
                  + 1i*TempMat(1:nrows_C,(1:nrows_C)+nrows_C);
          else
              Lyr(nlayer).MatKC = TempMat;
          end    
          %{
          if eq(size(TempMat,2),2*size(TempMat,1))
              Lyr(nlayer).MatKC = ...
                  TempMat(1:size(TempMat,1),1:size(TempMat,1)) ...
                  + 1i*TempMat(1:size(TempMat,1),(1:size(TempMat,1))+size(TempMat,1));
          else
              Lyr(nlayer).MatKC = TempMat;
          end    
          %}
          % === KR submatrix for layer 'nlayer' ===
          KR_filename = sprintf('Center_KR%d.agf',nlayer);

          TempMat = importdata(KR_filename,' ',0);
          if eq(size(TempMat,2),2*nrows_R)
              Lyr(nlayer).MatKR = ...
                  TempMat(1:nrows_C,1:nrows_R) ...
                  + 1i*TempMat(1:nrows_C,(1:nrows_R)+nrows_R);
          else
              Lyr(nlayer).MatKR = TempMat;
          end    
          %{
          if eq(size(TempMat,2),2*size(TempMat,1))
              Lyr(nlayer).MatKR = ...
                  TempMat(1:size(TempMat,1),1:size(TempMat,1)) ...
                  + 1i*TempMat(1:size(TempMat,1),(1:size(TempMat,1))+size(TempMat,1));
          else
              Lyr(nlayer).MatKR = TempMat;
          end    
          %}

          % === M submatrix for layer 'nlayer' ===
          M_filename = sprintf('Center_M%d.agf',nlayer);

          TempMat = importdata(M_filename,' ',0);
          if eq(size(TempMat,2),2*size(TempMat,1))
              Lyr(nlayer).MatM = ...
                  TempMat(1:size(TempMat,1),1:size(TempMat,1)) ...
                  + 1i*TempMat(1:size(TempMat,1),(1:size(TempMat,1))+size(TempMat,1));
          else
              Lyr(nlayer).MatM = TempMat;
          end    

          fprintf(1,'\t  <%s>', KL_filename);
          fprintf(1,'   <%s>', KC_filename);
          fprintf(1,'   <%s>', KR_filename);
          fprintf(1,'   <%s> \n', M_filename);            
        
          fclose(fid);        
          nlayer = nlayer + 1;
        
          filename = sprintf('Center_KC%d.agf',nlayer);
          fid = fopen(filename);
      end
      % nlayer = nlayer - 1;
  end

  for nlayer = 1:numel(Lyr) % DEBUG
      % disp(nlayer)
      % disp(size(Lyr(nlayer).MatKL));
      % disp(size(Lyr(nlayer).MatKC));
      % disp(size(Lyr(nlayer).MatKR));
  end

  % -------------------------------------------------------------------
  % Padding of additional layers to the left and right edges of the 
  % scattering region. 
  % -------------------------------------------------------------------

  nlayer= numel(Lyr); 

  for n = (nlayer+1):-1:2 % move layer 1 to 2, 2 to 3, ..., nlayer to nlayer+1
      Lyr(n) = Lyr(n-1);
  end 

  % Set layer 1 (left edge) as layer from left lead
  Lyr(1).MatM  = LeftLeadParam.MatM;
  Lyr(1).MatKL = LeftLeadParam.MatKL;
  Lyr(1).MatKC = LeftLeadParam.MatKC;
  Lyr(1).MatKR = (Lyr(2).MatKL)';

  % Set layer nlayer+2 (righte edge) as layer from right lead
  Lyr(nlayer+2).MatM  = RightLeadParam.MatM;
  Lyr(nlayer+2).MatKL = (Lyr(nlayer+1).MatKR)';
  Lyr(nlayer+2).MatKC = RightLeadParam.MatKC;
  Lyr(nlayer+2).MatKR = RightLeadParam.MatKR;

  nlayer= numel(Lyr); 

  CenterParam.nlayers = length(Lyr);
  CenterParam.Lyr = Lyr;
  CenterParam.MatKCL = Lyr(1).MatKL;
  CenterParam.MatKLC = (CenterParam.MatKCL)';
  CenterParam.MatKCR = Lyr(nlayer).MatKR;
  CenterParam.MatKRC = (CenterParam.MatKCR)';

  % -------------------------------------------------------------------
  % Check acoustic sum rule
  % -------------------------------------------------------------------

  % LeftLeadParam.MatKC = LeftLeadParam.MatKC - ...
  %     diag(sum([LeftLeadParam.MatKL LeftLeadParam.MatKC LeftLeadParam.MatKR],2));
  err_check = sum(abs(sum([LeftLeadParam.MatKL LeftLeadParam.MatKC LeftLeadParam.MatKR],2)));
  if gt(err_check,0)
      warning('++ Possible significant acoustic sum rule violation (Left lead bulk).');
  end

  % RightLeadParam.MatKC = RightLeadParam.MatKC - ...
  %     diag(sum([RightLeadParam.MatKL RightLeadParam.MatKC RightLeadParam.MatKR],2));
  err_check = sum(abs(sum([RightLeadParam.MatKL RightLeadParam.MatKC RightLeadParam.MatKR],2)));
  if gt(err_check,0)
      warning('++ Possible significant acoustic sum rule violation (Right lead bulk).');
  end

  err_check = sum(abs(sum([LeftLeadParam.MatKL LeftLeadParam.MatKC CenterParam.MatKLC],2)));
  if gt(err_check,0)
      warning('++ Possible significant acoustic sum rule violation (Left lead to channel).');
  end

  err_check = sum(abs(sum([CenterParam.MatKRC RightLeadParam.MatKC RightLeadParam.MatKR],2)));
  if gt(err_check,0)
      warning('++ Possible significant acoustic sum rule violation. (Right lead to channel).');
  end

  for nlayer = 1:length(CenterParam.Lyr)
      % CenterParam.Lyr(nlayer).MatKC = CenterParam.Lyr(nlayer).MatKC - ...
      %     diag(sum([CenterParam.Lyr(nlayer).MatKL CenterParam.Lyr(nlayer).MatKC CenterParam.Lyr(nlayer).MatKR],2));    
      err_check = sum(abs(sum([CenterParam.Lyr(nlayer).MatKL CenterParam.Lyr(nlayer).MatKC CenterParam.Lyr(nlayer).MatKR],2)));    
      if gt(err_check,0)
          warning('++ Possible significant acoustic sum rule violation (Channel layer %d).',nlayer);
      end    
  end

  fprintf(1,'\t  Acoustic Sum Rule check done. \n');

  % -------------------------------------------------------------------
  % Check matrix symmetry
  % -------------------------------------------------------------------

  err_check = sum(sum(abs(LeftLeadParam.MatKL - LeftLeadParam.MatKR'),2));
  if gt(err_check,0)
      warning('++ Possible significant matrix symmetry violation (Left lead bulk).');
  end

  err_check = sum(sum(abs(LeftLeadParam.MatKC - LeftLeadParam.MatKC'),2));
  if gt(err_check,0)
      warning('++ Possible significant matrix symmetry violation (Left lead bulk).');
  end

  err_check = sum(sum(abs(RightLeadParam.MatKL - RightLeadParam.MatKR'),2));
  if gt(err_check,0)
      warning('++ Possible significant matrix symmetry violation (Right lead bulk).');
  end

  err_check = sum(sum(abs(RightLeadParam.MatKC - RightLeadParam.MatKC'),2));
  if gt(err_check,0)
      warning('++ Possible significant matrix symmetry violation (Right lead bulk).');
  end

  for nlayer = 2:(length(CenterParam.Lyr)-1)
      err_check = sum(sum(abs(CenterParam.Lyr(nlayer).MatKC-CenterParam.Lyr(nlayer).MatKC'),2));    
      if gt(err_check,0)
          warning('++ Possible significant matrix symmetry violation (Channel layer %d).',nlayer);
      end    

      err_check = sum(sum(abs(CenterParam.Lyr(nlayer).MatKL-CenterParam.Lyr(nlayer-1).MatKR'),2));    
      if gt(err_check,0)
          warning('++ Possible significant matrix symmetry violation (Channel layer %d).',nlayer);
      end    
    
      % disp(size(CenterParam.Lyr(nlayer).MatKR));
      % disp(size(CenterParam.Lyr(nlayer+1).MatKL'));
      err_check = sum(sum(abs(CenterParam.Lyr(nlayer).MatKR-CenterParam.Lyr(nlayer+1).MatKL'),2));    
      if gt(err_check,0)
          warning('++ Possible significant matrix symmetry violation (Channel layer %d).',nlayer);
      end        
  end

  for nlayer = [1 length(CenterParam.Lyr)]
      err_check = sum(sum(abs(CenterParam.Lyr(nlayer).MatKC-CenterParam.Lyr(nlayer).MatKC'),2));    
      if gt(err_check,0)
        warning('++ Possible significant matrix symmetry violation (Channel layer %d).',nlayer);
      end    
  end

  fprintf(1,'\t  Matrix Symmetry check done. \n');

  % -------------------------------------------------------------------
  % Check size of layers in leads and channel
  % -------------------------------------------------------------------


  % -------------------------------------------------------------------
  % Check inter-layer coupling matrices between leads and channel
  % -------------------------------------------------------------------

  % if and(eq(size(LeftLeadParam.MatKC,1),size(CenterParam.Lyr(1).MatKC,1)),eq(size(LeftLeadParam.MatKC,2),size(CenterParam.Lyr(1).MatKC,2)))
  if eq(sum(abs(size(LeftLeadParam.MatKC) - size(CenterParam.Lyr(1).MatKC))),0)
      err_check = sum(sum(abs(LeftLeadParam.MatKL - CenterParam.Lyr(1).MatKL),2));
      if gt(err_check,0)
          warning('++ Possible significant inter-layer coupling matrix violation (Left lead to left contact).');
      end
  else
      warning('++ Left lead and contact matrix size mismatch.');     
  end

  if eq(sum(abs(size(RightLeadParam.MatKC) - size(CenterParam.Lyr(length(Lyr)).MatKC))),0)
      err_check = sum(sum(abs(RightLeadParam.MatKR - CenterParam.Lyr(length(Lyr)).MatKR),2));
      if gt(err_check,0)
          warning('++ Possible significant inter-layer coupling matrix violation (Right lead to right contact).');
      end
  else
      warning('++ Right lead and contact matrix size mismatch.');         
  end

  %{
  K_temp = CenterParam.MatKLC( 1:length(LeftLeadParam.MatKC), 1:length(LeftLeadParam.MatKC) );
  err_check = sum(sum(abs(LeftLeadParam.MatKR - K_temp),2));
  if gt(err_check,0)
      % fprintf(1,'\t ++ Possible significant inter-layer coupling matrix violation (Left lead to channel).\n');
      warning('++ Possible significant inter-layer coupling matrix violation (Left lead to channel).');     
  end

  n2 = size(CenterParam.MatKRC,2);
  K_temp = CenterParam.MatKRC( 1:length(RightLeadParam.MatKC), (1:length(RightLeadParam.MatKC))+n2-length(RightLeadParam.MatKC) );
  err_check = sum(sum(abs(RightLeadParam.MatKL - K_temp),2));
  if gt(err_check,0)
      % fprintf(1,'\t ++ Possible significant inter-layer coupling matrix violation (Right lead to channel).\n');
      warning('++ Possible significant inter-layer coupling matrix violation (Right lead to channel).');
  end
  %}

  fprintf(1,'\t  Coupling matrices check done. \n');

  cd(CurrDir);
end
