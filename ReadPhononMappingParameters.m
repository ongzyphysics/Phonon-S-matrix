function [LeftPhonDisp, RightPhonDisp] = ReadPhononMappingParameters(DataFilesDir,LeftPhonDisp_in,RightPhonDisp_in)
% Read parameters for mapping of AGF lattice to primitive bulk lattice
  Angstrom = 1E-10;

  CurrDir = pwd;
  cd(DataFilesDir);

  LeftPhonDisp = LeftPhonDisp_in;
  RightPhonDisp = RightPhonDisp_in;

  LeftPhonDisp.nbasis_bulk = [];
  LeftPhonDisp.BasisAtomRxvec = [];
  LeftPhonDisp.BasisAtomRyvec = [];
  LeftPhonDisp.BasisAtomRzvec = [];

  RightPhonDisp.nbasis_bulk = [];
  RightPhonDisp.BasisAtomRxvec = [];
  RightPhonDisp.BasisAtomRyvec = [];
  RightPhonDisp.BasisAtomRzvec = [];

  % ---------------------------------------------------------------------  
  % Left phonons
  % ---------------------------------------------------------------------
  filename = 'Left_Phonon_Parameters.agf';

  fid = fopen(filename,'r');
  if gt(fid,0)
      fprintf(1,'          <%s> \n',filename);
    
      while 1
          InputLine = fgetl(fid); % read text line

          if ~ischar(InputLine) % end of file
              break;
          end

          CheckBeginLine = [InputLine repmat(' ',1,30)]; % padded input text line

          % === Read force constant matrix information ===
          if strcmp(CheckBeginLine(1:numel('BEGIN_AGF2BULK_LATTICE')),'BEGIN_AGF2BULK_LATTICE')
              while 1
                  InputLine = fgetl(fid); % read text line
                  CheckEndLine = [InputLine repmat(' ',1,30)];
  
                  if strcmp(CheckEndLine(1:numel('END_AGF2BULK_LATTICE')),'END_AGF2BULK_LATTICE')
                      break;
                  else 
                      InputText = textscan(InputLine,'%d %d %f %f %f ','CommentStyle','#');
                      nbasis_agf = InputText{1};  % AGF unit-cell basis atom index
                      nbasis_bulk = InputText{2}; % Bulk lattice primitive unit-cell basis atom index
                      rfrac1 = InputText{3};
                      rfrac2 = InputText{4};
                      rfrac3 = InputText{5};

                      Rvec = rfrac1*LeftPhonDisp.Rvec1 + rfrac2*LeftPhonDisp.Rvec2 + rfrac3*LeftPhonDisp.Rvec3;
                      LeftPhonDisp.nbasis_bulk(nbasis_agf) = nbasis_bulk;
                      LeftPhonDisp.BasisAtomRxvec(nbasis_agf) = Rvec(1);
                      LeftPhonDisp.BasisAtomRyvec(nbasis_agf) = Rvec(2);
                      LeftPhonDisp.BasisAtomRzvec(nbasis_agf) = Rvec(3);
                      % fprintf(1,'<DEBUG> %4d %4d %12.6g %12.6g %12.6g \n',...
                      %   nbasis_agf,nbasis_bulk,Rvec(1),Rvec(2),Rvec(3));
                  end
              end
          end 
      end
      fclose(fid);
  else
      error(sprintf('<!!!> %s not found.',filename));
  end

  % ---------------------------------------------------------------------  
  % Right phonons
  % ---------------------------------------------------------------------
  filename = 'Right_Phonon_Parameters.agf';

  fid = fopen(filename,'r');
  if gt(fid,0)
      fprintf(1,'          <%s> \n',filename);
    
      while 1
          InputLine = fgetl(fid); % read text line

          if ~ischar(InputLine) % end of file
              break;
          end

          CheckBeginLine = [InputLine repmat(' ',1,30)]; % padded input text line

          % === Read force constant matrix information ===
          if strcmp(CheckBeginLine(1:numel('BEGIN_AGF2BULK_LATTICE')),'BEGIN_AGF2BULK_LATTICE')
              while 1
                  InputLine = fgetl(fid); % read text line
                  CheckEndLine = [InputLine repmat(' ',1,30)];
  
                  if strcmp(CheckEndLine(1:numel('END_AGF2BULK_LATTICE')),'END_AGF2BULK_LATTICE')
                      break;
                  else 
                      InputText = textscan(InputLine,'%d %d %f %f %f ','CommentStyle','#');
                      nbasis_agf = InputText{1};  % AGF unit-cell basis atom index
                      nbasis_bulk = InputText{2}; % Bulk lattice primitive unit-cell basis atom index
                      rfrac1 = InputText{3};
                      rfrac2 = InputText{4};
                      rfrac3 = InputText{5};

                      Rvec = rfrac1*RightPhonDisp.Rvec1 + rfrac2*RightPhonDisp.Rvec2 + rfrac3*RightPhonDisp.Rvec3;
                      RightPhonDisp.nbasis_bulk(nbasis_agf) = nbasis_bulk;
                      RightPhonDisp.BasisAtomRxvec(nbasis_agf) = Rvec(1);
                      RightPhonDisp.BasisAtomRyvec(nbasis_agf) = Rvec(2);
                      RightPhonDisp.BasisAtomRzvec(nbasis_agf) = Rvec(3);
                      % fprintf(1,'<DEBUG> %4d %4d %12.6g %12.6g %12.6g \n', ...
                      %   nbasis_agf,nbasis_bulk,Rvec(1),Rvec(2),Rvec(3));
                  end
              end
          end 
      end
      fclose(fid);
  else
      error(sprintf('<!!!> %s not found.',filename));
  end

  cd(CurrDir);
end
