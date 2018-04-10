function [LeftPhonDisp, RightPhonDisp] = ReadPhononDispersionParameters(DataFilesDir)

CurrDir = pwd;
Angstrom = 1E-10;

cd(DataFilesDir);

% === Read left lead phonon parameters ===
clear Rvec1 Rvec2 Rvec3 Gvec1 Gvec2 Gvec3 PrimCell 

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

        % === Read primitive lattice vectors ===
        if strcmp(CheckBeginLine(1:numel('BEGIN_BZ_LATTICE')),'BEGIN_BZ_LATTICE')
            InputLine = fgetl(fid); % read text line
            InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
            Rvec1 = double([InputText{1} InputText{2} InputText{3}])*Angstrom;

            InputLine = fgetl(fid); % read text line
            InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
            Rvec2 = double([InputText{1} InputText{2} InputText{3}])*Angstrom;

            InputLine = fgetl(fid); % read text line
            InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
            Rvec3 = double([InputText{1} InputText{2} InputText{3}])*Angstrom;
        end

        % === Read force constant matrix information ===
        if strcmp(CheckBeginLine(1:numel('BEGIN_FC_MATRIX_INDEX')),'BEGIN_FC_MATRIX_INDEX')
            while 1
                InputLine = fgetl(fid); % read text line
                CheckEndLine = [InputLine repmat(' ',1,30)];

                if strcmp(CheckEndLine(1:numel('END_FC_MATRIX_INDEX')),'END_FC_MATRIX_INDEX')
                    break;
                else 
                    InputText = textscan(InputLine,'%d %f %f %f ','CommentStyle','#');
                    ncell = InputText{1};  % transverse cell index
                    rfrac1 = InputText{2};
                    rfrac2 = InputText{3};
                    rfrac3 = InputText{4};
                    PrimCell(ncell).rvec = rfrac1*Rvec1 ...
                      + rfrac2*Rvec2 + rfrac3*Rvec3;
                end
            end
        end 
    end
    fclose(fid);

    % === Read mass-normalized force constant matrices ===
    M = importdata('Left_Phonon_M.agf',' ',0);
    InvSqrtM = diag(1./sqrt(diag(M)));
    for ncell = 1:1:numel(PrimCell)
        fname = sprintf('Left_Phonon_K%d.agf',ncell);
        TempK = importdata(fname,' ',0);
        if eq(size(TempK,2),2*size(TempK,1))
            K = TempK(1:size(TempK,1),1:size(TempK,1)) ...
                + 1i*TempK(1:size(TempK,1),(1:size(TempK,1))+size(TempK,1));
        else
            K = TempK;
        end
        PrimCell(ncell).H = InvSqrtM*K*InvSqrtM; 
    end

    Gvec1 = cross(Rvec2,Rvec3)/(Rvec1*cross(Rvec2,Rvec3)');
    Gvec2 = cross(Rvec3,Rvec1)/(Rvec2*cross(Rvec3,Rvec1)');
    Gvec3 = cross(Rvec1,Rvec2)/(Rvec3*cross(Rvec1,Rvec2)'); 

    LeftPhonDisp.PrimCell = PrimCell;
    LeftPhonDisp.Rvec1 = Rvec1;
    LeftPhonDisp.Rvec2 = Rvec2;
    LeftPhonDisp.Rvec3 = Rvec3;
    LeftPhonDisp.Gvec1 = Gvec1;
    LeftPhonDisp.Gvec2 = Gvec2;
    LeftPhonDisp.Gvec3 = Gvec3;
else
    error(sprintf('<!!!> %s not found.',filename));
end

% === Read right lead phonon parameters ===
clear Rvec1 Rvec2 Rvec3 Gvec1 Gvec2 Gvec3 PrimCell 

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

        % === Read primitive lattice vectors ===
        if strcmp(CheckBeginLine(1:numel('BEGIN_BZ_LATTICE')),'BEGIN_BZ_LATTICE')
            InputLine = fgetl(fid); % read text line
            InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
            Rvec1 = double([InputText{1} InputText{2} InputText{3}])*Angstrom;

            InputLine = fgetl(fid); % read text line
            InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
            Rvec2 = double([InputText{1} InputText{2} InputText{3}])*Angstrom;

            InputLine = fgetl(fid); % read text line
            InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
            Rvec3 = double([InputText{1} InputText{2} InputText{3}])*Angstrom;
        end

        % === Read force constant matrix information ===
        if strcmp(CheckBeginLine(1:numel('BEGIN_FC_MATRIX_INDEX')),'BEGIN_FC_MATRIX_INDEX')
            while 1
                InputLine = fgetl(fid); % read text line
                CheckEndLine = [InputLine repmat(' ',1,30)];

                if strcmp(CheckEndLine(1:numel('END_FC_MATRIX_INDEX')),'END_FC_MATRIX_INDEX')
                    break;
                else 
                    InputText = textscan(InputLine,'%d %f %f %f ','CommentStyle','#');
                    ncell = InputText{1};  % transverse cell index
                    rfrac1 = InputText{2};
                    rfrac2 = InputText{3};
                    rfrac3 = InputText{4};
                    PrimCell(ncell).rvec = rfrac1*Rvec1 ...
                      + rfrac2*Rvec2 + rfrac3*Rvec3;
                end
            end
        end 
    end
    fclose(fid);

    % === Read mass-normalized force constant matrices ===
    M = importdata('Right_Phonon_M.agf',' ',0);
    InvSqrtM = diag(1./sqrt(diag(M)));
    for ncell = 1:1:numel(PrimCell)
        fname = sprintf('Right_Phonon_K%d.agf',ncell);
        % K = importdata(fname,' ',0);
        TempK = importdata(fname,' ',0);
        if eq(size(TempK,2),2*size(TempK,1))
            K = TempK(1:size(TempK,1),1:size(TempK,1)) ...
                + 1i*TempK(1:size(TempK,1),(1:size(TempK,1))+size(TempK,1));
        else
            K = TempK;
        end

        PrimCell(ncell).H = InvSqrtM*K*InvSqrtM; 
    end

    Gvec1 = cross(Rvec2,Rvec3)/(Rvec1*cross(Rvec2,Rvec3)');
    Gvec2 = cross(Rvec3,Rvec1)/(Rvec2*cross(Rvec3,Rvec1)');
    Gvec3 = cross(Rvec1,Rvec2)/(Rvec3*cross(Rvec1,Rvec2)'); 

    RightPhonDisp.PrimCell = PrimCell;
    RightPhonDisp.Rvec1 = Rvec1;
    RightPhonDisp.Rvec2 = Rvec2;
    RightPhonDisp.Rvec3 = Rvec3;
    RightPhonDisp.Gvec1 = Gvec1;
    RightPhonDisp.Gvec2 = Gvec2;
    RightPhonDisp.Gvec3 = Gvec3;
else
    error(sprintf('<!!!> %s not found.',filename));
end

cd(CurrDir); 

end
