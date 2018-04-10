function ExtendedAtomisticGreensFunctionTransmission(InputFilesDir,OutputFilesDir)
% NOTE: InputFilesDir is the folder containing the input *.agf files
% NOTE: InputFilesDir is the folder containing the output *.dat files 

CurrDir = pwd; % current directory

tic;

% -----------------------------------------------------------------
% Read input frequency data points
% -----------------------------------------------------------------
[wvec, wwvec] = ReadTransmissionFrequency(InputFilesDir);
fprintf(1,'++ Input frequency range file <Input_Frequency.agf> read in %f seconds. \n\n',toc);
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% Read scattering results frequency data points
% -----------------------------------------------------------------
wvec_scatter = ReadScatteringFrequency(InputFilesDir);
fprintf(1,'++ Scattering frequency range file <Scattering_Output_Frequency.agf> read in %f seconds. \n\n',toc);
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% Read force constant matrices, mass matrices and lattice constants 
% -----------------------------------------------------------------
tic;
fprintf(1,'++ Input matrix files: \n');
[LeftLeadParam, RightLeadParam, CenterParam] = ReadSimulationParameters(InputFilesDir); 
fprintf(1,'   read in %f seconds. \n\n',toc);
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% Read phonon dispersion parameters and matrices 
% -----------------------------------------------------------------
tic;
fprintf(1,'++ Input phonon dispersion files: \n');
[LeftPhonDisp, RightPhonDisp] = ReadPhononDispersionParameters(InputFilesDir); 
[LeftPhonDisp, RightPhonDisp] = ReadPhononMappingParameters(InputFilesDir,LeftPhonDisp,RightPhonDisp);
fprintf(1,'   read in %f seconds. \n\n',toc);

% -----------------------------------------------------------------
% Convert force constants and atomic masses to H matrices (mass normalized matrices)
% -----------------------------------------------------------------
tic;
[Left, Center, Right] = ConvertMassNormalizedMatrices(LeftLeadParam,CenterParam,RightLeadParam);
fprintf(1,'++ Conversion of parameter took %f seconds. \n\n',toc);
% -----------------------------------------------------------------

for nw = 1:numel(wvec)
    % disp(nw);
    fprintf(1,'COMPUTATION: %4d of %4d \n\n',nw,length(wvec));
    % -----------------------------------------------------------------
    % Calculate surface green's functions for left lead
    % -----------------------------------------------------------------
    tic;
    LeftPhonon(nw) = ComputeSurfaceGreensFunction(wvec(nw),wwvec(nw),Left);
    LeftPhonon(nw) = MapBulkPhononModes(wvec(nw),LeftPhonon(nw),Left,LeftPhonDisp);
    LeftPhonon(nw) = RebuildExtendedModeEigenvector(wvec(nw),LeftPhonon(nw),Left,LeftPhonDisp);
    fprintf(1,'++ Left surface Greens function computation took %f seconds. \n\n',toc);
    % -----------------------------------------------------------------

    % -----------------------------------------------------------------
    % Calculate surface green's functions for right lead
    % -----------------------------------------------------------------
    tic;
    RightPhonon(nw) = ComputeSurfaceGreensFunction(wvec(nw),wwvec(nw),Right);
    RightPhonon(nw) = MapBulkPhononModes(wvec(nw),RightPhonon(nw),Right,RightPhonDisp);
    RightPhonon(nw) = RebuildExtendedModeEigenvector(wvec(nw),RightPhonon(nw),Right,RightPhonDisp);
    fprintf(1,'++ Right surface Greens function computation took %f seconds. \n\n',toc);
    % -----------------------------------------------------------------
    
    % -----------------------------------------------------------------
    % Compute frequency-dependent transmission
    % -----------------------------------------------------------------
    tic;
    PhononData(nw) = ComputeCenterPhononTransmission(wvec(nw),wwvec(nw),...
                     LeftPhonon(nw),RightPhonon(nw),Center,Left,Right);
    fprintf(1,'++ Phonon transmission computation took %f seconds. \n\n',toc);
    % -----------------------------------------------------------------

    % -----------------------------------------------------------------
    % Remove data used for frequency-dependent transmission to save memory 
    % -----------------------------------------------------------------
    LeftPhonon(nw).MatSurfGL = [];
    LeftPhonon(nw).MatSurfGR = [];
    LeftPhonon(nw).MatBulkG  = [];
    LeftPhonon(nw).U_plus      = [];
    LeftPhonon(nw).U_plus_adv  = [];
    LeftPhonon(nw).U_minus     = [];
    LeftPhonon(nw).U_minus_adv = [];
    LeftPhonon(nw).V_plus      = [];
    LeftPhonon(nw).V_plus_adv  = [];
    LeftPhonon(nw).V_minus     = [];
    LeftPhonon(nw).V_minus_adv = [];

    RightPhonon(nw).MatSurfGL = [];
    RightPhonon(nw).MatSurfGR = [];
    RightPhonon(nw).MatBulkG  = [];
    RightPhonon(nw).U_plus      = [];
    RightPhonon(nw).U_plus_adv  = [];
    RightPhonon(nw).U_minus     = [];
    RightPhonon(nw).U_minus_adv = [];
    RightPhonon(nw).V_plus      = [];
    RightPhonon(nw).V_plus_adv  = [];
    RightPhonon(nw).V_minus     = [];
    RightPhonon(nw).V_minus_adv = [];
end

% save ('TempData.mat','PhononData','LeftPhonon','RightPhonon','Left','Right','LeftPhonDisp','RightPhonDisp','wvec');

LeftPhonon = rmfield(LeftPhonon,{'MatSurfGL','MatSurfGR','MatBulkG', ...      
             'U_plus','U_plus_adv','U_minus','U_minus_adv', ... 
             'V_plus','V_plus_adv','V_minus','V_minus_adv'});
RightPhonon = rmfield(RightPhonon,{'MatSurfGL','MatSurfGR','MatBulkG', ...      
             'U_plus','U_plus_adv','U_minus','U_minus_adv', ... 
             'V_plus','V_plus_adv','V_minus','V_minus_adv'});

% -----------------------------------------------------------------
% Postprocessing transmission data for left and right lead mode
% -----------------------------------------------------------------
fprintf(1,'POSTPROCESSING DATA: \n\n',nw,length(wvec));

% -----------------------------------------------------------------
% Write text files containing overall and individual phonon mode transmission results
% Write text files containing phonon dispersion of left and right leads
% -----------------------------------------------------------------
tic;
fprintf(1,'++ Output files: \n');
WritePhononTransmissionOutput(OutputFilesDir,wvec,LeftPhonon,RightPhonon,PhononData);
WritePhononDispersionOutput(OutputFilesDir,Left,Right);
WritePhononScatteringOutput(OutputFilesDir,wvec_scatter,wvec,LeftPhonon,RightPhonon,PhononData);

fprintf(1,'   saved in %f seconds. \n\n',toc);
% -----------------------------------------------------------------
disp(' ');

return; % DEBUG LINE

end

