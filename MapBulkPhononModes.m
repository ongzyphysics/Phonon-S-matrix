function LeadPhonon = MapBulkPhononModes(w_in,LeadPhonon_in,Lead,LeadPhonDisp)
% This function maps the AGF modes to their bulk lattice phonon modes 
% and determines their primitive BZ wave vector and group velocity
% w_in = input frequency
% LeadPhonon_in = data structure storing lead phonon eigenvectors, velocities, wave vectors, etc
% Lead = data structure storing force constant matrices and transverse cell info for lead
% LeadPhonDisp = data structure for bulk lattice phonon dispersion

  w = w_in;
  LeadPhonon = LeadPhonon_in;

  % return; % DEBUG
  epsilon_cutoff_w = 1E-4; % matching cutoff parameter for frequency
  epsilon_cutoff_q = 1E-4; % matching cutoff parameter for wave vector
  % nwmax = numel(wvec);

  % if and(eq(Lead.n_tran1,1),eq(Lead.n_tran2,1))
  %   fprintf(1,'<DEBUG> System is quasi-1D \n');
  % end

  % ------------------------------------------------------------------------------
  % Find unfolded wave vectors 
  % ------------------------------------------------------------------------------

  Gvec1 = Lead.gvec_long;  % reciprocal lattice vector (AGF, non-primitive)
  Gvec2 = Lead.gvec_tran1; % reciprocal lattice vector (AGF, non-primitive)
  Gvec3 = Lead.gvec_tran2; % reciprocal lattice vector (AGF, non-primitive)
  PrimGvec1 = LeadPhonDisp.Gvec1; % reciprocal lattice vector (Brillouin Zone, primitive)
  PrimGvec2 = LeadPhonDisp.Gvec2; % reciprocal lattice vector (Brillouin Zone, primitive)
  PrimGvec3 = LeadPhonDisp.Gvec3; % reciprocal lattice vector (Brillouin Zone, primitive)
  
  PrimCell = LeadPhonDisp.PrimCell; % contains parameters for computing primitive phonon dispersion 
  
  ScaleFactor = abs(PrimGvec1*cross(PrimGvec2,PrimGvec3)')/abs(Gvec1*cross(Gvec2,Gvec3)');

  % ------------------------------------------------------------------------------
  % ------------------------------------------------------------------------------
  % === transverse wave vector(s) ===
  qfrac2_temp = [LeadPhonon.VecQ_tran1]/(2*pi); % fractional wave vector (transv. 1) (b/w 0 and 1)
  qfrac2_temp = mod(qfrac2_temp+0.5,1.0)-0.5;   % set value between -0.5 and 0.4999...
  qfrac3_temp = [LeadPhonon.VecQ_tran2]/(2*pi); % fractional wave vector (transv. 2) (b/w 0 and 1)
  qfrac3_temp = mod(qfrac3_temp+0.5,1.0)-0.5;   % set value between -0.5 and 0.4999...

  % === longitudinal wave vector ===
  vlong_temp  = [LeadPhonon.VecV_plus]; % AGF phonon longitudinal group velocity 
  qfrac1_temp = [LeadPhonon.VecQ_plus]/(2*pi); % fractional wave vector (long.) (b/w 0 and 1)  
  qfrac1_temp = mod(qfrac1_temp+0.5,1.0)-0.5;  % set value between -0.5 and 0.4999...
  p_temp      = round([LeadPhonon.VecT_plus]); % maximum transmission of propagating phonons

  [qfrac1_temp,qfrac2_temp,qfrac3_temp] = ...
    RedistributeDegenerateModes(qfrac1_temp,qfrac2_temp,qfrac3_temp,p_temp,epsilon_cutoff_q);
    % redistribute degenerate modes especially at BZ edge 

  nqmax = numel(p_temp); % number of lead channels (transmitting and nontransmitting) at frequency w
  qxvec = zeros(size(p_temp)); % allocate data for 'unfolded' wave vector x-component
  qyvec = zeros(size(p_temp)); % allocate data for 'unfolded' wave vector y-component
  qzvec = zeros(size(p_temp)); % allocate data for 'unfolded' wave vector z-component
  vxvec = zeros(size(p_temp)); % allocate data for 'unfolded' velocity x-component
  vyvec = zeros(size(p_temp)); % allocate data for 'unfolded' velocity y-component
  vzvec = zeros(size(p_temp)); % allocate data for 'unfolded' velocity z-component
  dwvec = -ones(size(p_temp)); % discrepancy between frequency of folded and unfolded k point
  phonbandvec = zeros(size(p_temp)); % phonon band index

  % U_temp = LeadPhonon.U_plus; % DEBUG

  for nq = 1:1:nqmax % loop over all channels (transmitting and nontransmitting) at frequency w
      qvec = qfrac1_temp(nq)*Gvec1 + qfrac2_temp(nq)*Gvec2 + qfrac3_temp(nq)*Gvec3; % 'folded' wave vector
      qfrac = [qfrac1_temp(nq) qfrac2_temp(nq) qfrac3_temp(nq)];
      if gt(p_temp(nq),0) % i.e. there exists a transmitting channel
          % === Match folded and unfolded k-points, identify band index of phonon mode ===
          [qvec_BZ, vel_BZ, ind_BZ, dw_BZ] = ...
            UnfoldWaveVector(w,qfrac,ScaleFactor,Lead,LeadPhonDisp,epsilon_cutoff_w);
  
          % === Match unfolded k-point to Brillouin Zone k-point
          qvec_BZ = BrillouinZoneMapping(qvec_BZ,ScaleFactor,PrimGvec1,PrimGvec2,PrimGvec3);  

          qxvec(nq) = qvec_BZ(1);
          qyvec(nq) = qvec_BZ(2);
          qzvec(nq) = qvec_BZ(3);
          vxvec(nq) = vel_BZ(1);
          vyvec(nq) = vel_BZ(2);
          vzvec(nq) = vel_BZ(3);
          dwvec(nq) = dw_BZ; 
          phonbandvec(nq) = ind_BZ;
      end
  end

  % LeadPhonon.U_plus = U_temp; % DEBUG
  % U_plus_before = LeadPhonon_in.U_plus; % DEBUG
  % U_plus_after = LeadPhonon.U_plus; % DEBUG
  % save('TempData','U_plus_before','U_plus_after','phonbandvec');

  LeadPhonon.VecQx_plus = qxvec;
  LeadPhonon.VecQy_plus = qyvec;
  LeadPhonon.VecQz_plus = qzvec;
  LeadPhonon.VecVelx_plus = vxvec;
  LeadPhonon.VecVely_plus = vyvec;
  LeadPhonon.VecVelz_plus = vzvec;
  LeadPhonon.VecDw_plus = dwvec;
  LeadPhonon.VecB_plus  = phonbandvec;

  % =================================================================================================
  % === transverse wave vector(s) ===
  qfrac2_temp = [LeadPhonon.VecQ_tran1]/(2*pi); % fractional wave vector (transv. 1) (b/w 0 and 1)
  qfrac2_temp = mod(qfrac2_temp+0.5,1.0)-0.5;       % set value between -0.5 and 0.4999...
  qfrac3_temp = [LeadPhonon.VecQ_tran2]/(2*pi); % fractional wave vector (transv. 2) (b/w 0 and 1)
  qfrac3_temp = mod(qfrac3_temp+0.5,1.0)-0.5;       % set value between -0.5 and 0.4999...

  % === longitudinal wave vector ===
  vlong_temp  = [LeadPhonon.VecV_plus_adv]; % AGF phonon longitudinal group velocity 
  qfrac1_temp = [LeadPhonon.VecQ_plus_adv]/(2*pi); % fractional wave vector (long.) (b/w 0 and 1)  
  qfrac1_temp = mod(qfrac1_temp+0.5,1.0)-0.5;           % set value between -0.5 and 0.4999...
  p_temp      = round([LeadPhonon.VecT_plus_adv]); % maximum transmission of propagating phonons

  [qfrac1_temp,qfrac2_temp,qfrac3_temp] = ... 
    RedistributeDegenerateModes(qfrac1_temp,qfrac2_temp,qfrac3_temp,p_temp,epsilon_cutoff_q);

  nqmax = numel(p_temp); % number of channels (transmitting and nontransmitting) at frequency w
  qxvec = zeros(size(p_temp)); % allocate data for 'unfolded' wave vector x-component
  qyvec = zeros(size(p_temp)); % allocate data for 'unfolded' wave vector y-component
  qzvec = zeros(size(p_temp)); % allocate data for 'unfolded' wave vector z-component
  vxvec = zeros(size(p_temp)); % allocate data for 'unfolded' velocity x-component
  vyvec = zeros(size(p_temp)); % allocate data for 'unfolded' velocity y-component
  vzvec = zeros(size(p_temp)); % allocate data for 'unfolded' velocity z-component
  dwvec = -ones(size(p_temp)); % discrepancy between frequency of folded and unfolded k point
  phonbandvec = zeros(size(p_temp)); % phonon band index

  for nq = 1:1:nqmax % loop over all channels (transmitting and nontransmitting) at frequency w
      qvec = qfrac1_temp(nq)*Gvec1 + qfrac2_temp(nq)*Gvec2 + qfrac3_temp(nq)*Gvec3; % 'folded' wave vector
      qfrac = [qfrac1_temp(nq) qfrac2_temp(nq) qfrac3_temp(nq)];
      if gt(p_temp(nq),0) % i.e. there exists a transmitting channel
          % === Match folded and unfolded k-points, identify band index of phonon mode ===
          [qvec_BZ, vel_BZ, ind_BZ, dw_BZ] = ...
            UnfoldWaveVector(w,qfrac,ScaleFactor,Lead,LeadPhonDisp,epsilon_cutoff_w);
  
          % === Match unfolded k-point to Brillouin Zone k-point
          qvec_BZ = BrillouinZoneMapping(qvec_BZ,ScaleFactor,PrimGvec1,PrimGvec2,PrimGvec3);

          qxvec(nq) = qvec_BZ(1);
          qyvec(nq) = qvec_BZ(2);
          qzvec(nq) = qvec_BZ(3);
          vxvec(nq) = vel_BZ(1);
          vyvec(nq) = vel_BZ(2);
          vzvec(nq) = vel_BZ(3);
          dwvec(nq) = dw_BZ; 
          phonbandvec(nq) = ind_BZ;
      end
  end

  LeadPhonon.VecQx_plus_adv = qxvec;
  LeadPhonon.VecQy_plus_adv = qyvec;
  LeadPhonon.VecQz_plus_adv = qzvec;
  LeadPhonon.VecVelx_plus_adv = vxvec;
  LeadPhonon.VecVely_plus_adv = vyvec;
  LeadPhonon.VecVelz_plus_adv = vzvec;
  LeadPhonon.VecDw_plus_adv = dwvec;
  LeadPhonon.VecB_plus_adv  = phonbandvec;

  % =================================================================================================
  % === transverse wave vector(s) ===
  qfrac2_temp = [LeadPhonon.VecQ_tran1]/(2*pi); % fractional wave vector (transv. 1) (b/w 0 and 1)
  qfrac2_temp = mod(qfrac2_temp+0.5,1.0)-0.5;       % set value between -0.5 and 0.4999...
  qfrac3_temp = [LeadPhonon.VecQ_tran2]/(2*pi); % fractional wave vector (transv. 2) (b/w 0 and 1)
  qfrac3_temp = mod(qfrac3_temp+0.5,1.0)-0.5;       % set value between -0.5 and 0.4999...

  % === longitudinal wave vector ===
  vlong_temp  = [LeadPhonon.VecV_minus]; % AGF phonon longitudinal group velocity 
  qfrac1_temp = [LeadPhonon.VecQ_minus]/(2*pi); % fractional wave vector (long.) (b/w 0 and 1)  
  qfrac1_temp = mod(qfrac1_temp+0.5,1.0)-0.5;           % set value between -0.5 and 0.4999...
  p_temp      = round([LeadPhonon.VecT_minus]); % maximum transmission of propagating phonons

  [qfrac1_temp,qfrac2_temp,qfrac3_temp] = ...
    RedistributeDegenerateModes(qfrac1_temp,qfrac2_temp,qfrac3_temp,p_temp,epsilon_cutoff_q);

  nqmax = numel(p_temp); % number of channels (transmitting and nontransmitting) at frequency w
  qxvec = zeros(size(p_temp)); % allocate data for 'unfolded' wave vector x-component
  qyvec = zeros(size(p_temp)); % allocate data for 'unfolded' wave vector y-component
  qzvec = zeros(size(p_temp)); % allocate data for 'unfolded' wave vector z-component
  vxvec = zeros(size(p_temp)); % allocate data for 'unfolded' velocity x-component
  vyvec = zeros(size(p_temp)); % allocate data for 'unfolded' velocity y-component
  vzvec = zeros(size(p_temp)); % allocate data for 'unfolded' velocity z-component
  dwvec = -ones(size(p_temp)); % discrepancy between frequency of folded and unfolded k point
  phonbandvec = zeros(size(p_temp)); % phonon band index

  for nq = 1:1:nqmax % loop over all channels (transmitting and nontransmitting) at frequency w
      qvec = qfrac1_temp(nq)*Gvec1 + qfrac2_temp(nq)*Gvec2 + qfrac3_temp(nq)*Gvec3; % 'folded' wave vector
      qfrac = [qfrac1_temp(nq) qfrac2_temp(nq) qfrac3_temp(nq)];
      if gt(p_temp(nq),0) % i.e. there exists a transmitting channel
          % === Match folded and unfolded k-points, identify band index of phonon mode ===
          [qvec_BZ, vel_BZ, ind_BZ, dw_BZ] = ...
            UnfoldWaveVector(w,qfrac,ScaleFactor,Lead,LeadPhonDisp,epsilon_cutoff_w);
  
          % === Match unfolded k-point to Brillouin Zone k-point
          qvec_BZ = BrillouinZoneMapping(qvec_BZ,ScaleFactor,PrimGvec1,PrimGvec2,PrimGvec3);

          qxvec(nq) = qvec_BZ(1);
          qyvec(nq) = qvec_BZ(2);
          qzvec(nq) = qvec_BZ(3);
          vxvec(nq) = vel_BZ(1);
          vyvec(nq) = vel_BZ(2);
          vzvec(nq) = vel_BZ(3);
          dwvec(nq) = dw_BZ; 
          phonbandvec(nq) = ind_BZ;
      end
  end

  LeadPhonon.VecQx_minus = qxvec;
  LeadPhonon.VecQy_minus = qyvec;
  LeadPhonon.VecQz_minus = qzvec;
  LeadPhonon.VecVelx_minus = vxvec;
  LeadPhonon.VecVely_minus = vyvec;
  LeadPhonon.VecVelz_minus = vzvec;
  LeadPhonon.VecDw_minus = dwvec;
  LeadPhonon.VecB_minus  = phonbandvec;

  % =================================================================================================
  % === transverse wave vector(s) ===
  qfrac2_temp = [LeadPhonon.VecQ_tran1]/(2*pi); % fractional wave vector (transv. 1) (b/w 0 and 1)
  qfrac2_temp = mod(qfrac2_temp+0.5,1.0)-0.5;       % set value between -0.5 and 0.4999...
  qfrac3_temp = [LeadPhonon.VecQ_tran2]/(2*pi); % fractional wave vector (transv. 2) (b/w 0 and 1)
  qfrac3_temp = mod(qfrac3_temp+0.5,1.0)-0.5;       % set value between -0.5 and 0.4999...

  % === longitudinal wave vector ===
  vlong_temp  = [LeadPhonon.VecV_minus_adv]; % AGF phonon longitudinal group velocity 
  qfrac1_temp = [LeadPhonon.VecQ_minus_adv]/(2*pi); % fractional wave vector (long.) (b/w 0 and 1)  
  qfrac1_temp = mod(qfrac1_temp+0.5,1.0)-0.5;           % set value between -0.5 and 0.4999...
  p_temp      = round([LeadPhonon.VecT_minus_adv]); % maximum transmission of propagating phonons

  [qfrac1_temp,qfrac2_temp,qfrac3_temp] = ... 
    RedistributeDegenerateModes(qfrac1_temp,qfrac2_temp,qfrac3_temp,p_temp,epsilon_cutoff_q);

  nqmax = numel(p_temp); % number of channels (transmitting and nontransmitting) at frequency w
  qxvec = zeros(size(p_temp)); % allocate data for 'unfolded' wave vector x-component
  qyvec = zeros(size(p_temp)); % allocate data for 'unfolded' wave vector y-component
  qzvec = zeros(size(p_temp)); % allocate data for 'unfolded' wave vector z-component
  vxvec = zeros(size(p_temp)); % allocate data for 'unfolded' group velocity x-component
  vyvec = zeros(size(p_temp)); % allocate data for 'unfolded' group velocity y-component
  vzvec = zeros(size(p_temp)); % allocate data for 'unfolded' group velocity z-component
  dwvec = -ones(size(p_temp)); % discrepancy between frequency of folded and unfolded k point
  phonbandvec = zeros(size(p_temp)); % phonon band index

  for nq = 1:1:nqmax % loop over all channels (transmitting and nontransmitting) at frequency w
      qvec = qfrac1_temp(nq)*Gvec1 + qfrac2_temp(nq)*Gvec2 + qfrac3_temp(nq)*Gvec3; % 'folded' wave vector
      qfrac = [qfrac1_temp(nq) qfrac2_temp(nq) qfrac3_temp(nq)];
      if gt(p_temp(nq),0) % i.e. there exists a transmitting channel
          % === Match folded and unfolded k-points, identify band index of phonon mode ===
          [qvec_BZ, vel_BZ, ind_BZ, dw_BZ] = ...
            UnfoldWaveVector(w,qfrac,ScaleFactor,Lead,LeadPhonDisp,epsilon_cutoff_w);
  
          % === Match unfolded k-point to Brillouin Zone k-point
          qvec_BZ = BrillouinZoneMapping(qvec_BZ,ScaleFactor,PrimGvec1,PrimGvec2,PrimGvec3);

          qxvec(nq) = qvec_BZ(1);
          qyvec(nq) = qvec_BZ(2);
          qzvec(nq) = qvec_BZ(3);
          vxvec(nq) = vel_BZ(1);
          vyvec(nq) = vel_BZ(2);
          vzvec(nq) = vel_BZ(3);
          dwvec(nq) = dw_BZ; 
          phonbandvec(nq) = ind_BZ;
      end
  end

  LeadPhonon.VecQx_minus_adv = qxvec;
  LeadPhonon.VecQy_minus_adv = qyvec;
  LeadPhonon.VecQz_minus_adv = qzvec;
  LeadPhonon.VecVelx_minus_adv = vxvec;
  LeadPhonon.VecVely_minus_adv = vyvec;
  LeadPhonon.VecVelz_minus_adv = vzvec;
  LeadPhonon.VecDw_minus_adv = dwvec;
  LeadPhonon.VecB_minus_adv  = phonbandvec;
  % ------------------------------------------------------------------------------

  % LeadPhonon = rmfield(LeadPhonon,'VecQ_plus');
  % LeadPhonon = rmfield(LeadPhonon,'VecQ_plus_adv');
  % LeadPhonon = rmfield(LeadPhonon,'VecQ_minus');
  % LeadPhonon = rmfield(LeadPhonon,'VecQ_minus_adv');
end


% ================================================================
% ================================================================
function [qfrac1_out,qfrac2_out,qfrac3_out] = ... 
  RedistributeDegenerateModes(qfrac1_in,qfrac2_in,qfrac3_in,p_in,epsilon_cutoff_q)
% qfrac1_in = 1D array with values between -0.5 to 0.4999... (for longitudinal)
% qfrac2_in = 1D array with values between -0.5 to 0.4999... (for transverse 1)
% qfrac3_in = 1D array with values between -0.5 to 0.4999... (for transverse 2)
% p_in = 1D array with values 0 (evanescent) or 1 (extended/propagating)

  for nq1 = 1:numel(p_in)
      qvec1 = [qfrac1_in(nq1); qfrac2_in(nq1); qfrac3_in(nq1)];
      nqvec1 = [nq1];

      if eq(p_in(nq1),1) % if this mode is extended
          for nq2 = (nq1+1):numel(p_in);
              qvec2 = [qfrac1_in(nq2); qfrac2_in(nq2); qfrac3_in(nq2)];        
              epsilon = norm(qvec1-qvec2)/norm(qvec1);
              if and(le(epsilon,epsilon_cutoff_q),eq(p_in(nq2),1)) % i.e. the other mode is degenerate and extended
                  nqvec1 = [nqvec1 nq2]; % store index of degenerate mode
              end
          end

         if eq(numel(nqvec1),2)
             % fprintf(1,'<DEBUG> Degeneracy = 2 \n');
             if eq(qvec1(2),-0.5)
                 qfrac2_in(nqvec1(2)) = 0.5;
                 % fprintf(1,'<DEBUG> %2d %2d %f %f \n',nqvec1(1),nqvec1(2),qfrac2_in(nqvec1(1)),qfrac2_in(nqvec1(2)));
             elseif eq(qvec1(3),-0.5)
                 qfrac3_in(nqvec1(2)) = 0.5;
                 % fprintf(1,'<DEBUG> %2d %2d %f %f \n',nqvec1(1),nqvec1(2),qfrac3_in(nqvec1(1)),qfrac3_in(nqvec1(2)));
             end
         elseif eq(numel(nqvec1),4)
             % fprintf(1,'<DEBUG> Degeneracy = 4 \n');
             qfrac2_in(nqvec1(2)) = -0.5;
             qfrac3_in(nqvec1(2)) = 0.5;
             qfrac2_in(nqvec1(3)) = 0.5;
             qfrac3_in(nqvec1(3)) = -0.5;
             qfrac2_in(nqvec1(4)) = 0.5;
             qfrac3_in(nqvec1(4)) = 0.5;
         else
             if gt(numel(nqvec1),1)
                 warning('<!!!> Degeneracy is neither 2 nor 4');
             end
         end
      end
  end

  qfrac1_out = qfrac1_in;  
  qfrac2_out = qfrac2_in;  
  qfrac3_out = qfrac3_in;  
end


% ================================================================
% ================================================================
function [qvec_out, vel_out, ind_out, dw_out] = ...
UnfoldWaveVector(w_in,qfrac_in,ScaleFactor_in,Lead_in,LeadPhonDisp_in,epsilon_cutoff_w)
% qvec_out == unfolded wave vector 
% vel_out  == corresponding 'unfolded' phonon group velocity
% ind_out  == corresponding 'unfolded' phonon band index
% dw_out   == frequency difference between unfolded and folded wave vector
  Gvec1 = Lead_in.gvec_long;  % reciprocal lattice vector (AGF)
  Gvec2 = Lead_in.gvec_tran1; % reciprocal lattice vector (AGF)
  Gvec3 = Lead_in.gvec_tran2; % reciprocal lattice vector (AGF)
  PrimCell = LeadPhonDisp_in.PrimCell;

  b1 = qfrac_in(1);
  b2 = qfrac_in(2);
  b3 = qfrac_in(3);

  sgn1 = -sign(b1 + 1E-6);
  sgn2 = -sign(b2 + 1E-6);
  sgn3 = -sign(b3 + 1E-6);
  qvec_in = b1*Gvec1 + b2*Gvec2 + b3*Gvec3; % original unfolded wave vector

  % === Set default values of return variables ===
  qvec_out = qvec_in;
  vel_out  = [0.0 0.0 0.0];
  ind_out  = 0;
  dw_out   = 0;

  % === Find unfolded wave vector in Brillouin Zone ===
  LoopExitFlag = 0;

  for ns1 = (0:ScaleFactor_in)*sgn1
      if eq(LoopExitFlag,1)
          break;
      end

      for ns2 = (0:ScaleFactor_in)*sgn2
          if eq(LoopExitFlag,1)
              break;
          end

          for ns3 = (0:ScaleFactor_in)*sgn3
              qvec_s = qvec_in + ns1*Gvec1 + ns2*Gvec2 + ns3*Gvec3; % 'unfolded'/shifted wave vector
              Hq_s = PrimCell(1).H*exp(-1i*2*pi*PrimCell(1).rvec*qvec_s');
              for np = 2:numel(PrimCell)
                  Hq_s = Hq_s + PrimCell(np).H*exp(-1i*2*pi*PrimCell(np).rvec*qvec_s');
              end
              Hq_s = 0.5*(Hq_s+Hq_s'); % make dynamical matrix completely hermitean
              wq_s = sort(real(sqrt(eig(Hq_s)))); % sorted frequency eigenvalues for 'unfolded' wave vector
  
              [dw_tmp, ind_tmp] = min(abs(wq_s-w_in));
                          
              if (eq(ns1,0)*eq(ns2,0)*eq(ns3,0)) % at original 'folded' k-point
                  dw_out = dw_tmp;
                  ind_out = ind_tmp;
                  qvec_out = qvec_s;                
              else                
                  if lt(dw_tmp,dw_out) 
                      dw_out = dw_tmp;
                      ind_out = ind_tmp;
                      qvec_out = qvec_s; % set 'unfolded' k-point                    
                  end
              end                                     

              epsilon = abs(dw_tmp/w_in);
              if lt(epsilon,epsilon_cutoff_w)
                  LoopExitFlag = 1;
                  break;
              end           
          end
      end      
  end

  % disp([w_in wq_s(ind_tmp)]); % DEBUG

  % === Find corresponding phonon group velocity of unfolded wave vector ==
  Hq = PrimCell(1).H*exp(-1i*2*pi*PrimCell(1).rvec*qvec_out');
  for np = 2:numel(PrimCell)
      Hq = Hq + PrimCell(np).H*exp(-1i*2*pi*PrimCell(np).rvec*qvec_out');
  end
  [Uq, wwq] = eig(Hq);
  [dw_tmp, ind_tmp] = min(abs(diag(wwq)-w_in^2));
  U_out = Uq(:,ind_tmp); % column vector describing eigenmode

  Hq_x = -PrimCell(1).H*exp(-1i*2*pi*PrimCell(1).rvec*qvec_out')*1i*PrimCell(1).rvec(1);
  Hq_y = -PrimCell(1).H*exp(-1i*2*pi*PrimCell(1).rvec*qvec_out')*1i*PrimCell(1).rvec(2);
  Hq_z = -PrimCell(1).H*exp(-1i*2*pi*PrimCell(1).rvec*qvec_out')*1i*PrimCell(1).rvec(3);
  for np = 2:numel(PrimCell)
      Hq_x = Hq_x - PrimCell(np).H*exp(-1i*2*pi*PrimCell(np).rvec*qvec_out') ...
             *1i*PrimCell(np).rvec(1);
      Hq_y = Hq_y - PrimCell(np).H*exp(-1i*2*pi*PrimCell(np).rvec*qvec_out') ...
             *1i*PrimCell(np).rvec(2);
      Hq_z = Hq_z - PrimCell(np).H*exp(-1i*2*pi*PrimCell(np).rvec*qvec_out') ...
             *1i*PrimCell(np).rvec(3);
  end

  vel_out(1,1) = real(1.0/(2.0*w_in)*U_out'*Hq_x*U_out);
  vel_out(1,2) = real(1.0/(2.0*w_in)*U_out'*Hq_y*U_out);
  vel_out(1,3) = real(1.0/(2.0*w_in)*U_out'*Hq_z*U_out);

  % dwvec(nq) = dw_s; 
  % phonbandvec(nq) = ind_s;
end


% ================================================================
% ================================================================
function qvec_out = BrillouinZoneMapping(qvec_in,ScaleFactor_in,PrimGvec1_in,PrimGvec2_in,PrimGvec3_in)
  qvec_s = qvec_in;

  % === Map wave vector to have positive fractional coefficients in reciprocal space ===
  fracqvec1_s = qvec_s*cross(PrimGvec2_in,PrimGvec3_in)'/(PrimGvec1_in*cross(PrimGvec2_in,PrimGvec3_in)');
  fracqvec2_s = qvec_s*cross(PrimGvec3_in,PrimGvec1_in)'/(PrimGvec2_in*cross(PrimGvec3_in,PrimGvec1_in)');
  fracqvec3_s = qvec_s*cross(PrimGvec1_in,PrimGvec2_in)'/(PrimGvec3_in*cross(PrimGvec1_in,PrimGvec2_in)');

  fracqvec1_s = mod(fracqvec1_s,1.0);
  fracqvec2_s = mod(fracqvec2_s,1.0);
  fracqvec3_s = mod(fracqvec3_s,1.0);

  qvec_s = fracqvec1_s*PrimGvec1_in + fracqvec2_s*PrimGvec2_in + fracqvec3_s*PrimGvec3_in;

  for ns1 = -1:1:0
      for ns2 = -1:1:0
          for ns3 = -1:1:0
              qvec_tmp = (fracqvec1_s+ns1)*PrimGvec1_in + (fracqvec2_s+ns2)*PrimGvec2_in ...
                         + (fracqvec3_s+ns3)*PrimGvec3_in; % shifted wave vector 
              if lt(norm(qvec_tmp),norm(qvec_s)) 
                  qvec_s = qvec_tmp;
              end
          end
      end
  end

  qvec_out = qvec_s;
end


