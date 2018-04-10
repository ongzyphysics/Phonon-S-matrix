function LeadPhonon = RebuildExtendedModeEigenvector(w_in,LeadPhonon_in,Lead_in,LeadPhonDisp_in)
  w = w_in;
  LeadPhonon = LeadPhonon_in;
  Lead = Lead_in;
  LeadPhonDisp = LeadPhonDisp_in;


  %{
  % save('TempData','w','U_before','U_after','p_temp','w_after',...
  %      'qxvec','qyvec','qzvec','qxvec_agf','qyvec_agf','qzvec_agf'); % DEBUG
  %}

  if and(gt(numel(LeadPhonDisp.nbasis_bulk),0),or(gt(Lead.n_tran1,1),gt(Lead.n_tran1,1)))
  % i.e. (1) AGF2BULK parameters are given and (2) there are multiple transverse momentum components 
      % === Parameters for U_plus ===
      p_temp = round([LeadPhonon.VecT_plus]); % maximum transmission of propagating phonons
      U_tmp = LeadPhonon.U_plus; 
      qxvec = LeadPhonon.VecQx_plus;
      qyvec = LeadPhonon.VecQy_plus;
      qzvec = LeadPhonon.VecQz_plus;
      nqmax = numel(qxvec);

      for nq = 1:1:nqmax % loop over all channels (transmitting and nontransmitting) at frequency w
          if gt(p_temp(nq),0) % i.e. there exists a transmitting channel
              qvec_BZ = [qxvec(nq) qyvec(nq) qzvec(nq)];
              [w_tmp, U_tmp(:,nq)] = FixEigenmode(w,U_tmp(:,nq),qvec_BZ,LeadPhonDisp);
          end
      end
      LeadPhonon.U_plus = U_tmp;

      % === Parameters for U_plus_adv ===
      p_temp = round([LeadPhonon.VecT_plus_adv]); % maximum transmission of propagating phonons
      U_tmp = LeadPhonon.U_plus_adv; 
      qxvec = LeadPhonon.VecQx_plus_adv;
      qyvec = LeadPhonon.VecQy_plus_adv;
      qzvec = LeadPhonon.VecQz_plus_adv;
      nqmax = numel(qxvec);

      for nq = 1:1:nqmax % loop over all channels (transmitting and nontransmitting) at frequency w
          if gt(p_temp(nq),0) % i.e. there exists a transmitting channel
              qvec_BZ = [qxvec(nq) qyvec(nq) qzvec(nq)];
              [w_tmp, U_tmp(:,nq)] = FixEigenmode(w,U_tmp(:,nq),qvec_BZ,LeadPhonDisp);
          end
      end
      LeadPhonon.U_plus_adv = U_tmp;

      % === Parameters for U_minus ===
      p_temp = round([LeadPhonon.VecT_minus]); % maximum transmission of propagating phonons
      U_tmp = LeadPhonon.U_minus; 
      qxvec = LeadPhonon.VecQx_minus;
      qyvec = LeadPhonon.VecQy_minus;
      qzvec = LeadPhonon.VecQz_minus;
      nqmax = numel(qxvec);
      for nq = 1:1:nqmax % loop over all channels (transmitting and nontransmitting) at frequency w
          if gt(p_temp(nq),0) % i.e. there exists a transmitting channel
              qvec_BZ = [qxvec(nq) qyvec(nq) qzvec(nq)];
              [w_tmp, U_tmp(:,nq)] = FixEigenmode(w,U_tmp(:,nq),qvec_BZ,LeadPhonDisp);
          end
      end
      LeadPhonon.U_minus = U_tmp;
    
      % === Parameters for U_minus_adv ===
      p_temp = round([LeadPhonon.VecT_minus_adv]); % maximum transmission of propagating phonons
      U_tmp = LeadPhonon.U_minus_adv; 
      qxvec = LeadPhonon.VecQx_minus_adv;
      qyvec = LeadPhonon.VecQy_minus_adv;
      qzvec = LeadPhonon.VecQz_minus_adv;
      nqmax = numel(qxvec);
      for nq = 1:1:nqmax % loop over all channels (transmitting and nontransmitting) at frequency w
          if gt(p_temp(nq),0) % i.e. there exists a transmitting channel
              qvec_BZ = [qxvec(nq) qyvec(nq) qzvec(nq)];
              [w_tmp, U_tmp(:,nq)] = FixEigenmode(w,U_tmp(:,nq),qvec_BZ,LeadPhonDisp);
          end
      end
      LeadPhonon.U_minus_adv = U_tmp;
  else 
  % i.e. no AGF2BULK parameters are given
  % in that case, we reverse the unfolding of the wave vectors
      Gvec1 = Lead.gvec_long;  % reciprocal lattice vector (AGF)
      Gvec2 = Lead.gvec_tran1; % reciprocal lattice vector (AGF)
      Gvec3 = Lead.gvec_tran2; % reciprocal lattice vector (AGF)

      % === Parameters for U_plus ===
      qfrac1 = LeadPhonon.VecQ_plus/(2*pi);  % longitudinal wave vector (-0.5 to 0.5)
      qfrac2 = LeadPhonon.VecQ_tran1/(2*pi); % transverse 1 wave vector (-0.5 to 0.5)
      qfrac3 = LeadPhonon.VecQ_tran2/(2*pi); % transverse 2 wave vector (-0.5 to 0.5)
      qfrac1 = mod(qfrac1+0.5,1.0)-0.5;
      qfrac2 = mod(qfrac2+0.5,1.0)-0.5;
      qfrac3 = mod(qfrac3+0.5,1.0)-0.5;

      nqmax = numel(qfrac1);
      for nq = 1:1:nqmax % loop over all channels (transmitting and nontransmitting) at frequency w
          qvec_agf = qfrac1(nq)*Gvec1 + qfrac2(nq)*Gvec2 + qfrac3(nq)*Gvec3;
          LeadPhonon.VecQx_plus(nq,1) = qvec_agf(1);
          LeadPhonon.VecQy_plus(nq,1) = qvec_agf(2);
          LeadPhonon.VecQz_plus(nq,1) = qvec_agf(3);
      end

      % === Parameters for U_plus_adv ===
      qfrac1 = LeadPhonon.VecQ_plus_adv/(2*pi);  % longitudinal wave vector (-0.5 to 0.5)
      qfrac2 = LeadPhonon.VecQ_tran1/(2*pi); % transverse 1 wave vector (-0.5 to 0.5)
      qfrac3 = LeadPhonon.VecQ_tran2/(2*pi); % transverse 2 wave vector (-0.5 to 0.5)
      qfrac1 = mod(qfrac1+0.5,1.0)-0.5;
      qfrac2 = mod(qfrac2+0.5,1.0)-0.5;
      qfrac3 = mod(qfrac3+0.5,1.0)-0.5;

      nqmax = numel(qfrac1);
      for nq = 1:1:nqmax % loop over all channels (transmitting and nontransmitting) at frequency w
          qvec_agf = qfrac1(nq)*Gvec1 + qfrac2(nq)*Gvec2 + qfrac3(nq)*Gvec3;
          LeadPhonon.VecQx_plus_adv(nq,1) = qvec_agf(1);
          LeadPhonon.VecQy_plus_adv(nq,1) = qvec_agf(2);
          LeadPhonon.VecQz_plus_adv(nq,1) = qvec_agf(3);
      end

      % === Parameters for U_minus ===
      qfrac1 = LeadPhonon.VecQ_minus/(2*pi);  % longitudinal wave vector (-0.5 to 0.5)
      qfrac2 = LeadPhonon.VecQ_tran1/(2*pi); % transverse 1 wave vector (-0.5 to 0.5)
      qfrac3 = LeadPhonon.VecQ_tran2/(2*pi); % transverse 2 wave vector (-0.5 to 0.5)
      qfrac1 = mod(qfrac1+0.5,1.0)-0.5;
      qfrac2 = mod(qfrac2+0.5,1.0)-0.5;
      qfrac3 = mod(qfrac3+0.5,1.0)-0.5;

      nqmax = numel(qfrac1);
      for nq = 1:1:nqmax % loop over all channels (transmitting and nontransmitting) at frequency w
          qvec_agf = qfrac1(nq)*Gvec1 + qfrac2(nq)*Gvec2 + qfrac3(nq)*Gvec3;
          LeadPhonon.VecQx_minus(nq,1) = qvec_agf(1);
          LeadPhonon.VecQy_minus(nq,1) = qvec_agf(2);
          LeadPhonon.VecQz_minus(nq,1) = qvec_agf(3);
      end

      % === Parameters for U_minus_adv ===
      qfrac1 = LeadPhonon.VecQ_minus_adv/(2*pi);  % longitudinal wave vector (-0.5 to 0.5)
      qfrac2 = LeadPhonon.VecQ_tran1/(2*pi); % transverse 1 wave vector (-0.5 to 0.5)
      qfrac3 = LeadPhonon.VecQ_tran2/(2*pi); % transverse 2 wave vector (-0.5 to 0.5)
      qfrac1 = mod(qfrac1+0.5,1.0)-0.5;
      qfrac2 = mod(qfrac2+0.5,1.0)-0.5;
      qfrac3 = mod(qfrac3+0.5,1.0)-0.5;

      nqmax = numel(qfrac1);
      for nq = 1:1:nqmax % loop over all channels (transmitting and nontransmitting) at frequency w
          qvec_agf = qfrac1(nq)*Gvec1 + qfrac2(nq)*Gvec2 + qfrac3(nq)*Gvec3;
          LeadPhonon.VecQx_minus_adv(nq,1) = qvec_agf(1);
          LeadPhonon.VecQy_minus_adv(nq,1) = qvec_agf(2);
          LeadPhonon.VecQz_minus_adv(nq,1) = qvec_agf(3);
      end

  end
end



% ==================================================================================
% ==================================================================================
function [w_out, U_out] = FixEigenmode(w_in,U_in,qvec_in,LeadPhonDisp_in)
% We rebuild the AGF mode eigenvector using the primitive lattice matrices and vectors
%   w_in    == eigenfrequency of mode
%   U_in    == eigenvector of mode
%   qvec_in == wave vector of mode
%   LeadPhonDisp_in == parameters for bulk lattice phonon dispersion

  % ==================================================================================
  % === Set variables ===
  PrimCell = LeadPhonDisp_in.PrimCell; % data structure with phonon dispersion information
  qvec = qvec_in; % unfolded wave vector of eigenmode
  ww = w_in*w_in; % frequency square of eigenmode 

  % == Compute all available eigenfrequencies for given wave vector ===
  Hq = PrimCell(1).H*exp(-1i*2*pi*PrimCell(1).rvec*qvec');
  for np = 2:numel(PrimCell)
      Hq = Hq + PrimCell(np).H*exp(-1i*2*pi*PrimCell(np).rvec*qvec');
  end
  Hq = 0.5*(Hq+Hq'); % make dynamical matrix completely hermitean
  [Uq,wwq] = eig(Hq);
  wwq = real(diag(wwq));

  % disp([ww sort(wwq)']);

  % === Determine corresponding primitive eigenvector for given unfolded wave vector ===
  [dww_tmp, ind_tmp] = min(abs(wwq-ww));
  U_tmp = Uq(:,ind_tmp); % reciprocal space eigenvector for primitive BZ
  w_out = sqrt(wwq(ind_tmp));

  % === Build folded eigenvector for AGF supercell lattice ===
  for n = 1:numel(LeadPhonDisp_in.nbasis_bulk) % loop over basis atoms in AGF cell 
      m = LeadPhonDisp_in.nbasis_bulk(n); % corresponding basis atom in primitive unit cell
      rvec(1,1) = LeadPhonDisp_in.BasisAtomRxvec(n); % x-component of lattice cell position
      rvec(1,2) = LeadPhonDisp_in.BasisAtomRyvec(n); % y-component of lattice cell position
      rvec(1,3) = LeadPhonDisp_in.BasisAtomRzvec(n); % z-component of lattice cell position
      U_agf((1:3)+3*(n-1),1) = U_tmp((1:3)+3*(m-1),1)*exp(1i*2*pi*qvec*rvec');
  end

  U_agf = U_agf/norm(U_agf); % eigenvector normalization

  U_out = U_in;

  n_tran = numel(U_out)/numel(U_agf); % number of transverse cells
  nsize_agf = numel(U_agf);

  for n = 1:n_tran
      SubU_i = U_in(1:nsize_agf,1);
      SubU_f = U_in((1:nsize_agf)+(n-1)*nsize_agf,1);
      exp_phase = SubU_i'*SubU_f/(SubU_f'*SubU_f);
      U_out((1:nsize_agf)+(n-1)*nsize_agf,1) = U_agf*exp_phase;
  end
  
  U_out = U_out/norm(U_out); % normalization
end

