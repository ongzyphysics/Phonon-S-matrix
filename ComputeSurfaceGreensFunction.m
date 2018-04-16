function LeadPhonon = ComputeSurfaceGreensFunction(w_in,ww_in,InputParam)
% NOTE: Compute the left and right-facing surface Greens function of given lead

  w = w_in;
  ww = ww_in;
  
  HC = InputParam.MatHC;
  HL = InputParam.MatHL;
  HR = InputParam.MatHR;
  D  = eye(size(HC));
  a_long = InputParam.a_long; % lattice constant
  n_tran1 = InputParam.n_tran1;
  n_tran2 = InputParam.n_tran2;
  n_tran = n_tran1*n_tran2; % number of transverse Fourier components

  if gt(n_tran,1) % i.e. if there are multiple transverse Fourier components
      SubInputParam = TransverseSubspaceTransform(InputParam); 
        % transform large matrices into Fourier-component submatrices
        % this step is somewhat inefficient as it is repeated for every frequency point 

      % === Compute the surface Green's function, eigenmodes and eigenvelocities for each Fourier component ===
      for nt1 = 1:1:n_tran1
          for nt2 = 1:1:n_tran2
              SubLeadPhonon(nt1,nt2) = ComputeSurfaceGreensFunction(w,ww,SubInputParam(nt1,nt2));
          end
      end

      LeadPhonon = CombineTransverseSubspaceComponents(SubLeadPhonon,InputParam);
  else 
      TauR = -HR;
      TauL = -HL;
    
      m_max = 40; % maximum number of recursive steps in self-consistent inner loop
            
      % === Compute surface Green's functions ===
      W = ww*D-HC; % initial surface hamiltonian      
      InvW = inv(W);
      WsR_old = W - TauR*InvW*TauL;
      WsL_old = W - TauL*InvW*TauR;
      Wb_old = W - TauR*InvW*TauL - TauL*InvW*TauR;
      TauR_old = -TauR*InvW*TauR;
      TauL_old = -TauL*InvW*TauL;
          
      for m = 1:m_max
          InvWb_old = inv(Wb_old);
          WsR_new = WsR_old - TauR_old*InvWb_old*TauL_old;
          WsL_new = WsL_old - TauL_old*InvWb_old*TauR_old;
          Wb_new = Wb_old - TauR_old*InvWb_old*TauL_old - TauL_old*InvWb_old*TauR_old;
          TauR_new = -TauR_old*InvWb_old*TauR_old;
          TauL_new = -TauL_old*InvWb_old*TauL_old;
          
          epsilon = norm(Wb_new-Wb_old)/norm(Wb_new);
            
          WsR_old = WsR_new;
          WsL_old = WsL_new;
          Wb_old = Wb_new;
          TauR_old = TauR_new;
          TauL_old = TauL_new;
              
          if or(ge(m,m_max),lt(epsilon,1E-8))
              WsR = WsR_old;
              WsL = WsL_old;
              InvWsR = inv(WsR_new);
              InvWsL = inv(WsL_new);
              InvWb = inv(Wb_new);
              break;
          end
      end
          
      % == Compute total transmission using Caroli formula and method === 
      Gamma_R = 1i*((TauR*InvWsR*TauL)-(TauR*InvWsR*TauL)'); 
      Gamma_L = 1i*((TauL*InvWsL*TauR)-(TauL*InvWsL*TauR)');
      Xi_negf = real(trace(Gamma_L*InvWb'*Gamma_R*InvWb)); % total transmittance 
        
      % === Compute transition amplitudes and transmission coeffs using mode analysis method === 
       
      % [tmp.UR,tmp.VR] = eig(a_long/(2*wvec(nw))*Gamma_R); 
        % this captures the extended modes but not the evanescent modes
      % [tmp.UR_adv,tmp.VR_adv] = eig(a_long/(2*wvec(nw))*Gamma_L); % ditto
      % tmp.dFR = 1E-3 * tmp.UR * tmp.VR * tmp.UR' / max(diag(abs(tmp.VR)));
      % tmp.dInvFR_adv = 1E-3 * tmp.UR_adv * tmp.VR_adv * tmp.UR_adv' / max(diag(abs(tmp.VR_adv)));
      % [tmp_U_plus,tmp_E_plus]         = eig(Gamma_R); 
        % this captures the extended modes but not the evanescent modes
      % [tmp_U_minus_adv,tmp_E_minus_adv] = eig(Gamma_L); % ditto
      % dF_plus = 0*1E-6 * Gamma_R / max(diag(abs(tmp_E_plus)));
      % dInvF_minus_adv = 0*1E-6 * Gamma_L / max(diag(abs(tmp_E_minus_adv)));
      % dInvF_minus = dF_plus;
      % dF_plus_adv = dInvF_minus_adv;
          
      F_plus = -InvWsR*TauL;            % Bloch matrix for right-moving, right-decaying retarded states
      F_plus_adv = -InvWsR'*TauL;       % Bloch matrix for left-moving, right-decaying retarded states
      InvF_minus = (-TauL*InvWsL')';    % Inverse Bloch matrix for left-moving, left-decaying advanced states
      InvF_minus_adv = (-TauL*InvWsL)'; % Inverse Bloch matrix for right-moving, left-decaying advanced states
  
      % F_plus = F_plus + dF_plus;
      % F_plus_adv = F_plus_adv + dF_plus_adv;
      % InvF_minus = InvF_minus + dInvF_minus;
      % InvF_minus_adv = InvF_minus_adv + dInvF_minus_adv;
  
      % D_plus = eye(size(F_plus));
      % D_plus_adv = eye(size(F_plus_adv));
      % D_minus = eye(size(InvF_minus));
      % D_minus_adv = eye(size(InvF_minus_adv));
  
  
      [U_plus,L_plus] = eig(F_plus);                      
        % Eigenvectors and values (phase factors) of retarded right-propagating Bloch matrix
      [U_plus_adv,L_plus_adv] = eig(F_plus_adv);          
        % Eigenvectors and values (phase factors) of advanced right-propagating Bloch matrix
      [U_minus,InvL_minus] = eig(InvF_minus);             
        % Eigenvectors and values (phase factors) of retarded left-propagating Bloch matrix
      [U_minus_adv,InvL_minus_adv] = eig(InvF_minus_adv); 
        % Eigenvectors and values (phase factors) of advanced left-propagating Bloch matrix
         
      Q_plus = imag(log(diag(L_plus))).* ...
               and(gt(diag(abs(L_plus)),0.95),lt(diag(abs(L_plus)),1.05));
        % Extract wave vectors of right-propagating retarded states (between 0 and 2 pi)
      Q_plus_adv = imag(log(diag(L_plus_adv))).* ... 
                   and(gt(diag(abs(L_plus_adv)),0.95),lt(diag(abs(L_plus_adv)),1.05));
        % Extract wave vectors of right-propagating advanced states (between 0 and 2 pi)
      Q_minus = -imag(log(diag(InvL_minus))).* ... 
                and(gt(diag(abs(InvL_minus)),0.95),lt(diag(abs(InvL_minus)),1.05));
        % Extract wave vectors of left-propagating retarded states (between 0 and 2 pi)
      Q_minus_adv = -imag(log(diag(InvL_minus_adv))).* ... 
                    and(gt(diag(abs(InvL_minus_adv)),0.95),lt(diag(abs(InvL_minus_adv)),1.05));
        % Extract wave vectors of left-propagating advanced states (between 0 and 2 pi)

      U_plus      = OrthonormEigenmodes(U_plus,Q_plus);
      U_plus_adv  = OrthonormEigenmodes(U_plus_adv,Q_plus_adv);
      U_minus     = OrthonormEigenmodes(U_minus,Q_minus);
      U_minus_adv = OrthonormEigenmodes(U_minus_adv,Q_minus_adv);
  
      V_plus = real(a_long/(2*w)*U_plus'*Gamma_R*U_plus); 
        % Group velocities of right-propagating retarded states       
      V_plus_adv = -real(a_long/(2*w)*U_plus_adv'*Gamma_R*U_plus_adv);
        % Group velocities of right-propagating advanced states
      V_minus = -real(a_long/(2*w)*U_minus'*Gamma_L*U_minus); 
        % Group velocities of left-propagating retarded states       
      V_minus_adv = real(a_long/(2*w)*U_minus_adv'*Gamma_L*U_minus_adv);
        % Group velocities of left-propagating advanced states

      % [S_plus,D_plus] = eig(V_plus);
      % [S_plus_adv,D_plus_adv] = eig(V_plus_adv);
      % [S_minus,D_minus] = eig(V_plus);
      % [S_minus_adv,D_minus_adv] = eig(V_plus_adv);

      %{
      V_plus = 0*V_plus;
      V_plus_adv = 0*V_plus_adv;
      V_minus = 0*V_minus;
      V_minus = 0*V_minus_adv;

      for n_mode = 1:numel(Q_plus)
          if and(gt(diag(abs(L_plus(n_mode,n_mode))),0.95),lt(diag(abs(L_plus(n_mode,n_mode))),1.05))
              VelxOp = -1i*a_long*HL*exp(-1i*Q_plus(n_mode)) + 1i*a_long*HR*exp(1i*Q_plus(n_mode));
              V_plus(n_mode,n_mode) = real(U_plus(:,n_mode)'*VelxOp*U_plus(:,n_mode)/(2*w));
          end

          if and(gt(diag(abs(L_plus_adv(n_mode,n_mode))),0.95),lt(diag(abs(L_plus_adv(n_mode,n_mode))),1.05))
              VelxOp = -1i*a_long*HL*exp(-1i*Q_plus_adv(n_mode)) + 1i*a_long*HR*exp(1i*Q_plus_adv(n_mode));
              V_plus_adv(n_mode,n_mode) = real(U_plus_adv(:,n_mode)'*VelxOp*U_plus_adv(:,n_mode)/(2*w));
          end

          if and(gt(diag(abs(InvL_minus(n_mode,n_mode))),0.95),lt(diag(abs(InvL_minus(n_mode,n_mode))),1.05))
              VelxOp = -1i*a_long*HL*exp(-1i*Q_minus(n_mode)) + 1i*a_long*HR*exp(1i*Q_minus(n_mode));
              V_minus(n_mode,n_mode) = real(U_minus(:,n_mode)'*VelxOp*U_minus(:,n_mode)/(2*w));
          end

          if and(gt(diag(abs(InvL_minus_adv(n_mode,n_mode))),0.95),lt(diag(abs(InvL_minus_adv(n_mode,n_mode))),1.05))
              VelxOp = -1i*a_long*HL*exp(-1i*Q_minus_adv(n_mode)) + 1i*a_long*HR*exp(1i*Q_minus_adv(n_mode));
              V_minus_adv(n_mode,n_mode) = real(U_minus_adv(:,n_mode)'*VelxOp*U_minus_adv(:,n_mode)/(2*w));
          end
      end
      %}

      %{
      save('-7','DebugData','F_plus','F_plus_adv','InvF_minus','InvF_minus_adv',...
           'U_plus','U_plus_adv','U_minus','U_minus_adv',...
           'Q_plus','Q_plus_adv','Q_minus','Q_minus_adv',...
           'V_plus','V_plus_adv','V_minus','V_minus_adv'...
      );
      %}

      BigG = InvWb;    
      t_L = 2*1i*w*sqrt(V_plus)*inv(U_plus)*BigG*inv(U_minus_adv')*sqrt(V_minus_adv)/a_long;
      T_plus = diag(t_L*t_L');      % Perfect left-to-right transmission (out-going) in lead 
      T_minus_adv = diag(t_L'*t_L); % Perfect left-to-right transmission (in-coming) in lead 
  
      t_R = 2*1i*w*sqrt(V_minus)*inv(U_minus)*BigG*inv(U_plus_adv')*sqrt(V_plus_adv)/a_long;
      T_minus = diag(t_R*t_R');     % Perfect right-to-left transmission (out-going) in lead 
      T_plus_adv = diag(t_R'*t_R);  % Perfect right-to-left transmission (in-coming) in lead 

      % === Store computed variables === 
      LeadPhonon.MatSurfGR = InvWsR; % Right surface Green's function
      LeadPhonon.MatSurfGL = InvWsL; % Left surface Green's function
      LeadPhonon.MatBulkG  = InvWb;  % Bulk Green's function

      LeadPhonon.U_plus      = U_plus;      % Right-propagating eigenmode matrix (retarded)
      LeadPhonon.U_plus_adv  = U_plus_adv;  % Right-propagating eigenmode matrix (advanced)
      LeadPhonon.U_minus     = U_minus;     % Left-propagating eigenmode matrix (retarded)
      LeadPhonon.U_minus_adv = U_minus_adv; % Left-propagating eigenmode matrix (advanced)

      LeadPhonon.V_plus      = V_plus;      % Velocity matrices of right-propagating eigenmodes (retarded)
      LeadPhonon.V_plus_adv  = V_plus_adv;  % Velocity matrices of right-propagating eigenmodes (advanced)
      LeadPhonon.V_minus     = V_minus;     % Velocity matrices of left-propagating eigenmodes (retarded)
      LeadPhonon.V_minus_adv = V_minus_adv; % Velocity matrices of left-propagating eigenmodes (advanced)
        
      LeadPhonon.VecV_plus      = diag(V_plus);      % 1D array of velocity of right-propagating eigenmodes (retarded)
      LeadPhonon.VecV_plus_adv  = diag(V_plus_adv);  % 1D array of velocity of right-propagating eigenmodes (advanced)
      LeadPhonon.VecV_minus     = diag(V_minus);     % 1D array of velocity of left-propagating eigenmodes (retarded)
      LeadPhonon.VecV_minus_adv = diag(V_minus_adv); % 1D array of velocity of left-propagating eigenmodes (advanced)
        
      LeadPhonon.VecQ_plus      = Q_plus;      % between 0 ans 2 pi
      LeadPhonon.VecQ_plus_adv  = Q_plus_adv;  % between 0 ans 2 pi
      LeadPhonon.VecQ_minus     = Q_minus;     % between 0 ans 2 pi
      LeadPhonon.VecQ_minus_adv = Q_minus_adv; % between 0 ans 2 pi

      LeadPhonon.VecT_plus      = T_plus;
      LeadPhonon.VecT_plus_adv  = T_plus_adv;
      LeadPhonon.VecT_minus     = T_minus;
      LeadPhonon.VecT_minus_adv = T_minus_adv;

      LeadPhonon.Xi_mode = real(sum(T_plus));
      LeadPhonon.Xi_negf = Xi_negf;
      LeadPhonon.VecQ_tran1 = 0.0*diag(V_plus);
      LeadPhonon.VecQ_tran2 = 0.0*diag(V_plus);
  end   

  % === Allocate variables for mapping to primitive lattice Brillouin zone === 
  LeadPhonon.VecQx_plus = [];
  LeadPhonon.VecQy_plus = [];
  LeadPhonon.VecQz_plus = [];
  LeadPhonon.VecVelx_plus = [];
  LeadPhonon.VecVely_plus = [];
  LeadPhonon.VecVelz_plus = [];
  LeadPhonon.VecDw_plus = [];
  LeadPhonon.VecB_plus  = [];

  LeadPhonon.VecQx_plus_adv = [];
  LeadPhonon.VecQy_plus_adv = [];
  LeadPhonon.VecQz_plus_adv = [];
  LeadPhonon.VecVelx_plus_adv = [];
  LeadPhonon.VecVely_plus_adv = [];
  LeadPhonon.VecVelz_plus_adv = [];
  LeadPhonon.VecDw_plus_adv = [];
  LeadPhonon.VecB_plus_adv  = [];

  LeadPhonon.VecQx_minus = [];
  LeadPhonon.VecQy_minus = [];
  LeadPhonon.VecQz_minus = [];
  LeadPhonon.VecVelx_minus = [];
  LeadPhonon.VecVely_minus = [];
  LeadPhonon.VecVelz_minus = [];
  LeadPhonon.VecDw_minus = [];
  LeadPhonon.VecB_minus  = [];

  LeadPhonon.VecQx_minus_adv = [];
  LeadPhonon.VecQy_minus_adv = [];
  LeadPhonon.VecQz_minus_adv = [];
  LeadPhonon.VecVelx_minus_adv = [];
  LeadPhonon.VecVely_minus_adv = [];
  LeadPhonon.VecVelz_minus_adv = [];
  LeadPhonon.VecDw_minus_adv = [];
  LeadPhonon.VecB_minus_adv  = [];
end 


% ================================================================
function U_out = OrthonormEigenmodes(U_in,Q_in)
% This function orthonormalizes the degenerate extended/propagating eigenmodes from the Bloch matrix
% using the Gram-Schmidt procedure
  U_out = U_in;

  nqmax = numel(Q_in); % number of modes

  for nq = 1:nqmax
      U_out(:,nq) = U_out(:,nq)/norm(U_out(:,nq));
  end


  FlagQ = zeros(size(Q_in));
  FlagQ(eq(abs(Q_in),0)) = 1; % set evanescent mode check flag to 1

  for nq1 = 1:nqmax      
      if eq(FlagQ(nq1),0) % i.e. mode has not been checked and is an extended mode
          FlagQ(nq1) = 1; % set check flag to 1
          q1_long = Q_in(nq1);  % longitudinal wave vector of mode 1
          U_long = U_in(:,nq1); % eigenvector of mode 1
          nqvec_long = nq1; % index of mode 1

          for nq2 = (nq1+1):nqmax % loop over remaining modes
              if eq(FlagQ(nq2),0)
                  q2_long = Q_in(nq2); % longitudinal wave vector of mode 2
                  epsilon = abs((q1_long - q2_long)/(0.5*q1_long+0.5*q2_long)); % relative diff. of long. wave vector

                  if lt(epsilon,1E-2) % i.e. mode 2 is numerically degenerate with mode 1
                      FlagQ(nq2) = 1;
                      U_long = [U_long U_in(:,nq2)]; % store eigenvector of mode 2
                      nqvec_long = [nqvec_long nq2]; % store index of mode 2
                  end
              end              
          end

          if gt(numel(nqvec_long),1)
              U_long = GramSchmidtOrthogonalize(U_long);
          end

          for n = 1:numel(nqvec_long)
              nq = nqvec_long(n);
              U_out(:,nq) = U_long(:,n);          
          end
      end 
  end  
end


% ================================================================
function U_out = GramSchmidtOrthogonalize(U_in)
  U_out = U_in;
  n_subspace = size(U_out,2);

  for n1 = 1:n_subspace
      v = U_out(:,n1);

      for n2 = 1:(n1-1)
          u = U_out(:,n2);
          u = u/norm(u);
          v = v - (u'*v)*u;
      end

      U_out(:,n1) = v/norm(v);
  end
end


