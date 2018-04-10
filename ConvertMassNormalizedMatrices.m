function [Left, Center, Right] = ConvertMassNormalizedMatrices(LeftLeadParam,CenterParam,RightLeadParam)

Lyr = CenterParam.Lyr;
nmax = length(Lyr);

for n = 1:1:nmax % we set up the inverse square root mass matrix for mass-normalizing the force constant matrices
    [V,D] = eig([Lyr(n).MatM]);
    Lyr(n).InvSqrtM = V*diag(1./sqrt(diag(D)))*V';    
end


[V,D] = eig([LeftLeadParam.MatM]);
Left.InvSqrtM = V*diag(1./sqrt(diag(D)))*V';
Left.MatHL = Left.InvSqrtM * LeftLeadParam.MatKL * Left.InvSqrtM;
Left.MatHC = Left.InvSqrtM * LeftLeadParam.MatKC * Left.InvSqrtM;
Left.MatHR = Left.InvSqrtM * LeftLeadParam.MatKR * Left.InvSqrtM;

Left.TransCell = LeftLeadParam.TransCell;
Left.a_long = LeftLeadParam.a_long; % longitudinal lattice spacing
Left.n_tran1 = LeftLeadParam.n_tran1; % number of Fourier components (transverse 1)
Left.n_tran2 = LeftLeadParam.n_tran2; % number of Fourier components (transverse 2)
Left.rvec_long = LeftLeadParam.rvec_long;   % lattice vector (longitudinal)
Left.rvec_tran1 = LeftLeadParam.rvec_tran1; % lattice vector (transverse 1)
Left.rvec_tran2 = LeftLeadParam.rvec_tran2; % lattice vector (transverse 2)
Left.gvec_long  = cross(Left.rvec_tran1,Left.rvec_tran2) ... % reciprocal lattice vector (longitudinal)
                  /(Left.rvec_long*cross(Left.rvec_tran1,Left.rvec_tran2)');
Left.gvec_tran1 = cross(Left.rvec_tran2,Left.rvec_long) ...  % reciprocal lattice vector (transverse 1)
                  /(Left.rvec_tran1*cross(Left.rvec_tran2,Left.rvec_long)');
Left.gvec_tran2 = cross(Left.rvec_long,Left.rvec_tran1) ...  % reciprocal lattice vector (transverse 2)
                  /(Left.rvec_tran2*cross(Left.rvec_long,Left.rvec_tran1)');


[V,D] = eig([RightLeadParam.MatM]);
Right.InvSqrtM = V*diag(1./sqrt(diag(D)))*V';
Right.MatHL = Right.InvSqrtM * RightLeadParam.MatKL * Right.InvSqrtM;
Right.MatHC = Right.InvSqrtM * RightLeadParam.MatKC * Right.InvSqrtM;
Right.MatHR = Right.InvSqrtM * RightLeadParam.MatKR * Right.InvSqrtM;

Right.TransCell = RightLeadParam.TransCell;
Right.a_long = RightLeadParam.a_long; % longitudinal lattice spacing
Right.n_tran1 = RightLeadParam.n_tran1; % number of Fourier components (transverse 1)
Right.n_tran2 = RightLeadParam.n_tran2; % number of Fourier components (transverse 2)
Right.rvec_long = RightLeadParam.rvec_long;   % lattice vector (longitudinal)
Right.rvec_tran1 = RightLeadParam.rvec_tran1; % lattice vector (transverse 1)
Right.rvec_tran2 = RightLeadParam.rvec_tran2; % lattice vector (transverse 2)
Right.gvec_long  = cross(Right.rvec_tran1,Right.rvec_tran2) ... % reciprocal lattice vector (longitudinal)
                   /(Right.rvec_long*cross(Right.rvec_tran1,Right.rvec_tran2)');
Right.gvec_tran1 = cross(Right.rvec_tran2,Right.rvec_long) ...  % reciprocal lattice vector (transverse 1) 
                   /(Right.rvec_tran1*cross(Right.rvec_tran2,Right.rvec_long)');
Right.gvec_tran2 = cross(Right.rvec_long,Right.rvec_tran1) ...  % reciprocal lattice vector (transverse 2) 
                   /(Right.rvec_tran2*cross(Right.rvec_long,Right.rvec_tran1)');


if eq(nmax,1)
    Lyr(1).MatHC = Lyr(1).InvSqrtM * Lyr(1).MatKC * Lyr(1).InvSqrtM;
else
    for n = 1:1:nmax
        Lyr(n).MatHC = Lyr(n).InvSqrtM * Lyr(n).MatKC * Lyr(n).InvSqrtM;
        if eq(n,1)
            Lyr(n).MatHR = Lyr(n).InvSqrtM * Lyr(n).MatKR * Lyr(n+1).InvSqrtM;
        elseif eq(n,nmax)
            Lyr(n).MatHL = Lyr(n).InvSqrtM * Lyr(n).MatKL * Lyr(n-1).InvSqrtM;
        else
            Lyr(n).MatHL = Lyr(n).InvSqrtM * Lyr(n).MatKL * Lyr(n-1).InvSqrtM;
            Lyr(n).MatHR = Lyr(n).InvSqrtM * Lyr(n).MatKR * Lyr(n+1).InvSqrtM;
        end
    end
end


HCL = Lyr(1).InvSqrtM * CenterParam.MatKCL * Left.InvSqrtM;
HLC = HCL';
HCR = Lyr(nmax).InvSqrtM * CenterParam.MatKCR * Right.InvSqrtM;
HRC = HCR';

Lyr = rmfield(Lyr,'MatM');
Lyr = rmfield(Lyr,'InvSqrtM');
Lyr = rmfield(Lyr,'MatKC');
Lyr = rmfield(Lyr,'MatKL');
Lyr = rmfield(Lyr,'MatKR');
Left  = rmfield(Left,'InvSqrtM');
Right = rmfield(Right,'InvSqrtM');

Center.Lyr = Lyr;
Center.HCL = HCL;
Center.HCR = HCR;
Center.HLC = HLC;
Center.HRC = HRC;



