function [scen_hat,report] = perturb_scenario_csi(scen_true,epsilon,seed)
%PERTURB_SCENARIO_CSI Generate a paired imperfect-CSI scenario.
% Model: Hhat=sqrt(1-epsilon^2)H+epsilon E, where E is circular complex
% Gaussian and is rescaled link-wise to match the energy of H. Reusing the
% same seed across epsilon values gives nested paired perturbation directions.

if nargin<3 || isempty(seed); seed=1; end
if ~isscalar(epsilon) || epsilon<0 || epsilon>=1
    error('epsilon must satisfy 0 <= epsilon < 1.');
end
scen_hat=scen_true;
if epsilon==0
    report=struct('epsilon',0,'nmse_Hd_user',0,'nmse_Hr_user',0, ...
        'nmse_Hd_tar',0,'nmse_Hr_tar',0,'nmse_G',0);
    return;
end
rng(seed,'twister');
[scen_hat.Hd_user,n1]=perturb_columns(scen_true.Hd_user,epsilon);
[scen_hat.Hr_user,n2]=perturb_columns(scen_true.Hr_user,epsilon);
[scen_hat.Hd_tar,n3]=perturb_columns(scen_true.Hd_tar,epsilon);
[scen_hat.Hr_tar,n4]=perturb_columns(scen_true.Hr_tar,epsilon);
[scen_hat.G,n5]=perturb_matrix(scen_true.G,epsilon);
report=struct('epsilon',epsilon,'nmse_Hd_user',n1,'nmse_Hr_user',n2, ...
    'nmse_Hd_tar',n3,'nmse_Hr_tar',n4,'nmse_G',n5);
end

function [Y,nmse]=perturb_columns(X,e)
Y=zeros(size(X),'like',X);
for jj=1:size(X,2)
    x=X(:,jj);
    z=(randn(size(x))+1i*randn(size(x)))/sqrt(2);
    nz=norm(z); nx=norm(x);
    if nz>0 && nx>0; z=z*(nx/nz); else; z=zeros(size(x),'like',x); end
    Y(:,jj)=sqrt(1-e^2)*x+e*z;
end
nmse=norm(Y-X,'fro')^2/max(norm(X,'fro')^2,realmin);
end

function [Y,nmse]=perturb_matrix(X,e)
z=(randn(size(X))+1i*randn(size(X)))/sqrt(2);
nz=norm(z,'fro'); nx=norm(X,'fro');
if nz>0 && nx>0; z=z*(nx/nz); else; z=zeros(size(X),'like',X); end
Y=sqrt(1-e^2)*X+e*z;
nmse=norm(Y-X,'fro')^2/max(norm(X,'fro')^2,realmin);
end
