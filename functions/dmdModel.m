
function [uf, vf, ug, vg] = dmd(uMat, vMat, dt, r, indexsetf, indexsetg, sortMethod, dataMethod)
% Function for computing dmd modes and reconstructions for f(x,t) and
% epsilon*g(x,t) for mode sensitivity. the reconstruction f(x,t) is
% considered the primary fluid structure while g(x,t) is the mode
% perturbation. There is much flexibility towards how these flow fields are
% constructed given the data and particular analysis.

% Required Inputs
% uMat: Matrix of x-component of velocity [m x n x p] Meshgrid style.
% vMat: Matrix of y-component of velocity [m x n x p] Meshgrid style.
% NOTE: m, n, and p are grid points in the y, x, and t-axes respectively.
% dt: time between time steps [scalar]
% r: rank appoximation/mode truncation number [scalar] (index)
% indexsetf: subset of modes used for the primary fluid structure f(x,t).
% NOTE: Index numbers must be less than r.
% indexsetg: subset of modes used for the mode perturbation epsilon*g(x,t). [scalar/vector]
% NOTE: Index numbers must be less than r.
% sortMethod:
% 'frequency' - mode hierarchy is ordered by increasing oscillation frequency
% 'amplitude' - mode hierarchy is ordered by decreasing amplitude 
% 'decay' - mode hierarchy is ordered by decreasing growth rate. 
% dataMethod:
% 'fourier' - data is subtracted by time-average flow before computing
% modes
% 'none' - standard dmd method is used and data is kept as is

% Outputs
% uf:  x-component of f(x,t) velocity [m x n x p]
% ug: x-component of epsilon*g(x,t) velocity [m x n x p]
% vf: y-component of f(x,t) velocity [m x n x p]
% vg: y-component of epsilon*g(x,t) velocity [m x n x p]

arguments
    uMat double
    vMat double
    dt(1,1) double
    r(1,1) double
    indexsetf(1,:) double
    indexsetg(1,:) double
    sortMethod(1,:) char ...
        {mustBeMember(sortMethod, {'frequency', 'amplitude', 'decay'})} = 'frequency'
    dataMethod(1,:) char ...
        {mustBeMember(dataMethod, {'fourier','dmd'})} = 'dmd'
end
% Reshape the velocity matrix for DMD
[ny,nx,nt] = size(uMat);
U = [reshape(uMat,[ny*nx,nt]); ...
     reshape(vMat,[ny*nx,nt])];
if strcmp(dataMethod, 'fourier')
    Umean = mean(U,2);
    U = U - Umean; % subtract data by the mean flow
end
tVec = 0:dt:nt*dt-dt;

% Extract subset of snapshots needed for DMD
U1 = U(:,1:end-1);
U2 = U(:,2:end);

% Compute SVD of U1
[Psi,Sigma,Phi]=svds(U1,min(size(U1)));
% take the r-rank approximation to project the pod modes onto the data
Psi = Psi(:,1:r);
Sigma = Sigma(1:r,1:r);
Phi = Phi(:,1:r);
% Compute S-matrix representing SVD-projected dynamics
S = Psi'*U2*Phi/Sigma;
% S is the linear map from one time step to another
% Eigenvalue decomposition of S
[Q,Lambda] = eig(S);

%Compute dynamic modes
Zeta = (U2*Phi/Sigma)*Q;
%Compute amplitude coefficients for dynamic modes
B = Zeta\U(:,1);
for ti = 2:length(tVec)
    B(:,ti) = Lambda*B(:,ti-1);
end

for i = 1:nt-1
    amp(i) = norm(B(i,:),2);
end
lambda = diag(Lambda);
omega = imag(log(lambda)/dt);

if strcmp(sortMethod, 'frequency')
    %oscillation frequency hierarchy
    [~,inds]=sort(abs(omega),'ascend');
elseif strcmp(sortMethod,'decay')
    %growth/decay rate hierarchy
    [~,inds]=sort(abs(lambda),'descend');
elseif strcmp(sortMethod, 'amplitude')
    %amplitude hierarchy
    [~,inds]=sort(amp,'descend');
end

figure
hold on
plot(omega(inds(indexsetf)),amp(inds(indexsetf)),'ro','LineWidth',1,'MarkerSize',4)
plot(omega(inds(indexsetg)),amp(inds(indexsetg)),'bo','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','b')
xlabel('imag($\omega_k$)','interpreter','latex','fontsize',14)
ylabel('$||b_k||_2$','interpreter','latex','fontsize',14)
ax=gca; ax.FontName = 'Times New Roman'; ax.FontSize = 14;
set(gcf,'color','white')
legend('modes of f(x,t)','modes of g(x,t)')

figure
plot(cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k--')
hold on
plot(real(lambda(inds(indexsetf))),imag(lambda(inds(indexsetf))),'ro','LineWidth',1,'MarkerSize',4)
plot(real(lambda(inds(indexsetg))),imag(lambda(inds(indexsetg))),'bo','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','b')
xlabel('real($\lambda_i$)','interpreter','latex','fontsize',14)
ylabel('imag($\lambda_i$)','interpreter','latex','fontsize',14)
ax=gca; ax.FontName = 'Times New Roman'; ax.FontSize = 14;
set(gcf,'color','white')
legend('unit circle','modes of f(x,t)','modes of g(x,t)')
axis equal

nu = 1:(ny*nx);
nv = (ny*nx)+1:2*(ny*nx);

%reconstruction of dataset f(x,t) and modes g(x,t)
Uf = real(Zeta(:,inds(indexsetf))*B(inds(indexsetf),:));
Ug = real(Zeta(:,inds(indexsetg))*B(inds(indexsetg),:));

if strcmp(dataMethod, 'fourier')
    umean=reshape(Umean(nu,:),[ny nx]);
    vmean=reshape(Umean(nv,:),[ny nx]);
    uf = reshape(Uf(nu,:),[ny nx nt]) + umean;
    vf = reshape(Uf(nv,:),[ny nx nt]) + vmean;
else
    uf = reshape(Uf(nu,:),[ny nx nt]);
    vf = reshape(Uf(nv,:),[ny nx nt]);
end

ug = reshape(Ug(nu,:),[ny nx nt]);
vg = reshape(Ug(nv,:),[ny nx nt]);

end