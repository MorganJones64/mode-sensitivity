%% Example for computing mode sensitivity quantities for turbulent channel flow from https://arxiv.org/abs/2410.20802
% First run the file DNSchannelformat.m to get the data, then run this script
clear all
addpath('./functions/')
load('channelflowData.mat') %data obtained from DNSchannelformat.m
%% Compute Modal Representations for mode sensitivity
r=nt-1;
% choose method of mode hierarchy
sortMethod = 'amplitude'; % or 'frequency'
% choose method for approximating dmd modes
dataMethod = 'fourier'; % or 'none'
% choose mode indices based on hierarchy
ind_g = [1:20]; %[1,2],  [11,12],  [13,14]
ind_f = setdiff(1:r, ind_g); % exclude the modes that are of g 
[uf, vf, ug, vg] = dmdModel(uMat, vMat, dt, r, ind_f, ind_g, sortMethod, dataMethod);
%% parameters and domain for passive tracers
clear sigmaS
xMinROI = 0;
xMaxROI = 2*pi;
yMinROI = -1;
yMaxROI = 1;
xstart = 256; %location of x=0
xend = 769;   %location of x=2*pi
ROInx = size(x(:,xstart:xend),2);
ROIny = size(xMat(:,xstart:xend),1);
intMethod = 'RK4';
umeanprof = mean(uMat(:,end,:),3)';
vmeanprof = mean(vMat(:,end,:),3)';
yVec=y';
xVec=x';
%% Set function handles for extrapolation of velocity outside ROI
%Forward FTLE
tVec = (0:1:nt-1)*dt;
uadv = umeanprof;

%aproximate perturbations of epsilon*g(x,t) using taylors frozen vortex
%hypothesis
ugtaylor = reshape(ug(:,xend,:),[ny,nt]);
vgtaylor = reshape(vg(:,xend,:),[ny,nt]);
ugF = interpolantG(ugtaylor,yVec,tVec);
vgF = interpolantG(vgtaylor,yVec,tVec);

%approximate an advection speed based on the time averaged flow at
%different y locations
uadvF = griddedInterpolant(yVec,umeanprof');

% Extrapolation functions for FTLE
ufExtrap = @(x, y, t, u) ufExtrapolate(x, y, t, u, umeanprof',yVec);
vfExtrap = @(x, y, t, v) vfExtrapolate(x, y, t, v, vmeanprof',yVec);
ugExtrap = @(x, y, t, u) extrapolateH(x, y, t, u, ugF,uadvF,xMaxROI);
vgExtrap = @(x, y, t, v) extrapolateH(x, y, t, v, vgF,uadvF,xMaxROI);
%% Set Integration Length and Number of Snapshots
%Forwards FTLE
%inteval of integration from t0 to t
tLength = 150; %round(6*pnt); 
tStep = 1;

%Backwards FTLE
%inteval of integration from t0 to t
%tLength = -100; 
%tStep = -1;

% First Frame to compute FTLE Field
fstart = 1;
% Last Frame to compute FTLE Field
fend = 1; 
% Frame Increment0
finc = 1;
% Frame loop vector
fLoop = fstart:finc:fend;
%% Preallocate space for Lagrangian Quantities
sigma_ftle = zeros([ROInx, ROIny, length(fLoop)]);
cseIntegral = zeros([ROInx, ROIny, length(fLoop)]);
deltaInfty = zeros([ROInx, ROIny, length(fLoop)]);
save('MSdata.mat','sigma_ftle','cseIntegral','deltaInfty','dt','tLength','xMaxROI','xMinROI','yMaxROI','yMinROI','fstart','fend','fLoop','intMethod')
%% Compute Mode Sensitivity Quantities. Runs for a Long time
for t0 = fstart
    [sigma_ftle(:, :, t0-fLoop(1)+1), cseIntegral(:, :, t0-fLoop(1)+1), deltaInfty(:, :, t0-fLoop(1)+1), xPos, yPos] = modeSensitivity(uf, vf, ...
        ug, vg,...
        xVec, yVec, ...
        t0, tLength, tStep, fstart, dt, ...
        xMinROI, xMaxROI, yMinROI, yMaxROI, ...
        ROInx, ROIny, intMethod, ...
        'extrap',true,'ufExtrap', ufExtrap, 'vfExtrap', vfExtrap, ...
        'ugExtrap', ugExtrap, 'vgExtrap', vgExtrap);
    save('MSdata.mat','sigma_ftle','cseIntegral','deltaInfty','xPos','yPos',"-append")
end
%% Compute Lagrangian Fields, FTLE, MS, and Zeta
FTLE = (1/abs(tLength*dt))*log(sigma_ftle);
MS = (deltaInfty.*cseIntegral).^2;
MS_scaled = log(MS)/abs(2*tLength*dt);
zeta = (1/abs(tLength*dt))*log(deltaInfty.*cseIntegral./sigma_ftle);
x=xPos(:,:,1);
y=yPos(:,:,1);
%% Plot the Finite-Time Lyapunov Exponent Field FTLE
ep=0.01;
width =1200; 
height = 400;
umin = 0;
umax = 1.5;
figure
set(gcf,'Position',[100 300 width height])
contourf(x,y,FTLE(:,:),200,'LineStyle','none')
yticks([-1, 0, 1])
colormap(jet)
clim([umin,umax])
c=colorbar;
c.Ticks =linspace(umin,umax,3);
c.FontSize = 20;
c.FontName = 'Times New Roman';
c.LineWidth = 1.2;
c.TickLength = .02;
xlabel('$x$','Interpreter','latex','fontsize',20);
ylabel('$y$','Interpreter','latex','fontsize',20);
set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
set(gcf,'color','w')
ax2 = gca; axis equal; box on; ax2.FontName = 'Times New Roman'; ax2.FontSize = 20;
%% Plot the Scaled Mode Sensitivity MS_scaled
ep=0.01;
width =1200; 
height = 400;
umin = -1.5;
umax = 0;
figure
set(gcf,'Position',[100 300 width height])
contourf(x,y,MS_scaled(:,:),200,'LineStyle','none')
yticks([-1, 0, 1])
colormap(jet)
clim([umin,umax])
c=colorbar;
c.Ticks =linspace(umin,umax,3);
c.FontSize =axislength;
c.FontName = 'Times New Roman';
c.LineWidth = 1.2;
c.TickLength = .02;
xlabel('$x$','Interpreter','latex','fontsize',20);
ylabel('$y$','Interpreter','latex','fontsize',20);
set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
set(gcf,'color','w')
ax2 = gca; axis equal; box on; ax2.FontName = 'Times New Roman'; ax2.FontSize = 20;
%% Plot The Lagrangian Response Zeta
ep=0.01;
width =1200; 
height = 400;
umin = -1.8;
umax = -1;
figure
set(gcf,'Position',[100 300 width height])
contourf(x,y,zeta(:,:),200,'LineStyle','none')
yticks([-1, 0, 1])
colormap(jet)
clim([umin,umax])
c=colorbar;
c.Ticks =linspace(umin,umax,3);
c.FontSize =axislength;
c.FontName = 'Times New Roman';
c.LineWidth = 1.2;
c.TickLength = .02;
xlabel('$x$','Interpreter','latex','fontsize',20);
ylabel('$y$','Interpreter','latex','fontsize',20);
set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
set(gcf,'color','w')
ax2 = gca; axis equal; box on; ax2.FontName = 'Times New Roman'; ax2.FontSize = 20;
%%
function u = ufExtrapolate(x, y, t, u, umean,yVec)
    F = griddedInterpolant(yVec, umean,'makima','nearest'); 
    u = F(y);
end

function v = vfExtrapolate(x, y, t, v, vmean,yVec)
    F = griddedInterpolant(yVec, vmean,'makima','nearest'); 
    u = F(y);
end

function F = interpolantG(uvec,yVec,tVec)
    [yi,ti] = ndgrid(yVec, tVec);
    F = griddedInterpolant(yi, ti, uvec,'makima','nearest'); 
end

function u = extrapolateH(x, y, t, u, F, uadvF,xwall)
    if isempty(x)
        uadv=nan;
    else
        uadv = uadvF(y);
    end
    u = F(y,t-(x-xwall)./uadv);
end
