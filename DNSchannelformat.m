%% Set File Path
clear all
close all
addpath('./turbmat-master/')
%% Set Data Parameters
authkey = 'edu.jhu.meneveau-hiSmxkae'; % request authkey by emailing turbulence@jhu.edu
dataset = 'channel';
spatialopt = 'Lag4'; % 4th order Lagrangian interpolation in space;
tempopt = 'PCHIP'; % Piecewise cubic Hermit interpolation in time;

dtDNS = 0.0065; % DNS timestep
dt = 3*dtDNS; %sub-sampled timestep for MS calculations
nt=565;
tstart = 0;
tend = dt*nt-dt;
tvec = tstart:dt:tend; %time vector for snapshots

dx = 8.0*pi/2048;
dy = 0.025;

% Set domain size and position
xmin = -pi;
xmax = 3*pi;
ymin = -1.0;
ymax = 1.0;
nx = int32( (xmax - xmin) / dx); %total number of streamwise points
ny = 2*int32( (ymax - ymin) / dy); %total number of wall-normal points
zoff = 1.5*pi; %plane location in the z-direction
npoints = nx*ny; %total number of data points per snapshot

% Create plane surface
x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);
[xMat, yMat] = meshgrid(x, y);
points(1,:) = xMat(:)';
points(2,:) = yMat(:)';
points(3,:) = zoff; 
%% Collect Snapshots. Runs for a Long Time.
i=1;
for time = tvec
    % Get the velocity at each point
    fprintf('\nRequesting velocity at (%ix%i) points for velocity contour plot t=%0.2f...\n',nx,ny,time);
    result3 = getVelocity(authkey, dataset, time, spatialopt, tempopt, npoints, points);
    % Calculate velocity
    u = result3(1,:);
    v = result3(2,:);
    uMat(:,:,i) = reshape(u, ny, nx);
    vMat(:,:,i) = reshape(v, ny, nx);
    i=i+1;
end

save("channelflowData.mat",'uMat','vMat','nx','ny','nt','x','y','xMat','yMat','tstart','tend','dt','dx','dy')