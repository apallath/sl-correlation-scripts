% Akash Pallath <akash.pallath@btech15.iitgn.ac.in>
% Created: January 2019
% Last modified: May 2019
% 
% Plot disp-disp correlations for topmost solid atom layer

clear;

%% Read input parameters
input_params = containers.Map('KeyType','char','ValueType','any');

input_file = fopen('all_inputs.input');
file_contents = textscan(input_file,'%s','delimiter','\n');
fclose(input_file);
paramlist = regexp(file_contents{1,1},'\s','split','once');
sz = size(paramlist);

for i = 1:sz(1)
    pair = paramlist{i,1};
    key = char(pair(1));
    val = char(pair(2));
    input_params(key) = val;
end
%% Parameters

% number of snapshots
nsteps = str2num(input_params('snapshots'));

% number of atoms in the region
natoms = str2num(input_params('nsolidtopatoms'));

% box limits
Lx_min = str2num(input_params('x_min'));
Lx_max = str2num(input_params('x_max'));
Lx = Lx_max - Lx_min;

Ly_min = str2num(input_params('y_min'));
Ly_max = str2num(input_params('y_max'));
Ly = Ly_max - Ly_min;

Lz_min = str2num(input_params('z_min'));
Lz_max = str2num(input_params('z_max'));
Lz = Lz_max - Lz_min;

ival = 0.005 * max([Lx, Ly, Lz]);

load('solid_top_disp_disp_corr_mean.mat');
load('solid_top_rvals.mat');
load('solid_top_initpos.mat');

%% read MSDs for normalization
msdfile = fopen(input_params('msd_file'));
fgetl(msdfile);
msds = cell2mat(textscan(msdfile, '%f %f %f %f\n', nsteps));
msdms = msds(:,4);
msdmean = mean(msdms);

% adjust binning
curg = circshift(curg,1);

%% plot figure
figure();
plot(rvals,curg);
hold on;
title('SC 100');
xlabel('r (LJ)');
ylabel('<dr(r)dr(0)> (arb)');

savefig('solid_top_disp_disp_corr_mean.fig');
saveas(gcf,'solid_top_disp_disp_corr_mean','epsc');
saveas(gcf,'solid_top_disp_disp_corr_mean','png');

%% find atom closest to center
mindist = Inf;
cdist = [15,15];
center = 0;
for j = 1:natoms
    rm = initlist(j,1:2);
    dist = norm(rm - cdist);
    if  dist < mindist
        center = j;
        mindist = dist;
    end        
end

%% plot surface

curlist = initlist;

for j = 1:natoms
    rm = curlist(j,1:2);
    rn = curlist(center,1:2);
    deltax = abs(rm(1) - rn(1));
    deltay = abs(rm(2) - rn(2));
    
    r = sqrt(min(deltax,Lx-deltax)^2 + min(deltay,Ly-deltay)^2);
 
    bin = floor(r/ival) + 1;
    curlist(j, 3) = curg(bin);
end

x = curlist(:,1);
y = curlist(:,2);
c = curlist(:,3);

c = c/nsteps;

figure();
colormap('default');
scatter(x,y,100,c,'filled');
hold on;
scatter(x(center),y(center),[],'red','filled');
title(sprintf('SC 100, k = %s',kval));

savefig('solid_top_disp_disp_corr_surface.fig');
saveas(gcf,'solid_top_disp_disp_corr_surface','epsc');
saveas(gcf,'solid_top_disp_disp_corr_surface','png');
