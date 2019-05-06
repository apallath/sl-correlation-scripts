% Akash Pallath <akash.pallath@iitgn.ac.in>
% Created: April 2019
% Last modified: May 2019
% 
% Calculate Z-Z correlations for the Zmin of liquid atom over each solid atom layer
%
% Important: Atoms in the dump file must be sorted by id, so that the
% sequence is the same at each timestep

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
% nsteps = str2num(input_params('snapshots'));
nsteps = 1;

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

ival = 0.005 * max([Lx, Ly]);

load(sprintf('zz_corr_mean.mat'));
load(sprintf('zz_corr_rvals.mat'));

%% Calculate init positions
inittopposfile = fopen(input_params('solid_init_top_positions_file'));
xylist = zeros(natoms,2);  %%initlist

for j = 1:9
    fgetl(inittopposfile); %read the line
end
vals = cell2mat(textscan(inittopposfile, '%f %f %f %f %f\n', natoms));
for j = 1:natoms
    xylist(j, :) = xylist(j, :) + [Lx_min + vals(j,3) * Lx, Ly_min + vals(j,4)*Ly];
end

%% Adjust binning
%curg = circshift(curg,1);
curg(1)=0;

%% plot figure
figure();
plot(rvals,curg);
hold on;
title(sprintf('SC 100'));
xlabel('r');
ylabel('<dz(0)dz(r)>');

savefig(sprintf('zz_corr_mean.fig'));
saveas(gcf,sprintf('zz_corr_mean'),'epsc');
saveas(gcf,sprintf('zz_corr_mean'),'png');

%% find atom closest to center
mindist = Inf;
cdist = [15,15];  
center = 0;
for j = 1:natoms
    rm = xylist(j,1:2);
    dist = norm(rm - cdist);
    if  dist < mindist
        center = j;
        mindist = dist;
    end        
end

%% plot surface

for j = 1:natoms
    rm = xylist(j,1:2);
    rn = xylist(center,1:2);
    deltax = abs(rm(1) - rn(1));
    deltay = abs(rm(2) - rn(2));
   
    r = sqrt(min(deltax,Lx-deltax)^2 + min(deltay,Ly-deltay)^2);
 
    bin = floor(r/ival) + 1;
    zlist(j,3) = curg(bin);
end

x = xylist(:,1);
y = xylist(:,2);
c = zlist(:,3);

figure();
colormap('default');
scatter(x,y,100,c,'filled');
hold on;
scatter(x(center),y(center),[],'red','filled');
title(sprintf('SC 100'));

savefig(sprintf('zz_corr_surface.fig'));
saveas(gcf,sprintf('zz_surface.eps'),'epsc');
saveas(gcf,sprintf('zz_surface.png'),'png');