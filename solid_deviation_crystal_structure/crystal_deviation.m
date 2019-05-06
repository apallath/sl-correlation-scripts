% Akash Pallath <akash.pallath@iitgn.ac.in>
% Created: February 2019
% Last modified: February 2019
% 
% Calculate crystal structure deviation from LAMMPS trajectories
%

% REQUIRED INPUT
kvals = [0.1, 1, 10, 100, 1000];

%% Read input parameters
input_params = containers.Map('KeyType','char','ValueType','any');

input_file = fopen('crystal_deviation.input');
file_contents = textscan(input_file,'%s','delimiter','\n');
fclose(input_file);
paramlist = regexp(file_contents{1,1},'\s','split','once');
sz = size(paramlist);

for i = 1:sz[1]
    pair = paramlist{i,1};
    key = char(pair(1));
    val = char(pair(2));
    input_params(key) = val;
end
%% Parameters

% number of steps
nsteps = str2num(input_params('snapshots'));

% number of atoms in the region
natoms = str2num(input_params('natoms'));

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

%% Calculate

mdevs = {};

kstrings = string(split(num2str(kvals)))';

for k = kstrings

%% Open the files containing the dumped data
eqfile = fopen(sprintf(input_params('equil_dump_file_pattern'), k));
prodfile = fopen(sprintf(input_params('prod_dump_file_pattern'), k));

% positions
origlist = zeros(natoms,3);
meanlist = zeros(natoms,1);

%% Calculate original positions
for j = 1:9
    fgetl(eqfile); %read the line
end
vals = cell2mat(textscan(eqfile, '%f %f %f %f %f\n', natoms));
for j = 1:natoms
    origlist(vals(j,1), :) = [Lx_min + vals(j,3) * Lx, Ly_min + vals(j,4)*Ly, Lz_min + vals(j,5)*Lz];
end

%% Calculate mean deviated positions, taking movement across periodic boundary conditions into account
for i = 1:nsteps    
    for j = 1:9
        fgetl(prodfile); %read the line
    end
    vals = cell2mat(textscan(prodfile, '%f %f %f %f %f\n', natoms));
    for j = 1:natoms
        rnew = [Lx_min + vals(j,3) * Lx, Ly_min + vals(j,4)*Ly, Lz_min + vals(j,5)*Lz];
        rorig = origlist(vals(j,1),:);
        deltax = abs(rnew(1) - rorig(1));
        deltay = abs(rnew(2) - rorig(2));
        deltaz = abs(rnew(3) - rorig(3));
        
        dev2 = min(deltax,Lx-deltax)^2 + min(deltay,Ly-deltay)^2 + min(deltaz,Lz-deltaz)^2;
        
        meanlist(vals(j,1),:) = meanlist(vals(j,1),:) + dev2;
    end
end

for j = 1:natoms
    meanlist(j) = meanlist(j) ./ nsteps;
end

meandev = mean(meanlist);
mdevs = [mdevs, meandev];
string(split(num2str(kvals)))';
end

scatter(kvals, cell2mat(mdevs),'filled');
hold on;
plot(kvals, cell2mat(mdevs));
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
labels = num2str(kvals');
text(kvals, cell2mat(mdevs), labels);
title('Deviation from Crystal Structure');
xlabel('Spring constant, k');
ylabel('\Sigma_i^{natoms}(r_i-r_{i,0})^2');
saveas(gcf,'deviation','epsc');
saveas(gcf,'deviation','png');
