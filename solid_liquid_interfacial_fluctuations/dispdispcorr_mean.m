% Akash Pallath <akash.pallath@btech15.iitgn.ac.in>
% Created: January 2019
% Last modified: May 2019
% 
% Calculate disp-disp correlations for topmost solid atom layer
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

for i = 1:sz[1]
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

%% Open the files containing the dumped data
posfile = fopen(input_params('solid_top_positions_file'));
meanposfile = fopen(input_params('solid_top_positions_file'));
initposfile = fopen(input_params('solid_init_top_positions_file'));

% positions
meanlist = zeros(natoms,3);
curlist = zeros(natoms,3);
initlist = zeros(natoms,3);

%% Variables
ival = 0.005 * max([Lx, Ly, Lz]);
rmax = 1.4 * max([Lx, Ly, Lz]);
rvals = 0:ival:rmax;
g = zeros(size(rvals));
gcount = zeros(size(rvals));

%% Calculate init positions
for j = 1:9
    fgetl(initposfile); %read the line
end
vals = cell2mat(textscan(initposfile, '%f %f %f %f %f\n', natoms));
for j = 1:natoms
    initlist(j, :) = initlist(j, :) + [Lx_min + vals(j,3) * Lx, Ly_min + vals(j,4)*Ly, Lz_min + vals(j,5)*Lz];
end

% SAVE
save('solid_top_initpos.mat','initlist');

%% Calculate mean positions
for i = 1:nsteps    
    for j = 1:9
        fgetl(meanposfile); %read the line
    end
    vals = cell2mat(textscan(meanposfile, '%f %f %f %f %f\n', natoms));
    for j = 1:natoms
        meanlist(j, :) = meanlist(j, :) + [Lx_min + vals(j,3) * Lx, Ly_min + vals(j,4)*Ly, Lz_min + vals(j,5)*Lz];
    end
end

for j = 1:natoms
    meanlist(j,:) = meanlist(j,:) ./ nsteps;
end

% SAVE
save('solid_top_meanpos.mat','meanlist');

%% Calculate correlations
for i = 1:nsteps
    disp(sprintf('Working on timestep %d/%d',i,nsteps));
    
    for j = 1:9
        fgetl(posfile); %read the line
    end
    vals = cell2mat(textscan(posfile, '%f %f %f %f %f\n', natoms));
    for j = 1:natoms
        curlist(j, :) = [Lx_min + vals(j,3) * Lx, Ly_min + vals(j,4)*Ly, Lz_min + vals(j,5)*Lz];
    end
    
    for m = 1:natoms-1
        for n = m+1:natoms
            if m ~= n
                rm = initlist(m, :);
                rn = initlist(n, :);
                rm_m = curlist(m, :);
                rn_m = curlist(n, :);
                um = rm_m - meanlist(m, :);
                un = rn_m - meanlist(n, :);
                
                deltax = abs(rm(1) - rn(1));
                deltay = abs(rm(2) - rn(2));
                deltaz = abs(rm(3) - rn(3));
                
                r = sqrt(min(deltax,Lx-deltax)^2 + min(deltay,Ly-deltay)^2 + deltaz^2);
                
                bin = floor(r/ival) + 1;
                
                g(bin) = g(bin) + norm(um)*norm(un);
                gcount(bin) = gcount(bin) + 1;
            end
        end
    end   

    %store correlations (till now) in a mat file
    curgcount = ones(size(rvals));
    for gci = 1:length(gcount)
        if gcount(gci) > 1
            curgcount(gci) = gcount(gci);
        end
    end
    curg = g./curgcount;    
    
    % SAVE
    save('solid_top_rvals.mat','rvals');    
    save('solid_top_disp_disp_corr_mean.mat','curg');
end
