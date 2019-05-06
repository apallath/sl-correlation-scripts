% Akash Pallath <akash.pallath@btech15.iitgn.ac.in>
% Created: April 2019
% Last modified: May 2019
% 
% Get z-values above each atom

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
natoms = str2num(input_params('natoms'));

% number of top solid atoms in the region
ntopatoms = str2num(input_params('nsolidtopatoms'));

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
posfile = fopen(input_params('positions_file'));
topposfile = fopen(input_params('solid_top_positions_file'));
inittopposfile = fopen(input_params('solid_init_top_positions_file'));

zvfile = fopen('Z_top.txt','w');

initlist = zeros(ntopatoms,3);
curlist = zeros(ntopatoms,4);

zlist = inf * ones(nsteps,ntopatoms);
%% Calculate min z positions
for i = 1:nsteps
    fprintf('Working on timestep %d/%d\n',i,nsteps);
    
    %%
    for j = 1:9
        fgetl(topposfile); %read the line
    end
    vals = cell2mat(textscan(topposfile, '%f %f %f %f %f\n', ntopatoms));
    for j = 1:ntopatoms
        curlist(j, :) = [Lx_min + vals(j,3) * Lx, Ly_min + vals(j,4)*Ly, Lz_min + vals(j,5)*Lz, vals(j,1)];
    end    
    
    %%    
    for j = 1:9
        fgetl(posfile); %read the line
    end
    
    vals = cell2mat(textscan(posfile, '%f %f %f %f %f\n', natoms));
    
    for j = 1:natoms        
        if vals(j,2) == 1
            xpos = Lx_min + vals(j,3) * Lx;
            ypos = Ly_min + vals(j,4) * Lx;
            zpos = Lz_min + vals(j,5) * Lz;
        
            for k = 1:ntopatoms
                rm = [xpos, ypos];
                rn = curlist(k,1:2);
                
                deltax = abs(rm(1) - rn(1));
                deltay = abs(rm(2) - rn(2));
                
                r = sqrt(min(deltax,Lx-deltax)^2 + min(deltay,Ly-deltay)^2);
                
                if r <= 1
                    zlist(i,k) = min(zlist(i,k), zpos);
                end
            end        
        end
    end
    
    fprintf(zvfile,"%d\n",i);
    
    for j=1:ntopatoms
        fprintf(zvfile,"%d %f\n",[curlist(j,4), zlist(i,j)]);
    end
end
