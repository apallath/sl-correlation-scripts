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

%% Open the files containing the dumped data
zfile = fopen('Z_top.txt');
zmfile = fopen('Z_top.txt');
inittopposfile = fopen(input_params('solid_init_top_positions_file'));
meantopposfile = fopen(input_params('solid_top_positions_file'));
curtopposfile = fopen(input_params('solid_top_positions_file'));

% positions
xylist = zeros(natoms,2);  %%initlist
surfzlist = zeros(natoms,1); %%surface z list
curzlist = zeros(natoms,1); %%current surface z list
zlist = zeros(natoms,1);   %%curlist
zmlist = zeros(natoms,1);  %%meanlist

%% Variables
ival = 0.005 * max([Lx, Ly]);
rmax = 2 * max([Lx, Ly]);
rvals = 0:ival:rmax;
g = zeros(size(rvals));
gcount = zeros(size(rvals));

%% Calculate init positions
for j = 1:9
    fgetl(inittopposfile); %read the line
end
vals = cell2mat(textscan(inittopposfile, '%f %f %f %f %f\n', natoms));
for j = 1:natoms
    xylist(j, :) = xylist(j, :) + [Lx_min + vals(j,3) * Lx, Ly_min + vals(j,4)*Ly];
end

%% calculate mean z values

for i = 1:nsteps 
    % get top surface positions
    for j = 1:9
        fgetl(meantopposfile); %read the line
    end
    vals = cell2mat(textscan(meantopposfile, '%f %f %f %f %f\n', natoms));
    for j = 1:natoms
        surfzlist(j,1) = Lz_min + vals(j,5) * Lz;
    end
    
    % get zmin values 
    fgetl(zmfile);
    vals = cell2mat(textscan(zmfile, '%f %f\n', natoms));
    for j = 1:natoms
        zmlist(j,1) = zmlist(j,1) + vals(j,2) - surfzlist(j,1); %distance from surface
    end
end  

zmlist(:,1) = zmlist(:,1)./ nsteps;

% SAVE
save('zz_corr_zmean.mat','zmlist');

%% Calculate correlations
for i = 1:nsteps
    disp(sprintf('Working on timestep %d/%d',i,nsteps));
    
    % get top surface positions
    for j = 1:9
        fgetl(curtopposfile); %read the line
    end
    vals = cell2mat(textscan(curtopposfile, '%f %f %f %f %f\n', natoms));
    for j = 1:natoms
        curzlist(j,1) = Lz_min + vals(j,5) * Lz;
    end
    
    %
    
    fgetl(zfile); %read the line

    vals = cell2mat(textscan(zfile, '%f %f %f\n', natoms));
    
    for j = 1:natoms
        zlist(j,:) = vals(j,2) - curzlist(j,1);
    end
    
    for m = 1:natoms-1
        for n = m+1:natoms
           if m ~= n
                rm = xylist(m,1:2);
                rn = xylist(n,1:2);
                rm_m = zlist(m,1);
                rn_m = zlist(n,1);
                um = rm_m - zmlist(m,1);
                un = rn_m - zmlist(n,1);
                
                deltax = abs(rm(1) - rn(1));
                deltay = abs(rm(2) - rn(2));             
                
                r = sqrt(min(deltax,Lx-deltax)^2 + min(deltay,Ly-deltay)^2);
               
                bin = floor(r/ival) + 1;
                g(bin) = g(bin) + um*un; % norm(um)*norm(un);
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
    save(sprintf('zz_corr_rvals.mat'),'rvals');    
    save(sprintf('zz_corr_mean.mat'),'curg');
end