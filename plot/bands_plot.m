clear variables
clc

filename = 'filename.dat';

%% ---------------------------------------------------------
%% Read file

file = fopen(filename);
if file == -1
    error('File cannot be opened')
end

header = fgets(file, 500);
fprintf(1, '%s\n\n', strtrim(header));

file_format_version = fread(file, 1, 'uint32');
if file_format_version ~= 103
    error(['Bad file format version ' num2str(file_format_version)]);
end

a = fread(file, [3 3], 'double')
b = fread(file, [3 3], 'double')

n_spins   = fread(file, 1, 'uint32')
n_kpoints = fread(file, 1, 'uint32')
n_bands   = fread(file, 1, 'uint32')
n_layers  = fread(file, 1, 'uint32')

supercell_height = fread(file, 1, 'double');
fermi_energy     = fread(file, 1, 'double');

energy_min = fread(file, 1, 'double');
energy_max = fread(file, 1, 'double');
cs_sq_max  = fread(file, 1, 'float');

ks = zeros(3, n_kpoints);
energies    = zeros(n_bands, n_kpoints);
occupations = zeros(n_bands, n_kpoints);
cs          = zeros(n_layers, n_bands, n_kpoints);

for ik = 1 : n_kpoints
    ks(:, ik)          = fread(file, [1 3], 'double');
    energies(:, ik)    = fread(file, [1 n_bands], 'double');
    occupations(:, ik) = fread(file, [1 n_bands], 'double');
    cs(:, :, ik)       = fread(file, [n_layers n_bands], 'float');
end

fclose(file);

%% ---------------------------------------------------------
%% Plot bands

figure
plot(energies' - fermi_energy, 'k');
