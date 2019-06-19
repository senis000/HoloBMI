%add_paths.m
[~, user_str]   = system('whoami');
user_str        = regexprep(user_str,'\r\n|\n|\r','');
%MAC MACHINE:
base_dir        =   fullfile('/', 'Users', user_str, 'Dropbox');
    
export_fig_path     = ...
    fullfile(base_dir, 'Code', 'export_fig', 'export_fig_11.19.15')
exist(export_fig_path)
addpath(genpath(export_fig_path)); 

%Add matlabexchange directories: 
matlabexchange_dir = fullfile(base_dir, 'Code', 'MatlabExchange'); 
addpath(genpath(matlabexchange_dir))