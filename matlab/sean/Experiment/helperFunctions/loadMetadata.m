metadata_file_path = '/Users/sean/Desktop/2018_06_14_SAY98_HD233/2018_06_14_SAY98_HD233.dat';

metadata = readtable(metadata_file_path);

if sum(isnan(metadata.Strain))
    disp('There is no Strain information in the metadata file provided. To continue do either:');
    disp('    (a) Enter strain information manually into JMP, resave as .dat and reload experiment in MATLAB');
    disp('    (b) If you have a metamorph log file containing the x/y coordinates of the frames, you can use that now');
    
    log_file_path = rename_logfile_if_necessary(input('Metamorph log file path: ', 's'));
    strains_str = input('Enter Column Strains (separated by comma, e.g. HD233, SAY98, HD233, SAY98): ', 's');
    nChannels = input('Enter # of frames / worm: ');
    strains = strsplit(strains_str, ',');
    disp('Clustering Columns. Ensure that in the plot, there is *1 COLOR* per column. If there are more than one colors, you will need to open JMP and assign strains manually');
    metadata.Strain = clusterColumns(log_file_path, strains, nChannels, 1);
end

function new_name = rename_logfile_if_necessary(old_name)
    if endsWith(old_name, '.LOG')
        cam_pos_fn_base = strsplit(old_name, '.');
        cam_pos_fn_base = cam_pos_fn_base(1);
        
        new_name = strcat(cam_pos_fn_base, '.txt');
        new_name = new_name{1};
        
        sprintf('Renaming Metamorph Log File from %s to %s.', old_name, new_name);
        movefile(old_name, new_name);
    else
        new_name = old_name;
    end
end