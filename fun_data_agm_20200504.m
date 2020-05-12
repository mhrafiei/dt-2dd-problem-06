function fun_data_agm_20200504(disl_num,res_val)
leglag = '%010d';

result_folders = dir(fullfile(pwd,[num2str(disl_num,leglag)],'result_*'));
result_masterfolder = fullfile(pwd,[num2str(disl_num,leglag)],'results');

if exist(result_masterfolder,'dir')==0
    mkdir(result_masterfolder);
end


c_val = 0;
for i0 = 1:length(result_folders)
   
    result_folder = fullfile(result_folders(i0).folder,result_folders(i0).name);
    
    dir_configs = dir(fullfile(result_folder,'configs_*'));
    
    dir_configs(end)=[];
    
    % determine mode of the size of data
    for i1 = 1:length(dir_configs)
        sourec_val      = fullfile(dir_configs(i1).folder,dir_configs(i1).name);
        s               = dir(sourec_val);
        b(i1,1)         = s.bytes/1000000;
    end
    min_val = mean(b)/2;
    
    % copy process
    for i1 = 1:length(dir_configs)
        c_val = c_val + 1;
        sourec_val      = fullfile(dir_configs(i1).folder,dir_configs(i1).name);
        destination_val = fullfile(result_masterfolder,['config_' num2str(c_val,leglag) '.mat']);
        
        s = dir(sourec_val);
        
        if s.bytes/1000000<min_val
            delete(sourec_val)
        else
            copyfile(sourec_val, destination_val)
        end
        
        disp(['folder# ' num2str(i0,leglag) ' out of ' num2str(length(result_folders),leglag) ' | file# ' num2str(i1,leglag) ' out of ' num2str(length(dir_configs),leglag) ])
    end
    
end


end
