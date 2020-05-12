clc
clear
close all
warning off

disl_num = [10,15,50,100,150,200];
res_val  = 80000;

leglag   = '%010d';

!rm master_0*
!rm ja_*

ja_sub_1to1 = [];
ja_sub_2to5 = [];

for i1 = 1:2
    
    for i0 = 1:length(disl_num)
        
        if i1 == 1
            script_val = [fileread('master_info_step1to1_marcc.txt'), ...
                newline, ['%%%%%%%%%%%%%%%%%%%%%%%%%'], ...
                newline, ['disl_num = ' num2str(disl_num(i0)) ';'], ...
                newline, ['res_val  = ' num2str(res_val) ';'], ...
                fileread('master_info_step1to1.txt')];
            
            Fun_Mcreate(string(script_val),['master_'  num2str(disl_num(i0),leglag) '_s1to1'])
            
            job_script = string('#!/bin/bash') + newline + ...
                string(['#SBATCH --job-name=s1to1_' num2str(disl_num(i0))]) + newline + ...
                string(fileread('job_info_step1to1_marcc.txt')) + ...
                newline + ["matlab -nodisplay -nosplash -nodesktop -r '" + ...
                ['master_'  num2str(disl_num(i0),leglag) '_s1to1'] + ";'"];
            
            Fun_Bcreate(job_script,['ja_'  num2str(disl_num(i0),leglag) '_s1to1'])
            ja_sub_1to1 = [ja_sub_1to1, newline, 'sbatch ', ['ja_'  num2str(disl_num(i0),leglag) '_s1to1.sh;']];
        else
            script_val = [fileread('master_info_step2to5_marcc.txt'), ...
                newline, ['%%%%%%%%%%%%%%%%%%%%%%%%%'], ...
                newline, ['disl_num = ' num2str(disl_num(i0)) ';'], ...
                newline, ['res_val  = ' num2str(res_val) ';'], ...
                fileread('master_info_step2to5.txt')];
            
            Fun_Mcreate(string(script_val),['master_'  num2str(disl_num(i0),leglag) '_s2to5'])
            
            job_script = string('#!/bin/bash') + newline + ...
                string(['#SBATCH --job-name=s2to5_' num2str(disl_num(i0))]) + newline + ...
                string(fileread('job_info_step2to5_marcc.txt')) + ...
                newline + ["matlab -nodisplay -nosplash -nodesktop -r '" + ...
                ['master_'  num2str(disl_num(i0),leglag) '_s2to5'] + ";'"];
            
            Fun_Bcreate(job_script,['ja_'  num2str(disl_num(i0),leglag) '_s2to5'])
            ja_sub_2to5 = [ja_sub_2to5, newline, 'sbatch ', ['ja_'  num2str(disl_num(i0),leglag) '_s2to5.sh;']];
        end
        
        
        
        
        
    end
end

Fun_Bcreate(ja_sub_1to1,'ja_sub_1to1')
Fun_Bcreate(ja_sub_2to5,'ja_sub_2to5')



% This function prints whatever is written in SCRIPT into  NAME.m and save
% it in the directory automaticaly

function Fun_Mcreate(SCRIPT,NAME)  %and use variable names that have meaning
NAME = sprintf('%s.m', NAME);
fid = fopen(NAME, 'wt'); %and use 't' with text files so eol are properly translated
%fprintf(fid,  SCRIPT);
fwrite(fid, SCRIPT);
fclose(fid);
% edit(NAME);

end

% SBATCH creator
function Fun_Bcreate(SCRIPT,NAME)
NAME = sprintf('%s.sh', NAME);
fid = fopen(NAME, 'wt'); %and use 't' with text files so eol are properly translated
%fprintf(fid,  SCRIPT);
fwrite(fid, SCRIPT);
fclose(fid);
% edit(NAME);
end
