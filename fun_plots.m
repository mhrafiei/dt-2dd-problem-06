function fun_plots(disl_num,res_val)

name_in = 'datain_scl';% 'datain_scl'
name_ou = 'dataoug_scl';% 'dataoul_scl'

file_ext    = ['_' num2str(disl_num,'%05d') '_' num2str(res_val,'%09d')];
folder_name = fullfile(pwd,['result' file_ext]);

dir_data = dir(fullfile(folder_name,'proc*'));


parfor i0 = 1:length(dir_data)
    [data(i0).in,data(i0).ou,data(i0).len] = fun_ass_nov232019(i0,dir_data,folder_name,name_in,name_ou);
    disp(['batch #: ' num2str(i0,'%07d')])
end


c = 0;
for i0 = 1:length(dir_data)
    c = c  + data(i0).len;
end

datain = nan(c,disl_num*3);
dataou = nan(c,disl_num);
c = 1;
for i0 = 1:length(dir_data)
    datain(c:c+data(i0).len-1,:) = data(i0).in;
    dataou(c:c+data(i0).len-1,:) = data(i0).ou;
    c = c + data(i0).len;
end

end

function [in_val,ou_val,len_val] = fun_ass_nov232019(i0,dir_data,folder_name,name_in,name_ou)

load(fullfile(folder_name,dir_data(i0).name),name_in,name_ou);

in_val = eval(name_in);
ou_val = eval(name_ou);
len_val = size(in_val,1);
end
