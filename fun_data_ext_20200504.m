function fun_data_ext_20200504(disl_num,res_val)
leglag     = '%010d';

%% folder info
% % % file_ext    = ['_' num2str(disl_num,'%05d') '_' num2str(res_val,'%09d')];
% % % folder_name = fullfile(pwd,['result' file_ext]);
folder_name = fullfile(pwd,[num2str(disl_num,leglag)],'results');

%% Read all configs 

dir_configs = dir(fullfile(folder_name,'config*'));
dir_configs(end) = [];
batch_num   = size(dir_configs,1);

% batch_num = 1000;
% load
tic
for i0 = 1:batch_num
    
    [out_val(i0).std_all,out_val(i0).std_len,out_val(i0).std_batch,out_val(i0).r,out_val(i0).a,out_val(i0).t,out_val(i0).f] = fun_ass_11172019_1(i0,folder_name,leglag);
    
    disp(['read stds | disl_num ' num2str(disl_num,leglag) ...
        ' res_val ' num2str(res_val,leglag) ...
        ' | batch # ' num2str(i0,leglag) ...
        ')' ' out of ' num2str(batch_num,leglag) ...
        ' | ' num2str(i0/batch_num*100,'%.f') '%'])
    
end
toc

res_val = 0;
for i0 = 1:batch_num
    res_val = res_val + length(out_val(i0).std_len);
    disp(['res_val value calculations '  num2str(i0,leglag)])
end

S = nan(res_val,1);
R = nan(res_val,disl_num);
A = nan(res_val,disl_num);    
T = nan(res_val,disl_num);
F = nan(res_val,disl_num);

% plug
c_val = 1;
for i0 = 1:batch_num
    range_val            = c_val : c_val + length(out_val(i0).std_len)-1;
    
    S(range_val,:) = out_val(i0).std_all;
    R(range_val,:) = out_val(i0).r;
    A(range_val,:) = out_val(i0).a;
    T(range_val,:) = out_val(i0).t;
    F(range_val,:) = out_val(i0).f;
    
    
    c_val = c_val + length(out_val(i0).std_len);
    
    disp(['plug stds | disl_num ' num2str(disl_num,leglag) ...
        ' res_val ' num2str(res_val,leglag) ...
        ' | batch # ' num2str(i0,leglag) ...
        ')' ' out of ' num2str(batch_num,leglag) ...
        ' | ' num2str(i0/batch_num*100,'%.f') '%'])
    
end

clear out_val

%% check uniqueness


[~,ind_unique] = unique([R,A,T],'rows');
len_old        = size(A,1);
len_new        = size(ind_unique,1);

S = S(ind_unique,:);
R = R(ind_unique,:);
A = A(ind_unique,:);
T = T(ind_unique,:);
F = F(ind_unique,:);

%% sort 
[S_sorted,ind_sorted]  = sort(S(:,1));
R_sorted               = R(ind_sorted,:);
A_sorted               = A(ind_sorted,:);
T_sorted               = T(ind_sorted,:);
F_sorted               = F(ind_sorted,:);

%% save
save(fullfile(folder_name,'data_raw.mat'),'S_sorted','R_sorted', ...
    'A_sorted','T_sorted','F_sorted','len_old','len_new','-v7.3')

%% save as dictionary in .txt for python
datain = [R_sorted, A_sorted, T_sorted];
dataou = F_sorted;
datain_list = fun_mat2list(datain);
dataou_list = fun_mat2list(dataou);

%dict_val = ["{'datain:' ",datain_list,', ',"'dataou:' ", dataou_list, '}'];
dict_val = append("{'datain': ",string(datain_list),", 'dataou': ", string(dataou_list), "}");
fid = fopen(fullfile(folder_name,'data_raw.txt'),'wt');
fprintf(fid, dict_val);
fclose(fid);

end

function [std_all,std_len,std_batch,r,a,t,f] = fun_ass_11172019_1(i0,folder_name,leglag)

load(fullfile(folder_name,['config_' num2str(i0,leglag)]),'std_all','r_li','a_li','t_li','frc_val_all');
std_len   = 1:length(std_all);
std_batch = i0;

std_all = (std_all);
r       = (r_li);
a       = (a_li);
t       = (t_li);

q           = (frc_val_all);
q(isnan(q)) = 0;
f           = squeeze(sum(q,1))';



end
