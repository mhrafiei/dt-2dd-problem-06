function fun_data_prp_20200504(disl_num,res_val)
leglag     = '%010d';
%% folder name
% % % file_ext    = ['_' num2str(disl_num,'%05d') '_' num2str(res_val,'%09d')];
% % % folder_name = fullfile(pwd,['result' file_ext]);
folder_name = fullfile(pwd,[num2str(disl_num,leglag)],'results');


%% delete 
% delete(fullfile(folder_name,'fig*'));
% delete(fullfile(folder_name,'data_ind*'));
% delete(fullfile(folder_name,'datafinali*'));
% delete(fullfile(folder_name,'mdl*'));

%% predefined parameters
d_star     = 3*.249e-9;
leglag     = '%010d';%['%0' num2str(ceil(log10(res_val))+1) 'd'];
rtt       =  [0.001:0.001:0.005];%[0.01:0.01:0.05];1-[0.125:-0.025:0.025];
rrs       =  1;
ver_rate  =  0.001;
data_rate =  0.95;
lb_val    = +0.05;
ub_val    = +0.95;
flb       = -4.652361404625608e+10;%-2.5e9;
fub       = +4.652361404625608e+10;%+2.5e9;

%% predefined functions
fun_scale       = @(x,lb,ub)   ((x-(ones(size(x,1),1)*min(x)))./((ones(size(x,1),1)*max(x))-(ones(size(x,1),1)*min(x))))*(ub-lb)+lb;

%% loading 
load(fullfile(folder_name,'data_raw.mat'),'R_sorted','A_sorted','T_sorted','F_sorted')

%% check repetition 
% already checked in fun_data_ext_*

%% histogram of all
close all

figure;hist(A_sorted(:),100);
xlabel('Range')
ylabel('Frequency')
set(gca,'fontname','times','fontsize',30)
pbaspect([1 1 1])

figure;hist(R_sorted(:),100);
xlabel('Range')
ylabel('Frequency')
set(gca,'fontname','times','fontsize',30)
pbaspect([1 1 1])

figure;hist(T_sorted(:),100);
xlabel('Range')
ylabel('Frequency')
set(gca,'fontname','times','fontsize',30)
pbaspect([1 1 1])

figure;hist(F_sorted(:),100);
xlabel('Range')
ylabel('Frequency')
set(gca,'fontname','times','fontsize',30)
pbaspect([1 1 1])

%% histogram of top data_rate %

top_val = floor(data_rate*size(A_sorted,1));

a_case = A_sorted(1:top_val,:);
figure;hist(a_case(:),100);
xlabel('Range')
ylabel('Frequency')
set(gca,'fontname','times','fontsize',30)
pbaspect([1 1 1])

r_case = R_sorted(1:top_val,:);
figure;hist(r_case(:),100);
xlabel('Range')
ylabel('Frequency')
set(gca,'fontname','times','fontsize',30)
pbaspect([1 1 1])

t_case = T_sorted(1:top_val,:);
figure;hist(t_case(:),100);
xlabel('Range')
ylabel('Frequency')
set(gca,'fontname','times','fontsize',30)
pbaspect([1 1 1])

f_case = F_sorted(1:top_val,:);
figure;hist(f_case(:),100);
xlabel('Range')
ylabel('Frequency')
set(gca,'fontname','times','fontsize',30)
pbaspect([1 1 1])


len_valold = size(A_sorted,1);
len_valnew = floor(len_valold * data_rate);

R = R_sorted(1:len_valnew,:); clear R_sorted;
A = A_sorted(1:len_valnew,:); clear A_sorted;
T = T_sorted(1:len_valnew,:); clear T_sorted;
F = F_sorted(1:len_valnew,:); clear F_sorted;

% % % %% plot histograms 
% % % for i0 = 1:disl_num
% % %     
% % %     fig = figure;
% % %     hist(R(:,i0),100);
% % %     savefig(fig,fullfile(folder_name,['fig_R_' num2str(i0,leglag)]))
% % %     
% % %     fig = figure;
% % %     hist(A(:,i0),100);
% % %     savefig(fig,fullfile(folder_name,['fig_A_' num2str(i0,leglag)]))
% % %     
% % %     fig = figure;
% % %     hist(T(:,i0),100);
% % %     savefig(fig,fullfile(folder_name,['fig_T_' num2str(i0,leglag)]))
% % %     
% % %     fig = figure;
% % %     hist(F(:,i0),100);
% % %     savefig(fig,fullfile(folder_name,['fig_F_' num2str(i0,leglag)]))
% % %     
% % %     close all
% % % end

datafinali_name = fullfile(folder_name,'data_inou.mat');


    
datain  = fun_scale([[zeros(1,disl_num),zeros(1,disl_num),zeros(1,disl_num)];[ones(1,disl_num),2*pi*ones(1,disl_num),2*pi*ones(1,disl_num)];[R,A,T]],lb_val,ub_val);
%dataou  = fun_scale([disl_num/3*flb*ones(1,disl_num);disl_num/3*fub*ones(1,disl_num);F],lb_val,ub_val);
dataou  = fun_scale([flb*ones(1,disl_num);fub*ones(1,disl_num);F],lb_val,ub_val);

datain(1:2,:) = [];
dataou(1:2,:) = [];

save(datafinali_name,'datain','dataou','-v7.3')

datain_list = fun_mat2list(datain);
dataou_list = fun_mat2list(dataou);

%dict_val = ["{'datain:' ",datain_list,', ',"'dataou:' ", dataou_list, '}'];
dict_val = append("{'datain': ",string(datain_list),", 'dataou': ", string(dataou_list), "}");
fid = fopen(fullfile(folder_name,'data_inou.txt'),'wt');
fprintf(fid, dict_val);
fclose(fid);

end

