function fun_data_cpt_20200504(disl_num,res_val)
leglag = '%010d';
range_val = 1:10; % plot the first 10 configs 

folder_name = fullfile(pwd,[num2str(disl_num,leglag)],'results');
load(fullfile(folder_name,'data_raw.mat'),'R_sorted','A_sorted','T_sorted','F_sorted')

for i0 = range_val
   
    max_val = max(abs(R_sorted(i0,:)));
    fun_figrot(A_sorted(i0,:),R_sorted(i0,:),T_sorted(i0,:),max_val)
end