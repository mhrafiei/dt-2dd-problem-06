function fun_data_raw_20200504(disl_num,res_val)
leglag     = '%010d';
batch_num  = 10000;
batch_size = floor(res_val./batch_num); % save every batch_size configurations and remove them from the memory

c_val = 1;
for i0 = 1:batch_num
    
    if i0 ~=batch_num
        ind_batch(i0).batch = [c_val:c_val+batch_size-1]';
    else
        ind_batch(i0).batch = [c_val:res_val]';
    end
    ind_batch(i0).len   = length(ind_batch(i0).batch);
    c_val = c_val + batch_size;
    
end

t000 = datetime;
t000.Format='yyyy_MM_dd_HH_mm_ss_SSS';

file_ext = ['_' num2str(disl_num,'%05d') '_' num2str(res_val,'%09d') '_' char(t000)];
folder_name = fullfile(pwd,[num2str(disl_num,leglag)],['result' file_ext]);

if ~exist(folder_name,'dir')
    mkdir(folder_name)
end

%% initial variables
d_star     = 3*.249e-9;
leglag     = ['%0' num2str(ceil(log10(res_val))+1) 'd'];
num_pbc    = 3;

save(fullfile(folder_name,'configs_info.mat'),'ind_batch','batch_size','batch_num','leglag','d_star')

% predefined genetic optimization variables
lb      = zeros(1,disl_num);
ub      = ones(1,disl_num)*2*pi;
A       = [];
b       = [];
Aeq     = [];
beq     = [];

bin_size       = 10;
Markov_val     = 1000; % it must be at least as much as the PopulationSize
PopulationSize = 20;
MaxGenerations = 50;

%% get the stresses
master_sigma;

%% optimized configurations and their pbcs
% compute res_val configurations of the disl_num dislocaitons such the
% standard deviation of bins frequency of forces become minimimized for
% both main configurations and also the pbc configurations
num_pbc_tot = ((2*num_pbc+1)^2-1)*disl_num;

% restart the matlab random function
rng('shuffle')    ;
s = rng;

% create close (d_star) configurations of dislocation positions other than
% their rotations

ind_rem     = round(((2*num_pbc+1)^2)/2);
parfor i0 = 1:batch_num
    
    %     fun_ass_11152019_1(i0,ind_batch,disl_num,num_pbc_tot,num_pbc,d_star,leglag,folder_name,ind_rem,batch_num);
    fun_ass_11162019(i0,ind_batch,disl_num,num_pbc_tot,num_pbc,d_star, ...
        leglag,folder_name,ind_rem,batch_num,Markov_val,MaxGenerations,PopulationSize,...
        bin_size,A,b,Aeq,beq,lb,ub,res_val)
    
end

disp('                                        ')
disp('========================================')
disp(['Dislocation # ' num2str(disl_num) ' & res_val ' num2str(res_val) ' has been completed'])
disp('========================================')
disp('                                        ')

% fig = figure; hist(frc_val_all(~isnan(frc_val_all(:))),100);hold on;
% h_all = hist(frc_val_all(~isnan(frc_val_all(:))),100);
% title(['Main & PBC Standard Deviation: ' num2str(std(h_all),'%0.0006f') ])
% savefig(fig,fullfile(folder_name,['fig_all']))
% 
% fig = figure; hist(frc_val_mai(~isnan(frc_val_mai(:))),100);hold on;
% h_mai = hist(frc_val_mai(~isnan(frc_val_mai(:))),100);
% title(['Main Standard Deviation: ' num2str(std(h_mai),'%0.0006f') ])
% savefig(fig,fullfile(folder_name,['fig_main']))
% 
% fig = figure; hist(frc_val_pbc(~isnan(frc_val_pbc(:))),100);hold on;
% h_pbc = hist(frc_val_pbc(~isnan(frc_val_pbc(:))),100);
% title(['pbcn Standard Deviation: ' num2str(std(h_pbc),'%0.0006f') ])
% savefig(fig,fullfile(folder_name,['fig_pbc']))

% close all
% save(fullfile(folder_name,['data_raw.mat']),'std_all','r_li','a_li','t_li','frc_val_all','-v7.3')


end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% functions

function  force_val = fun_ass_090119(B,r_li,a_li,t_li,disl_num)

fun_fg_loc_glb = @(a,r,g,t) cos(g).*((sin(conj(g)).*(cos(3.*a - 2.*t) + ...
    cos(a - 2.*t)))./(2.*r) - (cos(conj(g)).*(sin(a - 2.*t) + ...
    sin(3.*a - 2.*t) + 2.*sin(a)))./(2.*r)) + ...
    sin(g).*((cos(conj(g)).*(cos(3.*a - 2.*t) + ...
    cos(a - 2.*t)))./(2.*r) + (sin(conj(g)).*(sin(a + 2.*t) + ...
    sin(3.*a + 2.*t) - 2.*sin(a)))./(2.*r));

i0 = 1;
force_val = nan(disl_num,disl_num);
for i1 = 1:disl_num % base
    [r_loc_trgd,a_loc_trgd]=fun_polar_glob2loc(r_li(i0,i1),a_li(i0,i1),r_li(i0,:),a_li(i0,:));
    a_loc_trgd = a_loc_trgd - wrapToPi(B);
    for i2 = 1: disl_num % target
        force_val(i1,i2) = fun_fg_loc_glb(a_loc_trgd(i2),r_loc_trgd(i2),t_li(i0,i2),t_li(i0,i1));
    end
end
force_val(logical(eye(disl_num))) = NaN;
end
function [force_val_main,force_val_pbc] = fun_ass_100319(r_li,a_li,r_lipbc,a_lipbc,t_li,disl_num,B,B3,j0)

fun_fg_loc_glb = @(a,r,g,t) cos(g).*((sin(conj(g)).*(cos(3.*a - 2.*t) + ...
    cos(a - 2.*t)))./(2.*r) - (cos(conj(g)).*(sin(a - 2.*t) + ...
    sin(3.*a - 2.*t) + 2.*sin(a)))./(2.*r)) + ...
    sin(g).*((cos(conj(g)).*(cos(3.*a - 2.*t) + ...
    cos(a - 2.*t)))./(2.*r) + (sin(conj(g)).*(sin(a + 2.*t) + ...
    sin(3.*a + 2.*t) - 2.*sin(a)))./(2.*r));

% force_val     = nan(disl_num,disl_num);
% numtot_base   = length(r_lipbc)+disl_num;
numpbc_config = length(r_lipbc)/disl_num;

r_lipbc       = reshape(r_lipbc,numpbc_config,disl_num);
a_lipbc       = reshape(a_lipbc,numpbc_config,disl_num);

force_val_pbc = zeros(disl_num,disl_num);
for i0 = 1:length(numpbc_config) % go over all available pbc configs
    for i1 = 1:disl_num
        
        %         [r_loc_trgd,a_loc_trgd] = fun_polar_glob2loc(r_lipbc(i0,i1),a_lipbc(i0,i1),r_li,a_li);
        % %         [r_loc_trgd,a_loc_trgd] = fun_polar_glob2loc(r_lipbc(i0,i1),a_lipbc(i0,i1),r_li,a_li);
        %         a_loc_trgd              = a_loc_trgd - wrapToPi(B3(i0,i1,j0));
        for i2 = 1: disl_num % target
            %             force_val_pbc(i1,i2)    = force_val_pbc(i1,i2) + fun_fg_loc_glb(a_loc_trgd(i2),r_loc_trgd(i2),t_li(i2),t_li(i1));
            r_glb_base = r_lipbc(i0,i1);
            a_glb_base = a_lipbc(i0,i1);
            t_base     = t_li(1,i1);
            r_glb_trgd = r_li(1,i2);
            a_glb_trgd = a_li(1,i2);
            t_trgd     = t_li(1,i2);
            beta3      = B3(i0,i1,j0);
            force_val_pbc(i1,i2) = fun_force(r_glb_base,a_glb_base,t_base,r_glb_trgd,a_glb_trgd,t_trgd,beta3);
        end
    end
    %     force_val_pbc(logical(eye(disl_num))) = NaN;
end
force_val_pbc(logical(eye(disl_num))) = NaN;

i0 = 1;
force_val_main = nan(disl_num,disl_num);
% force_val2 = nan(disl_num,disl_num);
for i1 = 1:disl_num % base
    %     [r_loc_trgd,a_loc_trgd]=fun_polar_glob2loc(r_li(i0,i1),a_li(i0,i1),r_li(i0,:),a_li(i0,:));
    %     a_loc_trgd = a_loc_trgd - wrapToPi(B(j0,1));
    %
    % %     [r_loc_trgd2,a_loc_trgd2]=fun_polar_glob2loc(r_li(i0,i1),a_li(i0,i1),r_li(i0,:),a_li(i0,:)-B(j0,1));
    %
    %     [r_loc_trgd,a_loc_trgd]=fun_polar_glob2loc(r_li(i0,i1),a_li(i0,i1),r_li(i0,:),a_li(i0,:));
    for i2 = 1: disl_num % target
        %         force_val_main(i1,i2) = fun_fg_loc_glb(a_loc_trgd(i2),r_loc_trgd(i2),t_li(i0,i2),t_li(i0,i1));
        % %         force_val2(i1,i2) = fun_fg_loc_glb(a_loc_trgd2(i2),r_loc_trgd2(i2),t_li(i0,i2),t_li(i0,i1));
        r_glb_base = r_li(1,i1);
        a_glb_base = a_li(1,i1);
        t_base     = t_li(1,i1);
        r_glb_trgd = r_li(1,i2);
        a_glb_trgd = a_li(1,i2);
        t_trgd     = t_li(1,i2);
        beta3      = B(j0,1);
        force_val_main(i1,i2) = fun_force(r_glb_base,a_glb_base,t_base,r_glb_trgd,a_glb_trgd,t_trgd,beta3);
    end
end
force_val_main(logical(eye(disl_num))) = NaN;
% force_val = force_val_main + force_val_pbc;
end
function fun_ass_11162019(i0,ind_batch,disl_num,num_pbc_tot,num_pbc,d_star, ...
    leglag,folder_name,ind_rem,batch_num,Markov_val,MaxGenerations,PopulationSize,...
    bin_size,A,b,Aeq,beq,lb,ub,res_val)

% predefined data
num_val = ind_batch(i0).len ;

a_li = (nan(num_val ,disl_num));
r_li = (nan(num_val ,disl_num));

a_lipbc     = (nan(num_val,num_pbc_tot));
r_lipbc     = (nan(num_val,num_pbc_tot));
B           = (nan(num_val,1));
B3          = (nan((2*num_pbc+1)^2,disl_num, num_val));

frc_val_mai = nan(disl_num,disl_num,num_val);
frc_val_pbc = nan(disl_num,disl_num,num_val);
frc_val_all = nan(disl_num,disl_num,num_val);

std_all = nan(num_val,1);

t_li = (nan(num_val ,disl_num));

rng('shuffle')    ;
s = rng;

for i1 = 1:num_val
    % this mesh function not only figure out the dislocation locations but
    % rotate them along the principle axis and also transform them over
    % their centroid for both main and pbc configurations.
    
    [A_final,R_final,A_PBC,R_PBC,B(i1,1),B3(:,:,i1)] = fun_mesh_pbc_12142019(disl_num,d_star,num_pbc);
    
    A_PBC(ind_rem,:) = [];
    R_PBC(ind_rem,:) = [];
    
    a_li(i1,:)    = A_final;
    r_li(i1,:)    = R_final;
    a_lipbc(i1,:) = A_PBC(:);
    r_lipbc(i1,:) = R_PBC(:);
    
    %     figure;polarplot(A_PBC,R_PBC,'r*'); hold on; polarplot(A_final,R_final,'b*')
    
    sigma_val = nan(Markov_val,1);
    
    t_li_case = nan(Markov_val,disl_num);
    
    for i2 = 1:Markov_val
        
        t_li_case(i2,:) = rand(1,disl_num)*2*pi;
        force_val_p     = fun_ass_090119(0,r_li(i1,:),a_li(i1,:),t_li_case(i2,:),disl_num);
        f_val           = [-2.5e9;force_val_p(~isnan(force_val_p));+2.5e9];
        h_val           = hist(f_val(:),bin_size);
        sigma_val(i2,1) = std(h_val);
    end
    
    [~,ind_sorted] = sort(sigma_val);
    x0             = t_li_case(ind_sorted(1:PopulationSize),:);
    
    options = optimoptions('ga','ConstraintTolerance',1e-128, ...
        'Display','off', ...
        'FunctionTolerance',1e-128,'PopulationSize',PopulationSize, ...
        'InitialPopulationMatrix',x0, ...
        'UseParallel',false, ...
        'MaxGenerations',MaxGenerations);
    %'PlotFcns', @gaplotbestf, ...
    %'Display','iter', ...
    
    fun            = @(x)fun_fitness_100719(x,r_li(i1,:),a_li(i1,:),r_lipbc(i1,:),a_lipbc(i1,:),disl_num,B,B3,i1,bin_size);
% % %     for i2 = 1:100
% % %         tic
% % %          options = optimoptions('ga','ConstraintTolerance',1e-128, ...
% % %         'Display','none', ...
% % %         'FunctionTolerance',1e-128,'PopulationSize',PopulationSize, ...
% % %         'InitialPopulationMatrix',x0, ...
% % %         'UseParallel',false, ...
% % %         'MaxGenerations',1);
% % %         [t_li(i1,:),fit_val(i2,1),fit_exit,fit_output,fit_pop,fit_scores] = ga(fun,disl_num,A,b,Aeq,beq,lb,ub,[],options);
% % %         x0 = fit_pop;
% % %         disp(['i2: ' num2str(i2,'%03d') ' | ' num2str(fit_val(i2,1),'%0.3d') ' | time: ' num2str(round(toc,2))])
% % %     end
    
    [t_li(i1,:),fit_val,fit_exit,fit_output,fit_pop,fit_scores] = ga(fun,disl_num,A,b,Aeq,beq,lb,ub,[],options);
    
    %% 6 slip + 6 partials
    round_vals = [[0,1,2,3,4,5,6]  * pi/3,[1,3,5,7,9,11] * pi/6];
    round_asgn = [[0,1,2,3,4,5,0]  * pi/3,[1,3,5,7,9,11] * pi/6];
    t_temp = nan(13,length(t_li(i1,:)));
    for i2 = 1:13
        t_temp(i2,:) = abs(t_li(i1,:) - round_vals(i2));
    end
    
    [min_temp, ind_temp] = min(t_temp);
    t_li(i1,:) = round_asgn(ind_temp);
    
    
% % %     fig = figure('Name','progress','units','normalized','outerposition',[0 0 1 1]);
% % %     plot(fit_val([1,10:100]),'-r','linewidth',10);
% % %     xticks([1,10:100])
% % %     xticklabel({['1','','','','','50','','','','','100']})
% % %     xlabel('Generation Number');
% % %     ylabel('Fitness Value');
% % %     set(gca,'fontname','times','fontsize',70)
% % %     pbaspect([1 1 1])
% % %     set(gcf,'color','white')
% % %     frame_h = get(handle(gcf),'JavaFrame');
% % %     set(frame_h,'Maximized',1);
% % %     grid on
% % %     print(gcf, '-dtiff', '-r600', ['fig3d_' num2str(length(A_final)) '.tif']);
    
% % %     load matlab3
% % %     fun_figrot(A_final,R_final,t_li(i1,:))
    
    [frc_val_mai(:,:,i1),frc_val_pbc(:,:,i1)]  = fun_ass_100319(r_li(i1,:),a_li(i1,:),r_lipbc(i1,:),a_lipbc(i1,:),t_li(i1,:),disl_num,B,B3,i1);
    frc_val_all(:,:,i1)                        = frc_val_mai(:,:,i1) + frc_val_pbc(:,:,i1);
    
    f_all               = frc_val_all(:,:,i1);
    f_all               = [-2.5e9;f_all(~isnan(f_all));+2.5e9];
    h_all               = hist(f_all,100);
    std_all(i1,1)       = std(h_all);
    
    disp(['disl_num ' num2str(disl_num,leglag) ' res_val ' num2str(res_val,leglag) ' | batch# ' num2str(i0,leglag) ' out of ' num2str(batch_num,leglag) ' | case# ' num2str(i1,'%09d') ' out of ' num2str(num_val,leglag) ' | ' num2str(i0/batch_num*100,'%.f') '%'])

end

batch_id = i0;

t = datetime;
t.Format='yyyy_MM_dd_HH_mm_ss_SSS';


save(fullfile(folder_name,['configs_' num2str(i0,leglag) '_' char(t) '.mat']),'a_li','r_li','a_lipbc','r_lipbc','B','B3','t_li','frc_val_mai','frc_val_pbc','frc_val_all','std_all','batch_id','-v7.3')

clear

end
