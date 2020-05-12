function [A_final,R_final,A_PBC,R_PBC,B,B3] = fun_mesh_pbc_12142019(disl_num,d_star,num_pbc)
disl_num0    = disl_num - 1;
l0           = d_star*sin(pi/3);

% find i0_total
i0_p        = 1:disl_num;
i0_total    = find(disl_num0<=3.*(i0_p).*(i0_p+1) & disl_num0>=3.*(i0_p).*(i0_p-1)+1);

% add a bigger hexagon
ind_hexagon = randperm(5)-1;
% i0_total = i0_total+ind_hexagon(1).^2;
i0_total = ceil((disl_num-1)/6);
cond_val = true;
i0 = 0;
while cond_val
    i0 = i0 + 1;
    h_val(i0,1) = i0;
    
    if sum(h_val)>=i0_total
        i0_total = i0;
        cond_val = false;
    end
    
    
end

R_all = [];
A_all = [];

for i0  = 1:i0_total
    if rem(i0,2) == 0 
        lambda = [];
        l      = [];
        
        for j0 = 1:i0/2
            lambda(1,j0) = pi/2 - (pi/3/i0)*j0;
            l(1,j0)      = l0/sin(lambda(1,j0));
        end
        r = i0*[l(end:-1:1),l0,l(1:end-1)];
    else   
        lambda = [];
        l      = [];
        for j0 = 1:(i0+1)/2
            lambda(1,j0) = pi/2 - ( pi/3/(2*i0) + (j0-1)*(pi/3/i0)  );
            l(1,j0)      = l0/sin(lambda(1,j0));
        end
        r = i0*[l(end:-1:1),l(1:end-1)];   
    end
    R_i0    = repmat(r,1,6);
    del_val = pi/3/i0;
    A_i0    = 0:del_val:2*pi-del_val;
    
    R_all = [R_all,R_i0];
    A_all = [A_all,A_i0];   
end

R_all = [0,R_all];
A_all = [0,A_all];

% [aaa,rrr] = pol2cart(A_all',R_all');
% 
% ind_perm = nchoosek(1:length(A_all),2);
% 
% figure;polarplot(A_all,R_all,'.b','markersize',30)
% set(gca,'fontname','times','fontsize',30)
% hold on;
% for i0 = 1:size(ind_perm,1)
%     polarplot([A_all(ind_perm(i0,1)),A_all(ind_perm(i0,2))],[R_all(ind_perm(i0,1)),R_all(ind_perm(i0,2))],'--r','markersize',5)
% end



rng('shuffle')    ;
s = rng;

ind_final = randperm(length(A_all));

A_final = A_all(ind_final(1:disl_num));
R_final = R_all(ind_final(1:disl_num));

t = rand(1,length(A_final))*2*pi;
% fun_figrot(A_final,R_final,t)

% before centralizing (bc)
a_bc = A_final;
r_bc = R_final;
t_bc = t;


% % % save('matlab3.mat','A_final','R_final','t')

% fun_plotconfig(A_final,R_final)

% get the centroids of the main config
[a_lc,r_lc]      = fun_polar_centroid(A_final,R_final);

% transform the main config to the centroid
[r_lic0,a_lic0i] = fun_polar_glob2loc(r_lc,a_lc,R_final,A_final);

% Control
[a_cver,r_cver] = fun_polar_centroid(a_lic0i,r_lic0);
% wrap the transformed angles of the main config
a_lic0                 = wrapTo2Pi(a_lic0i);

% % compute the angle between principle axis and the local axis
B       = fun_productmoment(a_lic0',r_lic0');

% after centralizing (ac)
a_ac = a_lic0;
r_ac = r_lic0;
t_ac = t;

% fun_figrot(a_lic0,r_lic0,t)
% shift
A_final    = wrapTo2Pi(a_lic0    + repmat(B,1,disl_num));
R_final    = r_lic0;


% after shift (as)
a_as = A_final;
r_as = R_final;
t_as = wrapTo2Pi(t    + repmat(B,1,disl_num));

max_val = max(abs([r_bc,r_ac,r_as]));

% % % fun_figrot(a_bc,r_bc,t_bc,max_val)
% % % fun_figrot(a_ac,r_ac,t_ac,max_val)
% % % fun_figrot(a_as,r_as,t_as,max_val)

% fun_figrot(A_final,R_final,wrapTo2Pi(t    + repmat(B,1,disl_num)))

% fun_plotconfig(A_final,R_final)

%% PBC Section
% L_x = 2*max(abs(R_final.*cos(A_final)));
% L_y = 2*max(abs(R_final.*sin(A_final)));

L_x = 2.05*max(abs(R_final));
L_y = 2.05*max(abs(R_final));

%  figure;polarplot(A_final,R_final,'.','markersize',5)
% x=[];
% y=[];

% pbc_totnum = (2*num_pbc+1)^2*disl_num;
x = nan((2*num_pbc+1)^2,disl_num);
y = x;

for j0=1:disl_num
    count_val = 0;
    for i1=-num_pbc:num_pbc
        for i2=-num_pbc:num_pbc
            count_val = count_val + 1;
            x(count_val,j0) = R_final(j0)*cos(A_final(j0))+L_x*i1;
            y(count_val,j0) = R_final(j0)*sin(A_final(j0))+L_y*i2;
        end
    end
end
[A_PBC,R_PBC]=cart2pol(x,y);
A_PBC = wrapTo2Pi(A_PBC);

% % % figure;
% % % polarplot(A_PBC,R_PBC,'.r','markersize',5);
% % % hold on
% % % polarplot(A_final,R_final,'.b','markersize',5);
% % % set(gca,'fontname','times','fontsize',30);



%% caluclating Beta3
% centroid of all configurations (including the main one) in the global pcs
[a_lc_glb,r_lc_glb]      = fun_polar_centroid(A_PBC,R_PBC);

% local centroid of global pcs in all configurations
ind_main            = round(size(R_PBC,1)/2);
[r_lc_loc,a_lc_loc] = fun_polar_glob2loc(r_lc_glb,a_lc_glb,r_lc_glb(ind_main,:),a_lc_glb(ind_main,:));
% figure;polarplot(a_lc_loc,r_lc_loc,'*')

% calculating A_PBC_glb (not rotated angle of A_PBC in global PCS)
A_PBC_loc               = a_lic0;
R_PBC_loc               = r_lic0;
[R_PBC_glb,A_PBC_glb]   = fun_polar_loc2glob(r_lc_loc,a_lc_loc,R_PBC_loc,A_PBC_loc);
A_PBC_glb               = wrapTo2Pi(A_PBC_glb);
B3                      =  wrapTo2Pi(A_PBC - A_PBC_glb);

main_ind = round(size(A_PBC_glb,1)/2);
% % % figure;polarplot(A_PBC_glb,R_PBC_glb,'r.','markersize',5);
% % % hold on;polarplot(A_PBC_glb(main_ind,:),R_PBC_glb(main_ind,:),'b.','markersize',5);
% % % title(['PBC # = ' num2str(num_pbc) ' | ' num2str('Original')])
% % % set(gca,'fontname','times','fontsize',20)
% % % 
% % % figure;polarplot(A_PBC,R_PBC,'r.','markersize',5);
% % % hold on;polarplot(A_PBC(main_ind,:),R_PBC(main_ind,:),'b.','markersize',5);
% % % title(['PBC # = ' num2str(num_pbc) ' | ' num2str('Rotated Along Principle Axis')])
% % % set(gca,'fontname','times','fontsize',20)

end