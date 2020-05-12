function fit_val = fun_fitness_100719(x,r_li,a_li,r_lipbc,a_lipbc,disl_num,B,B3,j0,bin_size)

t_li = x;

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

force_val = force_val_main + force_val_pbc;

force_val(logical(eye(disl_num))) = NaN;
f_val = [-2.5e9;force_val(~isnan(force_val));+2.5e9;];

h_val = hist(f_val(:),bin_size);

fit_val = std(h_val);

end

