% compute the glide component of forces in a target dislocation comming
% from a base dislocation; beta 3 is if we have any ptoduct moment rotation

function f = fun_force(r_glb_base,a_glb_base,t_base,r_glb_trgd,a_glb_trgd,t_trgd,beta3)
% assumption: we are in 2st local coordinate system (figure 11) and a and r are the
% coordinate of the second dislocation and g is the rotation of the second
% dislocation with respect to the horizontal axis of the 2st local
% coordinate system and t is the rotation of the base dislocation with
% respect to the 2nd local coordinate system

fun_fg_loc_glb  = @(a,r,g,t) cos(g).*((sin(conj(g)).*(cos(3.*a - 2.*t) + ...
                             cos(a - 2.*t)))./(2.*r) - (cos(conj(g)).*(sin(a - 2.*t) + ...
                             sin(3.*a - 2.*t) + 2.*sin(a)))./(2.*r)) + ...
                             sin(g).*((cos(conj(g)).*(cos(3.*a - 2.*t) + ...
                             cos(a - 2.*t)))./(2.*r) + (sin(conj(g)).*(sin(a + 2.*t) + ...
                             sin(3.*a + 2.*t) - 2.*sin(a)))./(2.*r));
                         
a_glb_base              = wrapTo2Pi(a_glb_base - beta3);
a_glb_trgd              = wrapTo2Pi(a_glb_trgd - beta3);

[r_loc_trgd,a_loc_trgd] = fun_polar_glob2loc(r_glb_base,a_glb_base,r_glb_trgd,a_glb_trgd);
a_loc_trgd              = wrapTo2Pi(a_loc_trgd);

f                       = fun_fg_loc_glb(a_loc_trgd,r_loc_trgd,t_trgd,t_base);

end

% updated 06/26/2019
%% Inputs
% r_glb_base : coordinate of the base in the global PCS
% a_glb_base : coordinate of the base in the global PCS
% r_glb_trgd : coordinate of the trgd in the global PCS
% a_glb_trgd : coordinate of the trgd in the global PCS

%% Outputs
% r_loc_trgd : coordinate of the trgd in the local  PCS
% a_loc_trgd : coordinate of the trgd in the local  PCS

%% Example
% r_glb_base = 2.8;
% a_glb_base = pi-pi/6;
% r_glb_trgd = 4.4786;
% a_glb_trgd = 3.5348;
% r_loc_trgd = 3.555277766926234;
% a_loc_trgd = 4.21;

function [r_loc_trgd,a_loc_trgd]=fun_polar_glob2loc(r_glb_base,a_glb_base,r_glb_trgd,a_glb_trgd)

[x_glb_trgd,y_glb_trgd] = pol2cart(a_glb_trgd,r_glb_trgd);
[x_glb_base,y_glb_base] = pol2cart(a_glb_base,r_glb_base);

x_loc_trgd = x_glb_trgd - x_glb_base;
y_loc_trgd = y_glb_trgd - y_glb_base;

[a_loc_trgd,r_loc_trgd] = cart2pol(x_loc_trgd,y_loc_trgd);



% a_loc_trgd      = wrapTo2Pi(a_loc_trgd);
ind             = r_loc_trgd<0;
a_loc_trgd(ind) = wrapTo2Pi(a_loc_trgd(ind) + pi);
r_loc_trgd      = abs(r_loc_trgd);

end

% RLI(:,:,1) = [-(cos(AGI)./2 + 1./2).*(cos(AGL)./2 + 1./2).*((RGI.^2 - 2.*cos(AGI - AGL).*RGI.*RGL + RGL.^2)./(cos(AGI./2).^4.*cos(AGL./2).^4)).^(1./2)]; 
% RLI(:,:,2) = [+(cos(AGI)./2 + 1./2).*(cos(AGL)./2 + 1./2).*((RGI.^2 - 2.*cos(AGI - AGL).*RGI.*RGL + RGL.^2)./(cos(AGI./2).^4.*cos(AGL./2).^4)).^(1./2)];
% 
% ALI(:,:,1) = [-2.*atan((RGL - RGI + RGI.*(cos(AGI) + 1) - RGL.*(cos(AGL) + 1) + (cos(AGI)./2 + 1./2).*(cos(AGL)./2 + 1./2).*((RGI.^2 - 2.*cos(AGI - AGL).*RGI.*RGL + RGL.^2)./(cos(AGI./2).^4.*cos(AGL./2).^4)).^(1./2))./(RGI.*sin(AGI) - RGL.*sin(AGL)))];
% ALI(:,:,2) = [+2.*atan((RGI - RGL - RGI.*(cos(AGI) + 1) + RGL.*(cos(AGL) + 1) + (cos(AGI)./2 + 1./2).*(cos(AGL)./2 + 1./2).*((RGI.^2 - 2.*cos(AGI - AGL).*RGI.*RGL + RGL.^2)./(cos(AGI./2).^4.*cos(AGL./2).^4)).^(1./2))./(RGI.*sin(AGI) - RGL.*sin(AGL)))];
% 
% ALI      = wrapTo2Pi(ALI);
% ind      = RLI<0;
% ALI(ind) = wrapTo2Pi(ALI(ind) + pi);
% RLI      = abs(RLI);
% 
% RLI      = RLI(:,:,1);
% ALI      = ALI(:,:,1);
