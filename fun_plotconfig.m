function fun_plotconfig(a,r)

[x,y] = pol2cart(a,r);
x = x';
y = y';
a = a';
r = r';

n = length(a);
R = max(r);
figure
polarplot(a,r,'.r','markersize',30)
set(gca,'fontname','times','fontsize',40)
pax = gca;
pax.ThetaAxisUnits = 'radians';
rlim([0,R])
rticks([0 R/2 R])
rticklabels({' ',' ',' '})



% figure;
% polarplot(a,r,'.r','markersize',30)
hold on

% set(gca,'fontname','times','fontsize',30)

ind = nchoosek(1:n,2);

d   = floor(fun_norm([x(ind(:,1))-x(ind(:,2)),y(ind(:,1))-y(ind(:,2))])*10^10);

ind_min = find(d <12);

ind = ind(ind_min,:);

for i0 = 1:size(ind,1)
    
    polarplot([a(ind(i0,1));a(ind(i0,2))],[r(ind(i0,1));r(ind(i0,2))],'--r','linewidth',1)
    
end



end

function y = fun_norm(x)

y = sqrt(sum(x.^2,2));

end