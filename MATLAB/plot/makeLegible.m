function makeLegible(fontsize,linewidth)
%make current figure more legible by:
%-increasing font size
%-increasing line width (Optional)
%-setting background color to white

set(get(gca,'XLabel'),'FontSize',fontsize)
set(get(gca,'YLabel'),'FontSize',fontsize)
set(get(gca,'ZLabel'),'FontSize',fontsize)
set(get(gca,'Title'),'FontSize',fontsize)
%slightly smaller font size for axis tick labels and legend
set(gca,'FontSize',fontsize-2)
set(legend,'FontSize',fontsize-2) %okay if no legend exists?

% set(gca,'yminortick','on','xminortick','on','zminortick','on');

if nargin > 1
    %set line width for children that have 'LineWidth' property
    h = get(gca,'Children');
    hasLineWidth = isprop(h,'LineWidth');
    set(h(hasLineWidth),'LineWidth',linewidth);
end

set(gcf,'Color','white')
