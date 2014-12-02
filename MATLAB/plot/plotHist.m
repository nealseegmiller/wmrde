function [N,bins,h_hist,h_bars] = plotHist(x,bins,conf,h_axis)
%plot histogram
%INPUTS
%x:         n x 1 vector
%bins:      bins, or number of bins
%conf:      confidence for error bars, [] if none
%h_axis:    handle to existing axes
%OUTPUTS

if nargin < 4
    h_axis = [];
    if nargin < 3
        conf = [];
    end
end

if isempty(h_axis)
    set(figure,'name','histogram')
    ylabel('count')
else
    axes(h_axis)
end




if numel(bins) == 1
    bins = linspace(min(x),max(x),bins);
end

N = hist(x,bins);
h_before = get(gca,'Children');
hist(x,bins);
h_after = get(gca,'Children');
h_hist = setdiff(h_after,h_before);

h_bars = [];

if ~isempty(conf)
    sig = sqrt(chi2inv(conf,1));
    stdev = sqrt(x'*x/numel(x)); %assumes mean is zero
    h_bars(1) = plot(-sig*stdev*[1 1], [0 max(N)]);
    h_bars(2) = plot( sig*stdev*[1 1], [0 max(N)]);
end

if isempty(h_axis)
    make_legible(14)
end



