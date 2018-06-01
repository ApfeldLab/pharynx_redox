function plotWithRegions(regions, varargin)
%UNTITLED Plot the given data plus the given regions
%   Pass in a regions struct (see Constants.regions), and then whatever
%   else you normally pass in to plot.
    plot(varargin{:});
    fn = fieldnames(regions);
    colors = ['k', 'r'];
    for k=1:numel(fn)
        if ~isequal(fn{k}, 'medial_axis')
            c = colors(mod(k,numel(colors))+1);
            label = fn{k};
            vline(regions.(label), c, label);
        end
    end    
end