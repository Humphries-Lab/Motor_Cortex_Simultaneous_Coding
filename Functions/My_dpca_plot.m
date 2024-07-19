function My_dpca_plot(Xfull, W, V, plotFunction, varargin)

% dpca_plot(X, W, V, plotFunction, ...) 
% produces a plot of the dPCA results. X is the data matrix, W and V
% are decoder and encoder matrices, plotFunction is a
% pointer to to the function that plots one component (see dpca_plot_default()
% for the template)

% dpca_plot(..., 'PARAM1',val1, 'PARAM2',val2, ...) 
% specifies optional parameter name/value pairs:
%
%  'whichMarg'              - which marginalization each component comes
%                             from. Is provided as an output of the dpca()
%                             function.
%
%  'time'                   - time axis
%
%  'timeEvents'             - time-points that should be marked on each subplot
%
%  'ylims'                  - array of y-axis spans for each
%                             marginalization or a single value to be used
%                             for each marginalization
%
%  'componentsSignif'       - time-periods of significant classification for each
%                             component. See dpca_signifComponents()
%
%  'timeMarginalization'    - if provided, it will be shown on top, and 
%                             irrespective of significance (because
%                             significant classification is not assessed for 
%                             time components)
%
%  'legendSubplot'          - number of the legend subplot
%
%  'marginalizationNames'   - names of each marginalization
%
%  'marginalizationColours' - colours for each marginalization
%
%  'explainedVar'           - structure returned by the dpca_explainedVariance
%
%  'numCompToShow'          - number of components to show on the explained
%                             variance plots (default = 15)
%
%  'X_extra'                - data array used for plotting that can be larger
%                             (i.e. have more conditions) than the one used
%                             for dpca computations
%  'showNonsignificantComponents'
%                           - display non-signficant components when there
%                             are fewer significant components than
%                             subplots

% default input parameters
options = struct('time',           [], ...   
                 'whichMarg',      [], ...
                 'timeEvents',     [], ...
                 'ylims',          [], ...
                 'componentsSignif', [], ...
                 'timeMarginalization', [], ...
                 'legendSubplot',  [], ...
                 'marginalizationNames', [], ...
                 'marginalizationColours', [], ...
                 'explainedVar',   [], ...
                 'numCompToShow',  15, ...
                 'X_extra',        [], ...
                 'showNonsignificantComponents', false);

% read input parameters
optionNames = fieldnames(options);
if mod(length(varargin),2) == 1
	error('Please provide propertyName/propertyValue pairs')
end
for pair = reshape(varargin,2,[])    % pair is {propName; propValue}
	if any(strcmp(pair{1}, optionNames))
        options.(pair{1}) = pair{2};
    else
        error('%s is not a recognized parameter name', pair{1})
	end
end

% can't show more than there is
numCompToShow = min(options.numCompToShow, size(W,2));

X = Xfull(:,:)';
Xcen = bsxfun(@minus, X, mean(X));
XfullCen = bsxfun(@minus, Xfull, mean(X)');
N = size(X, 1);
dataDim = size(Xfull);
Z = Xcen * W;
%!!
%Z = bsxfun(@times, Z, 1./std(Z, [], 1));
%!!

toDisplayMargNames = 0;

% if there are 4 or less marginalizations, split them into rows
if ~isempty(options.whichMarg) && ...
   length(unique(options.whichMarg)) <= 4 && length(unique(options.whichMarg)) > 1

    % time marginalization, if specified, goes on top
    if ~isempty(options.timeMarginalization)
        margRowSeq = [options.timeMarginalization setdiff(1:max(options.whichMarg), options.timeMarginalization)];
    else
        margRowSeq = 1:max(options.whichMarg);
    end
    
    componentsToPlot = [];
    subplots = [];
    for i=1:length(margRowSeq)
        if ~isempty(options.componentsSignif) && margRowSeq(i) ~= options.timeMarginalization
            % selecting only significant components
            minL = min(length(options.whichMarg), size(options.componentsSignif,1));
            moreComponents = find(options.whichMarg(1:minL) == margRowSeq(i) & ...
                sum(options.componentsSignif(1:minL,:), 2)'~=0, 3);
            if options.showNonsignificantComponents && (length(moreComponents) < 3)
                % Optionally add non-significant components to fill subplots
                moreComponents = [moreComponents setdiff(find(options.whichMarg == margRowSeq(i), 3), moreComponents)];
            end
        else
            moreComponents = find(options.whichMarg == margRowSeq(i), 3);
        end
        componentsToPlot = [componentsToPlot moreComponents];
        subplots = [subplots (i-1)*4+2:(i-1)*4+2 + length(moreComponents) - 1];
    end
else
    % if there are more than 4 marginalizatons
    
    if isempty(options.whichMarg)
        % if there is no info about marginaliations
        componentsToPlot = 1:12;
    else
        % if there is info about marginaliations, select first 3 in each
        uni = unique(options.whichMarg);
        componentsToPlot = [];
        for u = 1:length(uni)
            componentsToPlot = [componentsToPlot find(options.whichMarg==uni(u), 2)];
        end
        componentsToPlot = sort(componentsToPlot);
        if length(componentsToPlot) > 12
            componentsToPlot = componentsToPlot(1:12);
        end
        
        toDisplayMargNames = 1;
    end
    subplots = [2 3 4 6 7 8 10 11 12 14 15 16];
    
    if numCompToShow < 12
        componentsToPlot = componentsToPlot(1:numCompToShow);
        subplots = subplots(1:numCompToShow);
    end
end
    
Zfull = reshape(Z(:,componentsToPlot)', [length(componentsToPlot) dataDim(2:end)]);

if ~isempty(options.X_extra)
    XF = options.X_extra(:,:)';
    XFcen = bsxfun(@minus, XF, mean(X));
    ZF = XFcen * W;
    %!!
    %ZF = bsxfun(@times, ZF, 1./std(ZF, [], 1));
    %!!
    dataDimFull = size(options.X_extra);
    Zfull = reshape(ZF(:,componentsToPlot)', [length(componentsToPlot) dataDimFull(2:end)]);
end

myFig = figure('Position', [0 0 1800 1000]);

% y-axis spans
if isempty(options.ylims)
    options.ylims = max(abs(Zfull(:))) * 1.1;
end
if length(options.ylims) == 1
    if ~isempty(options.whichMarg)
        options.ylims = repmat(options.ylims, [1 max(options.whichMarg)]);
    end
end

% plotting all components as subplots
for c = 1:length(componentsToPlot)
    cc = componentsToPlot(c);
    subplot(4, 4, subplots(c))
    
    if ~isempty(options.componentsSignif)
        signifTrace = options.componentsSignif(cc,:);
    else
        signifTrace = [];
    end
    
    if ~isempty(options.explainedVar)
        thisVar = options.explainedVar.componentVar(cc);
    else
        thisVar = [];
    end
    
    if ~isempty(options.whichMarg)
        thisYlim = options.ylims(options.whichMarg(cc));
        thisMarg = options.whichMarg(cc);
    else
        thisYlim = options.ylims;
        thisMarg = [];
    end
        
    dim = size(Xfull);
    cln = {c};
    for i=2:length(dim)
        cln{i} = ':';
    end

    %thisYlim = 5;
    % plot individual components using provided function
    plotFunction(Zfull(cln{:}), options.time, [-thisYlim thisYlim], ...
        thisVar, cc, options.timeEvents, ...
        signifTrace, thisMarg)
    
    if ismember(subplots(c), [2 6 10 14])
        if subplots(c) == 2 || subplots(c) == 14
            xlabel('Time (s)')
        else
            set(gca, 'XTickLabel', [])
        end
        ylabel('Normalized firing rate (Hz)')
    elseif ismember(subplots(c), [13 14 15 16])
        xlabel('Time (s)')
        set(gca, 'YTickLabel', [])
    else
        set(gca, 'XTickLabel', [])
        set(gca, 'YTickLabel', [])
    end
    
    if toDisplayMargNames && ~isempty(options.marginalizationNames)
        xx = xlim;
        yy = ylim;
        text(xx(1)+(xx(2)-xx(1))*0.1, yy(2)-(yy(2)-yy(1))*0.1, options.marginalizationNames(thisMarg))
    end
end 
%% plot IC vs speed


% colours for marginalizations
if isempty(options.marginalizationColours)
    if ~isempty(options.explainedVar)
        L = length(options.explainedVar.totalMarginalizedVar);
        options.marginalizationColours = lines(L);
    elseif ~isempty(options.whichMarg)
        L = length(unique(options.whichMarg));
        options.marginalizationColours = lines(L);
    else
        options.marginalizationColours = [];
    end
end

% red-to-blue colormap
r = [5 48 97]/256;       %# end
w = [.95 .95 .95];       %# middle
b = [103 0 31]/256;      %# start
c1 = zeros(128,3);
c2 = zeros(128,3);
for i=1:3
    c1(:,i) = linspace(r(i), w(i), 128);
    c2(:,i) = linspace(w(i), b(i), 128);
end
redBlue256 = [c1;c2];

colormap([options.marginalizationColours; redBlue256])

% if there are four marginalizations or less, display labels
if ~isempty(options.whichMarg) && ...
   length(unique(options.whichMarg)) <= 4 && length(unique(options.whichMarg)) > 1 ...
   && ~isempty(options.marginalizationNames)
   
    offsetX = 0.31;
    yposs = [0.9 0.65 0.45 0.25];

    for m = intersect(1:4, unique(options.whichMarg(componentsToPlot)))        
        row = find(margRowSeq == m, 1);
        subplot(4,4,(row-1)*4+2)
        pos = get(gca, 'Position');
        
        annotation('rectangle', [offsetX-0.005 pos(2) 0.015 pos(4)], ...
            'EdgeColor', 'none', 'FaceColor', options.marginalizationColours(m,:));
        
        annotation('textarrow', offsetX*[1 1], yposs(row)*[1 1], ...
            'string', options.marginalizationNames{m}, ...
            'HeadStyle', 'none', 'LineStyle', 'none', ...
            'TextRotation', 90);
    end
end


% cumulative explained variance
if ~isempty(options.explainedVar)
    axCum = subplot(4,4,1);
    hold on

%     % show signal variance if it's provided
%     if isfield(options.explainedVar, 'cumulativePCA_signal')
%         plot(1:numCompToShow, options.explainedVar.cumulativePCA_signal(1:numCompToShow), ...
%             '--k', 'LineWidth', 1)
%         plot(1:numCompToShow, options.explainedVar.cumulativeDPCA_signal(1:numCompToShow), ...
%             '--r', 'LineWidth', 1)
%         yy = [options.explainedVar.cumulativePCA_signal(1:numCompToShow) ...
%               options.explainedVar.cumulativeDPCA_signal(1:numCompToShow)];
%     end
    
    plot(1:numCompToShow, options.explainedVar.cumulativePCA(1:numCompToShow), ...
        '.-k', 'LineWidth', 1, 'MarkerSize', 15);
    if ~isempty(options.explainedVar.cumulativeDPCA)
        plot(1:numCompToShow, options.explainedVar.cumulativeDPCA(1:numCompToShow), ...
             '.-r', 'LineWidth', 1, 'MarkerSize', 15);
    end
    %yy = [options.explainedVar.cumulativePCA(1:numCompToShow) ...
    %    options.explainedVar.cumulativeDPCA(1:numCompToShow)];
    ylabel({'Explained variance (%)'})
        
    if isfield(options.explainedVar, 'totalVar_signal')
        plot([0 numCompToShow+1], options.explainedVar.totalVar_signal/options.explainedVar.totalVar*100*[1 1], 'k--')
    end
           
    %axis([0 numCompToShow+1 floor(min(yy-5)/10)*10 min(ceil(max(yy+5)/10)*10, 100)])
    axis([0 numCompToShow+1 0 100])
    xlabel('Component')
    if ~isempty(options.explainedVar.cumulativeDPCA)
        legend({'PCA', 'dPCA'}, 'Location', 'SouthEast');
    else
        legend({'PCA'}, 'Location', 'SouthEast');
    end
    legend boxoff
end




