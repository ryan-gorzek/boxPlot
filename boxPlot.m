function [xCoordinates, lgdObject] = boxPlot(inputData, NameValueArgs)
% BOXPLOT Plot customizable boxplot from vector or matrix input.
%
% MIT License
% Copyright (c) 2025 Ryan Gorzek
% https://github.com/ryan-gorzek/boxPlot/blob/main/LICENSE
% https://opensource.org/licenses/MIT
%
% Dependencies: none
% 
% Syntax:
%   xCoordinates = boxPlot(inputData)
%   [xCoordinates, lgdObject] = boxPlot(inputData, Name, Value)
%
% Description:
%   boxPlot creates highly customizable boxplots with support for grouping,
%   custom styling, outlier handling, and data point visualization.
%
% Input Arguments:
%
%   inputData -- Numeric vector or matrix of data for boxplot. 
%                If inputData is a matrix, each column will be plotted as a box.
%                If inputData is a vector, specify categorical name-value 
%                argument inputLabels to plot multiple boxes. 
%
% Name-Value Arguments:
%
%   inputLabels (default: []) 
%       Vector of categorical labels for vector input.
%
%   groupSize (default: 1) 
%       Scalar specifying the number of boxes per group.
%
%   labelGroups (default: false) 
%       Logical specifying whether to shrink x-axis labels to one per group.
%
%   boxLabels (default: [1:nBoxes]) 
%       Cell array of strings specifying x-axis labels.
%
%   boxColors (default: {[0.7, 0.7, 0.7]}) 
%       Cell array of RGB vectors specifying box fill colors.
%
%   boxAlpha (default: 1) 
%       Scalar specifying box face transparency (0-1).
%
%   boxEdgeColors (default: {[0.0, 0.0, 0.0]}) 
%       Cell array of RGB vectors specifying box edge colors.
%
%   boxEdgeWidth (default: 1) 
%       Scalar specifying box edge line width.
%
%   boxEdgeStyle (default: "-") 
%       String specifying box edge line style.
%
%   boxEdgeAlpha (default: 1) 
%       Scalar specifying box edge transparency (0-1).
%
%   boxSpacing (default: 0.75 for ungrouped, 1 for grouped) 
%       Scalar specifying spacing between boxes.
%
%   boxCurvature (default: [0, 0]) 
%       Two-element vector specifying box corner curvature.
%
%   boxCoordinates (default: [])
%       Vector specifying custom x-axis coordinates for each box.
%
%   boxOrientation (default: "vertical")
%       String specifying box orientation ("vertical" or "horizontal").
%
%   medianColors (default: {[0.0, 0.0, 0.0]}) 
%       Cell array of RGB vectors specifying median line colors.
%
%   medianWidth (default: 2) 
%       Scalar specifying median line width.
%
%   medianStyle (default: "-") 
%       String specifying median line style.
%
%   medianAlpha (default: 1) 
%       Scalar specifying median line transparency (0-1).
%
%   whiskerColors (default: {[0.0, 0.0, 0.0]}) 
%       Cell array of RGB vectors specifying whisker colors.
%
%   whiskerWidth (default: 1) 
%       Scalar specifying whisker line width.
%
%   whiskerStyle (default: "-") 
%       String specifying whisker line style.
%
%   whiskerAlpha (default: 1) 
%       Scalar specifying whisker transparency (0-1).
%
%   outlierColors (default: matches boxColors) 
%       Cell array of RGB vectors specifying outlier marker colors.
%
%   outlierSize (default: 25) 
%       Scalar specifying outlier marker size.
%
%   outlierStyle (default: "o") 
%       String specifying outlier marker style.
%
%   outlierAlpha (default: 1) 
%       Scalar specifying outlier transparency (0-1).
%
%   outlierJitter (default: "none") 
%       String specifying outlier jitter method ("none", "rand", "randn").
%
%   plotPoints (default: false) 
%       Logical or numeric array specifying whether to plot individual points.
%
%   pointColors (default: {[0.0, 0.0, 0.0]}) 
%       Cell array of RGB vectors specifying point colors.
%
%   pointSize (default: 50) 
%       Scalar specifying point marker size.
%
%   pointStyle (default: ".") 
%       String specifying point marker style.
%
%   pointAlpha (default: 1) 
%       Scalar specifying point transparency (0-1).
%
%   pointJitter (default: "none") 
%       String specifying point jitter method ("none", "rand", "randn").
%
%   jitterBound (default: 0.75) 
%       Scalar specifying maximum jitter displacement.
%
%   plotLines (default: false) 
%       Logical or numeric array specifying whether to connect points with lines.
%
%   lineColors (default: {[0.0, 0.0, 0.0]}) 
%       Cell array of RGB vectors specifying line colors.
%
%   lineWidth (default: 0.5) 
%       Scalar specifying line width.
%
%   lineStyle (default: "-") 
%       String specifying line style.
%
%   lineAlpha (default: 1) 
%       Scalar specifying line transparency (0-1).
%
%   plotLegend (default: false) 
%       Logical specifying whether to plot a legend.
%
%   lgdLabels (default: {}) 
%       Cell array of strings specifying legend labels.
%
%   lgdColors (default: {}) 
%       Cell array of RGB vectors specifying legend colors.
%
%   lgdColumns (default: 1) 
%       Scalar specifying number of legend columns.
%
%   lgdOrientation (default: "horizontal") 
%       String specifying legend orientation ("horizontal" or "vertical").
%
%   lgdBox (default: "off") 
%       String specifying legend box visibility ("on" or "off").
%
%   lgdFontSize (default: 12) 
%       Scalar specifying legend font size.
%
%   lgdLineHeight (default: 1) 
%       Scalar specifying legend marker line height.
%
%   lgdLineWidth (default: 1) 
%       Scalar specifying legend marker line width.
%
%   lgdLocation (default: "northeast") 
%       String specifying legend location.
%
%   lgdPosition (default: [0, 0, 0, 0]) 
%       Four-element vector specifying custom legend position.
%
% Output Arguments:
%
%   xCoordinates -- Vector of x-axis coordinates for each plotted box.
%
%   lgdObject -- Legend object (empty if no legend plotted).
%
% Example:
%   data = randn(100, 3);
%   boxPlot(data, 'boxColors', {[0.8 0.2 0.2], [0.2 0.8 0.2], [0.2 0.2 0.8]});
%

arguments
    inputData (:,:) {mustBeNumeric}
    
    NameValueArgs.inputLabels {mustBeVector, mustBeNumeric} = ...
        reshape(repmat(1:size(inputData,2), [size(inputData,1), 1]), [], 1)
    
    NameValueArgs.groupSize (1,1) {mustBeNumeric} = 1
    NameValueArgs.labelGroups (1,1) logical = false
    
    NameValueArgs.boxLabels {mustBeA(NameValueArgs.boxLabels, "cell")} = {}
    NameValueArgs.boxColors {mustBeA(NameValueArgs.boxColors, "cell")} = {}
    NameValueArgs.boxAlpha (1,1) {mustBeInRange(NameValueArgs.boxAlpha, 0, 1)} = 1
    NameValueArgs.boxEdgeColors {mustBeA(NameValueArgs.boxEdgeColors, "cell")} = {}
    NameValueArgs.boxEdgeWidth (1,1) {mustBeNumeric} = 1
    NameValueArgs.boxEdgeStyle (1,1) string = "-"
    NameValueArgs.boxEdgeAlpha (1,1) {mustBeInRange(NameValueArgs.boxEdgeAlpha, 0, 1)} = 1
    NameValueArgs.boxSpacing (1,1) {mustBeNumeric} = -1
    NameValueArgs.boxCurvature (1,2) double = [0, 0]
    NameValueArgs.boxCoordinates (:,1) {mustBeNumeric} = []
    NameValueArgs.boxOrientation {mustBeMember(NameValueArgs.boxOrientation, ["vertical", "horizontal"])} = "vertical"
    
    NameValueArgs.medianColors (1,:) {mustBeA(NameValueArgs.medianColors, "cell")} = {}
    NameValueArgs.medianWidth (1,1) {mustBeNumeric} = 2
    NameValueArgs.medianStyle (1,1) string = "-"
    NameValueArgs.medianAlpha (1,1) {mustBeInRange(NameValueArgs.medianAlpha, 0, 1)} = 1
    
    NameValueArgs.whiskerColors (1,:) {mustBeA(NameValueArgs.whiskerColors, "cell")} = {}
    NameValueArgs.whiskerWidth (1,1) {mustBeNumeric} = 1
    NameValueArgs.whiskerStyle (1,1) string = "-"
    NameValueArgs.whiskerAlpha (1,1) {mustBeInRange(NameValueArgs.whiskerAlpha, 0, 1)} = 1
    
    NameValueArgs.outlierColors (1,:) {mustBeA(NameValueArgs.outlierColors, "cell")} = {}
    NameValueArgs.outlierSize (1,1) {mustBeNumeric} = 25
    NameValueArgs.outlierStyle (1,1) string = "o"
    NameValueArgs.outlierAlpha (1,1) {mustBeInRange(NameValueArgs.outlierAlpha, 0, 1)} = 1
    NameValueArgs.outlierJitter {mustBeMember(NameValueArgs.outlierJitter, ["none", "rand", "randn"])} = "none"
    
    NameValueArgs.plotPoints {mustBeNumericOrLogical} = false
    NameValueArgs.pointColors {mustBeA(NameValueArgs.pointColors, "cell")} = {}
    NameValueArgs.pointSize (1,1) {mustBeNumeric} = 50
    NameValueArgs.pointStyle (1,1) string = "."
    NameValueArgs.pointAlpha (1,1) {mustBeInRange(NameValueArgs.pointAlpha, 0, 1)} = 1
    NameValueArgs.pointJitter {mustBeMember(NameValueArgs.pointJitter, ["none", "rand", "randn"])} = "none"
    
    NameValueArgs.jitterBound (1,1) {mustBeNumeric} = 0.75
    
    NameValueArgs.plotLines {mustBeNumericOrLogical} = false
    NameValueArgs.lineColors {mustBeA(NameValueArgs.lineColors, "cell")} = {}
    NameValueArgs.lineWidth (1,1) {mustBeNumeric} = 0.5
    NameValueArgs.lineStyle (1,1) string = "-"
    NameValueArgs.lineAlpha (1,1) {mustBeInRange(NameValueArgs.lineAlpha, 0, 1)} = 1
    
    NameValueArgs.plotLegend (1,1) logical = false
    NameValueArgs.lgdLabels (1,:) {mustBeA(NameValueArgs.lgdLabels, "cell")} = {}
    NameValueArgs.lgdColors (1,:) {mustBeA(NameValueArgs.lgdColors, "cell")} = {}
    NameValueArgs.lgdColumns (1,1) {mustBeNumeric} = 1
    NameValueArgs.lgdOrientation {mustBeMember(NameValueArgs.lgdOrientation, ["vertical", "horizontal"])} = "horizontal"
    NameValueArgs.lgdBox {mustBeMember(NameValueArgs.lgdBox, ["on", "off"])} = "off"
    NameValueArgs.lgdFontSize (1,1) {mustBeNumeric} = 12
    NameValueArgs.lgdLineHeight (1,1) {mustBeNumeric} = 1
    NameValueArgs.lgdLineWidth (1,1) {mustBeNumeric} = 1
    NameValueArgs.lgdLocation (1,1) string = "northeast"
    NameValueArgs.lgdPosition (1,4) {mustBeInRange(NameValueArgs.lgdPosition, 0, 1)} = [0, 0, 0, 0]
end

%% Get dataset dimensions
nBoxes = numel(unique(NameValueArgs.inputLabels));
nGroups = nBoxes / NameValueArgs.groupSize;
nSamples = size(inputData, 1);
nMissing = nnz(all(isnan(inputData), 2));

% Throw error if number of groups is not an integer
if rem(nGroups, 1) ~= 0
    error('Number of input boxes is not divisible by number of groups.'); 
end

% Set default boxSpacing based on grouping
if NameValueArgs.boxSpacing == -1
    if NameValueArgs.groupSize == 1
        NameValueArgs.boxSpacing = 0.75;
    else
        NameValueArgs.boxSpacing = 1;
    end
end

% Reshape input data into a column vector and get unique labels
inputData = reshape(inputData, [], 1);
uniqueLabels = sort(unique(NameValueArgs.inputLabels), 'ascend')';

%% Set default box labels
if isempty(NameValueArgs.boxLabels) && NameValueArgs.groupSize == 1
    NameValueArgs.boxLabels = cellstr(string(1:nBoxes));
    
elseif isempty(NameValueArgs.boxLabels) && ...
       NameValueArgs.groupSize > 1 && ...
       NameValueArgs.labelGroups == false
    NameValueArgs.boxLabels = cellstr(string(1:NameValueArgs.groupSize));
    
elseif isempty(NameValueArgs.boxLabels) && NameValueArgs.labelGroups == true
    NameValueArgs.boxLabels = cellstr(string(1:nGroups));
end

%% Set default colors
% Box fill colors
if isempty(NameValueArgs.boxColors)
    NameValueArgs.boxColors = repmat({[0.7, 0.7, 0.7]}, [1, nBoxes]);
elseif numel(NameValueArgs.boxColors) == 1
    NameValueArgs.boxColors = repmat(NameValueArgs.boxColors, [1, nBoxes]);
elseif numel(NameValueArgs.boxColors) == NameValueArgs.groupSize
    NameValueArgs.boxColors = repmat(NameValueArgs.boxColors, [1, nGroups]);
end

% Box edge colors
if isempty(NameValueArgs.boxEdgeColors)
    NameValueArgs.boxEdgeColors = repmat({[0.0, 0.0, 0.0]}, [1, nBoxes]);
elseif numel(NameValueArgs.boxEdgeColors) == 1
    NameValueArgs.boxEdgeColors = repmat(NameValueArgs.boxEdgeColors, [1, nBoxes]);
elseif numel(NameValueArgs.boxEdgeColors) == NameValueArgs.groupSize
    NameValueArgs.boxEdgeColors = repmat(NameValueArgs.boxEdgeColors, [1, nGroups]);
end

% Whisker colors
if isempty(NameValueArgs.whiskerColors)
    NameValueArgs.whiskerColors = repmat({[0.0, 0.0, 0.0]}, [1, nBoxes]);
elseif numel(NameValueArgs.whiskerColors) == 1
    NameValueArgs.whiskerColors = repmat(NameValueArgs.whiskerColors, [1, nBoxes]);
elseif numel(NameValueArgs.whiskerColors) == NameValueArgs.groupSize
    NameValueArgs.whiskerColors = repmat(NameValueArgs.whiskerColors, [1, nGroups]);
end

% Median line colors
if isempty(NameValueArgs.medianColors)
    NameValueArgs.medianColors = repmat({[0.0, 0.0, 0.0]}, [1, nBoxes]);
elseif numel(NameValueArgs.medianColors) == 1
    NameValueArgs.medianColors = repmat(NameValueArgs.medianColors, [1, nBoxes]);
elseif numel(NameValueArgs.medianColors) == NameValueArgs.groupSize
    NameValueArgs.medianColors = repmat(NameValueArgs.medianColors, [1, nGroups]);
end

% Outlier colors (default to box colors)
if isempty(NameValueArgs.outlierColors)
    NameValueArgs.outlierColors = NameValueArgs.boxColors;
elseif numel(NameValueArgs.outlierColors) == 1
    NameValueArgs.outlierColors = repmat(NameValueArgs.outlierColors, [1, nBoxes]);
elseif numel(NameValueArgs.outlierColors) == NameValueArgs.groupSize
    NameValueArgs.outlierColors = repmat(NameValueArgs.outlierColors, [1, nGroups]);
end

%% Set default point colors
if isempty(NameValueArgs.pointColors)
    NameValueArgs.pointColors = repmat({[0.0, 0.0, 0.0]}, [nSamples, nBoxes]);
elseif numel(NameValueArgs.pointColors) == 1
    NameValueArgs.pointColors = repmat(NameValueArgs.pointColors, [nSamples, nBoxes]);
elseif isvector(NameValueArgs.pointColors) && ...
       any(size(NameValueArgs.pointColors) == NameValueArgs.groupSize)
    NameValueArgs.pointColors = repmat(NameValueArgs.pointColors, [nSamples, nGroups]);
elseif isvector(NameValueArgs.pointColors) && ...
       any(size(NameValueArgs.pointColors) == nBoxes)
    NameValueArgs.pointColors = repmat(NameValueArgs.pointColors, [nSamples, 1]);
elseif all(ismember(size(NameValueArgs.pointColors), [nSamples - nMissing, NameValueArgs.groupSize]))
    NameValueArgs.pointColors = repmat(NameValueArgs.pointColors, [1, nGroups]);
end

% Handle "default" color stream specification for points
point_defaultIdx = cell2mat(cellfun(@(x) strcmp(x, "default"), ...
                                    NameValueArgs.pointColors, 'UniformOutput', false));
if any(point_defaultIdx, 'all')
    defaultColors = repmat(num2cell(colororder, 2), ...
                          [ceil(size(NameValueArgs.pointColors, 1) / 7), ...
                           size(NameValueArgs.pointColors, 2)]);
    defaultColors = defaultColors(1:size(NameValueArgs.pointColors, 1), :);
    NameValueArgs.pointColors(point_defaultIdx) = defaultColors(point_defaultIdx);
end

%% Set default line colors
if nGroups == nBoxes
    nLines = nBoxes - 1;
else
    nLines = nGroups * (NameValueArgs.groupSize - 1);
end

if isempty(NameValueArgs.lineColors)
    NameValueArgs.lineColors = repmat({[0.0, 0.0, 0.0]}, [nSamples, nLines]);
elseif numel(NameValueArgs.lineColors) == 1
    NameValueArgs.lineColors = repmat(NameValueArgs.lineColors, [nSamples, nLines]);
elseif isvector(NameValueArgs.lineColors) && ...
       any(size(NameValueArgs.lineColors) == nGroups) && ...
       nGroups ~= nBoxes
    NameValueArgs.lineColors = NameValueArgs.lineColors(kron(1:nGroups, ...
                                                        ones(nSamples, NameValueArgs.groupSize - 1)));
elseif isvector(NameValueArgs.lineColors) && ...
       any(size(NameValueArgs.lineColors) == nLines)
    NameValueArgs.lineColors = repmat(NameValueArgs.lineColors, [nSamples, 1]);
elseif all(ismember(size(NameValueArgs.lineColors), [nSamples - nMissing, nGroups])) && ...
       nGroups ~= nBoxes
    NameValueArgs.lineColors = NameValueArgs.lineColors(:, kron(1:nGroups, ...
                                                           ones(1, NameValueArgs.groupSize - 1)));
end

% Handle "default" color stream specification for lines
line_defaultIdx = cell2mat(cellfun(@(x) strcmp(x, "default"), ...
                                   NameValueArgs.lineColors, 'UniformOutput', false));
if any(line_defaultIdx, 'all')
    defaultColors = repmat(num2cell(colororder, 2), ...
                          [ceil(size(NameValueArgs.lineColors, 1) / 7), ...
                           size(NameValueArgs.lineColors, 2)]);
    defaultColors = defaultColors(1:size(NameValueArgs.lineColors, 1), :);
    NameValueArgs.lineColors(line_defaultIdx) = defaultColors(line_defaultIdx);
end

%% Validate box labels against groups
if NameValueArgs.labelGroups == true && ...
   numel(NameValueArgs.boxLabels) > nGroups
    error(['Number of box labels exceeds number of groups. ' ...
           'Did you mean to specify labelGroups = false?']);
elseif NameValueArgs.labelGroups == true && ...
       numel(NameValueArgs.boxLabels) < nGroups
    error('Insufficient number of box labels for number of groups.');
elseif NameValueArgs.labelGroups == false && ...
       NameValueArgs.groupSize ~= 1 && ...
       numel(NameValueArgs.boxLabels) == nGroups && ...
       numel(NameValueArgs.boxColors) == NameValueArgs.groupSize
    error(['Insufficient number of box labels for number of groups. ' ...
           'Did you mean to specify labelGroups = true?']);
end

%% Generate x-coordinates for boxes
if ~isempty(NameValueArgs.boxCoordinates)
    % Use custom coordinates if provided
    if numel(NameValueArgs.boxCoordinates) ~= nBoxes
        error('Number of custom box coordinates must equal number of boxes.');
    end
    xCoordinates = NameValueArgs.boxCoordinates';
elseif nGroups == 1
    xCoordinates = (0.70 : 0.65 : 0.70 + (0.65 * nBoxes) - 0.65) .* NameValueArgs.boxSpacing;
else
    xCoordinates = 0.70 : 0.55 : 0.70 + (0.55 * (nBoxes / nGroups) - 0.55); 
    initCoors = 0.70 : 0.55 : 0.70 + (0.55 * (nBoxes / nGroups) - 0.55);
    for grp = 2:nGroups
        xCoordinates = horzcat(xCoordinates, initCoors + xCoordinates(end) + 0.4);
    end
    xCoordinates = xCoordinates .* NameValueArgs.boxSpacing;
end

%% Generate jitter matrix
if strcmp(NameValueArgs.outlierJitter, "rand") || ...
   strcmp(NameValueArgs.pointJitter, "rand")
    jitterMat = rand(size(reshape(inputData, [], nBoxes))) .* (0.5 * NameValueArgs.jitterBound);
    jitterMat = jitterMat - mean(jitterMat, 1);
    
elseif strcmp(NameValueArgs.outlierJitter, "randn") || ...
       strcmp(NameValueArgs.pointJitter, "randn")
    randnMat = randn(size(reshape(inputData, [], nBoxes)));
    jitterMat = (randnMat ./ max(abs(randnMat), [], 1)) .* (0.5 * NameValueArgs.jitterBound);
    jitterMat = jitterMat - mean(jitterMat, 1);
    
else
    jitterMat = zeros(size(reshape(inputData, [], nBoxes)));
end

%% Initialize matrices for axis limits
maxMat = zeros(nBoxes, 3);
minMat = zeros(nBoxes, 3);

%% Determine if horizontal orientation
isHorizontal = strcmp(NameValueArgs.boxOrientation, "horizontal");

%% Plot boxes
for box = uniqueLabels
    % Get current box location and data
    boxNum = find(uniqueLabels == box);
    currData = inputData(NameValueArgs.inputLabels == box, 1);
    
    % Plot box if there are more than 4 data points
    if nnz(~isnan(currData)) > 4
        boxMedian = median(currData, 1, 'omitnan');
        
        upperQuantile = quantile(currData, 0.75);
        lowerQuantile = quantile(currData, 0.25);
        
        maxWhisker = upperQuantile + 1.5 * (upperQuantile - lowerQuantile); 
        upperWhisker = max(currData(currData < maxWhisker & currData >= upperQuantile));
        if isempty(upperWhisker), upperWhisker = maxWhisker; end
        
        minWhisker = lowerQuantile - 1.5 * (upperQuantile - lowerQuantile);
        lowerWhisker = min(currData(currData > minWhisker & currData <= lowerQuantile));
        if isempty(lowerWhisker), lowerWhisker = minWhisker; end
        
        hold on;
        
        if isHorizontal
            % Plot horizontal box
            rectangle('Position',  [lowerQuantile, xCoordinates(boxNum) - 0.25, upperQuantile - lowerQuantile, 0.5], ...
                      'FaceColor', [NameValueArgs.boxColors{boxNum}, NameValueArgs.boxAlpha], ...
                      'EdgeColor', [NameValueArgs.boxEdgeColors{boxNum}, NameValueArgs.boxEdgeAlpha], ...
                      'LineWidth',  NameValueArgs.boxEdgeWidth, ...
                      'LineStyle',  NameValueArgs.boxEdgeStyle, ...
                      'Curvature',  NameValueArgs.boxCurvature, ...
                      'Tag',       'Box');
            
            % Plot median line
            line([boxMedian, boxMedian], ...
                 [xCoordinates(boxNum) - 0.25, xCoordinates(boxNum) + 0.25], ...
                 'Color',     [NameValueArgs.medianColors{boxNum}, NameValueArgs.medianAlpha], ...
                 'LineWidth', NameValueArgs.medianWidth, ...
                 'LineStyle', NameValueArgs.medianStyle, ...
                 'Tag',       'Median');
            
            % Plot lower whisker
            line([lowerQuantile, lowerWhisker], ...
                 [xCoordinates(boxNum), xCoordinates(boxNum)], ...
                 'Color',     [NameValueArgs.whiskerColors{boxNum}, NameValueArgs.whiskerAlpha], ...
                 'LineWidth', NameValueArgs.whiskerWidth, ...
                 'LineStyle', NameValueArgs.whiskerStyle, ...
                 'Tag',       'Lower Whisker');
            
            % Plot upper whisker
            line([upperQuantile, upperWhisker], ...
                 [xCoordinates(boxNum), xCoordinates(boxNum)], ...
                 'Color',     [NameValueArgs.whiskerColors{boxNum}, NameValueArgs.whiskerAlpha], ...
                 'LineWidth', NameValueArgs.whiskerWidth, ...
                 'LineStyle', NameValueArgs.whiskerStyle, ...
                 'Tag',       'Upper Whisker');
            
            % Plot lower whisker bar
            line([lowerWhisker, lowerWhisker], ...
                 [xCoordinates(boxNum) - 0.1, xCoordinates(boxNum) + 0.1], ...
                 'Color',     [NameValueArgs.whiskerColors{boxNum}, NameValueArgs.whiskerAlpha], ...
                 'LineWidth', NameValueArgs.whiskerWidth, ...
                 'LineStyle', NameValueArgs.whiskerStyle, ...
                 'Tag',       'Lower Whisker Bar');
            
            % Plot upper whisker bar
            line([upperWhisker, upperWhisker], ...
                 [xCoordinates(boxNum) - 0.1, xCoordinates(boxNum) + 0.1], ...
                 'Color',     [NameValueArgs.whiskerColors{boxNum}, NameValueArgs.whiskerAlpha], ...
                 'LineWidth', NameValueArgs.whiskerWidth, ...
                 'LineStyle', NameValueArgs.whiskerStyle, ...
                 'Tag',       'Upper Whisker Bar');
        else
            % Plot vertical box (original orientation)
            rectangle('Position',  [xCoordinates(boxNum) - 0.25, lowerQuantile, 0.5, upperQuantile - lowerQuantile], ...
                      'FaceColor', [NameValueArgs.boxColors{boxNum}, NameValueArgs.boxAlpha], ...
                      'EdgeColor', [NameValueArgs.boxEdgeColors{boxNum}, NameValueArgs.boxEdgeAlpha], ...
                      'LineWidth',  NameValueArgs.boxEdgeWidth, ...
                      'LineStyle',  NameValueArgs.boxEdgeStyle, ...
                      'Curvature',  NameValueArgs.boxCurvature, ...
                      'Tag',       'Box');
            
            % Plot median line
            line([xCoordinates(boxNum) - 0.25, xCoordinates(boxNum) + 0.25], ...
                 [boxMedian, boxMedian], ...
                 'Color',     [NameValueArgs.medianColors{boxNum}, NameValueArgs.medianAlpha], ...
                 'LineWidth', NameValueArgs.medianWidth, ...
                 'LineStyle', NameValueArgs.medianStyle, ...
                 'Tag',       'Median');
            
            % Plot lower whisker
            line([xCoordinates(boxNum), xCoordinates(boxNum)], ...
                 [lowerQuantile, lowerWhisker], ...
                 'Color',     [NameValueArgs.whiskerColors{boxNum}, NameValueArgs.whiskerAlpha], ...
                 'LineWidth', NameValueArgs.whiskerWidth, ...
                 'LineStyle', NameValueArgs.whiskerStyle, ...
                 'Tag',       'Lower Whisker');
            
            % Plot upper whisker
            line([xCoordinates(boxNum), xCoordinates(boxNum)], ...
                 [upperQuantile, upperWhisker], ...
                 'Color',     [NameValueArgs.whiskerColors{boxNum}, NameValueArgs.whiskerAlpha], ...
                 'LineWidth', NameValueArgs.whiskerWidth, ...
                 'LineStyle', NameValueArgs.whiskerStyle, ...
                 'Tag',       'Upper Whisker');
            
            % Plot lower whisker bar
            line([xCoordinates(boxNum) - 0.1, xCoordinates(boxNum) + 0.1], ...
                 [lowerWhisker, lowerWhisker], ...
                 'Color',     [NameValueArgs.whiskerColors{boxNum}, NameValueArgs.whiskerAlpha], ...
                 'LineWidth', NameValueArgs.whiskerWidth, ...
                 'LineStyle', NameValueArgs.whiskerStyle, ...
                 'Tag',       'Lower Whisker Bar');
            
            % Plot upper whisker bar
            line([xCoordinates(boxNum) - 0.1, xCoordinates(boxNum) + 0.1], ...
                 [upperWhisker, upperWhisker], ...
                 'Color',     [NameValueArgs.whiskerColors{boxNum}, NameValueArgs.whiskerAlpha], ...
                 'LineWidth', NameValueArgs.whiskerWidth, ...
                 'LineStyle', NameValueArgs.whiskerStyle, ...
                 'Tag',       'Upper Whisker Bar');
        end
        
        % Plot upper outliers
        if any(currData > upperWhisker)
            if isHorizontal
                scatter(currData(currData > upperWhisker), ...
                        xCoordinates(boxNum) - jitterMat(currData > upperWhisker, boxNum), ...
                        NameValueArgs.outlierSize, ...
                        'MarkerFaceColor', NameValueArgs.outlierColors{boxNum}, ...
                        'MarkerEdgeColor', NameValueArgs.outlierColors{boxNum}, ...
                        'Marker',          NameValueArgs.outlierStyle, ...
                        'MarkerFaceAlpha', NameValueArgs.outlierAlpha, ...
                        'MarkerEdgeAlpha', NameValueArgs.outlierAlpha, ...
                        'Tag',             'Outlier');
            else
                scatter(xCoordinates(boxNum) - jitterMat(currData > upperWhisker, boxNum), ...
                        currData(currData > upperWhisker), ...
                        NameValueArgs.outlierSize, ...
                        'MarkerFaceColor', NameValueArgs.outlierColors{boxNum}, ...
                        'MarkerEdgeColor', NameValueArgs.outlierColors{boxNum}, ...
                        'Marker',          NameValueArgs.outlierStyle, ...
                        'MarkerFaceAlpha', NameValueArgs.outlierAlpha, ...
                        'MarkerEdgeAlpha', NameValueArgs.outlierAlpha, ...
                        'Tag',             'Outlier');
            end
            upperOutliers = currData(currData > upperWhisker);
        else
            upperOutliers = nan;
        end
        
        % Plot lower outliers
        if any(currData < lowerWhisker)
            if isHorizontal
                scatter(currData(currData < lowerWhisker), ...
                        xCoordinates(boxNum) - jitterMat(currData < lowerWhisker, boxNum), ...
                        NameValueArgs.outlierSize, ...
                        'MarkerFaceColor', NameValueArgs.outlierColors{boxNum}, ...
                        'MarkerEdgeColor', NameValueArgs.outlierColors{boxNum}, ...
                        'Marker',          NameValueArgs.outlierStyle, ...
                        'MarkerFaceAlpha', NameValueArgs.outlierAlpha, ...
                        'MarkerEdgeAlpha', NameValueArgs.outlierAlpha, ...
                        'Tag',             'Outlier');
            else
                scatter(xCoordinates(boxNum) - jitterMat(currData < lowerWhisker, boxNum), ...
                        currData(currData < lowerWhisker), ...
                        NameValueArgs.outlierSize, ...
                        'MarkerFaceColor', NameValueArgs.outlierColors{boxNum}, ...
                        'MarkerEdgeColor', NameValueArgs.outlierColors{boxNum}, ...
                        'Marker',          NameValueArgs.outlierStyle, ...
                        'MarkerFaceAlpha', NameValueArgs.outlierAlpha, ...
                        'MarkerEdgeAlpha', NameValueArgs.outlierAlpha, ...
                        'Tag',             'Outlier');
            end
            lowerOutliers = currData(currData < lowerWhisker);
        else
            lowerOutliers = nan;
        end
        
        % Store min/max from each box to set axis limits
        maxMat(boxNum, :) = [upperQuantile, upperWhisker, max(upperOutliers)];
        minMat(boxNum, :) = [lowerQuantile, lowerWhisker, min(lowerOutliers)];
        
    elseif nnz(~isnan(currData)) == 0
        % Store min/max from each box to set axis limits
        maxMat(boxNum, :) = repmat(max(currData), [1, 3]);
        minMat(boxNum, :) = repmat(min(currData), [1, 3]);
    
    elseif nnz(~isnan(currData)) <= 4
        boxMedian = median(currData, 1, 'omitnan');
        
        hold on;
        
        if isHorizontal
            % Plot median line
            line([boxMedian, boxMedian], ...
                 [xCoordinates(boxNum) - 0.25, xCoordinates(boxNum) + 0.25], ...
                 'Color',     [NameValueArgs.medianColors{boxNum}, NameValueArgs.medianAlpha], ...
                 'LineWidth', NameValueArgs.medianWidth, ...
                 'LineStyle', NameValueArgs.medianStyle, ...
                 'Tag',       'Median');
            
            % Scatter outliers
            scatter(currData, ...
                    xCoordinates(boxNum) - jitterMat(:, boxNum), ...
                    NameValueArgs.outlierSize, ...
                    'MarkerFaceColor', NameValueArgs.outlierColors{boxNum}, ...
                    'MarkerEdgeColor', NameValueArgs.outlierColors{boxNum}, ...
                    'Marker',          NameValueArgs.outlierStyle, ...
                    'MarkerFaceAlpha', NameValueArgs.outlierAlpha, ...
                    'MarkerEdgeAlpha', NameValueArgs.outlierAlpha, ...
                    'Tag',             'Outlier');
        else
            % Plot median line
            line([xCoordinates(boxNum) - 0.25, xCoordinates(boxNum) + 0.25], ...
                 [boxMedian, boxMedian], ...
                 'Color',     [NameValueArgs.medianColors{boxNum}, NameValueArgs.medianAlpha], ...
                 'LineWidth', NameValueArgs.medianWidth, ...
                 'LineStyle', NameValueArgs.medianStyle, ...
                 'Tag',       'Median');
            
            % Scatter outliers
            scatter(xCoordinates(boxNum) - jitterMat(:, boxNum), ...
                    currData, ...
                    NameValueArgs.outlierSize, ...
                    'MarkerFaceColor', NameValueArgs.outlierColors{boxNum}, ...
                    'MarkerEdgeColor', NameValueArgs.outlierColors{boxNum}, ...
                    'Marker',          NameValueArgs.outlierStyle, ...
                    'MarkerFaceAlpha', NameValueArgs.outlierAlpha, ...
                    'MarkerEdgeAlpha', NameValueArgs.outlierAlpha, ...
                    'Tag',             'Outlier');
        end
        
        % Store min/max from each box to set axis limits
        maxMat(boxNum, :) = repmat(max(currData), [1, 3]);
        minMat(boxNum, :) = repmat(min(currData), [1, 3]);
    end
end

%% Connect groups with lines if specified
if (isnumeric(NameValueArgs.plotLines) || NameValueArgs.plotLines == true) && ...
   nGroups == nBoxes
    currData = reshape(inputData, [], nBoxes);
    
    if isnumeric(NameValueArgs.plotLines)
        lineIdx = NameValueArgs.plotLines; 
    else
        lineIdx = 1:nBoxes; 
    end
    
    connectionNum = 1;
    for connection = 1:numel(lineIdx) - 1
        for sample = find(~all(isnan(currData), 2))'
            sampleIdx = find(~all(isnan(currData), 2)) == sample;
            
            if isHorizontal
                plot(currData(sample, lineIdx(connection:connection + 1)), ...
                     xCoordinates(lineIdx(connection:connection + 1)) - ...
                     jitterMat(sample, lineIdx(connection:connection + 1)), ...
                     'Color',     [NameValueArgs.lineColors{sampleIdx, connectionNum}, ...
                                  NameValueArgs.lineAlpha], ...
                     'LineWidth', NameValueArgs.lineWidth, ...
                     'LineStyle', NameValueArgs.lineStyle, ...
                     'Tag',       'Line');
            else
                plot(xCoordinates(lineIdx(connection:connection + 1)) - ...
                     jitterMat(sample, lineIdx(connection:connection + 1)), ...
                     currData(sample, lineIdx(connection:connection + 1)), ...
                     'Color',     [NameValueArgs.lineColors{sampleIdx, connectionNum}, ...
                                  NameValueArgs.lineAlpha], ...
                     'LineWidth', NameValueArgs.lineWidth, ...
                     'LineStyle', NameValueArgs.lineStyle, ...
                     'Tag',       'Line');
            end
        end
        connectionNum = connectionNum + 1;
    end
    
elseif (isnumeric(NameValueArgs.plotLines) || NameValueArgs.plotLines == true) && ...
       nGroups ~= nBoxes
    currData = reshape(inputData, [], nBoxes);
    
    if isnumeric(NameValueArgs.plotLines)
        lineIdx = NameValueArgs.plotLines;
    else
        lineIdx = 1:nBoxes / nGroups;
    end
    
    connectionNum = 1;
    for group = 1:nGroups
        boxIdx = reshape(1:nBoxes, [], nGroups)'; 
        boxIdx = boxIdx(:, lineIdx);
        
        for connection = 1:numel(boxIdx(group, :)) - 1
            for sample = find(~all(isnan(currData), 2))'
                if isHorizontal
                    plot(currData(sample, boxIdx(group, connection:connection + 1)), ...
                         xCoordinates(boxIdx(group, connection:connection + 1)) - ...
                         jitterMat(sample, boxIdx(group, connection:connection + 1)), ...
                         'Color',     [NameValueArgs.lineColors{find(~all(isnan(currData), 2)) == sample, connectionNum}, ...
                                      NameValueArgs.lineAlpha], ...
                         'LineWidth', NameValueArgs.lineWidth, ...
                         'LineStyle', NameValueArgs.lineStyle, ...
                         'Tag',       'Line');
                else
                    plot(xCoordinates(boxIdx(group, connection:connection + 1)) - ...
                         jitterMat(sample, boxIdx(group, connection:connection + 1)), ...
                         currData(sample, boxIdx(group, connection:connection + 1)), ...
                         'Color',     [NameValueArgs.lineColors{find(~all(isnan(currData), 2)) == sample, connectionNum}, ...
                                      NameValueArgs.lineAlpha], ...
                         'LineWidth', NameValueArgs.lineWidth, ...
                         'LineStyle', NameValueArgs.lineStyle, ...
                         'Tag',       'Line');
                end
            end
            connectionNum = connectionNum + 1;
        end
    end
end

%% Plot individual data points if specified
if (isnumeric(NameValueArgs.plotPoints) || NameValueArgs.plotPoints == true)
    currData = reshape(inputData, [], nBoxes);
    
    if isnumeric(NameValueArgs.plotPoints)
        pointIdx = NameValueArgs.plotPoints;
    else
        pointIdx = 1:nBoxes;
    end
    
    for sample = find(~all(isnan(currData), 2))'
        if isHorizontal
            scatter(currData(sample, pointIdx), ...
                    xCoordinates(pointIdx) - jitterMat(sample, pointIdx), ...
                    NameValueArgs.pointSize, ...
                    vertcat(NameValueArgs.pointColors{find(~all(isnan(currData), 2)) == sample, pointIdx}), ...
                    'Marker',          NameValueArgs.pointStyle, ...
                    'MarkerFaceAlpha', NameValueArgs.pointAlpha, ...
                    'MarkerEdgeAlpha', NameValueArgs.pointAlpha, ...
                    'Tag',             'Point');
        else
            scatter(xCoordinates(pointIdx) - jitterMat(sample, pointIdx), ...
                    currData(sample, pointIdx), ...
                    NameValueArgs.pointSize, ...
                    vertcat(NameValueArgs.pointColors{find(~all(isnan(currData), 2)) == sample, pointIdx}), ...
                    'Marker',          NameValueArgs.pointStyle, ...
                    'MarkerFaceAlpha', NameValueArgs.pointAlpha, ...
                    'MarkerEdgeAlpha', NameValueArgs.pointAlpha, ...
                    'Tag',             'Point');
        end
    end
end

%% Set axis limits and labels
if isHorizontal
    ylim([0.2 * NameValueArgs.boxSpacing, xCoordinates(end) + (0.5 * NameValueArgs.boxSpacing)]);
    
    % Set y-axis tick labels
    xLabPos = [];
    if NameValueArgs.labelGroups == true
        for group = 1:(nBoxes / nGroups):nBoxes
            xLabPos = horzcat(xLabPos, median(xCoordinates(group : group + (nBoxes / nGroups) - 1)));
        end
        set(gca, 'ytick',      xLabPos, ...
                 'yticklabel', NameValueArgs.boxLabels);
    else 
        set(gca, 'ytick',      xCoordinates, ...
                 'yticklabel', NameValueArgs.boxLabels);
    end
else
    xlim([0.2 * NameValueArgs.boxSpacing, xCoordinates(end) + (0.5 * NameValueArgs.boxSpacing)]);
    
    % Set x-axis tick labels
    xLabPos = [];
    if NameValueArgs.labelGroups == true
        for group = 1:(nBoxes / nGroups):nBoxes
            xLabPos = horzcat(xLabPos, median(xCoordinates(group : group + (nBoxes / nGroups) - 1)));
        end
        set(gca, 'xtick',      xLabPos, ...
                 'xticklabel', NameValueArgs.boxLabels);
    else 
        set(gca, 'xtick',      xCoordinates, ...
                 'xticklabel', NameValueArgs.boxLabels);
    end
end

%% Set figure and axes appearance
set(gcf, 'color', 'w');
set(gca, 'box',        'off', ...
         'XColor',     'k', ...
         'YColor',     'k', ...
         'TickDir',    'out', ...
         'TickLength', [0.01, 0.01], ...
         'FontSize',   13, ...
         'LineWidth',  1);

%% Set data axis limits
if ~all(isnan(maxMat), 'all')
    yUpper = max(maxMat, [], 'all');
    yLower = min(minMat, [], 'all');
    yExt = (yUpper - yLower) * 0.2;
    if isHorizontal
        xlim([yLower - yExt, yUpper + yExt]);
    else
        ylim([yLower - yExt, yUpper + yExt]);
    end
end

%% Plot legend if specified
if NameValueArgs.plotLegend == true && ...
   ~isempty(NameValueArgs.lgdLabels) && ...
   ~isempty(NameValueArgs.lgdColors)
    
    if NameValueArgs.lgdFontSize > 2
        NameValueArgs.lgdLineHeight = (NameValueArgs.lgdFontSize - 2) * NameValueArgs.lgdLineHeight;
    end
    
    for lgdEntry = 1:numel(NameValueArgs.lgdColors)
        currData = inputData(NameValueArgs.inputLabels == lgdEntry, 1);
        med = median(currData, 1, 'omitnan');
        
        if isHorizontal
            % Plot line with box color for legend
            line([med, med], ...
                 [xCoordinates(lgdEntry) - 0.25, xCoordinates(lgdEntry) + 0.25], ...
                 'Color',     NameValueArgs.lgdColors{lgdEntry}, ...
                 'LineWidth', NameValueArgs.lgdLineHeight, ...
                 'Tag',       'Legend Line');
            
            % Plot white line over line with box color for legend
            line([med, med], ...
                 [xCoordinates(lgdEntry) - 0.25, xCoordinates(lgdEntry) + 0.25], ...
                 'Color',     [1.0, 1.0, 1.0], ...
                 'LineWidth', NameValueArgs.lgdLineHeight, ...
                 'Tag',       'Legend Line Cover');
        else
            % Plot line with box color for legend
            line([xCoordinates(lgdEntry) - 0.25, xCoordinates(lgdEntry) + 0.25], ...
                 [med, med], ...
                 'Color',     NameValueArgs.lgdColors{lgdEntry}, ...
                 'LineWidth', NameValueArgs.lgdLineHeight, ...
                 'Tag',       'Legend Line');
            
            % Plot white line over line with box color for legend
            line([xCoordinates(lgdEntry) - 0.25, xCoordinates(lgdEntry) + 0.25], ...
                 [med, med], ...
                 'Color',     [1.0, 1.0, 1.0], ...
                 'LineWidth', NameValueArgs.lgdLineHeight, ...
                 'Tag',       'Legend Line Cover');
        end
        
        set(gca, 'Children', circshift(gca().Children, -2, 1));
    end
    
    warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
    
    lgdObject = legend(findobj(gca, 'Tag', 'Legend Line'), ...
                       strcat('\fontsize{', num2str(NameValueArgs.lgdFontSize), '}', NameValueArgs.lgdLabels), ...
                       'AutoUpdate',  'off', ...
                       'NumColumns',  NameValueArgs.lgdColumns, ...
                       'Orientation', NameValueArgs.lgdOrientation, ...
                       'Box',         NameValueArgs.lgdBox, ...
                       'Location',    NameValueArgs.lgdLocation);
    
    lgdObject.ItemTokenSize = [30 * NameValueArgs.lgdLineWidth, 9];
    
    if any(NameValueArgs.lgdPosition)
        lgdObject.Position = NameValueArgs.lgdPosition;
    end
else
    lgdObject = [];
end

end