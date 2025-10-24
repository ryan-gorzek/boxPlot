# boxPlot - A Better Boxplot for MATLAB

A highly customizable boxplot function for MATLAB that provides extensive control over appearance, supports grouped data, individual data point visualization, and paired data connections.

![plot](https://github.com/ryan-gorzek/boxPlot/blob/main/plots/boxPlot_Examples.jpeg)
## Why boxPlot?

MATLAB's built-in `boxplot` and `boxchart` functions have limitations in customization and can be cumbersome to style. **boxPlot** offers:

- **Full color control** - Customize box fills, edges, medians, whiskers, and outliers independently
- **Transparency support** - Set alpha values for all visual elements
- **Grouped data** - Easy handling of factorial designs with automatic spacing
- **Individual points** - Overlay raw data with customizable jitter
- **Paired connections** - Connect related observations with lines (ideal for repeated measures)
- **Publication-ready** - Clean, professional look with minimal code
- **Intuitive syntax** - Modern name-value pair arguments with sensible defaults
- **Zero dependencies** - Works with base MATLAB, no toolboxes required

## Installation

1. Download `boxPlot.m`
2. Add it to your MATLAB path or working directory
3. Start plotting!

```matlab
addpath('/path/to/boxPlot');
```

## Quick Start

### Basic Usage

```matlab
% Generate sample data
data = randn(100, 3) + [5, 7, 6];

% Create boxplot with custom colors
boxPlot(data, ...
    'boxLabels', {'Group A', 'Group B', 'Group C'}, ...
    'boxColors', {[0.8, 0.3, 0.3], [0.3, 0.6, 0.8], [0.4, 0.8, 0.4]});
```

### Show Individual Data Points

```matlab
boxPlot(data, ...
    'boxColors', {[0.7, 0.7, 0.7], [0.7, 0.7, 0.7], [0.7, 0.7, 0.7]}, ...
    'plotPoints', true, ...
    'pointJitter', 'rand', ...
    'pointAlpha', 0.4);
```

### Paired Data with Connecting Lines

```matlab
% Paired measurements (e.g., before/after)
baseline = randn(30, 1) * 2 + 20;
followup = baseline + randn(30, 1) * 1.5 + 3;
pairedData = [baseline, followup];

boxPlot(pairedData, ...
    'boxLabels', {'Baseline', 'Follow-up'}, ...
    'plotLines', true, ...
    'plotPoints', true, ...
    'pointJitter', 'rand');
```

### Grouped Data with Legend

```matlab
% 2x2 design: 2 timepoints, 2 conditions
groupedData = randn(50, 4) + [10, 12, 11, 13];

boxPlot(groupedData, ...
    'groupSize', 2, ...
    'labelGroups', true, ...
    'boxLabels', {'Time 1', 'Time 2'}, ...
    'boxColors', {[0.4, 0.6, 0.9], [0.9, 0.5, 0.4], ...
                  [0.4, 0.6, 0.9], [0.9, 0.5, 0.4]}, ...
    'plotLegend', true, ...
    'lgdLabels', {'Condition A', 'Condition B'}, ...
    'lgdColors', {[0.4, 0.6, 0.9], [0.9, 0.5, 0.4]});
```

### Custom Box Coordinates

```matlab
% Non-uniform spacing (e.g., for time course data)
data = randn(50, 4) + [10, 15, 20, 30];
customCoords = [1, 2, 4, 7];  % Specify exact x-positions

boxPlot(data, ...
    'boxCoordinates', customCoords, ...
    'boxLabels', {'0 hrs', '24 hrs', '72 hrs', '168 hrs'});
```

## Key Features

### Comprehensive Customization

Control every aspect of your boxplot:
- **Boxes**: Fill color, edge color, transparency, line width, corner curvature
- **Medians**: Color, line width, line style, transparency
- **Whiskers**: Color, line width, line style, transparency
- **Outliers**: Color, size, marker style, transparency, jitter
- **Points**: Individual data point overlay with full styling control
- **Lines**: Connect paired observations across groups

### Grouping Support

Handle complex experimental designs with ease:
- Specify `groupSize` to define factorial structures
- Use `labelGroups` to show one label per group
- Automatically handles spacing and alignment

### Smart Defaults

- Sensible default colors (gray boxes, black medians/whiskers)
- Automatic axis scaling with appropriate padding
- Automatic missing data (NaN) handling
- Publication-ready appearance out of the box

### Flexibility

- Matrix input: each column becomes a box
- Vector input with labels: specify categorical grouping
- Handles datasets with varying numbers of observations

## Syntax

```matlab
xCoordinates = boxPlot(inputData)
[xCoordinates, lgdObject] = boxPlot(inputData, Name, Value)
```

### Input Arguments

- `inputData` - Numeric matrix or vector containing data to plot

### Name-Value Arguments (Selected)

| Argument | Default | Description |
|----------|---------|-------------|
| `boxLabels` | `{1:N}` | Cell array of x-axis labels |
| `boxColors` | Gray | Cell array of RGB colors for box fills |
| `boxAlpha` | `1` | Box fill transparency (0-1) |
| `boxEdgeWidth` | `1` | Width of box edges |
| `medianWidth` | `2` | Width of median line |
| `outlierSize` | `30` | Size of outlier markers |
| `plotPoints` | `false` | Overlay individual data points |
| `pointJitter` | `'none'` | Jitter method: `'none'`, `'rand'`, `'randn'` |
| `plotLines` | `false` | Connect observations with lines |
| `groupSize` | `1` | Number of boxes per group |
| `labelGroups` | `false` | Use one label per group |
| `plotLegend` | `false` | Display legend |

See full documentation in the function header for all 40+ customization options.

### Output Arguments

- `xCoordinates` - X-axis positions of boxes (useful for adding annotations)
- `lgdObject` - Legend object (if legend is displayed)

## Examples

The included `boxPlot_Examples.m` script demonstrates:

1. Basic boxplot with custom colors
2. Boxplot with individual data points
3. Grouped boxplots with legend
4. Paired data with connecting lines
5. Customized appearance with thick lines and styling
6. Multiple groups with color-coded subject tracking
7. Minimal clean style for publications
8. Wide data with custom spacing
9. Transparent boxes with visible data distribution
10. Publication-ready multi-panel comparison

## Advanced Usage

### Custom Color Schemes

```matlab
% Define a custom color palette
colors = {[0.85, 0.33, 0.10], ...  % Orange
          [0.93, 0.69, 0.13], ...  % Yellow
          [0.49, 0.18, 0.56], ...  % Purple
          [0.47, 0.67, 0.19]};     % Green

boxPlot(data, 'boxColors', colors);
```

### Transparent Boxes with Data

```matlab
% Great for showing full data distribution
boxPlot(data, ...
    'boxAlpha', 0.3, ...           % Transparent boxes
    'plotPoints', true, ...
    'pointSize', 25, ...
    'pointAlpha', 0.6, ...
    'pointJitter', 'randn');       % Normal jitter
```

### Repeated Measures Design

```matlab
% Track individual subjects across conditions
nSubjects = 20;
conditions = 4;
data = randn(nSubjects, conditions);

% Create unique color for each subject
subjectColors = cell(nSubjects, conditions);
cmap = lines(nSubjects);
for i = 1:nSubjects
    for j = 1:conditions
        subjectColors{i,j} = cmap(i,:);
    end
end

boxPlot(data, ...
    'plotLines', true, ...
    'lineColors', subjectColors(:, 1:conditions-1), ...
    'plotPoints', true, ...
    'pointColors', subjectColors);
```

### Styling for Publications

```matlab
% Minimal, clean style
boxPlot(data, ...
    'boxColors', {[0.95, 0.95, 0.95]}, ...  % Light gray
    'boxEdgeWidth', 1.5, ...
    'medianWidth', 2.5, ...
    'whiskerWidth', 1.5);

% Add grid for readability
grid on;
set(gca, 'GridAlpha', 0.15);
```

## Tips and Tricks

1. **Jitter control**: Use `jitterBound` to adjust the spread of jittered points
2. **Axis limits**: Automatically calculated with 20% padding, or set manually with `ylim` and `xlim`
3. **Box spacing**: Adjust to your preferred look with `boxSpacing`

## Requirements

- MATLAB R2019b or later (for arguments block syntax)
- No additional toolboxes required

## License

MIT License - see LICENSE file for details

Copyright (c) 2025 Ryan Gorzek
