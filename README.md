# boxPlot - Enhanced Boxplot for MATLAB

A highly customizable boxplot function for MATLAB that provides extensive control over appearance, supports grouped data, individual data point visualization, and paired data connections.

## Why boxPlot?

MATLAB's built-in `boxplot` function has limitations in customization and can be cumbersome to style. **boxPlot** offers:

- **Full color control** - Customize box fills, edges, medians, whiskers, and outliers independently
- **Transparency support** - Set alpha values for all visual elements
- **Grouped data** - Easy handling of factorial designs with automatic spacing
- **Individual points** - Overlay raw data with customizable jitter
- **Paired connections** - Connect related observations with lines (ideal for repeated measures)
- **Publication-ready** - Clean, professional output with minimal code
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
- Handles missing data (NaN) gracefully
- Publication-ready appearance out of the box

### Flexible Input

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

The included `boxPlot_examples.m` script demonstrates:

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

Run the examples:
```matlab
boxPlot_examples
```

## Comparison with MATLAB's boxplot

| Feature | MATLAB boxplot | boxPlot |
|---------|----------------|---------|
| Custom box colors | Limited | Full RGB control per box |
| Transparency | Not supported | All elements support alpha |
| Individual points | Requires separate code | Built-in with jitter |
| Paired connections | Not supported | Built-in line plotting |
| Grouped data | Complex syntax | Simple `groupSize` parameter |
| Legend | Manual creation | Automatic with styling |
| Dependencies | Statistics Toolbox | None - base MATLAB only |
| Modern syntax | Older style | Name-value pairs |

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

1. **Color consistency**: Define your color scheme once and reuse it across figures
2. **Jitter control**: Use `jitterBound` to adjust the spread of jittered points
3. **Missing data**: The function handles NaN values automatically
4. **Axis limits**: Automatically calculated with 20% padding, or set manually with `ylim`
5. **Export**: Use MATLAB's `exportgraphics` for high-quality figure export

```matlab
exportgraphics(gcf, 'myBoxplot.png', 'Resolution', 300);
```

## Requirements

- MATLAB R2019b or later (for arguments block syntax)
- No additional toolboxes required

## License

MIT License - see LICENSE file for details

Copyright (c) 2022 Ryan Gorzek

## Contributing

Contributions, issues, and feature requests are welcome! Feel free to check the [issues page](../../issues).

## Citation

If you use this function in your research, please cite:

```
Gorzek, R. (2022). boxPlot: Enhanced boxplot visualization for MATLAB. 
GitHub repository, https://github.com/gorzek-ryan/matlab_viz
```

## Changelog

### Version 1.0.0
- Initial release
- Full customization of box appearance
- Support for grouped data and legends
- Individual point plotting with jitter
- Paired data connections
- Comprehensive documentation and examples

## Support

For questions, bug reports, or feature requests:
- Open an issue on GitHub
- Check existing documentation and examples
- Review the function header for detailed parameter descriptions

## Related Functions

This is part of a larger collection of MATLAB visualization tools. Check out the [matlab_viz repository](https://github.com/gorzek-ryan/matlab_viz) for more enhanced plotting functions.

---

**Made with ❤️ for better scientific visualization in MATLAB**