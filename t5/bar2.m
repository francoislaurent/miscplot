function varargout=bar2(y, varargin)
%BAR2   Bar graph with styling of individual columns.
%   [...] = BAR2(Y, KEY1, VALUE1, ...)
%   [...] = BAR2(X, Y, ...)
%   [PATCHES, LINES] = BAR2(Y)
%   [PATCHES, LINES, TEXTS] = BAR2(Y, 'Legend', LEGEND)
%   [PATCHES, TEXTS] = BAR2(Y, 'EdgeColor', 'none', 'Legend', LEGEND)
%   [LINES] = BAR2(Y, 'FaceColor', 'none')
%   with optional abscissa matrix X and bar heights Y.
%   X and Y should have the same size, with as many rows as
%   bars per group and as many columns as groups of bars.
%   PATCHES, LINES and TEXTS are graphical element handles
%   and are individually returned if non-empty.
%   For example, if 'EdgeColor' is set to 'none', no lines
%   are drawn and LINES is not returned.
%   Optional key-value pairs are (keys are case-insensitive):
%   'FaceColor': cell array of colors
%   'EdgeColor': cell array of colors
%   'FaceAlpha': array of floats in the [0,1] range
%   'EdgeAlpha': array of floats in the [0,1] range; ignored but if 0;
%       if 0, no lines are drawn and LINES is not returned
%   'LineWidth': array of floats
%   'Legend': cell array of strings (as many entries as bars in a group)
%   'LegendLocation': any of 'TopLeft', 'TopMiddle', 'TopRight'
%   'LegendWidth': relative width of the squares in the legend
%   'LegendTextWidth': relative width of the space dedicated to text
%       in the legend
%   'LegendFontSize': font size of text elements in the legend

% parse args
firstarg = 1;
if ischar(varargin{firstarg})
    x = [];
else
    x = y;
    y = varargin{firstarg};
    firstarg += 1;
end%if
dn = 2;
facecolor = 'b';
edgecolor = 'k';
facealpha = 1;
edgealpha = 1;
linewidth = .5;
legend = {};
legendlocation = 'topright';
lg_relative_width = .5;
lg_relative_text_width = 3;
lg_font_size = [];
for a = firstarg:2:nargin-1
    switch lower(varargin{a})
    case 'interval'
        dn = varargin{a+1};
    case 'facecolor'
        facecolor = varargin{a+1};
    case 'edgecolor'
        edgecolor = varargin{a+1};
    case 'facealpha'
        facealpha = varargin{a+1};
    case 'edgealpha'
        edgealpha = varargin{a+1};
    case 'linewidth'
        linewidth = varargin{a+1};
    case 'legend'
        legend = varargin{a+1};
    case 'legendlocation'
        legendlocation = varargin{a+1};
    case 'legendwidth'
        lg_relative_width = varargin{a+1};
    case 'legendtextwidth'
        lg_relative_text_width = varargin{a+1};
    case 'legendfontsize'
        lg_font_size = varargin{a+1};
    end%switch
end%for
if ~iscell(facecolor)
    facecolor = {facecolor};
end%if
if ~iscell(edgecolor)
    edgecolor = {edgecolor};
end%if
text_prop = {'VerticalAlignment', 'baseline'};
if ~isempty(lg_font_size)
    text_prop = [text_prop, {'FontSize', lg_font_size}];
end%if

sz = size(y);
sz_fc = size(facecolor);
if any(sz_fc < sz)
    facecolor = repmat(facecolor, ceil(sz ./ sz_fc));
end%if
sz_ec = size(edgecolor);
if any(sz_ec < sz)
    edgecolor = repmat(edgecolor, ceil(sz ./ sz_ec));
end%if
sz_fa = size(facealpha);
if any(sz_fa < sz)
    facealpha = repmat(facealpha, ceil(sz ./ sz_fa));
end%if
sz_ea = size(edgealpha);
if any(sz_ea < sz)
    edgealpha = repmat(edgealpha, ceil(sz ./ sz_ea));
end%if
sz_lw = size(linewidth);
if any(sz_lw < sz)
    linewidth = repmat(linewidth, ceil(sz ./ sz_lw));
end%if

if isempty(x)
    m = sz(2); n = sz(1);
    dx = .5;
    width = 1 / ((1+dx) * (n+dn));
    iid = width * (1+dx);
    xgrid = 0:iid:iid*(n-1)+eps;
    xgrid -= mean(xgrid);
    x = repmat(1:m, n, 1) + repmat(xgrid(:), 1, m);
end%if

draw_area = false(n,m);
draw_line = false(n,m);
for i = 1:n
    for j = 1:m
        draw_area(i,j) = ~strcmpi(facecolor{i,j}, 'none') && 0 < facealpha(i,j);
        draw_line(i,j) = 0 < linewidth(i,j) && ~strcmpi(edgecolor{i,j}, 'none') && 0 < edgealpha(i,j);
    end%for
end%for

area_h = [];
line_h = [];
text_h = [];
hold on
for i = 1:n
    for j = 1:m
        if draw_area(i,j)
            x_ = x(i,j) + [1, 1, -1, -1] * (width / 2);
            y_ = [0, y(i,j), y(i,j), 0];
            area_h(end+1) = patch(x_, y_, 'k', 'FaceColor', facecolor{i,j}, 'EdgeColor', 'none', ...
                'FaceAlpha', facealpha(i,j));
        end%if
        if draw_line(i,j)
            line_h(end+1) = plot(x_, y_, '-', 'Color', edgecolor{i,j}, 'LineWidth', linewidth(i,j));%, 'Alpha', edgealpha(i,j));
        end%if
    end%for
end%for

set(gca, 'XTick', 1:m)

if ~isempty(legend)
    xl = xlim(); yl = ylim();
    aspect_ratio = (yl(2)-yl(1)) / (xl(2)-xl(1));
    lg_text_width = width * lg_relative_text_width;
    lg_width = width * lg_relative_width;
    lg_height = lg_width * aspect_ratio;
    x_ = [0:lg_width:lg_width*(m-1)+eps, lg_width*(m+.2)];
    w_ = lg_width * m + lg_text_width;
    x_ = repmat(x_, 1, n) + reshape(repmat(0:w_:w_*(n-1)+eps, m+1, 1), 1, (m+1)*n);
    lg_y = yl(1) + .9 * (yl(2)-yl(1));
    switch lower(legendlocation)
    case {'topleft' 'left'}
        lg_x = xl(1)+(x_+.5*lg_text_width);%min(x(:))
    case {'topmiddle' 'middle'}
        lg_x = mean(x(:))+(x_ - mean(x_));
    case {'topright' 'right'}
        lg_x = xl(2)+(x_-x_(end)-lg_text_width);%max(x(:))
    otherwise
        error(['unexpected legend location: ', legendlocation])
    end%switch
    y_ = lg_y + lg_height * [0 0 1 1];
    k = 0;
    for i = 1:n
        for j = 1:m
            k = k + 1;
            x_ = lg_x(k);
            if draw_area
                area_h(end+1) = patch(x_ + lg_width * [0 1 1 0], y_, 'k', ...
                    'FaceColor', facecolor{i,j}, 'EdgeColor', 'none', 'FaceAlpha', facealpha(i,j));
            end%if
            if j == 1
                x__ = x_ + lg_width * [1 0 0 1];
                y__ = y_;
            elseif j == m
                x__ = x_ + lg_width * [0 1 1 0];
                y__ = y_;
            else
                x__ = [x_+lg_width, x_, nan, x_, x_+lg_width];
                y__ = [lg_y, lg_y, nan, lg_y+lg_height, lg_y+lg_height];
            end%if
            if draw_line
                line_h(end+1) = plot(x__, y__, '-', 'Color', edgecolor{i,j}, 'LineWidth', linewidth(i,j));
            end%if
        end%for
        k = k + 1;
        text_h(end+1) = text(lg_x(k), lg_y, legend{i}, text_prop{:});
    end%for
end%if

nargout_ = 0;
if nargout_ < nargout && ~isempty(area_h)
    nargout_ = nargout_ + 1;
    varargout{nargout_} = area_h;
end%if
if nargout_ < nargout && ~isempty(line_h)
    nargout_ = nargout_ + 1;
    varargout{nargout_} = line_h;
end%if
if nargout_ < nargout && ~isempty(line_h)
    nargout_ = nargout_ + 1;
    varargout{nargout_} = line_h;
end%if

end%function
