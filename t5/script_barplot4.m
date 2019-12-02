
Fig = 5

% first column is [n] in field name behavior_..._window_[n], second column is window size in seconds:
window_sizes = {1 1};%; 2 5; 3 15};

% first column is index in array behavior_proba_..., second column is color:
behaviors = { ...
%    1,  [0,0,0], [.41,.41,.41], [.65,.65,.65]; ...
    2,  [1,0,0], [1,.41,.41],  [1,.65,.65]; ...
%    3,  [0,1,0], [.57,1,.57],  [.76,1,.76]; ...
    4,  [0,0,1], [.54,.54,1],  [.72,.72,1]; ...
    5,  [0,1,1], [.62,1,1],   [.8,1,1]};

adjust = cell(0,3);
wd = pwd();
cd('../data')
%dirs = ls('-d T5_plots/bar_plots/*/')
%adjust = {1, 'LegendLocation', 'topleft'};

if Fig == 4
dirs = ['T5_plots/bar_plots/13A08/'; ...
        'T5_plots/bar_plots/32E04/'; ...
        'T5_plots/bar_plots/11A07/'; ...
        'T5_plots/bar_plots/39H05/'; ...
        'T5_plots/bar_plots/88C11/'; ...
        '66A12/                   '];
adjust = {3, 'LegendLayout', 'vertical';
          3, 'LegendLocation', 'topleft';
          3, 'LegendHorizontalMargin', .05};

elseif Fig == 5
dirs = ['T5_plots/bar_plots/11F06/'; ...
        'T5_plots/bar_plots/23A05/'; ...
        'T5_plots/bar_plots/23D10/'; ...
        'T5_plots/bar_plots/27B12/'; ...
        'T5_plots/bar_plots/29d11/'; ...
        'T5_plots/bar_plots/30E06/'; ...
        'T5_plots/bar_plots/61A01/'; ...
        'T5_plots/bar_plots/65B04/'; ...
        'T5_plots/bar_plots/77H11/'; ...
        'T5_plots/bar_plots/82B06/'; ...
        'T5_plots/bar_plots/91H01/'];
adjust = {1, 'LegendLayout', 'vertical';
          1, 'LegendHorizontalMargin', .9;
          4, 'LegendLayout', 'vertical';
          4, 'LegendHorizontalMargin', .9};
end%if
cd(wd)

% parse the field names
fields = {};
for f = 1:size(window_sizes,1)
    fields{end+1} = ['behaviour_proba_after_start_window_', num2str(window_sizes{f,1})];
end%for

ms = {load('../data/probabilities.mat')};

for d = 1:size(dirs,1)
    dir_ = strtrim(dirs(d,:));
    ms{2} = load(['../data/',dir_,'probabilities.mat']);

% arrange the data and colors
table = zeros(numel(fields), size(behaviors,1));
colors = cell(size(table));
%for f = 1:numel(fields)
for t = 1:numel(ms)
    m = ms{t};
    _table = m.probabilities.(fields{f});
    for b = 1:size(behaviors,1)
        table(t,b) = _table(behaviors{b,1});
        colors{t,b} = behaviors{b,2*t};
    end%for
end%for

larva = dir_(1:end-1);
larva = strsplit(larva, '/'){end};

optional_args = {};
lw = .9;
ltw = 5.5;
lfs = 14;
for a = 1:size(adjust,1)
    if d == adjust{a,1}
        switch adjust{a,2}
        case 'LegendWidth'
            lw = adjust{a,3};
        case 'LegendTextWidth'
            ltw = adjust{a,3};
        case 'LegendFontSize'
            lfs = adjust{a,3};
        otherwise
            optional_args = [optional_args adjust(a,[2,3])];
        end%switch
    end%if
end%for
% plot
figure();
set(gca, 'FontSize', 18)
ylim([0,50]); % BEFORE bar2
bar2(100 * table, ...
    'FaceColor', colors, 'EdgeColor', repmat(colors(1,:),size(colors,1),1), ...
    'LineWidth', 1, ...
    'Legend', {'control', 'silencing'}, ...
    'LegendWidth', lw, 'LegendTextWidth', ltw, 'LegendFontSize', lfs, optional_args{:});

set(gca, 'XTickLabel', {'bend', 'hunch', 'back-up'})
ylabel('Behavioral probability (%)')

pause(1)
%print('-depsc', [larva,'_3probabilities.eps'])
file = [larva,'_3probabilities.png'];
print('-dpng', '-r300', '-F:30', file)
close
system(['convert ', file, ' -transparent white ', file]);

end%for
