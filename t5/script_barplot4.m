
% first column is [n] in field name behavior_..._window_[n], second column is window size in seconds:
window_sizes = {1 1};%; 2 5; 3 15};

% first column is index in array behavior_proba_..., second column is color:
behaviors = { ...
%    1,  [0,0,0], [.41,.41,.41], [.65,.65,.65]; ...
    2,  [1,0,0], [1,.41,.41],  [1,.65,.65]; ...
%    3,  [0,1,0], [.57,1,.57],  [.76,1,.76]; ...
    4,  [0,0,1], [.54,.54,1],  [.72,.72,1]; ...
    5,  [0,1,1], [.62,1,1],   [.8,1,1]};

wd = pwd();
cd('../data')
dirs = ls('-d T5_plots/bar_plots/*/')
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
if d == 1
    optional_args = {'LegendLocation', 'topleft'};
end%if
% plot
figure();
set(gca, 'FontSize', 18)
ylim([0,50]); % BEFORE bar2
bar2(100 * table, ...
    'FaceColor', colors, 'EdgeColor', repmat(colors(1,:),size(colors,1),1), ...
    'LineWidth', 1, ...
    'Legend', {'control', 'silencing'}, ...
    'LegendWidth', .9, 'LegendTextWidth', 5.5, 'LegendFontSize', 14, optional_args{:});

set(gca, 'XTickLabel', {'bend', 'hunch', 'back-up'})
ylabel('Behavioral probability (%)')

pause(1)
print('-depsc', [larva,'_3probabilities.eps'])
close

end%for
