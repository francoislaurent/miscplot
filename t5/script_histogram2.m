
% load file (adapt the filepath)
m1 = load('../data/cumulative.mat');
m2 = load('../data/cumulative_p3.mat');
m = {m1, m2};
lgd = {'high', 'low'};

% first column is [n] in field name behavior_..._window_[n], second column is window size in seconds:
window_sizes = {1 1; 2 5; 3 15};

% first column is index in array behavior_proba_..., second column is color:
behaviors = { ...
    1,  [0,0,0], [.65,.65,.65]; ...
    2,  [1,0,0], [1,.65,.65]; ...
    3,  [0,1,0], [.76,1,.76]; ...
    4,  [0,0,1], [.72,.72,1]; ...
    5,  [0,1,1], [.8,1,1]};

% parse the field names
fields = {};
for f = 1:size(window_sizes,1)
    fields{end+1} = ['behaviour_proba_after_start_window_', num2str(window_sizes{f,1})];
end%for

for f = 1:numel(fields)

% arrange the data and colors
table = zeros(2, size(behaviors,1));
colors = cell(size(table));
for t = 1:numel(m)
    table_ = m{t}.cumulative.(fields{f});
    for b = 1:size(behaviors,1)
        table(t,b) = table_(behaviors{b,1});
        colors{t,b} = behaviors{b,t+1};
    end%for
end%for

% plot
figure();
set(gca, 'FontSize', 18)
bar2(100 * table, ...
    'Interval', 1, ...
    'FaceColor', colors, 'EdgeColor', repmat(colors(1,:),size(colors,1),1), ...
    'LineWidth', 1, ...
    'Legend', lgd, ...
    'LegendWidth', .9, 'LegendTextWidth', 3, 'LegendFontSize', 14);

set(gca, 'XTickLabel', {'crawl', 'bend', 'stop', 'hunch', 'back-up'})
ylabel('Behavioral probability (%)')

pause(1)
print('-depsc', ['histogram2_window_',num2str(window_sizes{f,2}),'.eps'])
close

end%for
