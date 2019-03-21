
m = load('../data/globale_probabilities_screen_scale.mat');

sample_attr = 'globale_probabilities_start_window_1';

quantile_low = .1;
quantile_high = .9;

ref = 'FCF_attP2_1500062_UAS_TNT_2_0003';

behavior = {'crawl',    [0,0,0]; ...
    'bend', [1,0,0]; ...
    'stop', [0,1,0]; ...
    'hunch',    [0,0,1]; ...
    'back-up',  [0,1,1]};

c = @(behavior_pair, alpha) blend(behavior_pair, behavior, alpha);
a1 = .5; % alpha for the bottom-left part
a2 = .75; % alpha for the top-right part
duel = {'hunch bend',   c('hunch bend',a1),     c('hunch bend',a2); ...
    'bend back-up',     c('bend back-up',a1),   c('bend back-up',a2); ...
    'back-up hunch',    c('back-up hunch',a1),  c('back-up hunch',a2)};

marker_size = 14;
font_size = 18;

%
function s_ = capitalize(s)
    s_ = [upper(s(1)), s(2:end)];
end%function

ref = find(strcmp(ref, m.neuron_names));

table = m.(sample_attr);
for d = 1:size(duel,1)
    space = find(duel{d,1} == ' ');
    x_behavior = duel{d,1}(1:space-1);
    y_behavior = duel{d,1}(space+1:end);
    x_col = find(strcmp(x_behavior, behavior(:,1)));
    y_col = find(strcmp(y_behavior, behavior(:,1)));
    x = 100 * table(:, x_col);
    y = 100 * table(:, y_col);

    x_ = sort(x);
    x_q_low = numel(x) * quantile_low;
    x_threshold_low = (x_(floor(x_q_low)) + x_(ceil(x_q_low))) / 2;
    x_q_high = numel(x) * quantile_high;
    x_threshold_high = (x_(floor(x_q_high)) + x_(ceil(x_q_high))) / 2;
    y_ = sort(y);
    y_q_low = numel(y) * quantile_low;
    y_threshold_low = (y_(floor(y_q_low)) + y_(ceil(y_q_low))) / 2;
    y_q_high = numel(y) * quantile_high;
    y_threshold_high = (y_(floor(y_q_high)) + y_(ceil(y_q_high))) / 2;

    figure();
    h = scatter(x, y, marker_size, 'k', 'LineWidth', 1.5);
    set(gca, 'FontSize', font_size)
    hold on
    xl = xlim();
    yl = ylim();
    patch_prop = {'EdgeColor', 'none'};
    x_ = [x_threshold_high, xl(2)];
    y_ = [yl(1), y_threshold_low];
    patch(x_([1,2,2,1]), y_([1,1,2,2]), 'FaceColor', behavior{x_col,2}, patch_prop{:});
    x_ = [xl(1), x_threshold_low];
    y_ = [y_threshold_high, yl(2)];
    patch(x_([1,2,2,1]), y_([1,1,2,2]), 'FaceColor', behavior{y_col,2}, patch_prop{:});
    x_ = [xl(1), x_threshold_low];
    y_ = [yl(1), y_threshold_high];
    patch(x_([1,2,2,1]), y_([1,1,2,2]), 'FaceColor', duel{d,2}, patch_prop{:});
    x_ = [x_threshold_high, xl(2)];
    y_ = [y_threshold_low, yl(2)];
    patch(x_([1,2,2,1]), y_([1,1,2,2]), 'FaceColor', duel{d,3}, patch_prop{:});

    if ~isempty(ref)
        scatter(x(ref), y(ref), marker_size * 4, 'r', 'LineWidth', 3);
    end%if

    xlabel([capitalize(x_behavior) ' probability (%)'])
    ylabel([capitalize(y_behavior) ' probability (%)'])

    pause(1)
    print('-depsc', ['duel_' x_behavior '_' y_behavior '.eps'])

end%for
