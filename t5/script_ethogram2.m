
file = '../data/ethogramme.mat';

behaviors = {'hunch_large_squares', 'blue'; ...
    'cast_large_squares',   'red'; ...
    'back_large_squares',   'cyan'; ...
    'stop_large_squares',   'green'; ...
    'run_large_squares',    'black'; ...
    'roll_large_squares',   'yellow'; ...
    'small_motion_squares', .1 * [1,1,1]};

ethogram1_prop = {file, ...
    'Behavior', behaviors, ...
    'StimulusOnset', 45, ...
    'MaxTimeFromOnset', 10, ...
    'OnsetTolerance', [.2,.5], ...
    'FillGaps', 1, ...
    'Level', 3, ...
    'NoResponse', 'hide', ...
    'LineWidth', 3};

fontsize = 18;

figure();

%subplot(3, 1, 1)
subplot('Position', [.14 .68 .8 .27])
set(gca, 'FontSize', fontsize)
all_but_hunch = behaviors(2:end,1);
all_but_cast = behaviors([1,3:end],1);
ethogram1(ethogram1_prop{:}, ...
    'Exclude', {1, all_but_hunch; 2, all_but_cast; 3, 'small_motion_squares'});
set(gca, 'XTickLabel', {})
yticklabel = get(gca, 'YTickLabel');
yticklabel{1} = '';
set(gca, 'YTickLabel', yticklabel)

%subplot(3, 1, 2)
subplot('Position', [.14 .38 .8 .27])
set(gca, 'FontSize', fontsize)
cast = {'cast_large_squares', 'small_motion_squares'};
ethogram1(ethogram1_prop{:}, ...
    'Exclude', {1, all_but_hunch; 2, cast; 3, 'small_motion_squares'});
set(gca, 'XTickLabel', {})
yticklabel = get(gca, 'YTickLabel');
yticklabel{1} = '';
set(gca, 'YTickLabel', yticklabel)
ylabel('Number of animals')

%subplot(3, 1, 3)
subplot('Position', [.14 .08 .8 .27])
set(gca, 'FontSize', fontsize)
hunch = {'hunch_large_squares', 'small_motion_squares'};
ethogram1(ethogram1_prop{:}, ...
    'Exclude', {1, hunch; 2, 'small_motion_squares'; 3, 'small_motion_squares'});
yticklabel = get(gca, 'YTickLabel');
yticklabel{1} = '';
set(gca, 'YTickLabel', yticklabel)
xlabel('Time from stimulus onset (s)')

pause(1)
print('-depsc', 'ethogram2.eps', '-S480,960')
