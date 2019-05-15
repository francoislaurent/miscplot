
behaviors = {'hunch_large_squares', 'blue'; ...
    'cast_large_squares',   'red'; ...
    'back_large_squares',   'cyan'; ...
    'stop_large_squares',   'green'; ...
    'run_large_squares',    'black'; ...
    'roll_large_squares',   'yellow'; ...
    'small_motion_squares', .5 * [1,1,1]};


ethogram = load('../data/ethogramme.mat');
args = {ethogram, ...
    [105, 115, 125], ...
    'Behavior', behaviors, ...
    'MaxTimeFromOnset', [7,8], ...
    'OnsetTolerance', [.3,.6], ...
    'FillGaps', 1, ...
    'LineWidth', 2};
fs = 18;


if 1
figure();
set(gca, 'FontSize', fs)

ethogram2(args{:}, 'Include', {1, 'hunch_large_squares'; 2, 'small_motion_squares'});

xlabel('Time from first stimulus onset (s)')
ylabel('Number of animals')

print('-depsc', 'train_hunch_small_any.eps')
close
end%if


if 1
figure();
set(gca, 'FontSize', fs)

ethogram2(args{:}, 'Include', {1, 'hunch_large_squares'; 2, 'cast_large_squares'});

xlabel('Time from first stimulus onset (s)')
ylabel('Number of animals')

print('-depsc', 'train_hunch_bend_any.eps')
close
end%if


if 1
figure();
set(gca, 'fontsize', fs)

ethogram2(args{:}, 'include', {1, 'cast_large_squares'; 2, 'hunch_large_squares'});

xlabel('Time from first stimulus onset (s)')
ylabel('Number of animals')

print('-depsc', 'train_bend_hunch_any.eps')
close
end%if


if 1
figure();
set(gca, 'FontSize', fs)

ethogram2(args{:}, 'Include', {1, 'small_motion_squares'; 2, 'hunch_large_squares'});

xlabel('Time from first stimulus onset (s)')
ylabel('Number of animals')

print('-depsc', 'train_small_hunch_any.eps')
close
end%if


if 1
figure();
set(gca, 'FontSize', fs)

ethogram2(args{:}, 'Include', {1, {'back_large_squares', 'roll_large_squares'}});

xlabel('Time from first stimulus onset (s)')
ylabel('Number of animals')

print('-depsc', 'train_back_any_any.eps')
close
end%if


if 1
figure();
set(gca, 'FontSize', fs)

ethogram2(args{:}, 'Include', 'small_motion_squares');

xlabel('Time from first stimulus onset (s)')
ylabel('Number of animals')

print('-depsc', 'train_small_small_small.eps')
close
end%if


if 1
figure();
set(gca, 'FontSize', fs)

ethogram2(args{:}, 'Include', 'cast_large_squares');

xlabel('Time from first stimulus onset (s)')
ylabel('Number of animals')

print('-depsc', 'train_bend_bend_bend.eps')
close
end%if


if 0
order1 = ethogram2(ethogram, ...
    [105, 115, 125], ...
    'Behavior', behaviors, ...
    'MaxTimeFromOnset', [7,8], ...
    'OnsetTolerance', [.3,.6], ...
    'FillGaps', 1, ...
    %'Include', 'hunch_large_squares', ...
    'Include', {2, 'hunch_large_squares'}, 'Exclude', {1, 'hunch_large_squares'}, ...
    %'Exclude', {1, 'small_motion_squares'; 2, 'small_motion_squares'; 3, 'small_motion_squares'}, ...
    'LineWidth', 2);

order2 = ethogram2(ethogram, ...
    [105, 115, 125], ...
    'Behavior', behaviors, ...
    'MaxTimeFromOnset', [7,8], ...
    'OnsetTolerance', [.3,.6], ...
    'FillGaps', 1, ...
    'Include', {1, 'hunch_large_squares'}, 'Exclude', {2, 'hunch_large_squares'}, ...
    'ymin', numel(order1), ...
    'LineWidth', 2);

xlabel('Time from first stimulus onset (s)')
ylabel('Number of animals')

print('-depsc', 'stim_train2.eps')
close
end%if
