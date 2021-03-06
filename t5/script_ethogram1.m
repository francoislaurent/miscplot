
Fig = 5

behaviors = {'hunch_large_squares', 'blue'; ...
    'cast_large_squares',   'red'; ...
    'back_large_squares',   'cyan'; ...
    'stop_large_squares',   'green'; ...
    'run_large_squares',    'black'; ...
    'roll_large_squares',   'yellow'; ...
    'small_motion_squares', .5 * [1,1,1]};

ntrials = 1000;

figure();
set(gca, 'FontSize', 18)

%seed = .1;
%while round(seed) != seed
%    seed = rand('seed');
%end%while
%printf('%.0f\n', seed)
rand('seed', 15698812726790480012274691666152165506220855954338243033166107128567733436938819560761317475076882758950135850311990626557494799892480)

ethogram = load('../data/ethogramme.mat');
ix = randperm(numel(ethogram.ethogramme))(1:ntrials);
ethogram.ethogramme = ethogram.ethogramme(ix);

m = ethogram1(ethogram, ...
    'Behavior', behaviors, ...
    'StimulusOnset', 45, ...
    'MaxTimeFromOnset', 3, ...
    'OnsetTolerance', [.2,.5], ...
    'FillGaps', 1, ...
    'Level', 2, ...
    %'NoResponse', 'hide', ...
    %'Exclude', {1, 'small_motion_squares'; 2, 'small_motion_squares'; 3, 'small_motion_squares'}, ...
    'LineWidth', 3);

xlabel('Time from stimulus onset (s)')
ylabel('Number of animals')

%print('-depsc', 'control1.eps')
if Fig == 4
print('-dpng', '-r600', '-F:24', 'control1.png')
close
elseif Fig == 5
file = 'control-small1.png';
print('-dpng', '-r300', '-F:30', file)
close
system(['convert ', file, ' -transparent white ', file]);
end%if

