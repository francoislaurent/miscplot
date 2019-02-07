
figure();
set(gca, 'FontSize', 18)

m = ethogram1('../data/ethogramme.mat',
    'StimulusOnset', 45,
    'MaxTimeFromOnset', 10,
    'OnsetTolerance', [.2,.5],
    'FillGaps', 1,
    'Level', 3,
    'NoResponse', 'hide',
    'Exclude', {1, 'small_motion_squares'; 2, 'small_motion_squares'; 3, 'small_motion_squares'},
    'LineWidth', 3);

xlabel('Time from stimulus onset (s)')
ylabel('Number of animals')
