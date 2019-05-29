
behaviors = { ...
    'hunch',    'hunch_large_squares',  'blue'; ...
    'bend',     'cast_large_squares',   'red'; ...
    'back',     'back_large_squares',   'cyan'; ...
    'stop',     'stop_large_squares',   'green'; ...
    'run',      'run_large_squares',    'black'; ...
    'roll',     'roll_large_squares',   'yellow'; ...
    'small',    'small_motion_squares', .5 * [1,1,1]};

ntrials = 1380;
rand('seed', 15698812726790480012274691666152165506220855954338243033166107128567733436938819560761317475076882758950135850311990626557494799892480)

ethogram = load('../data/ethogramme.mat');
ix = randperm(numel(ethogram.ethogramme))(1:ntrials);
%ethogram.ethogramme = ethogram.ethogramme(ix);

figure();
set(gca, 'FontSize', 18)

[order, counts] = ethogram3(ethogram, ...
    'Behavior', behaviors, ...
    'EventMaxResponseTime', .3, ...
    'EventOnsetTolerance', .2, ...
    'FirstOnsetAsOrigin', true, ...
    'EventOnsetLineWidth', 3, ...
    %'MaxWindow', [42,48], ...
    'EventOnsetTime', 45, ...
    'Grouping', {'transition=1,state(-1)=hunch,state;transition=0,state(0)=hunch;transition=1,state(-1)=bend,state;transition=0,state(0)=bend;transition=1,state(-1)=back,state;transition=0,state(0)=back;transition=1,state(-1)=stop,state;transition=0,state(0)=stop;transition=1,state(-1)=run,state;transition=0,state(0)=run;transition=1,state(-1)=roll,state;transition=0,state(0)=roll;transition=1,state(-1)=small,state;transition=0,state(0)=small'}, ...
    'ReturnGroupCounts', true, ...
    'FillGaps', 1);

xlim([-3,3])
%ylim([0,405])
xlabel('Time from stimulus onset (s)')

ylabel('Number of animals')

%print('-depsc', ['by_previous_state.eps'])
%close
