
stimulus_onset = [45, 83];

behaviours = {'hunch_large', 'blue'; ...
    'cast_large',   'red'; ...
    'back_large',   'cyan'; ...
    'stop_large',   'green'; ...
    'run_large',    'black'; ...
    'roll_large',   'yellow'; ...
    'small_motion', .5 * [1,1,1]};
feature = 'mean';

lw = 3;

m = load('../data/behaviour_probability.mat');
m = m.behaviour_proba;
t = m.t - stimulus_onset(1);

figure();
set(gca, 'LineWidth', 2)
set(gca, 'FontSize', 18)
hold on

for b = 1:size(behaviours, 1)%:-1:1
    y = 100 * m.([behaviours{b,1},'_',feature]);
    plot(t, y, '-', 'Color', behaviours{b,2}, 'LineWidth', lw)
end%for
yl = ylim();
yl = [0,80];
plot([0,0], yl, 'm-', 'LineWidth', lw+1)
for s = 2:numel(stimulus_onset)
    plot((stimulus_onset(s)-stimulus_onset(1))*[1,1], yl, 'm-', 'LineWidth', lw-1)
end%for
xlim(t([1,end]))
ylim(yl)
box on
set(gca, 'XTick', [0, 20, 40])

xlabel('Time from stimulus onset (s)')
ylabel('Behavioral probability (%)')

print('-depsc', 'timeproba1.eps')
close

