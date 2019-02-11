
function m = ethogram1(m, varargin)
%ETHOGRAM1   Plot ethogram.
%   DATA_STRUCT=ETHOGRAM1(DATA_FILE_OR_STRUCT, KEY1, VALUE1, ...)
%   DATA_STRUCT is a structure with an 'ethogramme' field and optionally a
%   'ethogramme_data' field.
%   The 'ethogramme' field is an array of structures with behavior names as fields.
%   Each behavior field should be a structure with an 'x' attribute that contains
%   a 4-element array.
%   The 'ethogramme_data' is a structure with a 'protocol' field.
%   Optional key-value pairs can be (keys are case-insensitive):
%   'StimulusOnset': stimulus onset time(s); vertical bars are plotted
%   'MaxTimeFromOnset': segment half-width from first stimulus onset,
%       or [before-onset, after-onset] time pair
%   'OnsetTolerance': after-stimulus time tolerance for response detection,
%       or [before-onset, after-onset] time pair
%   'Behavior': behavior description (cell array); the first column contains
%       behavior names as available in DATA_STRUCT.ethogramme; the second column
%       contains colors; row order matters
%   'Level': number of successive behaviors after stimulus onset for grouping traces
%       by behavior ordering
%   'NoResponse': if is 'hide', then trials which response to stimulus happens more
%       then MaxTimeFromOnset after onset are not plotted
%   'Exclude': cell array of response behaviors not to be plotted;
%       the first column is response level (from 1 to 'Level'), the second column is
%       a behavior or cell array of behaviors; undefined behavior is 'none'
%   'FillGaps': fill gaps between successive states
%   'LineColor': color of the stimulus onset vertical line(s)
%   'LineWidth': width of the stimulus onset vertical line(s)

if ischar(m)
    m = load(m);
end%if

stimulus_onset = []; % seconds
max_time_from_onset = []; % seconds
onset_tol = [0, 0]; % seconds
dt = .05; % seconds
end_tol = 0; % seconds
fill_gaps = 0; % 0 or 1
n_levels = 1;
no_response = 'show';
exclude = {};
o_color = 'm';
o_linewidth = 2;
behavior = {};
%    'hunch_large',  'blue';
%    'cast_large',   'red';
%    'back_large',   'cyan';
%    'stop_large',   'green';
%    'run_large',    'black';
%    'roll_large',   'yellow';
%    'small_motion', [26, 26, 26] / 256};
suffix = '';%'squares';

for v = 1:2:nargin-1
    switch lower(varargin{v})
        case 'stimulusonset'
            stimulus_onset = varargin{v+1};
        case 'maxtimefromonset'
            max_time_from_onset = varargin{v+1};
        case 'onsettolerance'
            onset_tol = varargin{v+1};
        case 'end'
            end_tol = varargin{v+1};
        case 'behavior'
            behavior = varargin{v+1};
        case 'fillgaps'
            fill_gaps = varargin{v+1};
        case 'level'
            n_levels = varargin{v+1};
        case 'noresponse'
            no_response = varargin{v+1};
        case 'exclude'
            exclude = varargin{v+1};
        case 'linecolor'
            o_color = varargin{v+1};
        case 'linewidth'
            o_linewidth = varargin{v+1};
    end%switch
end%for

% parse user arguments
if strcmpi(stimulus_onset, 'none')
    stimulus_onset = [];
elseif isempty(stimulus_onset)
    p = m.ethogramme_data.protocol;
    stims = cut(p, '#');
    stimulus_onset = [];
    for s_ = 1:numel(stims)
        if strcmp(stims{s_}, 'n')
            break
        end%if
        specs = cut(stims{s_}, '_');
        specs = specs{end};
        onset = cut(specs, 's');
        onset = onset{1};
        stimulus_onset(end+1) = str2num(onset);
    end%for
end%if
if isempty(behavior)
    trial = m.ethogramme(1);
    behaviors = fieldnames(trial);
    matched = false(size(behaviors));
    larva_behavior = { ...
        'hunch',  [0,0,1]; ...
        'cast',   [1,0,0]; ...
        'back',   [0,1,1]; ...
        'stop',   [0,1,0]; ...
        'run',    [0,0,0]; ...
        'roll',   [1,1,0]};
    behavior = cell(0, 2);
    for b = 1:size(larva_behavior,1)
        matches = find(strncmp(larva_behavior{b,1}, behaviors, numel(larva_behavior{b,1})));
        for m_ = 1:numel(matches)
            behavior(end+1,:) = ...
                {behaviors{matches(m_)}, larva_behavior{b,2}/m_};
        end%for
        matched(matches) = 1;
    end%for
    unmatched = find(~matched);
    for m_ = 1:numel(unmatched)
        extra_color = rand(1, 3);
        extra_color = extra_color .* (.9 / sqrt(sum(extra_color .* extra_color)));
        behavior(end+1,:) = ...
            {behaviors{unmatched(m_)}, extra_color};
    end%for
end%if
if isscalar(max_time_from_onset)
    max_time_from_onset = [max_time_from_onset, max_time_from_onset];
elseif ~isempty(max_time_from_onset)
    max_time_from_onset(1) = abs(max_time_from_onset(1));
end%if
if isscalar(onset_tol)
    onset_tol = [0, onset_tol];
end%if
no_response = ~any(strcmpi(no_response, {'hide', 'none'}));
if isempty(exclude)
    exclude = cell(0,2); % just make sure first dimension is null
else
    if ~iscell(exclude)
        exclude = {1, {exclude}};
    elseif size(exclude,2) == 2 && isnumeric(exclude{1})
        for e = 1:size(exclude,1)
            if ~iscell(exclude{e,2})
                exclude{e,2} = exclude(e,2);
            end%if
        end%for
    else
        exclude = {1, exclude};
    end%if
    for e1 = 1:size(exclude,1)
        exclude_ = exclude{e1,2};
        for e2 = 1:numel(exclude_)
            if strcmpi(exclude_{e2}, 'none')
                exclude_{e2} = -1;
            else
                e_ = find(strcmp(exclude_{e2}, behavior(:,1)));
                if isempty(e_)
                    error(['no such behavior: ', exclude_{e2}])
                end%if
                exclude_{e2} = e_;
            end%if
        end%for
        exclude{e1,2} = [exclude_{:}];
    end%for
end%if

nbehaviors = size(behavior,1);
no_state = size(behavior,1) + 1;
trials = {};
for e = 1:numel(m.ethogramme)
    trial = m.ethogramme(e);
    states = fieldnames(trial);
    for b = 1:nbehaviors
        empty = isempty(trial.(behavior{b,1}));
        if ~empty
            break
        end%if
    end%while
    if empty
        continue
    end%if
    N = zeros(1, nbehaviors);
    Npost = zeros(1, nbehaviors);
    X = {};
    for b = 1:nbehaviors
        ts = trial.(behavior{b,1});
        if isempty(stimulus_onset)
            N(b) = numel(ts);
            for t = 1:numel(ts)
                x = ts(t).x([2,3]);
                X{end+1} = x;
            end%for
            Npost = 1; % so that Npost is ignored
        elseif isempty(max_time_from_onset)
            N(b) = numel(ts);
            for t = 1:numel(ts)
                x = ts(t).x([2,3]);
                x = x - stimulus_onset(1);
                if -onset_tol(1)<=x(1)
                    Npost(b) = Npost(b) + 1;
                end%if
                X{end+1} = x;
            end%for
        else
            for t = 1:numel(ts)
                x = ts(t).x([2,3]);
                x = x - stimulus_onset(1);
                if -max_time_from_onset(1)<x(2) && x(1)<max_time_from_onset(2)
                    N(b) = N(b) + 1;
                    X{end+1} = x;
                    if -onset_tol(1)<=x(1)
                        Npost(b) = Npost(b) + 1;
                    end%if
                end%if
            end%for
        end%if
    end%for
    X = [X{:}];
    if isempty(X) || sum(Npost) < 1
        continue
    end%if
    if ~(isempty(stimulus_onset) || isempty(max_time_from_onset) || end_tol<0)
        if ~(min(X(1,:))<=-max_time_from_onset(1)+end_tol && max_time_from_onset(2)-end_tol<=max(X(2,:)))
            continue
        end%if
    end%if
    trials{1,end+1} = N;
    trials{2,end} = X;
end%for

ntrials = size(trials,2);
if n_levels == 0 || isempty(stimulus_onset)
    order = 1:ntrials;
else
    group = cell([2, nbehaviors, (nbehaviors+1)*ones(1,n_levels-1)]);
    stops = cell(size(group));
    for t = 1:ntrials
        N = cumsum(trials{1,t});
        s_onsets = trials{2,t}(1,:);
        [s_onsets, order] = sort(s_onsets);
        s_stops = trials{2,t}(2,order);
        % fill the gaps if any
        if fill_gaps
            s_isis = s_onsets(2:end) - s_stops(1:end-1); % inter-state intervals
            gap = find(dt < s_isis);
            for g = 1:numel(gap)
                s_half_isi = s_isis(gap(g)) * .5;
                s_onset = s_onsets(gap(g)+1) - s_half_isi;
                s_stop = s_stops(gap(g)) + s_half_isi;
                %s_onsets(gap(g)+1) = s_onset; % should gap filling affect the ordering?
                trials{2,t}(1,order(gap(g)+1)) = s_onset;
                trials{2,t}(2,order(gap(g))) = s_stop;
            end%for
        end%if
        % first (new) state after onset and its end time
        candidate_first_at_onset = find(-onset_tol(1) <= s_onsets);
        [i_, best] = min(abs(s_onsets(candidate_first_at_onset)));
        first_post_onset = candidate_first_at_onset(best);
        if ~isempty(first_post_onset)
            s_onset = s_onsets(first_post_onset);
            if s_onset <= onset_tol(2)
                state = first_post_onset;
                col = 1;
            else
                % or baseline state at onset (and its end time)
                baseline = first_post_onset - 1; % rather baseline than first
                state = baseline;
                col = 2;
            end%if
        else
            % baseline state at onset (and its end time)
            last = numel(order);
            state = last;
            col = 2;
        end%if
        if state < 1
            continue
        end%if
        s_stop = s_stops(state);
        states = order(state:end);
        if n_levels < numel(states)
            states = states(1:n_levels);
        end%if
        for s = 1:numel(states)
            states(s) = find(states(s) <= N, 1);
        end%for
        continue_ = false;
        for e1 = 1:size(exclude,1)
            level = exclude{e1,1};
            exclude_ = exclude{e1,2};
            if level <= numel(states)
                if any(exclude_==states(level))
                    continue_ = true; break
                end%if
            elseif any(exclude_==-1)
                continue_ = true; break
            end%if
        end%for
        if continue_
            continue
        end%if
        if numel(states) < n_levels
            states = [states, no_state * ones(1, n_levels - numel(states))];
        end%if
        ix = mat2cell([col, states], 1, ones(1, numel(states)+1));
        group{ix{:}}(end+1)= t;
        stops{ix{:}}(end+1) = s_stop;
    end%for
    order = {};
    if n_levels == 1
        for g = 1:size(group,2)
            [stops_, order_] = sort(stops{1,g});
            order{end+1} = group{1,g}(order_(end:-1:1));
        end%for
        if no_response
            for g = size(group,2):-1:1
                [stops_, order_] = sort(stops{2,g});
                order{end+1} = group{2,g}(order_);
            end%for
        end%if
    else
        for g = 1:size(group,2)
            subgroup_stops = subgroup(stops, 1, g);
            subgroup_group = subgroup(group, 1, g);
            for h = 1:numel(subgroup_group)
                if ~isempty(subgroup_group{h})
                    [stops_, order_] = sort(subgroup_stops{h});
                    order{end+1} = subgroup_group{h}(order_(end:-1:1));
                end%if
            end%for
        end%for
        if no_response
            for g = size(group,2):-1:1
                subgroup_stops = subgroup(stops, 2, g);
                subgroup_group = subgroup(group, 2, g);
                for h = numel(subgroup_group):-1:1
                    if ~isempty(subgroup_group{h})
                        [stops_, order_] = sort(subgroup_stops{h});
                        order{end+1} = subgroup_group{h}(order_);
                    end%if
                end%for
            end%for
        end%if
    end%if
    order = [order{:}];
end%if

ntrials = numel(order);
y0 = ntrials; dy = 1;
y = y0;
f = zeros(1, nbehaviors);
faces = cell(1, nbehaviors);
vertices = cell(1, nbehaviors);
[vertices{:}] = deal({});
[faces{:}] = deal({});
for t = 1:ntrials
    n = trials{1,order(t)};
    x = trials{2,order(t)};
    k = 1;
    for b = 1:nbehaviors
        for s = 1:n(b)
            vertices{b}{end+1,1} = x([1,1,2,2], k);
            vertices{b}{end,2} = [y-dy; y; y; y-dy];
            faces{b}{end+1,1} = f(b)+1:f(b)+4;
            k = k + 1;
            f(b) = f(b) + 4;
        end%for
    end%for
    y = y - dy;
end%for

for b = 1:nbehaviors
    vertices{b} = cell2mat(vertices{b});
    faces{b} = cell2mat(faces{b});
end%for

gca();
for b = 1:nbehaviors
    patch('Faces', faces{b}, 'Vertices', vertices{b}, 'FaceColor', behavior{b,2}, 'EdgeColor', 'none')
    hold on
end%for

for o = 1:numel(stimulus_onset)
    x0 = stimulus_onset(o) - stimulus_onset(1);
    plot([x0, x0], [0, y0], '-', 'Color', o_color, 'LineWidth', o_linewidth);
end%for

if ~isempty(max_time_from_onset)
    xlim([-max_time_from_onset(1), max_time_from_onset(2)])
end%if
ylim([0, y0])

% nested functions
function ss=cut(s, sep)
    seps_ = find(s == sep);
    if isempty(seps_)
        ss = {s};
    else
        ss = {s(1:seps_(1)-1)};
        for s_ = 1:numel(seps_)-1
            ss{end+1} = s(seps_(s_)+1:seps_(s_+1)-1);
        end%for
        ss{end+1} = s(seps_(end)+1:end);
    end%if
end%function
function m_=subgroup(m_, c, g)
    s_ = size(m_);
    m_ = m_(c,g,:);
    if 3 < numel(s_)
        s_ = s_(3:end);
        m_ = permute(reshape(m_, s_), numel(s_):-1:1);
    end%if
end%function

end%function
