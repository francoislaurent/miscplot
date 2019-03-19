
function m = ethogram2(m, stimulus_onset, varargin)
%ETHOGRAM2   Plot ethogram.
%   DATA_STRUCT=ETHOGRAM2(DATA_FILE_OR_STRUCT, STIMULUS_ONSET, KEY1, VALUE1, ...)
%   DATA_STRUCT is a structure with an 'ethogramme' field and optionally a
%   'ethogramme_data' field.
%   STIMULUS_ONSET is a scalar or array of stimulus onset times in seconds.
%   The 'ethogramme' field is an array of structures with behavior names as fields.
%   Each behavior field should be a structure with an 'x' attribute that contains
%   a 4-element array.
%   The 'ethogramme_data' is a structure with a 'protocol' field.
%   Optional key-value pairs can be (keys are case-insensitive):
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

max_time_from_onset = []; % seconds
onset_duration = []; % seconds
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

for v = 1:2:nargin-2
    switch lower(varargin{v})
        case 'maxtimefromonset'
            max_time_from_onset = varargin{v+1};
        case 'onsettolerance'
            onset_tol = varargin{v+1};
        case 'onsetduration'
            onset_duration = varargin{v+1};
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

if isempty(stimulus_onset)
    error('no stimulus onset defined')
end%if
stimulus_onset = sort(stimulus_onset);
relative_onset = stimulus_onset - stimulus_onset(1);
stimulus_max_interval = relative_onset(end);

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
        if isempty(max_time_from_onset)
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
                if -max_time_from_onset(1)<x(2) && x(1)<stimulus_max_interval+max_time_from_onset(2)
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
    if ~(isempty(max_time_from_onset) || end_tol<0)
        if ~(min(X(1,:))<=-max_time_from_onset(1)+end_tol && ...
            stimulus_max_interval+max_time_from_onset(2)-end_tol<=max(X(2,:)))
            continue
        end%if
    end%if
    trials{1,end+1} = N;
    trials{2,end} = X;
end%for

n_stims = numel(stimulus_onset);
ntrials = size(trials,2);
if n_levels == 0 || isempty(stimulus_onset)
    order = 1:ntrials;
else
    if no_response
        group = cell([2*ones(1,n_stims),nbehaviors*ones(1,n_stims),(nbehaviors+1)*ones(1,(n_levels-1)*n_stims)]);
    else
        group = cell([nbehaviors*ones(1,n_stims),(nbehaviors+1)*ones(1,(n_levels-1)*n_stims)]);
    end%if
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
        states = no_state * ones(1,numel(size(group)));
        continue_ = false;
        for n = 1:n_stims
            sn_onsets = s_onsets - relative_onset(n);
            sn_stops = s_stops - relative_onset(n);
            % first (new) state after onset and its end time
            candidate_first_at_onset = find(-onset_tol(1) <= sn_onsets);
            [i_, best] = min(abs(sn_onsets(candidate_first_at_onset)));
            first_post_onset = candidate_first_at_onset(best);
            if ~isempty(first_post_onset)
                sn_onset = sn_onsets(first_post_onset);
                if sn_onset <= onset_tol(2)
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
            if ~no_response && col == 2
                continue_ = true; break
            end%if
            if state < 1
                continue_ = true; break % TODO: check this
                continue
            end%if
            if n == 1
                s1_stop = sn_stops(state);
            end%if
            sn_states = order(state:end);
            if n_levels < numel(sn_states)
                sn_states = sn_states(1:n_levels);
            end%if
            for s = 1:numel(sn_states)
                sn_states(s) = find(sn_states(s) <= N, 1);
            end%for
            for e1 = 1:size(exclude,1)
                level = exclude{e1,1};
                exclude_ = exclude{e1,2};
                if level <= numel(sn_states)
                    if any(exclude_==sn_states(level))
                        continue_ = true; break
                    end%if
                elseif any(exclude_==-1)
                    continue_ = true; break
                end%if
            end%for
            if continue_
                break
            end%if
            if no_response
                states(n) = col;
                states(n_stims+n:n_stims:numel(sn_states)*n_stims+n) = sn_states;
            else
                states(n:n_stims:(numel(sn_states)-1)*n_stims+n) = sn_states;
            end%if
        end%for
        if continue_
            continue
        end%if
        ix = mat2cell(states, 1, ones(1, numel(states)));
        group{ix{:}}(end+1) = t;
        stops{ix{:}}(end+1) = s1_stop;
    end%for
    order = {};
    if true%no_response
        group = permute(group, numel(size(group)):-1:1);
        stops = permute(stops, numel(size(stops)):-1:1);
        for g = 1:numel(group)
            [stops_, order_] = sort(stops{g});
            order{end+1} = group{g}(order_(end:-1:1));
        end%for
    else
        rep_to_every_stim = cell(1,numel(size(group)));
        [rep_to_every_stim{:}] = deal(:);
        [rep_to_every_stim{1:n_stims}] = deal(1);
        sub_group = group(rep_to_every_stim{:});
        sub_stops = stops(rep_to_every_stim{:});
        sub_group = permute(sub_group, numel(size(sub_group)):-1:1);
        sub_stops = permute(sub_stops, numel(size(sub_stops)):-1:1);
        for g = 1:numel(sub_group)
            [stops_, order_] = sort(sub_stops{g});
            order{end+1} = sub_group{g}(order_(end:-1:1));
        end%for
    end%if
    order = [order{:}];
end%if
if isempty(order)
    error('not any trial was selected')
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
    if ~isempty(faces{b})
        patch('Faces', faces{b}, 'Vertices', vertices{b}, 'FaceColor', behavior{b,2}, 'EdgeColor', 'none')
        hold on
    end%if
end%for

for o = 1:numel(stimulus_onset)
    x0 = relative_onset(o);
    plot([x0, x0], [0, y0], '-', 'Color', o_color, 'LineWidth', o_linewidth);
end%for

if ~isempty(max_time_from_onset)
    xlim([-max_time_from_onset(1), stimulus_max_interval+max_time_from_onset(2)])
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
