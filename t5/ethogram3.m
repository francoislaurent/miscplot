function order=ethogram3(m, varargin)
%ETHOGRAM3   Plot ethogram.
%   ORDER=ETHOGRAM3(DATA_FILE_OR_STRUCT, KEY1, VALUE1, ...)

if ischar(m)
    m = load(m);
end%if
try
    m = m.ethogramme;
end%try

behavior_description = {};
max_window = [];
event = [];
event_type = [];
event_onset_time = [];
event_onset_tolerance = 0;
event_min_response_time = 0;
event_max_response_time = 0;
event_transitions_only = false;
grouping = {};
sort_by = {'onset=1,end'};
fill_gaps = 0; % boolean
ymin = []; ymax = [];

for v = 1:2:nargin-2
    switch lower(varargin{v})
    case 'behavior'
        behavior_description = varargin{v+1};
    case 'maxwindow'
        max_window = varargin{v+1};
    case 'eventonsettime'
        event_onset_time = varargin{v+1};
    case 'eventonsettolerance'
        event_onset_tolerance = varargin{v+1};
    case 'eventminresponsetime'
        event_min_response_time = varargin{v+1};
    case 'eventmaxresponsetime'
        event_max_response_time = varargin{v+1};
    case 'eventtype'
        event_type = varargin{v+1};
    case 'eventtime'
        event_time = varargin{v+1};
    case 'event'
        event = varargin{v+1};
    case 'transitionsonlyevent'
        event_transitions_only = varargin{v+1};
    case 'grouping'
        grouping = varargin{v+1};
    case 'sortby'
        sort_by = varargin{v+1};
    case 'fillgaps'
        fill_gaps = varargin{v+1};
    case 'ymin'
        ymin = varargin{v+1};
    case 'ymax'
        ymax = varargin{v+1};
    otherwise
        warning(['unknown argument: ''' varargin{v} ''''])
    end%switch
end%for

% behavior description
if isempty(behavior_description)
    behavior_description = { ...
        'hunch',    'hunch_large_squares',  'blue'; ...
        'bend',     'cast_large_squares',   'red'; ...
        'back',     'back_large_squares',   'cyan'; ...
        'stop',     'stop_large_squares',   'green'; ...
        'run',      'run_large_squares',    'black'; ...
        'roll',     'roll_large_squares',   'yellow'; ...
        'small',    'small_motion_squares', .5 * [1,1,1]; ...
        };
end%if
specified_behaviors = behavior_description(:,2);

behaviors = fieldnames(m(1));
if isempty(behavior_description)
    behavior_description = cell(0, 3);
end%if
specified_behaviors = behavior_description(:,2);
for b = 1:numel(behaviors)
    behavior = behaviors{b};
    if ~any(strcmp(behavior, specified_behaviors))
        extra_color = rand(1, 3);
        extra_color = extra_color .* (.9 / sqrt(sum(extra_color .* extra_color)));
        behavior_description(end+1,:) = {behavior, behavior, extra_color};
    end%if
end%for

% data structuring
ntrials = numel(m);
nbehaviors = size(behavior_description,1);
trials = cell(ntrials, nbehaviors);
trial_defined = true(ntrials,1);
for t=1:ntrials
    for b=1:nbehaviors
        t0t1 = cell(0,1);
        N = m(t).(behavior_description{b,2});
        for n=1:numel(N)
            t0t1{end+1,1} = reshape(N(n).x([2,3]),1,2);
        end%for
        t0t1 = cell2mat(t0t1);
        if ~isempty(t0t1) && ~isempty(max_window)
            t0t1 = t0t1(t0t1(:,1)<max_window(2) & max_window(1)<t0t1(:,2),:);
        end%if
        trials{t,b} = t0t1;
    end%for
end%for
clear m N t0t1
trials = trials(trial_defined,:);
ntrials = size(trials,1)
original_trial_index = find(trial_defined);

% event selection
%events = cell(size(trials));
%for t=1:ntrials
%    for b=1:nbehaviors
%        events{t,b} = {[], (1:size(trials{t,b},1))'};
%    end%for
%end%for

criteria = cell(0,2);

% parse event type
if isempty(event_type)
    event_type = {};
elseif ischar(event_type)
    event_type = {event_type};
end%if
event_type_cumulated_filter = empty(size(event_type));
event_type_global_individual_filter = empty(size(event_type));
event_type_individual_filters = empty(size(event_type));
event_type_onset = zeros(size(event_type));
for i=1:numel(event_type)
    if isempty(event_type{i})
        event_type{i} = {};
        continue
    end%if
    % parse event_type
    event_type{i} = strsplit(event_type{i}, ';');
    for e=1:numel(event_type{i})
        e_ = event_type{i}{e};
        event_type{i}{e} = cell(1,0);
        % central event marker
        e__ = '^';
        if ~isempty(e_) && e_(1) == e__
            if event_type_onset(i)
                error('multiple event type onsets')
            end%if
            event_type_onset(i) = e;
            %event_type{i}{e}{1,end+1} = e__;
            e_ = e_(2:end);
        end%if
        % joker marker
        e__ = '.';
        all_ = ~isempty(e_) && e_(1) == e__;
        if all_
            e_ = e_(2:end);
        end%if
        % exclusion marker
        e__ = '~';
        neg_ = ~isempty(e_) && e_(1) == e__;
        if neg_
            if all_
                error('both `.` and `~` found')
            end%if
            %event_type{i}{e}{1,end+1} = e__;
            e_ = e_(2:end);
        end%if
        % post cardinal filters
        pre_filters = {};
        post_filters = {};
        if ~isempty(e_) && e_(end) == ']'
            bracket = strfind(e_, '[')(end);
            filtstr = e_(bracket+1:end-1);
            post_filters = parse_filters(filtstr);
            e_ = e_(1:bracket-1);
        end%if
        % cardinal
        card = 1;
        if ~isempty(e_) && e_(end) == '}'
            bracket = strfind(e_, '{');
            cardstr = e_(bracket+1:end-1);
            card = strsplit(cardstr, ':');
            e_ = e_(1:bracket-1);
            if ~any(numel(card) == [1,2])
                error(['wrong cardinal: ''' cardstr ''''])
            end%if
            for c=1:numel(card)
                if isempty(card{c})
                    if c == 1
                        card{c} = 1;
                    else
                        card{c} = inf;
                    end%if
                else
                    card{c} = str2num(card{c});
                end%if
            end%for
            card = [card{:}];
        else
            % post_filters is actually pre_filters
            pre_filters = post_filters;
            post_filters = {};
        end%if
        event_type{i}{e}{1,end+1} = card;
        % pre cardinal filters
        if isempty(pre_filters) && ~isempty(e_) && e_(end) == ']'
            bracket = strfind(e_, '[');
            filtstr = e_(bracket+1:end-1);
            pre_filters = parse_filters(filtstr);
            e_ = e_(1:bracket-1);
        end%if
        % remove parentheses
        block_ = ~isempty(e_) && e_(1) == '(' && e_(end) == ')';
        if block_
            e_ = e_(2:end-1);
        end%if
        % behaviors
        if isempty(e_)
            if ~all_
                warning('missing joker `.`')
            end%if
            e_ = 1:nbehaviors;
            individual_filters = cell(size(e_));
        else
            if all_
                error('both `.` and behavior(s) found')
            end%if
            e_ = strsplit(e_, ',');
            individual_filters = cell(size(e_));
            for b=1:numel(e_)
                e__ = e_{b};
                if e__(end) == ']'
                    bracket = strfind(e__, '[');
                    filtstr = e__(bracket+1:end-1);
                    individual_filters{b} = parse_filters(filtstr);
                    e__ = e__(1:bracket-1);
                end%if
                e_{b} = find(strcmp(e__, behavior_description(:,1)));
                if ~isscalar(e_{b})
                    error(['wrong behavior: ''' e__ ''''])
                end%if
            end%for
            e_ = [e_{:}];
            if neg_
                if ~all(cellfun(@isempty, individual_filters, 'UniformOutput', true))
                    error('individual filters on complemented behaviors not supported')
                end%if
                e__ = true(1,nbehaviors);
                e__(e_) = false;
                e_ = find(e__);
            elseif ~block_ && ~isempty(pre_filters) && 1 < numel(e_) && isempty(individual_filters{end})
                individual_filters{end} = pre_filters;
                pre_filters = {};
            end%if
        end%if
        if ~isempty(post_filters)
            error('cumulated filters not implemented yet')
        end%if
        if ~all(cellfun(@isempty, individual_filters, 'UniformOutput', true))
            error('individual filters not implemented yet')
        end%if
        event_type{i}{e}{1,end+1} = e_;
        event_type_cumulated_filter{i}{e} = post_filters;
        event_type_global_individual_filter{i}{e} = pre_filters;
        event_type_individual_filters{i}{e} = individual_filters;
    end%for
    if event_type_onset(i) == 0
        event_type_onset(i) = 1;
    end%if
    % filter onset type
    %onset = false(1,nbehaviors);
    %onset(event_type{i}{event_type_onset(i)}{end}) = true;
    %[events{:,find(~onset)}] = deal(cell(0,2));
end%for

% temporal reconstruction of the trials
bhr_col = 1;
begin_col = 2;
end_col = 3;
trials_ = cell(ntrials,1);
bs_ = 1:nbehaviors;
for t=1:size(trials,1)
    xs = cell(size(trials,2),1);
    ns = zeros(1,size(trials,2));
    for b=1:size(trials,2)
        xs{b} = trials{t,b};
        ns(b) = size(xs{b},1);
    end%for
    if sum(ns)==0
        continue
    end%if
    xs = [repelems(bs_, [bs_;ns])', cell2mat(xs)];
    [js,js] = sort(xs(:,begin_col));
    trials_{t,1} = xs(js,:);
end%for
trials = trials_;
clear xs ns is js trials_

% onset time
if isempty(event_onset_time)
    error('not implemented: EventOnsetTime not defined')
else
    if isscalar(event_transitions_only) && ~isscalar(event_onset_time)
        event_transitions_only = repmat(event_transitions_only, 1, numel(event_onset_time));
    end%if
    criteria = cell(1,numel(event_onset_time));
    events = cell(ntrials,numel(event_onset_time));
    [events{:}] = deal(cell(0,2));
    for t=1:size(trials,1)
            x = trials{t};
            if isempty(x)
                %error('missing states')
                continue
            end%if
            j = 1:size(x,1);
            for e=1:numel(event_onset_time)
                x0 = x(j,begin_col) - (event_onset_time(e)+event_min_response_time);
                if event_onset_tolerance
                    candidate = (0 <= x0) & (x0 <= event_onset_tolerance);
                    if ~any(candidate)
                        candidate = (-event_onset_tolerance <= x0) & (x0 <= (event_onset_tolerance+event_max_response_time-event_min_response_time));
                    end%if
                else
                    candidate = (0 <= x0) & (x0 <= (event_max_response_time-event_min_response_time));
                end%if
                if any(candidate)
                    candidate = j(find(candidate, 1));
                    b = x(candidate,bhr_col);
                    if event_transitions_only(e)
                        events{t,e}(end+1,:) = {b, candidate};
                    else
                        events{t,e}(end+1,:) = {[b, 1], candidate};
                    end%if
                elseif ~event_transitions_only(e)
                    x1 = x(j,end_col) - (event_onset_time(e)+event_min_response_time);
                    candidate = j(find(-event_onset_tolerance < x1, 1)); % x0s are ordered, so as x1s
                    b = x(candidate,bhr_col);
                    events{t,e}(end+1,:) = {[b, 0], candidate};
                end%if
            end%for
    end%for
    for e=1:numel(criteria)
        criteria{e}{end+1,1} = 'state';%'onset'
        criteria{e}{end,2} = 1:nbehaviors;%1:numel(event_onset_time);
        if ~event_transitions_only(e)
            criteria{e}{end+1,1} = 'transition';
            criteria{e}{end,2} = [0,1];
        end%if
    end%for
    clear x x0 j js
    [events{find(any(cellfun(@isempty,events),2)),:}] = deal(cell(0,2));
end%if

% event type
if ~isempty(event_type)
    if size(events,2) ~= size(event_type,2)
        if size(events,2)==1 && size(event_type,1)==1 && 1<size(event_type,2)
            event_type = permute(event_type, [2,1]);
        else
            error('wrong number of event types')
        end%if
    end%if
    for t=1:size(events,1)
        trial = trials{t};
        events_t = events(t,:);
        [events{t,:}] = deal(cell(0,2));
        for i=1:size(event_type,1)
            ks = cell(1,size(events,2));
            [ks{:}] = deal(cell(0,2));
            for e=1:size(events,2)
                for k=1:size(events_t{e},1)
                    js = events_t{e}{k,2};
                    if isempty(js)
                        continue
                    end%if
                    % assert size(js,2)==1
                    if isempty(event_type{i,e})
                        is = [js;js;js]';
                        ks{e}(end+1,:) = {[events_t{e}{k,1}, i], is};
                        continue
                    end%if
                    is = zeros(0,3);
                    for j=1:numel(js)
                        d0 = js(j);
                        discard = false;
                        p = event_type_onset(i,e);
                        if p == 0
                            error('event_type_onset==0')
                        end%if
                        dpre = d0 - 1;
                        dpost = d0; % if ~discard, dpost will be dmax+1, just like dpre is dmin-1
                        for q=p:numel(event_type{i,e})
                            [card,behaviors] = event_type{i,e}{q}{:};
                            filter = event_type_global_individual_filter{i,e}{q};
                            c = 0;
                            while ~discard && c < card(1)
                                if dpost <= size(trial,1)
                                    x = trial(dpost,:);
                                    behavior = x(bhr_col);
                                    discard = ~any(behaviors == behavior);
                                    if ~discard && ~isempty(filter)
                                        discard = apply_filter(filter, x);
                                    end%if
                                else
                                    discard = true;
                                end%if
                                c = c + 1;
                                dpost = dpost + 1;
                            end%while
                            if discard
                                break
                            end%if
                            if ~isscalar(card)
                                while c < card(2) && dpost <= size(trial,1) && any(behaviors == trial(dpost,bhr_col))
                                    if isempty(filter)
                                        x = trial(dpost,:);
                                        stop = apply_filter(filter, x);
                                        if stop
                                            dpost = dpost + 1;
                                            break
                                        end%if
                                    end%if
                                    c = c + 1;
                                    dpost = dpost + 1;
                                end%while
                            end%if
                        end%for
                        if discard
                            continue
                        end%if
                        for q=p-1:-1:1
                            [card,behaviors] = event_type{i,e}{q}{:};
                            filter = event_type_global_individual_filter{i,e}{q};
                            c = 0;
                            while ~discard && c < card(1)
                                if d == 0
                                    discard = true;
                                else
                                    x = trial(dpre,:);
                                    behavior = x(bhr_col);
                                    discard = ~any(behaviors == behavior);
                                    if ~discard && ~isempty(filter)
                                        discard = apply_filter(filter, x);
                                    end%if
                                end%if
                                c = c + 1;
                                dpre = dpre - 1;
                            end%while
                            if discard
                                break
                            end%if
                            if isscalar(card)
                                while c < card(2) && 0 < dpre && any(behaviors == trial(dpre,bhr_col))
                                    if isempty(filter)
                                        x = trial(dpre,:);
                                        stop = apply_filter(filter, x);
                                        if stop
                                            dpre = dpre - 1;
                                            break
                                        end%if
                                    end%if
                                    c = c + 1;
                                    dpre = dpre - 1;
                                end%while
                            end%if
                        end%for
                        if ~discard
                            is(end+1,:) = [dpre+1, d0, dpost-1];
                        end%if
                    end%for
                    if 0 < size(is,1)
                        ks{e}(end+1,:) = {[events_t{e}{k,1}, i], is};
                    end%if
                end%for
            end%for
            if ~any(cellfun(@isempty, ks))
                for e=1:size(events,2)
                    events{t,e} = [events{t,e}; ks{e}];
                end%for
            end%if
        end%for
    end%for
    for e=1:numel(criteria)
        criteria{e}{end+1,1} = 'type';
        criteria{e}{end,2} = 1:size(event_type,1);
    end%for
    %[events{find(any(cellfun(@isempty,events),2)),:}] = deal(cell(0,2));
end%if


group = {reshape(1:ntrials,ntrials,1)};
if ~isempty(grouping)
    if ~iscell(grouping)
        grouping = strsplit(grouping, ';');
    end%if
    G = numel(grouping);
    for g=1:G
        new_group = cell(0,1);
        if ~iscell(grouping{g})
            grouping{g} = strsplit(grouping{g}, ';');
            for h=1:numel(grouping{g})
            	grouping{g}{h} = strsplit(grouping{g}{h}, ',');
            end%for
        end%if
        for h=1:numel(grouping{g})
        events_ = events;
        if 1 < numel(grouping{g}{h})
            events_ = apply_filters(grouping{g}{h}(1:end-1), events_, trials, group, behavior_description, bhr_col, criteria);
        end%if
        gh = grouping{g}{h}{end};
        if strncmp(gh, 'state', 5)
            cursor = [];
            if 5 < numel(gh)
                if ~( gh(6) == '(' && gh(end) == ')' )
                    warning('parentheses expected after ''state''')
                end%if
                cursor = str2num(gh(7:end-1));
            end%if
            if isempty(cursor)
                cursor = 0;
            end%if
            if false%cursor == 0
                for i=1:numel(group)
                    for e=1:size(events_,2)
                        [i,e]
                        new_group(end+1,1) = group{i}(~cellfun(@isempty, events_(group{i},e), 'UniformOutput', true));
                        new_group{end}
                    end%for
                end%for
            else
                for i=1:numel(group)
                    %group{i} = reshape(group{i},numel(group{i}),1);
                    group__ = cell(nbehaviors,1);
                    for t=group{i}'
                        for e=1:size(events_,2)
                            evts_ = events_{t,e};
                            for k=1:size(evts_,1)
                                s_ = evts_{k,2}(2)+cursor;
                                if 0<s_ && s_<=size(trials{t},1)
                                    b_ = trials{t}(s_,bhr_col);
                                    group__{b_}(end+1,1) = t;
                                end%if
                            end%for
                        end%for
                    end%for
                    new_group = [new_group; group__];
                end%for
            end%if
        else
            if strcmp(gh,'type')
                events_ = events_(:,1);
            end%if
            feature = find(strcmp(gh, criteria{1}(:,1)));
            F = criteria{1}{feature,2}; % workaround
            for i=1:numel(group)
                %group{i} = reshape(group{i},numel(group{i}),1);
                group__ = cell(numel(F),1);
                for e=1:size(events_,2)
                    F = criteria{e}{feature,2};
                    for t=group{i}'
                        evts_ = events_{t,e};
                        if isempty(evts_)
                            continue
                        end%if
                        K = size(evts_,1);
                        f = zeros(K,1);
                        for k=1:K
                            f(k) = evts_{k,1}(feature);
                        end%for
                        for f_=1:numel(F)
                            if any(F(f_)==f)
                                group__{f_}(end+1,1) = t;
                            end%if
                        end%for
                    end%for
                end%for
                new_group = [new_group; group__];
            end%for
        end%if
        end%for
        new_group(cellfun(@isempty, new_group)) = [];
        group = new_group;
    end%for
end%if
if ~isempty(sort_by)
    sort_by = strsplit(sort_by, ';');
    for s=1:numel(sort_by)
    sort_by_ = strsplit(sort_by{s}, ',');
    if 1 < numel(sort_by_)
        events_ = apply_filters(sort_by_(1:end-1), events, trials, group, behavior_description, bhr_col, criteria);
    end%if
    if 1 < size(events_,2)
        warning('SortBy should select a single event onset')
    end%if
    sort_by_ = sort_by_{end};
    sort_by_ = strsplit(sort_by_, '.');
    s_ = sort_by_{1};
    cursor = [];
    if strncmp(s_, 'state', 5)
        if 5 < numel(s_)
            if ~( s_(6) == '(' && s_(end) == ')' )
                warning('parentheses expected after ''state''')
            end%if
            cursor = str2num(s_(7:end-1));
        end%if
        sort_by_ = sort_by_(2:end);
    end%if
    if isempty(cursor)
        cursor = 0;
    end%if
    switch sort_by_{1}
    case 'begin'
        score = @(a)a(begin_col);
    case 'end'
        score = @(a)a(end_col);
    case 'duration'
        score = @(a)a(end_col)-a(begin_col);
    otherwise
        error('not supported')
    end%switch
    if isscalar(sort_by_)
        sorting_order = 'descend';
    else
        sorting_order = sort_by_{2};
    end%if
    for g=1:numel(group)
        scores = nan(numel(group{g}),1);
        for i=1:numel(group{g})
            t = group{g}(i);
            evt_ = zeros(0,1);
            for e=1:size(events_,2)
                if ~isempty(events_{t,e})
                    evt_ = [evt_; cell2mat(events_{t,e}(:,2))(:,2)];
                end%if
            end%for
            if isempty(evt_)
                %warning(['no event left for trial ' num2str(original_trial_index(t))]);
                continue
            elseif ~isscalar(evt_)
                %warning(['multiple events for trial ' num2str(original_trial_index(t)) ' for sorting'])
                evt_ = min(evt_);
            end%if
            scores(i) = score(trials{t}(evt_,:));
        end%for
        [scores, order] = sort(scores, sorting_order);
        order = order(~isnan(scores));
        if ~isempty(order)
            group{g} = group{g}(order);
        end%if
    end%for
    end%for
end%if
order = cell2mat(group);
if any(diff(sort(order))==0)
    warning('duplicate trials')
end%if
ntrials = size(order,1)

if fill_gaps
    for t=1:ntrials
        stop = trials{order(t)}(1:end-1,end_col);
        start = trials{order(t)}(2:end,begin_col);
        hisi = (start - stop) * .5;
        if any(hisi < 0)
            error('overlap')
        end%if
        trials{order(t)}(1:end-1,end_col) = stop + hisi;
        trials{order(t)}(2:end,begin_col) = start - hisi;
    end%for
end%if

if isempty(ymax)
    if isempty(ymin)
        ymin = 0;
        ymax = ntrials;
    else
        ymax = ymin + ntrials;
    end%if
elseif isempty(ymin)
    ymin = ymax - ntrials;
end%if
dy = 1;
y = ymax;
f = zeros(1, nbehaviors);
faces = cell(1, nbehaviors);
vertices = cell(1, nbehaviors);
[vertices{:}] = deal({});
[faces{:}] = deal({});
for t = 1:ntrials
    z = trials{order(t)};
    for k=1:size(z,1)
        b = z(k,bhr_col);
        x = z(k,[begin_col,end_col])';
        vertices{b}{end+1,1} = x([1 1 2 2]);
        vertices{b}{end,2} = [y-dy; y; y; y-dy];
        faces{b}{end+1,1} = f(b)+1:f(b)+4;
        f(b) = f(b) + 4;
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
        patch('Faces', faces{b}, 'Vertices', vertices{b}, 'FaceColor', behavior_description{b,3}, 'EdgeColor', 'none')
        hold on
    end%if
end%for

function filters=parse_filters(str_)
    filters = strsplit(str_, ';');
    for f_=1:numel(filters)
        filters{f_} = strsplit(filters{f_}, ',');
    end%for
    error_=@(f__,g__)error(['malformed expression: `' filters{f__}{g__} '`']); % error message
    ops_ = {'<=','>=','==','<','>'};
    LEFT_ = 0;
    RIGHT_ = 1;
    for f_=1:numel(filters)
        for g_=1:numel(filters{f_})
            members_ = filters{f_}(g_);
            for o_=1:numel(ops_)
                members__ = {};
                for m_=1:2:numel(members_)
                    submembers_ = strsplit(members_{m_}, ops_{o_});
                    if 3 < numel(submembers_)
                        error_(f_,g_)
                    end%if
                    members__ = [members__, submembers_(1)];
                    for s_=2:numel(submembers_)
                        members__ = [members__, ops_(o_), submembers_(s_)];
                    end%for
                end%for
                members_ = members__;
            end%for
            for m_=1:2:numel(members_)
                members_{m_} = strtrim(members_{m_});
            end%for
            switch numel(members_)
            case 1
                filters{f_}{g_} = members_;
                continue
            case 3
                op_ = members_{2};
                if isletter(members_{1}(1))
                    variable_ = members_{1};
                    value_ = {members_{3}, RIGHT_, op_};
                else
                    value_ = {members_{1}, LEFT_, op_};
                    variable_ = members_{3};
                end%if
            case 5
                variable_ = members_{3};
                value_ = {members_{1}, LEFT_, members_{2}; members_{5}, RIGHT_, members_{4}};
            otherwise
                % cannot be
                error_(f_,g_)
            end%switch
            if ~isletter(variable_(1))
                error_(f_,g_)
            end%if
            filter_ = cell(1,size(value_,1));
            for v_=1:size(value_,1)
                val_ = str2num(value_{v_,1});
                if isempty(val_)
                    error_(f_,g_)
                end%if
                pos_ = value_{v_,2};
                op_ = value_{v_,3};
                switch op_
                case '<='
                    if pos_ == LEFT_
                        filter_{v_} = @(a)val_<=a;
                    else%if pos_ == RIGHT_
                        filter_{v_} = @(a)a<=val_;
                    end%if
                case '>='
                    if pos_ == LEFT_
                        filter_{v_} = @(a)a<=val_;
                    else%if pos_ == RIGHT_
                        filter_{v_} = @(a)val_<=a;
                    end%if
                case '=='
                    if size(value_,1)==1
                        filter_{v_} = @(a)a==val_;
                    else
                        error_(f_,g_)
                    end%if
                case '<'
                    if pos_ == LEFT_
                        filter_{v_} = @(a)val_<a;
                    else%if pos_ == RIGHT_
                        filter_{v_} = @(a)a<val_;
                    end%if
                case '>'
                    if pos_ == LEFT_
                        filter_{v_} = @(a)a<val_;
                    else%if pos_ == RIGHT_
                        filter_{v_} = @(a)val_<a;
                    end%if
                end%switch
            end%for
            filters{f_}{g_} = [{variable_}, filter_];
        end%for
    end%for
end%function

function cells=empty(siz)
    cells = cell(siz);
    [cells{:}] = deal({});
end%function

function discard=apply_filter(filter, x)
    duration = [];
    U = numel(filter);
    union = false;
    u = 1;
    while ~union && u<=U
        V = numel(filter{u});
        intersection = true;
        v = 1;
        while intersection && v<=V
            filt = filter{u}{v};
            variable = filt{1};
            switch variable
            case 'duration'
                if isempty(duration)
                    duration = x(end_col)-x(begin_col);
                end%if
                val = duration;
            otherwise
                error(['unsupported feature: `' variable '`'])
            end%switch
            W = numel(filt);
            w = 2;
            while intersection && w<=W
                intersection = filt{w}(val);
                w = w + 1;
            end%while
            v = v + 1;
        end%while
        union = intersection;
        u = u + 1;
    end%while
    discard = ~union;
end%function

function events_=apply_filters(filters_, events_, trials_, groups_, behavior_description, bhr_col, criteria)
    % only for Grouping and SortBy
    ntrials_ = size(trials_,1);
    [ntrials__, nevts_] = size(events_);
    if ntrials_ ~= ntrials__
        error('wrong size for events')
    end%if
    if isempty(groups_)
        groups_ = {1:ntrials_};
    end%if
    %if 1 < nevts_ % only for EventType
    %    if size(filters_,2) ~= nevts_
    %        if size(filters_,1) == nevts_
    %            filters_ = filters_';
    %        elseif ~isscalar(filters_)
    %            error('wrong filter size')
    %        end%if
    %    end%if
    %end%if
    for f_=1:numel(filters_)
        filter_ = filters_{f_};
        equalsign_ = strfind(filter_, '=');
        if ~isscalar(equalsign_)
            error('malformed expression')
        end%if
        feature_ = filter_(1:equalsign_-1);
        value_ = str2num(filter_(equalsign_+1:end));
        if strcmp(feature_, 'onset')
            events_ = events_(:, value_); % select column or make the other columns empty?
        elseif strncmp(feature_, 'state', 5)
            value_ = find(strcmp(value_,behavior_description(:,1)));
            cursor_ = [];
            if 5 < numel(feature_)
                if ~( feature_(6) == '(' && feature_(end) == ')' )
                    warning('parentheses expected after ''state''')
                end%if
                cursor_ = str2num(feature_(7:end-1));
            end%if
            if isempty(cursor_)
                cursor_ = 0;
            end%if
            if false%cursor_ == 0
                for i_=1:numel(groups_)
                    for b_=[1:value_-1,value_+1,nbhvrs_]
                        for t_=groups_{i_}
                            events_{t_,e_} = cell(0,2);
                        end%for
                    end%for
                end%for
            else
                for i_=1:numel(groups_)
                    for t_=groups_{i_}
                        for e_=1:size(events_,2)
                            evts_ = events_{t_,e_};
                            select_ = false(size(evts_,1),1);
                            for k_=1:size(evts_,1)
                                s_ = evts_{k_,2}(2)+cursor_;
                                if 0<s_ && s_<=size(trials_{t_})
                                    b_ = trials_{t_}(s_,bhr_col);
                                    select_(k_) = b_ == value_;
                                end%if
                            end%for
                            events_{t_,e_} = evts_(select_,:);
                        end%for
                    end%for
                end%for
            end%if
        else
            feature_ = find(strcmp(feature_, criteria{1}(:,1))); % workaround
            for i_=1:numel(groups_)
                for t_=groups_{i_}'
                    for e_=1:size(events_,2)
                        evts_ = events_{t_,e_};
                        select_ = true(size(evts_,1),1);
                        for k_=1:size(evts_,1)
                            select_(k_) = evts_{k_,1}(feature_) == value_;
                        end%for
                        events_{t_,e_} = evts_(select_,:);
                    end%for
                end%for
            end%for
        end%if
    end%for
end%function

end%function
