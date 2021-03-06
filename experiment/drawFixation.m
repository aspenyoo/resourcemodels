function drawFixation(windowPtr)

nPriorities = length(priorityset);

% get experimental settings
settings = getExperimentalSettings();
cmap = settings.cmap;   % colormap
nColors = size(cmap,1);     % number of unique colors in this continuum


% screen settings
center = settings.screenCenter;
ppd = settings.ppd;
bgColor = settings.bgColor;
fgColor = settings.fgColor;
minpriority = settings.minpriority;
maxpriority = settings.maxpriority;
fixationsize = settings.fixationsize; % size of fixation (dva)

% fixation settings
fixationsize = fixationsize*ppd;                    % fixation size (pixels)
quadrantdirection = [1 1; -1 1; -1 -1; 1 -1];       % moves counterclockwise from 
smoothline = 1;                                     % smooth lines setting 1
fixationlinewidth = 4;                              % thickness of priority barss (pixels)
bordertopthickness = 1;                             % thickness of black border on priority bars (pixels)
dotradius = 4;                                      % radius of center dot of fixation (pixels)
dotrect = [center-dotradius center+dotradius];      % rect for center dot of fixation (pixels)
% fixationsize = 30; % radius of fixation in pixels

% % open screen and get screen settings
% screenNumber = max(Screen('Screens'));
% topPriorityLevel = MaxPriority(screenNumber); %topPriorityLevel1 = MaxPriority(window);
% Priority(topPriorityLevel);
% [windowPtr, screenRect] = Screen('OpenWindow', screenNumber, [0 0 0]);
Screen('BlendFunction', windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% [~, h] = Screen('WindowSize', windowPtr);


% get colors for bars
pVec = linspace(minpriority,maxpriority,nColors); % colors will correspond to these proportions

colors = nan(3,2*nPriorities);
xy = zeros(2,nPriorities*2);
for ipriority = 1:nPriorities
    priority = priorityset(ipriority);  % current priority
    
    [~,idx] = min(abs(priority-pVec));
    colors(:,2*ipriority-1) = cmap(idx,:);
    colors(:,2*ipriority) = cmap(idx,:);
    
    xy(:,2*ipriority) = round(sqrt(2)/2*(fixationsize-dotradius)*priority.*quadrantdirection(ipriority,:));
end
colors = colors*255;
% xyb = xy+(sqrt(2)/2*(dotradius+bordertopthickness))*sign(xy);
xy = xy + dotradius*sign(xy);                 % end coordinates of priority lines   
xyb = xy + bordertopthickness*sign(xy);       % end coordintes for the black border
% xy = xy+(sqrt(2)/2*dotradius)*sign(xy);

circlerect = [center-fixationsize center+fixationsize];
% Screen('fillRect',windowPtr,bgColor);       % background color
% Screen('FillOval', windowPtr, 255, circlerect);    % white background circle fixation
Screen('DrawLines', windowPtr, xyb, fixationlinewidth+2.5, 0, center, smoothline);          % priority lines
Screen('DrawLines', windowPtr, xy, fixationlinewidth ,colors, center, smoothline);          % priority lines
Screen('FillOval', windowPtr, 0, dotrect);
Screen('FrameOval', windowPtr, 0, circlerect, 2);    % circle fixation
% Screen('Flip', windowPtr);
% 
% pause;
% sca