function drawPrecue(windowPtr,priorityset)

nPriorities = length(priorityset);

% get experimental settings
settings = getExperimentalSettings();

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
quadrantdirection = [1 1; -1 1; -1 -1; 1 -1];       % moves clockwise from bottom right
smoothline = 1;                                     % smooth lines setting 1
fixationlinewidth = 4;                              % thickness of priority barss (pixels)
bordertopthickness = 1;                             % thickness of black border on priority bars (pixels)
dotradius = 4;                                      % radius of center dot of fixation (pixels)
dotrect = [center-dotradius center+dotradius];      % rect for center dot of fixation (pixels)
circlerect = [center-fixationsize center+fixationsize];

% color settings and variables
cmap = settings.cmap;   % L*a*b starting and ending colors
cslope = 1/(maxpriority-minpriority);
cintercept = -cslope*minpriority;
% nColors = size(cmap,1);     % number of unique colors in this continuum
% pVec = linspace(minpriority,maxpriority,nColors); % colors will correspond to these proportions

% plot precue
Screen('BlendFunction', windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % line to allow blending
Screen('FillOval', windowPtr, 0, dotrect);              % dot at center of screen
Screen('FrameOval', windowPtr, 0, circlerect, 2);       % circle fixation

xy = [fixationsize -fixationsize 0 0; 0 0 -fixationsize fixationsize];
Screen('DrawLines', windowPtr, xy, 1, 0, center);

% get colors and sizes for each quadrant
% colors = nan(3,2*nPriorities);
% xy = zeros(2,nPriorities*2);
for ipriority = 1:nPriorities
    priority = priorityset(ipriority);  % current priority
    
%     [~,idx] = min(abs(priority-pVec));
    %     colors(:,2*ipriority-1) = cmap(idx,:);
    %     colors(:,2*ipriority) = cmap(idx,:);
    p = cslope*priority+cintercept;
    color = lab2rgb((1-p).*cmap(1,:)+p.*cmap(2,:))*255;
    arclength = priority*90;
    Screen('FrameArc', windowPtr, color, circlerect, 90*ipriority+45-arclength/2, arclength, 4);
    %     xy(:,2*ipriority) = round(sqrt(2)/2*(fixationsize-dotradius)*priority.*quadrantdirection(ipriority,:));
end
% colors = colors*255;                          % color coordinates from 0 to 255
% xy = xy + dotradius*sign(xy);                 % end coordinates of priority lines   
% xyb = xy + bordertopthickness*sign(xy);       % end coordintes for the black border




% Screen('DrawLines', windowPtr, xyb, fixationlinewidth+2.5, 0, center, smoothline);          % priority lines
% Screen('DrawLines', windowPtr, xy, fixationlinewidth ,colors, center, smoothline);          % priority lines
% Screen('FillOval', windowPtr, 0, dotrect);
% Screen('FrameOval', windowPtr, 0, circlerect, 2);    % circle fixation
