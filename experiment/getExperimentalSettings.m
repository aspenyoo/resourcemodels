function settings = getExperimentalSettings()

h = figure(99);
cmap = colormap('parula'); % get a rough colormap
close(h)

prioritySets = [0.6 0.3  0.1   0.0;
                0.5 0.25 0.125 0.125];
      
% experiment info
settings.nTrials = 120;
settings.nRuns = 10;

% screen info
screenDistance = 60;                      % distance between observer and Screen (in cm)
screenHeight = 28.5;                        % height of screen (cm)
ScreenNumber=max(Screen('Screens'));       % use external Screen if exists
[w, h]=Screen('WindowSize', ScreenNumber);  % Screen resolution
screenResolution = [w h];                 % Screen resolution
screenCenter = screenResolution/2;       % Screen center
screenAngle = 2*(atand((screenHeight/2) / screenDistance)) ; % total visual angle of Screen
ppd = screenResolution(2) / screenAngle;  % pixels per degree
settings.apertureSize = 12;                % radius of aperature (dva)            
settings.bgColor = 125;
settings.fgColor = 200;
settings.minpriority = min(prioritySets(:));
settings.maxpriority = max(prioritySets(:));

% stimulus info 
fixationsize = 1;                     % fixation size (dva)
dotradius = 4;                                      % radius of center dot of fixation (pixels)
stimulusecc = 7;            % annulus stimulus eccentricity is on (dva)
jitter = 0.5;               % SD of uniform jitter distribution (dva)

% adding variables into the struct
settings.fixationsize = fixationsize;
settings.dotradius = dotradius;
settings.stimulusecc = stimulusecc;
settings.jitter = jitter;
settings.screenResolution = screenResolution;
settings.screenCenter = screenCenter;
settings.screenDistance = screenDistance;
settings.screenAngle = screenAngle;
settings.ppd = ppd;
settings.prioritySets = prioritySets;
settings.cmap = cmap;

end