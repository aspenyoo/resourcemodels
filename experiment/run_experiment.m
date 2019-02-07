% RUN_EXPERIMENT runs an experiment with different
%
% code substantially modified from Afredo's file pri_4targ_wm_LngTR.m


function run_experiment(subjectID, run, iset, do_el)
if nargin<4; do_el=0; end

% Screen Preferences
% Screen('Preference', 'SkipSyncTests', 1); %comment out when actually running a participant
Screen('Preference', 'VisualDebugLevel', 1); %turns PTB welcome screen black
AssertOpenGL;

laSemina=sum(100*clock);
rng(laSemina);

try
    % Open onscreen window with default settings
    screenNumber = max(Screen('Screens'));
    topPriorityLevel = MaxPriority(screenNumber); %topPriorityLevel1 = MaxPriority(window);
    Priority(topPriorityLevel);
    [windowPtr] = Screen('OpenWindow', screenNumber, [0 0 0]);
    Screen('BlendFunction', windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % define dot/fixation/pizza slices field parameter
%     dot_w = 0.2;           % width of dot (deg)
    
    % timing (seconds)
    time_start = 1.9;
    time_fixExpand = 0.1;
    time_precue = 0.7;
    time_interval = 0.1;
    time_present = 0.1;
    time_ISI = 1;
    time_responseCue = 0.8;
    time_feedback = 0.8;
    time_ITI = 1;
    
    % relevant response keys
    KbName('UnifyKeyNames');
    if ismac == 1 %IsOSX
        ESCa_key = KbName('escape'); % press this key to abort
    else
        ESCa_key = KbName('escape');
    end
    SPC_key = KbName('space');
    
    % define colors
    white = WhiteIndex(screenNumber) ;
    black = BlackIndex(screenNumber);
    darkgray = white/2.2;
    
    
    % =======================================================
    % ASPEN STUFF BELOW
    %====================================
    settings = getExperimentalSettings;
    
    % experimental settings
    nTrials = settings.nTrials;
    nRuns = settings.nRuns;
    nTrialsPerRun = nTrials/nRuns;
    prioritySets = settings.prioritySets;
    prioritySet = prioritySets(iset,:);
    
    % screen settings
    center = settings.screenCenter;
    ppd = settings.ppd;
    apertureSize= settings.apertureSize*ppd ;          % (pixels)
    apertureRect = [center-apertureSize center+apertureSize];
    
    % fixation settings
    fixationsize = settings.fixationsize;       % size of fixation (dva)
    fixationsize = fixationsize*ppd;            % size of fixation (pixels)
    dotradius = settings.dotradius;             % radius of fixation dot (pixels)
    circlerect = [center-fixationsize center+fixationsize];
    dotrect = [center-dotradius center+dotradius];      % rect for center dot of fixation (pixels)
    
    % stimulus settings
    nItems = 4;
    stimulusecc = settings.stimulusecc*ppd;       % annulus stimulus eccentricity is on (pixels)
    jitter = settings.jitter*ppd;                 % uniform jitter distribution (pixels)
    itemrad = 0.2;                                % radius of item (dva)
    itemsize = itemrad * ppd;                     % size of items (pixels)
    
    % create folder for subjectID if it doesn't exist
    folderDir = sprintf('output/%s/',subjectID);
    if ~exist(folderDir,'dir')
        mkdir(folderDir)
    end
    filename = sprintf('%s%s_pricond%d_designMat.mat',folderDir,subjectID,iset);
    if (run == 1) % if first run
        
        % ====== MAKE DESIGN MATRIX ======
        % column meanings for designMat for nItems = 4
        %   1-4: priorities of items in quadrant 1-4
        %     5: priority of target
        %   6-9: polar angle of items (radians)
        % 10-13: x coordinates of items relative to center of screen (pixels)
        % 14-17: x coordinates of items relative to center of screen (pixels)
        
        nonzeroPrioritySet = prioritySet(prioritySet~=0);
        nNonzeroPriorities = length(nonzeroPrioritySet);
        condmat = cell(1,nNonzeroPriorities);
        for ipriority = 1:nNonzeroPriorities
            priority = nonzeroPrioritySet(ipriority);
            tempmat = perms(prioritySet);       % all permutations of configurations
            tempmat = [tempmat priority*ones(size(tempmat,1),1)]; % current priority is the target in all these trials
            condmat{ipriority} = tempmat;
        end
        
        trial = 0;
        % multiply the condmats according to priority
        multiplier = roundn(nonzeroPrioritySet./min(nonzeroPrioritySet),-4);
        % do something for the case that multiplier does not result in integers
        %     if (sum(round(multiplier) == multiplier) ~= nItems)
        %         sprintf('fuck')
        %     end
        designMat = [];
        for ipriority = 1:nNonzeroPriorities
            tempmat = repmat(condmat{ipriority},multiplier(ipriority),1);
            designMat = [designMat; tempmat];
        end
        nrowsDesignMat = size(designMat,1);
        
        % locations of items
        possibleAngles = (10:10:80)/180*pi;
        randangles = possibleAngles(randi(length(possibleAngles),nrowsDesignMat,nItems));
        randangles(:,2) = randangles(:,2)+pi/2;
        randangles(:,3) = randangles(:,3)+pi;
        randangles(:,4) = randangles(:,4)+3*pi/2;
        designMat = [designMat randangles];
        
        % getting the correct number of trials
        designMat = repmat(designMat,ceil(nTrials/nrowsDesignMat),1);
        designMat(randperm(nrowsDesignMat),:) = designMat;
        designMat = designMat(1:nTrials,:);
        
        %     designMat = [.6 .3 .1 0 .6 pi/180*[20 110 200 290]];
        % change to pixel locations (relative to center)
        [X,Y] = pol2cart(designMat(:,end-nItems+1:end), stimulusecc);
        X = X + (rand(size(X))-0.5)*jitter;
        Y = Y + (rand(size(Y))-0.5)*jitter;
        designMat = [designMat round(X) round(Y)];
        
        save(filename,'designMat','trial')
    else
        load(filename,'designMat','trial')
    end
    
    % Get EyeTracker and Provide info
    if do_el
        el=EyelinkInitDefaults(windowPtr);
        Eyelink('Initialize','PsychEyelinkDispatchCallback') % initialises the eyetracker
        status=Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,AREA');
        Eyelink('command', 'sample_rate = 500');
        Eyelink('command','calibration_type=HV9');% updating number of callibration dots
        Eyelink('command', 'enable_automatic_calibration = YES');
        Eyelink('command','file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
        Eyelink('command','file_sample_data  = LEFT,RIGHT,GAZE,AREA,GAZERES,STATUS');
        % s=Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
        % s=Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
        %Calibrate the eye tracker
        EyelinkDoTrackerSetup(el);
        if status~=0
            error('link_sample_data error, status: ',status)
        end
        %get gaze data from EyeTracker
        edfFile=sprintf('%seyedata_%s_%02d.edf', folderDir, subjectID, run);
        Eyelink('openfile',edfFile);
    end
    
    
    % ============== FIRST SCREEN =================
    Screen('FillRect', windowPtr, black);
    Screen('FillOval',windowPtr, darkgray, apertureRect);
    runText=sprintf('This is Run %01.f \n', run);
    DrawFormattedText(windowPtr, [runText, 'Press the spacebar when you are ready to continue.'],'center', 'center', white);
    Screen('Flip',windowPtr);
    
    %Wait for keyboard press
    resp = 0;
    while resp == 0
        [resp, ~] = checkForResp(SPC_key, ESCa_key);
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            if do_el
                Eyelink('ShutDown');
            end
            return;
        end
    end
    
    clear resp;
    %     expStart(1)=GetSecs;
    if do_el
        eyeTstrt=Eyelink('TrackerTime');
        Eyelink('startrecording');
        Eyelink('Message','xDAT %i', 10);
    end
    %Draw Aperture
    Screen('FillOval',windowPtr, darkgray, apertureRect);
    %Draw fixations
    Screen('FillOval', windowPtr, 0, dotrect);
    Screen('FrameOval', windowPtr, 0, circlerect, 2);    % circle fixation
    % Show it:
    Screen('Flip', windowPtr);
    WaitSecs(time_start);
    
    
    %% ==============================================================
    %                    BEGIN TRIALS
    % ================================================================
    for itrial=trial+1:nTrialsPerRun*run
        
        % =========== 1: fixation cross and circle expand ================
        % write to eyelink
        %         trialStart(itrial)=GetSecs;
        if do_el
            Eyelink('Message','TarX %s', num2str(0));
            Eyelink('Message','TarY %s', num2str(0));
            Eyelink('Message','xDAT %i',1);
        end
        %Draw Aperture
        Screen('FillOval',windowPtr, darkgray, apertureRect);
        %Draw fixation
        Screen('FillOval', windowPtr, black, dotrect);
        Screen('FrameOval', windowPtr, black, [center-fixationsize-5 center+fixationsize+5], 2);    % circle fixation
        
        Screen('Flip', windowPtr);
        WaitSecs(time_fixExpand);
        
        
        % =============== 2: precue ===============
        %write to eyetrck data
        if do_el
            Eyelink('Message','TarX %s', num2str(0));
            Eyelink('Message','TarY %s', num2str(0));
            Eyelink('Message','xDAT %i',2);
        end
        %Draw Aperture
        Screen('FillOval',windowPtr, darkgray, apertureRect);
        %Draw fixation
        drawPrecue(windowPtr,designMat(itrial,1:nItems))
        % Show it:
        Screen('Flip', windowPtr);
        % Wait
        WaitSecs(time_precue);
        
        
        % ============= 3: fixation circle ==============
        %write to eyetrck data
        if do_el
            Eyelink('Message','TarX %s', num2str(0));
            Eyelink('Message','TarY %s', num2str(0));
            Eyelink('Message','xDAT %i',3);
        end
        %Draw Aperture
        Screen('FillOval',windowPtr, darkgray, apertureRect);
        %Draw fixation
        Screen('FillOval', windowPtr, 0, dotrect);
        Screen('FrameOval', windowPtr, 0, circlerect, 2);    % circle fixation
        Screen('Flip', windowPtr);
        % Wait
        WaitSecs(time_interval);
        
        
        % ============== 4: targets appear ================
        %write to eyetrck data
        if do_el
            Eyelink('Message','TarX %s', num2str(0));
            Eyelink('Message','TarY %s', num2str(0));
            Eyelink('Message','xDAT %i',4);
        end
        
        % item info
        items_xy = [designMat(itrial,2*nItems+1+(1:nItems)); designMat(itrial,3*nItems+1+(1:nItems))];
        
        %Draw Aperture
        Screen('FillOval',windowPtr, darkgray, apertureRect);
        %Draw fixation
        Screen('FillOval', windowPtr, 0, dotrect);
        Screen('FrameOval', windowPtr, 0, circlerect, 2);    % circle fixation
        %Draw Dots
        Screen('DrawDots', windowPtr, items_xy, repmat(itemsize,1,nItems), repmat(white,1,nItems), center,1);
        % Show it:
        Screen('Flip', windowPtr);
        % Wait
        WaitSecs(time_present);
        
        
        % ============== 5: delay =============
        %write to eyetrck data
        if do_el
            Eyelink('Message','TarX %s', num2str(0));
            Eyelink('Message','TarY %s', num2str(0));
            Eyelink('Message','xDAT %i',5);
        end
        %Draw Aperture
        Screen('FillOval',windowPtr, darkgray, apertureRect);
        %Draw fixation
        Screen('FillOval', windowPtr, 0, dotrect);
        Screen('FrameOval', windowPtr, 0, circlerect, 2);    % circle fixation
        
        Screen('Flip', windowPtr);
        %         delayStart(itrial)=GetSecs;
        % Delay
        WaitSecs(time_ISI);
        
        
        % =========== 6: test cue appears (response) ==========
        
        % get target info
        targetQuad = find(designMat(itrial,1:nItems)==designMat(itrial,nItems+1));
        if (length(targetQuad)>1) % in cases where there are two targets with same priority
            targetQuad = targetQuad(ceil(rand*length(targetQuad)));
        end
        
        % get target info
        target_xy = [designMat(itrial,2*nItems+1+targetQuad); designMat(itrial,3*nItems+1+targetQuad)];
        
        %write target location to eyetrck data
        if do_el
            Eyelink('Message','TarX %s', num2str(target_xy(1)/ppd));
            Eyelink('Message','TarY %s', num2str(target_xy(2)/ppd));
            Eyelink('Message','xDAT %i',6);
        end
        
        %Draw Aperture
        Screen('FillOval',windowPtr, darkgray, apertureRect);
        %Draw fixation
        Screen('FillOval', windowPtr, 0, dotrect);
        Screen('FrameOval', windowPtr, 0, circlerect, 2);    % circle fixation
        %Quadrant test
        Screen('FrameArc',windowPtr,white, circlerect,90*targetQuad,90,4);
        % Show it:
        Screen('Flip', windowPtr);
        % Wait
        WaitSecs(time_responseCue);
        
        
        % =========== 7: feedback (test target actual location) =========
        %write target location to eyetrck data
        if do_el
            Eyelink('Message','TarX %s', num2str(target_xy(1)/ppd));
            Eyelink('Message','TarY %s', num2str(target_xy(2)/ppd));
            Eyelink('Message','xDAT %i',7);
        end
        
        %Draw Aperture
        Screen('FillOval',windowPtr, darkgray, apertureRect);
        %Draw fixation
        Screen('FillOval', windowPtr, 0, dotrect);
        Screen('FrameOval', windowPtr, 0, circlerect, 2);    % circle fixation
        %Draw Test Dot
        Screen('DrawDots', windowPtr, target_xy, itemsize, white, center,1);
        % Show it:
        Screen('Flip', windowPtr);
        %wait
        WaitSecs(time_feedback);
        
        
        % ============== 8: ITI ================
        %write to eyetrck data
        if do_el
            Eyelink('Message','TarX %s', num2str(0));
            Eyelink('Message','TarY %s', num2str(0));
            Eyelink('Message','xDAT %i',8);
        end
        %Draw Aperture
        Screen('FillOval',windowPtr, darkgray, apertureRect);
        %Draw fixation
        Screen('FillOval', windowPtr, 0, dotrect);
        Screen('FrameOval', windowPtr, 0, circlerect, 2);    % circle fixation
        
        Screen('Flip', windowPtr);
        WaitSecs(time_ITI);
        %         trialEnd(itrial) = GetSecs;
        %         trialTime=trialEnd(itrial)-trialStart(itrial);
        if do_el
            eyeTend=Eyelink('TrackerTime');
            trackerTime=eyeTend-eyeTstrt;
        end
         
        %save output matrix
        save(filename,'designMat','trial')
        trial = trial+1;
    end
    
    % Conclusion Screen
    Screen('FillRect', windowPtr, black);
    Screen('FillOval',windowPtr, darkgray, apertureRect);
    conclText=sprintf('Run %01.f is over.', run);
    DrawFormattedText(windowPtr, conclText,'center', 'center', white);
    Screen('Flip',windowPtr);
    
    % Stop ET Recording
    if do_el
        Eyelink('stoprecording');
        %Eyelink('closefile');
        Eyelink('ReceiveFile',edfFile,pwd,1);
        Eyelink('Shutdown')
    end
    
    %Wait
    KbStrokeWait;
    sca;
    
catch
    
    % This section is executed in case an error happens in the
    % experiment code implemented between try and catch...
    ShowCursor;
    Screen('CloseAll'); % AKA (sca)
    psychrethrow(psychlasterror);
end