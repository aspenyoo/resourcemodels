%% SAUCE 12/7/17
% Not perfect but will run and output what you need
% may crash if run>10 ... still need to check for this 
% OG=old gangster (aka: from old script/old data collected
%% ASTHETICS remove hard coded numbers for trials... make the calculations based off of one variable do not put hard numbers 
% for the trials such as the start and end trials and so forth... this will have to be looked at more when you compute the 
% saving variables needs to be more streamline... maybe change the output matrix to one huge file 1_10runs and include all thins in here
% Do we want to have the esc option during the whole experiment??? it will catch shit I guess {can do and just make sure you stop ET and save ET file}
% could concatenate all trial data into one big file once run 20 is done or just create a new script to run to do so 
% add time marks for all screens and also add a begin experiment time that
%% Lines of Code we need for coomunicating with TMS
% look more into the siwm_OG.m script to figure shit out
% was/is located in /Volumes/data/SIWM_TMS/s01_sess01/ses01/140801_stim01.mat
% % TMS (piece of code we need)
%                 daq=DaqDeviceIndex([],0);
% if DaqDevice Index does not work use 
%                 daq=[2];
%                 err=DaqAOut(daq,0,0);
%                 err2=DaqAOut(daq,0,1);
%                 err=DaqAOut(daq,0,0);
%% Trial Duration
%New Scanner TR= 1.3s {1300 ms}
  %%%% Total Duration for each Run= 280.8s [4.68 min] {216 TRs without the
  %%%% two we wait for after scanner outputs first 5}
  %%% Total TRs to input into scanner == <-218->
%% Screen durations for NEW RT (1300ms) 
% For New Scanner WITH (1.3s TR [1300ms])
        % 0: NO HUGE Waittime after bak tick we wait for two TRs then start recording             
        % 1: Fixation_cross Expands .1s (100ms)           
        % 2: Probability Slice Cues .7s (700ms)     %%% OG .5s (500ms) OG still used for subj1          
        % 3: Fixation-cross         .1s (100ms)
        % 4: Saccade Cues           .7s (700ms)     %%% OG .5s (500ms) OG still used for subj1                 
        % 5: Delay                10.1s (10100ms) 
        % 6: Response Cue           .8s (800ms)     %%% OG .6s (600ms) OG
        % 7: Feedback               .8s (800ms)     %%% OG .6s (600ms) OG
        % 8: ITI          [8.8 10.1 11.4] (8800ms, 10100ms, 11400ms) %%%% OG [9.6, 10.8, 12.0] (9600ms, 10800ms, 12000ms) OG
        %%%%% ???s for 12 trials %%%%%%                    
%%
function pri_4targ_wm_LngTR(subjectID,run,do_el,scanner,TMS) 
%% Screen Preferences 
%Screen('Preference', 'SkipSyncTests', 1);%comment out when actually running a participant 
Screen('Preference', 'VisualDebugLevel', 1); %turns PTB welcome screen black   
AssertOpenGL;
rng('default')
laSemina=sum(100*clock);
rng(laSemina);
if nargin<=2;
    do_el=0; 
    scanner=0;    
    TMS=0;
elseif nargin<=3;
    scanner=0;
    TMS=0;
elseif nargin<=4 && scanner==9;
    TMS=0;    
end
try
    %% Open onscreen window with default settings:
    screenNumber = max(Screen('Screens'));
    topPriorityLevel = MaxPriority(screenNumber); %topPriorityLevel1 = MaxPriority(window);
    Priority(topPriorityLevel); 
    if scanner==1
        [window,screenRect] = Screen('OpenWindow', screenNumber, [0 0 0]);
        screenRect=[0 0 1280 1024];
    else
        [window,screenRect] = Screen('OpenWindow', screenNumber, [0 0 0]);
    end
    Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);
    [xCenter, yCenter] = RectCenter(screenRect);
    [center(1), center(2)] = RectCenter(screenRect); 
    [width, height]=Screen('DisplaySize', screenNumber);    
%    HideCursor; 
    if scanner==1
        % &&& CBI &&& settings are for MRI screen
        width=1024;
        height=768;
        %Prisma
        theDisplaySize = [42.5 33];% [53 33]; %[53.3 35.5]; % height is 33
        distBtwScrAndPart = 63 +9.5; %58;     
    elseif TMS==1
        % &&& TMS Room &&& settings are for behavioral testing rooms
        distBtwScrAndPart = 56; %distance between monitor and participant
        theDisplaySize = [37.5 30]; %height is 30
    else
        % &&& KINGPIN &&& settings are for behavioral testing rooms
        distBtwScrAndPart = 56; %distance between monitor and participant
        theDisplaySize = [37.5 30]; %height is 30        
    end
    %Pixels per Degree
    ppd = pi * (screenRect(3)-screenRect(1)) / atan((theDisplaySize(1))/distBtwScrAndPart/2) / 360; 
%OG     ppd1 = pi * screenRect(3)/ atan(scrnWidth/distBtwScrAndPart/2) / 360; 
%OG     ppd2 =screenRect(4)/ degScrnHeight; % pixels per degree 
    %% define dot/fixation/pizza slices field parameter
    n_dots       = 1;            % number of dots
    max_d       = 10;           % maximum radius of  annulus (in degrees)[boundary for dots] 
    min_d       = 9;            % minumum 
    dot_w       = .2;           % width of dot (deg)  
    apertureSize= 12 ;          % in degrees 
    centerCircle= .6 ;          % size in degrees for center circle % line thickness: .02
    stretchCenterCircle=.6*1.2; %size in degrees
    outrTstWedgThick=.07;       %thickness of test cue
    fxationThick=.07* ppd;      %thickness of fixation
    pizzaSliceThick=.07* ppd;   %thickness of pizza cuts
    degPizzaStart=10;
    degPizzaEnd=20;
    deg1stSmllPizzaStrt=13;
    deg1stSmllPizzaEnd=26;
    deg2ndSmllPizzaStrt=15;
    deg2ndSmllPizzaEnd=30;
    %% Timing 
    if scanner==1
        initlWaitOnStrtRecord = 1.9;
        fixExpnd_time=.1;
        probzCue_time=.7;
        fixIntrvn_time=.1;
        targts_time=.7;
        theDelay_time=10.1;%10.8
        goCue_time=.8;
        feedBack_time=.8;
        iti_times_time=[8.8 10.1 11.4];%[8.8 10 11.2];
    elseif scanner ==9
        initlWaitOnStrtRecord = 1.9;
        fixExpnd_time=.1;
        probzCue_time=.7;
        fixIntrvn_time=.1;
        targts_time=.7;
        theDelay_time=10.1;%10.8
        goCue_time=.8;
        feedBack_time=.8;
        iti_times_time=[8.8 10.1 11.4];
    else 
        initlWaitOnStrtRecord = 1.9;
        fixExpnd_time=.1;
        probzCue_time=.7;
        fixIntrvn_time=.1;
        targts_time=.7;
        theDelay_time=3;%10.8
        goCue_time=.8;
        feedBack_time=.8;
        iti_times_time=[3 4.5 6];%[8.8 10 11.2];     
    end
    %% define keys to listen to 
    KbName('UnifyKeyNames');
    if ismac == 1 %IsOSX
        ESCa_key = KbName('escape'); % press this key to abort
    else
        %ESCa_key = KbName('esc');
        ESCa_key = KbName('escape');
    end
    STRT_key = [KbName('5%') KbName('5')];  %[resp=34 on mac] should be lower-case t at prisma? (or %5, which is top-row 5, or 5, which is numpad 5)
    SPC_key = KbName('space');
    %% define colors
    white = WhiteIndex(screenNumber) ;
    black = BlackIndex(screenNumber);
    semiBlack= black+85;
    darkgray = white/2.2;  
       lightgray = white/1.4;
    %% define max/min dot radius present
    dot_s = dot_w * ppd;
    rmax = max_d * ppd;	% maximum radius of annulus (pixels from center)
    rmin = min_d * ppd; % minimum
    radius_to_present = rmax * sqrt(rand(n_dots,1));
    radius_to_present(radius_to_present<rmin) = rmin;
    %% define positions and angles
    positionOfMainCircle = [(xCenter-(centerCircle*ppd)) (yCenter-(centerCircle*ppd)) (xCenter+(centerCircle*ppd)) ...
        (yCenter+(centerCircle*ppd))] ; 
    positionOfExpandedCircle = [(xCenter-(stretchCenterCircle*ppd)) (yCenter-(stretchCenterCircle*ppd)) ...
        (xCenter+(stretchCenterCircle*ppd)) (yCenter+(stretchCenterCircle*ppd))] ; 
    positionOfCircleForOuterWedge= [(xCenter-((centerCircle+outrTstWedgThick)*ppd)) ...
        (yCenter-((centerCircle+outrTstWedgThick)*ppd)) (xCenter+((centerCircle+outrTstWedgThick)*ppd)) ...
        (yCenter+((centerCircle+outrTstWedgThick)*ppd))]; 
    theAperture = [(xCenter-(apertureSize*ppd)) (yCenter-(apertureSize*ppd)) (xCenter+(apertureSize*ppd)) ...
          (yCenter+(apertureSize*ppd))] ;
    %start of draw and end of draw
    sizeAngle  = [90 90 90 90 360] ;%%%%%%%%%%%%%%%%%%%%%%% because of how FrameArc & FillArc read angles our starting point
    startAngle = [90 180 270 360] ; %%%%%% WARNING %%%%%%%% will be 90(90==0 degree starting point and increases clockwise)
    %% creating wifi arc parameters
    wifi_arc_1=.17;
    wifi_arc_2=.27;
    wifi_arc_3=.38;
    wifi_arc_4=.49;
    wifi_InnerArcPosition= [(xCenter-(wifi_arc_1*ppd)) (yCenter-(wifi_arc_1*ppd)) (xCenter+(wifi_arc_1*ppd)) (yCenter+(wifi_arc_1*ppd))] ;
    wifi_secndArcPosition=[(xCenter-(wifi_arc_2*ppd)) (yCenter-(wifi_arc_2*ppd)) (xCenter+(wifi_arc_2*ppd)) (yCenter+(wifi_arc_2*ppd))] ;
    wifi_thirdArcPosition=[(xCenter-(wifi_arc_3*ppd)) (yCenter-(wifi_arc_3*ppd)) (xCenter+(wifi_arc_3*ppd)) (yCenter+(wifi_arc_3*ppd))] ;
    wifi_OuterArcPosition=[(xCenter-(wifi_arc_4*ppd)) (yCenter-(wifi_arc_4*ppd)) (xCenter+(wifi_arc_4*ppd)) (yCenter+(wifi_arc_4*ppd))] ;  
    %% number of trials %%  
    numbTrialsPerRun = 12; % could this be tweeked with in order for it to teach itself what trial number is best ???
    outputMatrix=zeros(numbTrialsPerRun,30);% Allocate Space for Output variables
    totTrials=120;% total numbr of trial to aim for in one scan session  
    expStart=zeros(1,1);
    trialStart=zeros(numbTrialsPerRun,1);
    delayStart=zeros(numbTrialsPerRun,1);
    trialEnd=zeros(numbTrialsPerRun,1);   
    %% what needs to be done every 10th run     
    if run==1 || run==11 || run==21 || run==31
        %% Creating .csv file
        if 1==numel(num2str(run))
            dataRunFile=sprintf('P%d_runs0%g_%g_ABprobz.csv', subjectID,run,(run+9));
        else
            dataRunFile=sprintf('P%d_runs%g_%g_ABprobz.csv', subjectID,run,(run+9));
        end    
        fid = fopen(dataRunFile, 'a'); %open file stream
        fprintf(fid, ...
            'Subject\tRunNum\tTrialNum\tTrialDur\tTrialITI\tQuadTested\tProbzTested\tLocTested_x\tLocTested_y\tAngleTested\tProbzQuad1\tLocDotQuad1_x\tLocDotQuad1_y\tAngleDot1\tProbzQuad2\tLocDotquad2_x\tLocDotQuad2_y\tAngleDot2\tProbzQuad3\tLocDotQuad3_x\tLocDotQuad3_y\tAngleDot3\tProbzQuad4\tLocDotQuad4_x\tLocDotQuad4_y\tAngleDot4\t%s\n',...
            datestr(now));
        %% randomizing probability
        % is this what is presented at each quadrant????? the probability
        % for each of the queadrants???? need to test and run this
        % SAUCE 12/7/17: I do not think so since we change the index of
        % these to account for the quadrant list to be tested
        allPrbz_prbz_quds_iti_anAngsTst=zeros(totTrials,4,5);
        v=[0.6 0.3 0.1 0];
        q=unique(perms(v),'rows'); % all possible permutations
        matrix1=zeros(totTrials, 4);
        %%%%% O why the fuck have a 201 matrix with a repetition of all of
        %%%%% these numbers when we just want to sample orm he original q
        %%%%% matrix
        matrix1(1:120,:)=repmat(q(19:24,:),20,1);
        matrix1(121:180,:)=repmat(q(13:18,:),10,1);
        matrix1(181:204,:)=repmat(q(7:12,:),4,1); % changed output from 108 to 100
        matrix1(201:end,:)=[];
        allPrbz_prbz_quds_iti_anAngsTst(:,:,1)=datasample(matrix1,totTrials);
        %% list of probability being tested      
        a=[0.6 0.3 0.1];
        repA=[repmat((a(1)),(120*.6),1);repmat((a(2)),(120*.3),1);repmat((a(3)),(120*.1),1)];
        theTestProbzShuff=Shuffle(repA);
        % Getting Rid of consecutive 10% probz
        while sum(diff(find(theTestProbzShuff==.1))==1)==1 %&& sum(diff(find(theTestProbzShuff==.1))==2)==1
            theTestProbzShuff=Shuffle(theTestProbzShuff);
        end
        theTestProbzShuff=[theTestProbzShuff,zeros(totTrials,3)];
        %% list of quadrants being test
        repQuads=repmat((1:4)',(totTrials/4),1);
        % save all into matrix
        allPrbz_prbz_quds_iti_anAngsTst(:,:,2)=theTestProbzShuff;
        allPrbz_prbz_quds_iti_anAngsTst(:,:,3)= [Shuffle(repQuads),zeros(totTrials,3)];
        %% list of ITIs
        itiDelay=zeros(totTrials,1);
        for yy=1:(totTrials/numbTrialsPerRun);
            if yy==1
                corrIndx1=1;
            else
                corrIndx1=((yy-1).*numbTrialsPerRun)+1;
            end
            corrIndx2=yy.*numbTrialsPerRun;
            itiDelay(corrIndx1:corrIndx2,1)=Shuffle([repmat(iti_times_time(1),3,1);repmat(iti_times_time(2),6,1);repmat(iti_times_time(3),3,1)]);   
        end
        clear('yy')
        allPrbz_prbz_quds_iti_anAngsTst(:,:,4)= [itiDelay,zeros(totTrials,3)]; 
        %% list of angles for dot locations
        angleList1=Shuffle([(repmat([10:10:80],1,15))']);%;((datasample([10:10:80],4)'))]);
        angleList2=Shuffle([(repmat([100:10:170],1,15))']);%;((datasample([100:10:170],4)'))]);
        angleList3=Shuffle([(repmat([190:10:260],1,15))']);%;((datasample([190:10:260],4)'))]);
        angleList4=Shuffle([(repmat([280:10:350],1,15))']);%;((datasample([280:10:350],4)'))]);
        allAngPerQuad=[angleList1,angleList2,angleList3,angleList4];
        allPrbz_prbz_quds_iti_anAngsTst(:,:,5)=allAngPerQuad;
        if 1==numel(num2str(run))
            datafile =sprintf('Prbz_qudsTst_iti_angsLst_P%g_runs0%g_%g.mat',subjectID,run,(run+9));
        else
            datafile =sprintf('Prbz_qudsTst_iti_angsLst_P%g_runs%g_%g.mat',subjectID,run,(run+9));
        end
        save(datafile,'allPrbz_prbz_quds_iti_anAngsTst');
        % dataRunFile=sprintf('P%d_runs%02.f_%g_ABprobz.csv', subjectID,updateRun,(updateRun+9))
        %Turn degrees to pixels        
        allCords=zeros(totTrials,8);
        for numOfAngList=1:min(size(allAngPerQuad))
            for rowsOfAng=1:max(size(allAngPerQuad))
                indexToAlloc=[1,3,5,7];
                strtColm=indexToAlloc(numOfAngList);
                allCords(rowsOfAng,strtColm)=(cosd(allAngPerQuad(rowsOfAng,numOfAngList))).*radius_to_present;
                allCords(rowsOfAng,(strtColm+1))=(sind(allAngPerQuad(rowsOfAng,numOfAngList))).*radius_to_present;
            end
        end
        pixLocOfdots=zeros(2,totTrials,4);
        pixLocOfdots(:,:,1)=transpose(allCords(:,1:2));
        pixLocOfdots(:,:,2)=transpose(allCords(:,3:4));
        pixLocOfdots(:,:,3)=transpose(allCords(:,5:6));
        pixLocOfdots(:,:,4)=transpose(allCords(:,7:8));
        if 1==numel(num2str(run))
            datafile2 =sprintf('pixLocOfdots_P%g_runs0%g_%g.mat',subjectID,run,(run+9));
        else
            datafile2 =sprintf('pixLocOfdots_P%g_runs%g_%g.mat',subjectID,run,(run+9));
        end
        % dataRunFile=sprintf('P%d_runs%02.f_%g_ABprobz.csv', subjectID,updateRun,(updateRun+9))
        save(datafile2,'pixLocOfdots');
        startTrial=1;
        numTrials=numbTrialsPerRun;      
    elseif 2<=run && run<=10 || 12<=run && run<=20 || 22<=run && run<=30 || 32<=run && run<=40
        if 2<=run && run<=10
            updatedRun=1;
            preRun=run;
        elseif 12<=run && run<=20
            updatedRun=11;
            preRun=run-10;
        elseif 22<=run && run<=30
            updatedRun=21;
            preRun=run-20;
        elseif 32<=run && run<=40
            updatedRun=31;
            preRun=run-30;
        end
        % Change what you load and then also load the PIX and location
%         if 1==numel(num2str(updatedRun))
%             dataRunFile=sprintf('P%d_runs0%g_%g_ABprobz.csv', subjectID,updatedRun,(updatedRun+9));
%         else 
%             dataRunFile=sprintf('P%d_runs%g_%g_ABprobz.csv', subjectID,updatedRun,(updatedRun+9));
%         end
        dataRunFile=sprintf('P%d_runs%02.f_%g_ABprobz.csv', subjectID,updatedRun,(updatedRun+9));
        fid = fopen(dataRunFile, 'a'); %open file stream
        fprintf(fid, ...
            'Subject\tRunNum\tTrialNum\tTrialDur\tTrialITI\tQuadTested\tProbzTested\tLocTested_x\tLocTested_y\tAngleTested\tProbzQuad1\tLocDotQuad1_x\tLocDotQuad1_y\tAngleDot1\tProbzQuad2\tLocDotquad2_x\tLocDotQuad2_y\tAngleDot2\tProbzQuad3\tLocDotQuad3_x\tLocDotQuad3_y\tAngleDot3\tProbzQuad4\tLocDotQuad4_x\tLocDotQuad4_y\tAngleDot4\t%s\n',...
            datestr(now));
        %%% this may have to be changed to account for the preRun variabel
        %%% instead of the run variabel in the if then statement
        % SAUCE 12/7/17: no it does not because it is just a naming convention 
%         if 1==numel(num2str(run))
%             matrxFile1=sprintf('Prbz_qudsTst_iti_angsLst_P%g_runs0%g_%g.mat', subjectID, updatedRun, (updatedRun+9));    
%         elseif run==10
%             matrxFile1=sprintf('Prbz_qudsTst_iti_angsLst_P%g_runs0%g_%g.mat', subjectID, updatedRun, (updatedRun+9));
%         else
%             matrxFile1=sprintf('Prbz_qudsTst_iti_angsLst_P%g_runs%g_%g.mat', subjectID, updatedRun, (updatedRun+9)); 
%         end
        % you can use this below to get rid of the if statements... 
        matrxFile1=sprintf('Prbz_qudsTst_iti_angsLst_P%g_runs%02.f_%g.mat', subjectID, updatedRun, (updatedRun+9));
        load (matrxFile1)
%         if 1==numel(num2str(run))
%             matrxFile2=sprintf('pixLocOfdots_P%g_runs0%g_%g.mat', subjectID, updatedRun, (updatedRun+9));
%         elseif run==10
%             matrxFile2=sprintf('pixLocOfdots_P%g_runs0%g_%g.mat', subjectID, updatedRun, (updatedRun+9));
%         else 
%             matrxFile2=sprintf('pixLocOfdots_P%g_runs%g_%g.mat', subjectID, updatedRun, (updatedRun+9));   
%         end
        matrxFile2=sprintf('pixLocOfdots_P%g_runs%02.f_%g.mat', subjectID, updatedRun, (updatedRun+9));   
        load (matrxFile2)
        startTrial=(((preRun-1).*numbTrialsPerRun)+1);
        numTrials=startTrial+11;
    end
    %% Flip the window before start to get a gray screen 
    Screen('FillRect', window, darkgray);
    %% Get EyeTracker and Provide info 
    if do_el
        if scanner==1
            Eyelink('SetAddress','192.168.1.5')
        end
        el=EyelinkInitDefaults(window);
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
        if 1==numel(num2str(run))
            edfFile=sprintf('ED_%g_0%g.edf', subjectID, run);
        else 
            edfFile=sprintf('ED_%g_%g.edf', subjectID, run);
        end
        Eyelink('openfile',edfFile);
    end
    %% Saving Screen/Stimulus/OtherShit parameters 
    if 1==numel(num2str(run))
        varsBefStim=(['all_Variables_bef_stim_Loop_P',num2str(subjectID),'_r0',num2str(run),'.mat']);
    else
        varsBefStim=(['all_Variables_bef_stim_Loop_P',num2str(subjectID),'_r',num2str(run),'.mat']);
    end
    save (varsBefStim)
    
%% #########################################################################%
%% ######################### start the loop of trials ######################%
%% #########################################################################%

for trial=startTrial:numTrials  
      trialNum=sprintf('$ Trial %s $', num2str(trial));
      disp(trialNum)
      %% Get Probabilities
      rowInd = allPrbz_prbz_quds_iti_anAngsTst(trial,:,1);% this is where you get the probabilities  
      %% Gettingcorrect Probability to test from matrixOfProbz
      probTested=allPrbz_prbz_quds_iti_anAngsTst(trial,1,2);
      indxTofindAng= find(rowInd==probTested);
      %Swaping index to match quadrant being tested 
      rowInd([indxTofindAng allPrbz_prbz_quds_iti_anAngsTst(trial,1,3)])=rowInd([allPrbz_prbz_quds_iti_anAngsTst(trial,1,3) indxTofindAng]);
      newArangInd=rowInd;
      if allPrbz_prbz_quds_iti_anAngsTst(trial,1,3)==1
          quadTested=90;
      elseif allPrbz_prbz_quds_iti_anAngsTst(trial,1,3)==2
          quadTested=180;
      elseif allPrbz_prbz_quds_iti_anAngsTst(trial,1,3)==3
          quadTested=270;
      elseif allPrbz_prbz_quds_iti_anAngsTst(trial,1,3)==4
          quadTested=360;
      end
      %% Test Target Location
      if quadTested == 90
          testDot= pixLocOfdots(:,trial,1);
      elseif quadTested == 180
          testDot= pixLocOfdots(:,trial,2);
      elseif quadTested == 270
          testDot= pixLocOfdots(:,trial,3);
      elseif quadTested == 360
          testDot= pixLocOfdots(:,trial,4);
      end           
      %get target location
      tarX=(testDot(1))./ppd; %want in degrees of visual angle   
      tarY=(testDot(2))./ppd;
      %% Intro Screen
      Screen('FillRect', window, black);
      Screen('TextFont',window, 'Helvetica');
      Screen('TextSize',window, 20);
      Screen('TextStyle', window, 1+2);
      if trial == startTrial
          Screen('FillOval',window, darkgray, theAperture);
          runText=sprintf('This is Run %01.f \n', run);
          if scanner==1
              DrawFormattedText(window, [runText, 'Waiting for Scanner'],'center', 'center', white);
              resp = 0;
              Screen('Flip',window);
              while resp == 0
                  [resp, ~] = checkForResp(STRT_key, ESCa_key);
                  if resp == -1
                      Screen('CloseAll'); ShowCursor;
                      if do_el
                          Eyelink('ShutDown');
                      end
                      return;
                  end
              end
              clear resp;
              expStart(trial)=GetSecs; % This is the Start time of the trial % WILL HAVE TO SAVE
              %After back tick we will start recording ET and draw first fixation cross
              if do_el
                  eyeTstrt=Eyelink('TrackerTime');
                  Eyelink('startrecording');
                  Eyelink('Message','xDAT %i', 10);
              end
              %Draw Aperture
              Screen('FillOval',window, darkgray, theAperture);
              %Draw Circle
              Screen('FrameArc',window, lightgray,positionOfMainCircle,startAngle(1),sizeAngle(5),(fxationThick), (fxationThick));
              %Draw fixation cross
              Screen('DrawLine', window, lightgray, xCenter, positionOfMainCircle(2), xCenter, positionOfMainCircle(4), (fxationThick));
              Screen('DrawLine', window, lightgray, positionOfMainCircle(1), yCenter, positionOfMainCircle(3), yCenter, fxationThick);
              %Draw Center Dart Board
              Screen('FrameArc',window, lightgray, wifi_InnerArcPosition, startAngle(1), sizeAngle(5), (fxationThick), (fxationThick));
              for drawDartBrd=1:length(startAngle);
                  Screen('FrameArc',window, lightgray, wifi_secndArcPosition, (startAngle(drawDartBrd)+deg2ndSmllPizzaStrt), (sizeAngle(1)-deg2ndSmllPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
                  Screen('FrameArc',window, lightgray, wifi_thirdArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
                  Screen('FrameArc',window, lightgray, wifi_OuterArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
              end
              % Show it:
              Screen('Flip', window);
              %Wait for 1.8 seconds to pass first 5 being spit out by
              %scanner because we want to skip first two TRs
              resp = 0;
              while (GetSecs-expStart(trial)) < initlWaitOnStrtRecord
                  [resp, ~] = checkForResp([], ESCa_key);
                  if resp == -1
                      sca;
                      fprintf('ESC pressed during pre-trial wait time\n');
                      ShowCursor;
                      return;
                  end
              end
              clear resp;
              resp = 0;
              % wait till next 5 is spit out that way we wait a total of
              % 2TR's before experiment starts
              while resp == 0
                  [resp, ~] = checkForResp(STRT_key, ESCa_key);
                  if resp == -1
                      Screen('CloseAll'); ShowCursor;
                      if do_el
                          Eyelink('ShutDown');
                      end
                      return;
                  end
              end
              clear resp;
          else
              %DrawFormattedText(window, conclText,'center', 'center', white);
              DrawFormattedText(window, [runText, 'Press the spacebar when you are ready to continue.'],'center', 'center', white);
              % Show it
              Screen('Flip',window);
              %Wait for keyboard press:
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
              expStart(trial)=GetSecs;
              if do_el
                  eyeTstrt=Eyelink('TrackerTime');
                  Eyelink('startrecording');
                  Eyelink('Message','xDAT %i', 10);
              end
              %Draw Aperture
              Screen('FillOval',window, darkgray, theAperture);
              %Draw Circle
              Screen('FrameArc',window, lightgray,positionOfMainCircle,startAngle(1),sizeAngle(5),(fxationThick), (fxationThick));
              %Draw fixation cross
              Screen('DrawLine', window, lightgray, xCenter, positionOfMainCircle(2), xCenter, positionOfMainCircle(4), (fxationThick));
              Screen('DrawLine', window, lightgray, positionOfMainCircle(1), yCenter, positionOfMainCircle(3), yCenter, fxationThick);
              %Draw Center Dart Board
              Screen('FrameArc',window, lightgray, wifi_InnerArcPosition, startAngle(1), sizeAngle(5), (fxationThick), (fxationThick));
              for drawDartBrd=1:length(startAngle);
                  Screen('FrameArc',window, lightgray, wifi_secndArcPosition, (startAngle(drawDartBrd)+deg2ndSmllPizzaStrt), (sizeAngle(1)-deg2ndSmllPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
                  Screen('FrameArc',window, lightgray, wifi_thirdArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
                  Screen('FrameArc',window, lightgray, wifi_OuterArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
              end
              % Show it:
              Screen('Flip', window);
              % Wait
              resp = 0;
              while (GetSecs-expStart(trial)) < initlWaitOnStrtRecord
                  [resp, ~] = checkForResp([], ESCa_key);
                  if resp == -1
                      sca;
                      fprintf('ESC pressed during pre-trial wait time\n');
                      ShowCursor;
                      return;
                  end
              end
              clear resp;
          end
      end 
 %$$$$$$ Stoped here with ESC  key implementation $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$     
      %% 1st Screen: Fixation cross and circle expand
      %write to eyetrck data     
      trialStart(trial)=GetSecs;      
      if do_el
          Eyelink('Message','TarX %s', num2str(0));
          Eyelink('Message','TarY %s', num2str(0));
          Eyelink('Message','xDAT %i',1);
      end
      %Draw Aperture
      Screen('FillOval',window, darkgray, theAperture);
      %Draw Circle
      Screen('FrameArc',window, lightgray,positionOfExpandedCircle,startAngle(1),sizeAngle(5),(fxationThick), (fxationThick));
      %Draw fixation cross
      Screen('DrawLine', window, lightgray, xCenter, positionOfExpandedCircle(2), xCenter, positionOfExpandedCircle(4), (fxationThick));
      Screen('DrawLine', window, lightgray, positionOfExpandedCircle(1), yCenter, positionOfExpandedCircle(3), yCenter, fxationThick);
      %Draw Center Dart Board
      Screen('FrameArc',window, lightgray, wifi_InnerArcPosition, startAngle(1), sizeAngle(5), (fxationThick), (fxationThick));
      for drawDartBrd=1:length(startAngle)
          Screen('FrameArc',window, lightgray, wifi_secndArcPosition, (startAngle(drawDartBrd)+deg2ndSmllPizzaStrt), (sizeAngle(1)-deg2ndSmllPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_thirdArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_OuterArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
      end
      % Show it:
      Screen('Flip', window);
      % Wait              
%       resp = 0;
%       while (GetSecs-trialStrt) < initlWaitOnStrtRecord
%           [resp, ~] = checkForResp([], ESCa_key);
%           if resp == -1
%               Screen('CloseAll'); ShowCursor;
%               Eyelink('StopRecording');
%               Eyelink('ShutDown');
%               return;
%           end
%       end
%       clear resp;    
      WaitSecs(fixExpnd_time);
      %% 2nd Screen: Fixation cross and probability wedges 
      %write to eyetrck data
      if do_el
          Eyelink('Message','TarX %s', num2str(0));
          Eyelink('Message','TarY %s', num2str(0));
          Eyelink('Message','xDAT %i',2);
      end
      %Draw Aperture
      Screen('FillOval',window, darkgray, theAperture);
      %Draw Circle
      Screen('FrameArc',window, lightgray,positionOfMainCircle,startAngle(1),sizeAngle(5),(fxationThick), (fxationThick));
      %Draw fixation cross
      Screen('DrawLine', window, lightgray, xCenter, positionOfMainCircle(2), xCenter, positionOfMainCircle(4), (fxationThick));
      Screen('DrawLine', window, lightgray, positionOfMainCircle(1), yCenter, positionOfMainCircle(3), yCenter, fxationThick);
      %Draw Center Dart Board
      Screen('FrameArc',window, semiBlack, wifi_InnerArcPosition, startAngle(1), sizeAngle(5), (fxationThick), (fxationThick));
      for drawDartBrd=1:length(startAngle)
          Screen('FrameArc',window, semiBlack, wifi_secndArcPosition, (startAngle(drawDartBrd)+deg2ndSmllPizzaStrt), (sizeAngle(1)-deg2ndSmllPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, semiBlack, wifi_thirdArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, semiBlack, wifi_OuterArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
      end
      for makeSlices=1:length(newArangInd)
          if newArangInd(makeSlices)==0 
          elseif newArangInd(makeSlices)==.10
          Screen('FrameArc',window, white, wifi_InnerArcPosition, (startAngle(makeSlices)+deg1stSmllPizzaStrt), (sizeAngle(1)-deg1stSmllPizzaEnd), (fxationThick), (fxationThick));
          Screen('FrameArc',window, white, wifi_secndArcPosition, (startAngle(makeSlices)+deg2ndSmllPizzaStrt), (sizeAngle(1)-deg2ndSmllPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));  
          elseif newArangInd(makeSlices)==.30
          Screen('FrameArc',window, white, wifi_InnerArcPosition, (startAngle(makeSlices)+deg1stSmllPizzaStrt), (sizeAngle(1)-deg1stSmllPizzaEnd), (fxationThick), (fxationThick));
          Screen('FrameArc',window, white, wifi_secndArcPosition, (startAngle(makeSlices)+deg2ndSmllPizzaStrt), (sizeAngle(1)-deg2ndSmllPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, white, wifi_thirdArcPosition, (startAngle(makeSlices)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          elseif newArangInd(makeSlices)==.60
          Screen('FrameArc',window, white, wifi_InnerArcPosition, (startAngle(makeSlices)+deg1stSmllPizzaStrt), (sizeAngle(1)-deg1stSmllPizzaEnd), (fxationThick), (fxationThick));
          Screen('FrameArc',window, white, wifi_secndArcPosition, (startAngle(makeSlices)+deg2ndSmllPizzaStrt), (sizeAngle(1)-deg2ndSmllPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, white, wifi_thirdArcPosition, (startAngle(makeSlices)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, white, wifi_OuterArcPosition, (startAngle(makeSlices)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          end
      end      
      % Show it:
      Screen('Flip', window);
      % Wait
      WaitSecs(probzCue_time); 
      %% 3rd Screen: fixation cross 
      %write to eyetrck data
      if do_el
          Eyelink('Message','TarX %s', num2str(0));
          Eyelink('Message','TarY %s', num2str(0));
          Eyelink('Message','xDAT %i',3);
      end
      %Draw Aperture
      Screen('FillOval',window, darkgray, theAperture);
      %Draw Circle
      Screen('FrameArc',window, lightgray,positionOfMainCircle,startAngle(1),sizeAngle(5),(fxationThick), (fxationThick));
      %Draw fixation cross
      Screen('DrawLine', window, lightgray, xCenter, positionOfMainCircle(2), xCenter, positionOfMainCircle(4), (fxationThick));
      Screen('DrawLine', window, lightgray, positionOfMainCircle(1), yCenter, positionOfMainCircle(3), yCenter, fxationThick);
      %Draw Center Dart Board
      Screen('FrameArc',window, lightgray, wifi_InnerArcPosition, startAngle(1), sizeAngle(5), (fxationThick), (fxationThick));
      for drawDartBrd=1:length(startAngle)
          Screen('FrameArc',window, lightgray, wifi_secndArcPosition, (startAngle(drawDartBrd)+deg2ndSmllPizzaStrt), (sizeAngle(1)-deg2ndSmllPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_thirdArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_OuterArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
      end      % Show it:
      Screen('Flip', window);
      % Wait
      WaitSecs(fixIntrvn_time);
      %% 4th Screen: fixation cross and targets
      %write to eyetrck data
      if do_el
          Eyelink('Message','TarX %s', num2str(0));
          Eyelink('Message','TarY %s', num2str(0));
          Eyelink('Message','xDAT %i',4);
      end
      %Draw Aperture
      Screen('FillOval',window, darkgray, theAperture);
      %Draw Circle
      Screen('FrameArc',window, lightgray,positionOfMainCircle,startAngle(1),sizeAngle(5),(fxationThick), (fxationThick));
      %Draw fixation cross
      Screen('DrawLine', window, lightgray, xCenter, positionOfMainCircle(2), xCenter, positionOfMainCircle(4), (fxationThick));
      Screen('DrawLine', window, lightgray, positionOfMainCircle(1), yCenter, positionOfMainCircle(3), yCenter, fxationThick);
      %Draw Center Dart Board
      Screen('FrameArc',window, lightgray, wifi_InnerArcPosition, startAngle(1), sizeAngle(5), (fxationThick), (fxationThick));
      for drawDartBrd=1:length(startAngle)
          Screen('FrameArc',window, lightgray, wifi_secndArcPosition, (startAngle(drawDartBrd)+deg2ndSmllPizzaStrt), (sizeAngle(1)-deg2ndSmllPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_thirdArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_OuterArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
      end
      %Draw Dots
      Screen('DrawDots', window, pixLocOfdots(:,trial,1), dot_s, white, center,1);
      Screen('DrawDots', window, pixLocOfdots(:,trial,2), dot_s, white, center,1);
      Screen('DrawDots', window, pixLocOfdots(:,trial,3), dot_s, white, center,1);
      Screen('DrawDots', window, pixLocOfdots(:,trial,4), dot_s, white, center,1);
      % Show it:
      Screen('Flip', window);
      % Wait
      WaitSecs(targts_time);
      %% 5th Screen: Delay 
      %write to eyetrck data
      if do_el
          Eyelink('Message','TarX %s', num2str(0));
          Eyelink('Message','TarY %s', num2str(0));
          Eyelink('Message','xDAT %i',5);
      end
      %Draw Aperture
      Screen('FillOval',window, darkgray, theAperture);
      %Draw Circle
      Screen('FrameArc',window, lightgray,positionOfMainCircle,startAngle(1),sizeAngle(5),(fxationThick), (fxationThick));
      %Draw fixation cross
      Screen('DrawLine', window, lightgray, xCenter, positionOfMainCircle(2), xCenter, positionOfMainCircle(4), (fxationThick));
      Screen('DrawLine', window, lightgray, positionOfMainCircle(1), yCenter, positionOfMainCircle(3), yCenter, fxationThick);
      %Draw Center Dart Board
      Screen('FrameArc',window, lightgray, wifi_InnerArcPosition, startAngle(1), sizeAngle(5), (fxationThick), (fxationThick));
      for drawDartBrd=1:length(startAngle)
          Screen('FrameArc',window, lightgray, wifi_secndArcPosition, (startAngle(drawDartBrd)+deg2ndSmllPizzaStrt), (sizeAngle(1)-deg2ndSmllPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_thirdArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_OuterArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
      end
      % Show it:
      Screen('Flip', window);
      delayStart(trial)=GetSecs;
      % Delay
      WaitSecs(theDelay_time);     
      %% 6th Screen: Test Cue (response)
      %write target location to eyetrck data
      if do_el
          Eyelink('Message','TarX %s', num2str(tarX));
          Eyelink('Message','TarY %s', num2str(tarY));
          Eyelink('Message','xDAT %i',6);
      end
      %Draw Aperture
      Screen('FillOval',window, darkgray, theAperture);
      %Draw Circle
      Screen('FrameArc',window, lightgray,positionOfMainCircle,startAngle(1),sizeAngle(5),(fxationThick), (fxationThick));
      %Draw fixation cross
      Screen('DrawLine', window, lightgray, xCenter, positionOfMainCircle(2), xCenter, positionOfMainCircle(4), (fxationThick));
      Screen('DrawLine', window, lightgray, positionOfMainCircle(1), yCenter, positionOfMainCircle(3), yCenter, fxationThick);
      %Draw Center Dart Board
      Screen('FrameArc',window, lightgray, wifi_InnerArcPosition, startAngle(1), sizeAngle(5), (fxationThick), (fxationThick));
      for drawDartBrd=1:length(startAngle)
          Screen('FrameArc',window, lightgray, wifi_secndArcPosition, (startAngle(drawDartBrd)+deg2ndSmllPizzaStrt), (sizeAngle(1)-deg2ndSmllPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_thirdArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_OuterArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
      end
      %Quadrant test   
      Screen('FrameArc',window, white,positionOfCircleForOuterWedge,quadTested,sizeAngle(1),(outrTstWedgThick*ppd),(outrTstWedgThick*ppd));
      % Show it:
      Screen('Flip', window); 
      % Wait 
      WaitSecs(goCue_time);
      %% 7th Screen: Feedback (test target actual location)
      %write target location to eyetrck data
      if do_el
          Eyelink('Message','TarX %s', num2str(tarX));
          Eyelink('Message','TarY %s', num2str(tarY));
          Eyelink('Message','xDAT %i',7);
      end
      %Draw Aperture
      Screen('FillOval',window, darkgray, theAperture);
      %Draw Circle
      Screen('FrameArc',window, lightgray,positionOfMainCircle,startAngle(1),sizeAngle(5),(fxationThick), (fxationThick));
      %Draw fixation cross
      Screen('DrawLine', window, lightgray, xCenter, positionOfMainCircle(2), xCenter, positionOfMainCircle(4), (fxationThick));
      Screen('DrawLine', window, lightgray, positionOfMainCircle(1), yCenter, positionOfMainCircle(3), yCenter, fxationThick);
      %Draw Center Dart Board
      Screen('FrameArc',window, lightgray, wifi_InnerArcPosition, startAngle(1), sizeAngle(5), (fxationThick), (fxationThick));
      for drawDartBrd=1:length(startAngle)
          Screen('FrameArc',window, lightgray, wifi_secndArcPosition, (startAngle(drawDartBrd)+deg2ndSmllPizzaStrt), (sizeAngle(1)-deg2ndSmllPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_thirdArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_OuterArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
      end
      %Draw Test Dot
      Screen('DrawDots', window, testDot, dot_s, white, center,1);
      % Show it:
      Screen('Flip', window);
      %wait
      WaitSecs(feedBack_time);
      %% 8th ITI
      %write to eyetrck data
      if do_el
          Eyelink('Message','TarX %s', num2str(0));
          Eyelink('Message','TarY %s', num2str(0));
          Eyelink('Message','xDAT %i',8);
      end
      %Draw Aperture
      Screen('FillOval',window, darkgray, theAperture);
      %Draw Circle
      Screen('FrameArc',window, lightgray,positionOfMainCircle,startAngle(1),sizeAngle(5),(fxationThick), (fxationThick));
      %Draw fixation cross
      Screen('DrawLine', window, lightgray, xCenter, positionOfMainCircle(2), xCenter, positionOfMainCircle(4), (fxationThick));
      Screen('DrawLine', window, lightgray, positionOfMainCircle(1), yCenter, positionOfMainCircle(3), yCenter, fxationThick);
      %Draw Center Dart Board
      Screen('FrameArc',window, lightgray, wifi_InnerArcPosition, startAngle(1), sizeAngle(5), (fxationThick), (fxationThick));
      for drawDartBrd=1:length(startAngle)
          Screen('FrameArc',window, lightgray, wifi_secndArcPosition, (startAngle(drawDartBrd)+deg2ndSmllPizzaStrt), (sizeAngle(1)-deg2ndSmllPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_thirdArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
          Screen('FrameArc',window, lightgray, wifi_OuterArcPosition, (startAngle(drawDartBrd)+degPizzaStart), (sizeAngle(1)-degPizzaEnd), (pizzaSliceThick), (pizzaSliceThick));
      end
      % Show it:
      Screen('Flip', window);
      % Wait
      WaitSecs(allPrbz_prbz_quds_iti_anAngsTst(trial,1,4));
      trialEnd(trial)=GetSecs;
      trialTime=trialEnd(trial)-trialStart(trial);
      if do_el
          eyeTend=Eyelink('TrackerTime');
          trackerTime=eyeTend-eyeTstrt;
      end
%% Save the data for Trial #  
      %location of tested dot in degrees
      outputMatrix(trial,1)=tarX;
      outputMatrix(trial,2)=tarY;
      %location of tested dot in Pixels
      outputMatrix(trial,3)=testDot(1);
      outputMatrix(trial,4)=testDot(2);
      % angle of tested dot
      if allPrbz_prbz_quds_iti_anAngsTst(trial,1,3)==4
          outputMatrix(trial,5)=round(360-(acosd((testDot(1))./radius_to_present)));%round(360-(acosd((testDot(1))./r)),2); % not compatible w/Matlab on Prisma
      elseif allPrbz_prbz_quds_iti_anAngsTst(trial,1,3)==3
          outputMatrix(trial,5)=round(360-(acosd((testDot(1))./radius_to_present)));%round(360-(acosd((testDot(1))./r)),2);
      else
          outputMatrix(trial,5)=round(acosd((testDot(1))./radius_to_present));%round(acosd((testDot(1))./r),2); %  
      end   
      %quadrant where tested dot was located
      outputMatrix(trial,6)=allPrbz_prbz_quds_iti_anAngsTst(trial,1,3);
      %probability tested for quadrant
      outputMatrix(trial,7)=probTested;
      %actual cartesian location of all dots
      outputMatrix(trial,8)=pixLocOfdots(1,trial,1)./ppd;
      outputMatrix(trial,9)=pixLocOfdots(2,trial,1)./ppd;
      outputMatrix(trial,10)=pixLocOfdots(1,trial,2)./ppd;
      outputMatrix(trial,11)=pixLocOfdots(2,trial,2)./ppd;
      outputMatrix(trial,12)=pixLocOfdots(1,trial,3)./ppd;
      outputMatrix(trial,13)=pixLocOfdots(2,trial,3)./ppd;
      outputMatrix(trial,14)=pixLocOfdots(1,trial,4)./ppd;
      outputMatrix(trial,15)=pixLocOfdots(2,trial,4)./ppd;
      %location of all dots in Pixels
      outputMatrix(trial,16)=pixLocOfdots(1,trial,1);
      outputMatrix(trial,17)=pixLocOfdots(2,trial,1);
      outputMatrix(trial,18)=pixLocOfdots(1,trial,2);
      outputMatrix(trial,19)=pixLocOfdots(2,trial,2);
      outputMatrix(trial,20)=pixLocOfdots(1,trial,3);
      outputMatrix(trial,21)=pixLocOfdots(2,trial,3);
      outputMatrix(trial,22)=pixLocOfdots(1,trial,4);
      outputMatrix(trial,23)=pixLocOfdots(2,trial,4);
      %Probability of all quadrants
       outputMatrix(trial,24)=newArangInd(1);
       outputMatrix(trial,25)=newArangInd(2);
       outputMatrix(trial,26)=newArangInd(3);
       outputMatrix(trial,27)=newArangInd(4);
      %% Trial Time Data
      %time for trial
      outputMatrix(trial,28)=trialTime;
      if do_el
          outputMatrix(trial,29)=trackerTime;
      else
          outputMatrix(trial,29)=nan;
      end
      %ITI  
      outputMatrix(trial,30)=allPrbz_prbz_quds_iti_anAngsTst(trial,1,4);
      %save output matrix
      if 1==numel(num2str(run))
        datafile2 =sprintf('outptMatrx_P%g_r0%g.mat',subjectID,run);
      else
        datafile2 =sprintf('outptMatrx_P%g_r%g.mat',subjectID,run);
      end
      save(datafile2,'outputMatrix'); 
     %% Save to .csv
     fprintf(fid , ...
         '%d\t%d\t%d\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%d\t%s\n',...
         subjectID, run,trial, trialTime, allPrbz_prbz_quds_iti_anAngsTst(trial,1,4),...
         allPrbz_prbz_quds_iti_anAngsTst(trial,1,3), probTested, tarX, tarY, outputMatrix(trial,5),...
         newArangInd(1), outputMatrix(trial,8), outputMatrix(trial,9), allPrbz_prbz_quds_iti_anAngsTst(trial,1,5),...
         newArangInd(2), outputMatrix(trial,10), outputMatrix(trial,11), allPrbz_prbz_quds_iti_anAngsTst(trial,2,5),...
         newArangInd(3), outputMatrix(trial,12), outputMatrix(trial,13), allPrbz_prbz_quds_iti_anAngsTst(trial,3,5),...
         newArangInd(4), outputMatrix(trial,14), outputMatrix(trial,15), allPrbz_prbz_quds_iti_anAngsTst(trial,4,5), datestr(now));
end
     %% Close .csv
     fclose(fid); %close file stream      
     %% Conclusion Screen
     Screen('FillRect', window, black);
     Screen('FillOval',window, darkgray, theAperture);
     conclText=sprintf('Run %01.f is over.', run);
     DrawFormattedText(window, conclText,'center', 'center', white);
     Screen('Flip',window);
     %% Stop ET Recording 
     if do_el
         Eyelink('stoprecording');
         %Eyelink('closefile');
         Eyelink('ReceiveFile',edfFile,pwd,1);
         Eyelink('Shutdown') 
     end
     %Wait
     KbStrokeWait;
     % Saving all Workspace 
%% So this will be confusing because it will save some variables that are just at the end because of the trial loop 
 % we want this in order to save the experiment start trial and all other
 % screen preferences .... so maybe have it before the trial loop
 % starts...????
     if 1==numel(num2str(run))
         theFileName=(['allVarbls_P',num2str(subjectID),'_r0',num2str(run),'.mat']);     
     else
         theFileName=(['allVarbls_P',num2str(subjectID),'_r',num2str(run),'.mat']); 
     end
     save(theFileName)
     %% Done. Show cursor and close window.
     ShowCursor;
     Screen('CloseAll');
     sca;
     if scanner
        disp('');
        disp('!!!!! Make sure to "clear all" and "clear mex" if on windows computer !!!!!!!');
        disp('');
     else 
        disp('');
        disp('If in lab experiment room: restart MATLAB for every new block/run');
        disp('');
     end
catch
    % This section is executed in case an error happens in the
    % experiment code implemented between try and catch...
    ShowCursor;
    Screen('CloseAll'); % AKA (sca)
    psychrethrow(psychlasterror);
end;