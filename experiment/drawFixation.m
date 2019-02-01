% drawFixation(priorityset)

screenNumber = max(Screen('Screens'));
topPriorityLevel = MaxPriority(screenNumber); %topPriorityLevel1 = MaxPriority(window);
Priority(topPriorityLevel);
[window, screenRect] = Screen('OpenWindow', screenNumber, [0 0 0]);
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
[~, h] = Screen('WindowSize', window);
[x_c, y_c] = RectCenter(screenRect);
center = [x_c y_c];