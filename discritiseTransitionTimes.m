function [ timeOnInFrame ] = discritiseTransitionTimes...
    (times, states, expT, nFrames)
% Calculate the time on in each frame based on the transition
% times and states
 
        stateCount = 1;
        timeOnInFrame = zeros(1,nFrames);
         
        % frame start times and end times
        frameStart = (0:nFrames-1) * expT;
        frameEnd = (1:nFrames) * expT;
                  
        for f = 1:nFrames %loop over frames
             
            % first case the state lasts for the entire frame and the
            % fluorophore was in an on state
            if frameStart(f) >= times(stateCount) && frameEnd(f) <= times(stateCount+1) && states(stateCount) ~= 0
 
                timeOnInFrame(f) = states(stateCount)*expT;
 
            % second case the state transition occurs within the frame
            elseif frameStart(f) >= times(stateCount) && frameEnd(f) > times(stateCount+1)
                % fluorophore on time is added to any previous on time from
                % incase there was a previous state in this frame
                timeOnInFrame(f) = timeOnInFrame(f) + (times(stateCount+1) - frameStart(f))...
                    * states(stateCount);
 
                % state count is stepped forward
                stateCount = stateCount + 1;
                times(stateCount) = times(stateCount);
                times(stateCount+1) = times(stateCount+1);
 
                % make sure the end of the frame is reached including and
                % adding any further transitions in the frame
                while frameEnd(f) > times(stateCount+1)
                    timeOnInFrame(f) = timeOnInFrame(f) + ...
                        (times(stateCount+1) - times(stateCount)) * states(stateCount);
                    stateCount = stateCount + 1;
                    times(stateCount) = times(stateCount);
                    times(stateCount+1) = times(stateCount+1);
                end
                 
                % add any contributions from the final state that ends this
                % frame
                timeOnInFrame(f) = timeOnInFrame(f) + ...
                    (frameEnd(f) - times(stateCount)) * states(stateCount);
%             else
%                 fprintf('Sanity Check'); % were any cases missed?
            end
             
        end
 
end