function direction=find_upcoming_direction(current_target,next_target)
%% find_upcoming_direction finds the angle of the upcoming arm reach based on the targets' position

% INPUTS
% current_target: [x,y] position of the current target
% next_target: [x,y] position of the following target

% OUTPUT
% direction=angle between vector direction and the vector [1 0] [rad]

% 
% 23/01/2023
% Andrea Colins Rodriguez

tmp_dir=next_target-current_target;
direction=atan2(tmp_dir(:,2),tmp_dir(:,1));

end