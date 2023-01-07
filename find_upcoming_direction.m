function direction=find_upcoming_direction(current_target,next_target)
%% current_target [x1,y1] position
%% next_target [x2,y2] position
%% direction=angle between vector direction and [1 0] (rad)

tmp_dir=next_target-current_target;
direction=atan2(tmp_dir(:,2,:),tmp_dir(:,1,:));

end