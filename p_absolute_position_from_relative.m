function [mic_pos, source_pos]  = p_absolute_position_from_relative(M, array_loc, spacing, rel_source_loc, rotation)
 
    array_rel = [-1*[-(M-1)/2:(M-1)/2]' * spacing,...
                  zeros(M, 2)];  % M * 3 array of Catreisan coordinates of mics
                                 %   relative to array center.
%     az = atan2(rel_source_loc(2), rel_source_loc(1)); % azimuth of source
    az = 0;
    ph = (az*180/pi + rotation)*pi/180;

    RM = [cos(ph) -sin(ph) 0;sin(ph) cos(ph) 0;0 0 1];
    mic_pos = array_rel * RM' + repmat(array_loc, M, 1);

    source_pos = rel_source_loc * RM' + repmat(array_loc, size(rel_source_loc, 1), 1);
    
end
