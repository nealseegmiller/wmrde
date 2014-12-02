function [steer_cmd,yawrate,wheelvel,steer_rate] = zoeArcController(speed,turnrad,steer)
    %convert arc commands into wheel velocity commands
    %INPUTS
    %speed:         (m/s) commanded speed
    %turnrad:       (m) commanded turn radius (positive -> left turn, negative -> right turn)
    %steer:         [front rear]' measured steer angles (rad)
    %OUTPUTS
    %steer_cmd:     commanded steer angles (rad)
    %yawrate:       angular vel about body z axis (=yaw rate on flat ground)
    %wheelvel:      4x1 commanded wheel velocities (rad/s)
    %steer_rate:    2x1 steering rates (rad/s) due to feedback

    L = 1.91;
    B = 1.64;
    rad = .325;
    
    min_rad = 3; %(m) minimum turn radius
    max_dvel = 1.5; %(m/s) 0.5*max difference btwn left/right wheel vel due to feedback

    Kp = 1; %proportional gain
    
    if abs(turnrad) >= 1000 || turnrad == 0
        %drive straight
        steer_cmd = zeros(2,1);
        yawrate = 0;
    else
        %enforce radius limit
        turnrad = sign(turnrad)*max(min_rad,abs(turnrad));

        steer_cmd = atan(L/2/turnrad)*[1 -1]';
        yawrate = speed/turnrad;
    end

    %feedforward velocity
    rcs = 1./cos(steer_cmd); %reciprocal of cosine of commanded steer angles
    
    H = [rcs(1),-B/2; %front left
         rcs(1), B/2; %front right
         rcs(2),-B/2; %rear left
         rcs(2), B/2]; %rear right
    wheelvel = H*[speed yawrate]';

    %add feedback control
    
    err = steer_cmd - steer;
    dvel = Kp*err;
    %enforce limits
    dvel(dvel > max_dvel) = max_dvel;
    dvel(dvel < -max_dvel) = -max_dvel;
    
    wheelvel = wheelvel + dvel([1 1 2 2]).*[-1 1 -1 1]';

    %convert from m/s to rad/s, divide by wheel radius
    wheelvel = wheelvel/rad;
    
    steer_rate = 2*dvel/B;

end

