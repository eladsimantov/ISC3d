function [R_thigh, R_shank, R_foot] = isc_joint2rotm(pelvis, hip, knee, ankle, side)
    % Inputs: 
    %   Nx3 joint angles in DEGREES! (Vicon anatomical convention),
    %   side "L" or "R" because a some angles need to be reversed.
    % Outputs: 
    %   3x3xN rotation matrices for each segment in global frame
    % Usage:
    %   [R_thigh, R_shank, R_fosot] = isc_joint2rotm(pelvis, hip, knee, ankle, "L")
    % 
    arguments (Input)
        pelvis (:,3) double
        hip (:,3) double
        knee (:,3) double
        ankle (:,3) double
        side (1,1) char {mustBeMember(side,{'L','R'})}
    end

    arguments (Output)
        R_thigh (3,3,:) double
        R_shank (3,3,:) double
        R_foot (3,3,:) double
    end

    % Unpack to meaningful angles (TO RADIANS). 
    pelvis = deg2rad(pelvis);
    hip    = deg2rad(hip);
    knee   = deg2rad(knee);
    ankle  = deg2rad(ankle);

    [pelvicTilt, pelvicObliquity, pelvicRotation] = deal(pelvis(:,1),pelvis(:,2),pelvis(:,3));
    [hipFlexion, hipAdduction, hipRotation] = deal(hip(:,1),hip(:,2),hip(:,3));
    [kneeFlexion, kneeAdduction, kneeRotation] = deal(knee(:,1),knee(:,2),knee(:,3));
    [ankleDorsiflexion, ankleInversion, ankleRotation] = deal(ankle(:,1),ankle(:,2),ankle(:,3));
    offsetAnkleDorsiflex = pi/2; % 90 degree offset for foot
    offsetAnkleDorsiflex = pi/2*ones(size(pelvicRotation)); % 90 degree offset for foot

    % Redefine our own Coordinate systems
    % We use the ISB recommendation Wu et al 1995.
    % Pelvis CS: X = Forward, Y = Up,   Z = Right. 
    % Thigh  CS: X = Forward, Y = Up,   Z = Right.
    % Shank  CS: X = Forward, Y = Up,   Z = Right.
    % Foot   CS: X = Forward, Y = Up,   Z = Right. (not rotated yet)
    
    % Euler Intrinsic rotations (first z then x then y = Rz*Rx*Ry)
    % We must adjust the rotation directions according to each
    % individual joint and differently for each side.
    % see https://help.vicon.com/space/Nexus216/11611411/Plug-in+Gait+kinematic+variables
    % see https://media.isbweb.org/images/documents/standards/Wu%20and%20Cavanagh%20J%20Biomech%2028%20(1995)%201258-1261.pdf

    % ----------- Pelvis ----------- 
    % Tilt/Z, Obliquity/X, Rotation/Y
    % Tilt/Z -> Add -ve (negative sign) because +ve pelvic tilt when PSIS
    %   is above the ASIS, which is facing downward and our Z axis
    %   definition (ISB) faces to the right of the body.
    % Obliquity/X -> +Left and -Right because obliquity is when the same side
    %   is tilted above about the X axis facing forward (Sagittal axis).
    % Rotation/Y -> +ve Y axis (I hope..)

    % ----------- Hip -----------
    % Flexion/Z, Adduction/X, Rotation/Y
    % Flexion/Z -> +Right +Left because both sides are +ve Z axis
    % Adduction/X -> +Right -Left because adduction is inward and X
    %   faces forward.
    % Rotation/Y -> +Right -Left because inward rotation (+ve) is about
    %   segment axis which is Y and left side rotating inward is
    %   about the -ve segment axis.

    % ----------- Knee -----------
    % Flexion/Z, Adduction/X, Rotation/Y
    % Flexion/Z -> -Right -Left because both sides are -ve Z axis
    % Adduction/X -> +ve for Varus (Knee bent out such that legs are bowed outward) so 
    %   +Right -Left because when X faces forward Varus is
    %   different for R/L.
    % Rotation/Y -> +ve for inward rotation so 
    %   +Right -Left by the Right hand rule (or left hand rule)

    % ----------- Ankle -----------
    % Dorsiflexion/Z, Rotation/X, Inversion/Y 
    % The ankle joint is weird in the sense that it has flipped the
    % rotations order ZYX rather than ZXY as other joints. In
    % addition it has a constant dorsiflexion offset of 90 degrees
    % because of anatomical position of the foot.
    % Dorsiflexion/Z -> +Right +Left because both sides are +ve Z
    %   axis (some define plantarflexion as +ve so notice both this and the offset!)
    % Rotation/X -> +ve for internal (point foot toward sagittal axis) so 
    %   +Right -Left because when X faces upward (because it was
    %   rotated by dorsiflexion offset about 90 degrees) and we are
    %   in intrinsic rotations.
    % Inversion/Y -> +ve for inward inversion of foot about segment axis so 
    %   -Right +Left about the rotated Y axis of the foot.
    if strcmp(side,'L')
        Rp = eul2rotm([-pelvicTilt, pelvicObliquity, pelvicRotation],'ZXY');  % Tilt/Z, Obliquity/X', Rotation/Y''
        Rh = eul2rotm([hipFlexion, -hipAdduction, -hipRotation],'ZXY');  % Flexion/Z, Adduction/X', Rotation/Y''
        Rk = eul2rotm([-kneeFlexion, -kneeAdduction, -kneeRotation],'ZXY');  % Flexion/Z, Adduction/X', Rotation/Y''
        
        ROffsetAnkle = eul2rotm([offsetAnkleDorsiflex,zeros(length(offsetAnkleDorsiflex),2)],"ZYX");
        Ra = eul2rotm([ankleDorsiflexion, ...
            -ankleRotation, ankleInversion],'ZXY');  % Dorsiflexion/Z, Rotation/X', Inversion/Y''
        Rtmp = pagemtimes(ROffsetAnkle,Ra);
        Ra = eul2rotm([offsetAnkleDorsiflex+ankleDorsiflexion, ...
            -ankleRotation, -ankleInversion],'ZXY');  % Dorsiflexion/Z, Rotation/X', Inversion/Y''
    elseif strcmp(side,'R')
        Rp = eul2rotm([-pelvicTilt, -pelvicObliquity, pelvicRotation],'ZXY');  % Tilt/Z, Obliquity/X', Rotation/Y''
        Rh = eul2rotm([hipFlexion, hipAdduction, hipRotation],'ZXY');  % Flexion/Z, Adduction/X', Rotation/Y''
        Rk = eul2rotm([-kneeFlexion, kneeAdduction, kneeRotation],'ZXY');  % Flexion/Z, Adduction/X', Rotation/Y''
        
        Ra = eul2rotm([offsetAnkleDorsiflex+ankleDorsiflexion, ...
            ankleRotation, -ankleInversion],'ZXY');  % Dorsiflexion/Z, Rotation/X', Inversion/Y''
    end
    % Calculate global segment orientations (postmultiply for
    % intrinsic rotations)
    R_thigh = pagemtimes(Rp, Rh);
    R_shank = pagemtimes(R_thigh, Rk);
    R_foot  = pagemtimes(R_shank, Ra);
end



