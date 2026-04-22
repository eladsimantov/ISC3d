function [alpha, beta, gamma] = isc_rotm2elev3D(RotmArray,side,is_foot)
%ISC_ROTM2ELEV3D  Computes elevation angles arrays for the thigh,
%shank and foot segments given the rotation matrices arrays. This 
%is extended to 3D elevation angles with a special interpertation for
%the foot elevation angle to avoid rapid jumps in the frontal plane elevation angle.
%This function also takes laterality into account (left or right) to correctly compute the angles.
% Input: 
%       RotmArray - 3x3xN rotation matrices (global←segment)
%       side - character indicating the side of the body ('L' for left, 'R' for right)
%       is_foot - boolean indicating whether the segment is the foot (true) or not (false)
% Output: 
%       alpha - sagittal elevation angles (XY plane) in DEGREE
%       beta  - frontal elevation angle (ZY plane) in DEGREE
%       gamma - transverse elevation angle (XZ plane) in DEGREE
% Usage:
%       [alpha_thigh, beta_thigh, gamma_thigh] = isc_rotm2elev3D(R_thigh, 'L', false)
%       [alpha_shank, beta_shank, gamma_shank] = isc_rotm2elev3D(R_shank, 'L', false)
%       [alpha_foot, beta_foot, gamma_foot] = isc_rotm2elev3D(R_foot, 'L', true)
%   -----------------------------------------------------------------------------
%   References:
%   [1] Siman Tov, E., & Krausz, N. E. (2026). "Extending the Law of 
%       Intersegmental Coordination: Implications for Powered 
%       Prosthetic Controls." arXiv:2602.02181.
%       PDF:  https://arxiv.org/pdf/2602.02181.pdf
%   [2] N. A. Borghese, L. Bianchi, and Francesco Lacquaniti, 
%       "Kinematic determinants of human locomotion," 
%       The Journal of Physiology, vol. 494, no. 3, pp. 863–879.
%   -----------------------------------------------------------------------------
%   
%                       Extended Documentation
% 
% --------------------- Borghese formulas --------------------- %
% The origin of ISC according to Borghese formulas:   
% alpha = atan2((distalCoords(:,1) - proximalCoords(:,1)), (proximalCoords(:,2) - distalCoords(:,2))); % Sagittal (XY)
% beta  = atan2((proximalCoords(:,3) - distalCoords(:,3)), (proximalCoords(:,2) - distalCoords(:,2))); % Frontal (YZ)
% 
% The interpretation of the angles using rotation matrices is as follows:
% alpha = atan2(-segmentAxis(:,1), segmentAxis(:,2)); % sagittal
% beta  = atan2(segmentAxis(:,3), segmentAxis(:,2)); % frontal
% 
% -------------------- Laterality --------------------- %
% on the right leg, beta is defined as position about the (-X) axis, and on the left (+X).
% on the right leg, gamma is defined as position about the (-Z) axis, and on the left (+Z).
% Our definition assumes elevation angles in 3D are positive for moving the leg:
%   alpha: forward in sagittal plane (flexion)          (+Z) for Right, (+Z) for Left
%   beta: outward in frontal plane (abduction)          (-X) for Right, (+X) for Left 
%   gamma: outward in transverse (external rotation).   (-Y) for Right, (+Y) for Left
% 
% -------------------- Foot special case --------------------- %
% We define a new interpertation specifically dealing with the foot segment
% We compute the elevation angles of the **normal to foot axis**, (x-axis). 
% If we think of thigh-foot coordination in the Frontal plane, 
% we need a better/new way to define the elevation angles for
% the foot segment because if we just use the standard approach, 
% we get rapid jumps in the elevation angle in frontal plane. 
% This definition does not interfere with the elevation angles in the 
% saggittal plane but rather it just adds a 90 degree offset.
% This is the way.
arguments (Input)
    RotmArray (3,3,:) double
    side (1,1) char {mustBeMember(side,{'L','R'})}
    is_foot (1,1) logical = false
end
arguments (Output)
    alpha (:,1) double
    beta (:,1) double
    gamma (:,1) double
end

% Extract the axes from the rotation matrix array as defined in the ISB frames.
normalAxis = squeeze(RotmArray(:,1,:))'; % N×3, forward direction (x-axis)
segmentAxis = squeeze(RotmArray(:,2,:))'; % N×3, segment up direction (y-axis)
zAxis = squeeze(RotmArray(:,3,:))'; % N×3, right direction (z-axis)

is_left = strcmp(side, 'L');
is_right = strcmp(side, 'R');
is_thigh_shank = ~is_foot;

if is_left && is_foot
    % For the foot segment on the left side.
    % Given atan2(towards axis, from axis)
    % for a (+alpha) on left side, (+Y) axis goes towards the (-X) axis.
    % for a (+beta) on left side, (+Y) axis goes towards the (+Z) axis.
    % for a (+gamma) on left side, (+Z) axis goes towards the (+X) axis.
    alpha = rad2deg(atan2(-normalAxis(:,1),  normalAxis(:,2)));     % sagittal (XY in ISB lab frame)
    beta  = rad2deg(atan2( normalAxis(:,3),  normalAxis(:,2)));     % frontal (ZY in ISB lab frame)
    gamma = rad2deg(atan2( zAxis(:,1),       zAxis(:,3)));          % transverse (XZ in ISB lab frame)
elseif is_left && is_thigh_shank
    % For thigh and shank segments on the left side
    % Given atan2(towards axis, from axis)
    % for a (+alpha) on left side, (+Y) axis goes towards the (-X) axis.
    % for a (+beta) on left side, (+Y) axis goes towards the (+Z) axis. 
    % for a (+gamma) on left side, (+Z) axis goes towards the (+X) axis.
    alpha = rad2deg(atan2(-segmentAxis(:,1), segmentAxis(:,2)));    % sagittal (XY in ISB lab frame)
    beta  = rad2deg(atan2( segmentAxis(:,3), segmentAxis(:,2)));    % frontal (ZY in ISB lab frame)
    gamma = rad2deg(atan2( zAxis(:,1),       zAxis(:,3)));          % transverse (XZ in ISB lab frame)
elseif is_right && is_thigh_shank
    % For thigh and shank segments on the right side
    % Given atan2(towards axis, from axis)
    % for a (+alpha) on right side, (+Y) axis goes towards the (-X) axis.
    % for a (+beta) on right side, (+Y) axis goes towards the (-Z) axis.
    % for a (+gamma) on right side, (+Z) axis goes towards the (-X) axis.
    alpha = rad2deg(atan2(-segmentAxis(:,1), segmentAxis(:,2)));    % sagittal (XY in ISB lab frame)
    beta  = rad2deg(atan2(-segmentAxis(:,3), segmentAxis(:,2)));    % frontal (ZY in ISB lab frame)
    gamma = rad2deg(atan2(-zAxis(:,1),       zAxis(:,3)));          % transverse (XZ in ISB lab frame)
elseif is_right && is_foot
    % For the foot segment on the right side
    % Given atan2(towards axis, from axis)
    % for a (+alpha) on right side, (+Y) axis goes towards the (-X) axis.
    % for a (+beta) on right side, (+Y) axis goes towards the (-Z) axis.
    % for a (+gamma) on right side, (+Z) axis goes towards the (-X) axis.
    alpha = rad2deg(atan2(-normalAxis(:,1),  normalAxis(:,2)));     % sagittal (XY in ISB lab frame)
    beta  = rad2deg(atan2(-normalAxis(:,3),  normalAxis(:,2)));     % frontal (ZY in ISB lab frame)
    gamma = rad2deg(atan2(-zAxis(:,1),       zAxis(:,3)));          % transverse (XZ in ISB lab frame)
end
end