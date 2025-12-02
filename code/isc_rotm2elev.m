function [alpha, beta] = isc_rotm2elev(RotmArray)
% Input: RotmArray - 3x3xN rotation matrices (global←segment)
% Output: alpha - sagittal elevation angles (XY plane) in DEGREE
%         beta  - frontal elevation angle (XZ plane) in DEGREE
% Usage:
%   [alpha_thigh, beta_thigh] = isc_rotm2elev(R_thigh)
%   [alpha_shank, beta_shank] = isc_rotm2elev(R_shank)
%   [alpha_foot, beta_foot] = isc_rotm2elev(R_foot)

arguments (Input)
    RotmArray (3,3,:) double
end
arguments (Output)
    alpha (:,1) double
    beta (:,1) double
end
segmentAxis = squeeze(RotmArray(:,2,:))'; % N×3, segment direction (y-axis)
% Borghese formula:   
% alpha = atan2((distalCoords(:,1) - proximalCoords(:,1)), (proximalCoords(:,2) - distalCoords(:,2))); % Sagittal (XY)
% beta  = atan2((proximalCoords(:,3) - distalCoords(:,3)), (proximalCoords(:,2) - distalCoords(:,2))); % Frontal (YZ)

% our segment Y axis is exactly the proximal minus distal (vector from distal to proximal) coords in 3D!
alpha = rad2deg(atan2(-segmentAxis(:,1),segmentAxis(:,2))); % sagittal (XY in ISB lab frame)
beta  = rad2deg(atan2(segmentAxis(:,3), segmentAxis(:,2))); % frontal (XZ in ISB lab frame)
end
