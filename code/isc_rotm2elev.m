function [alpha, beta] = isc_rotm2elev(RotmArray)
%ISC_ROTM2ELEV  Computes elevation angles arrays for the thigh,
%shank and foot segments given the rotation matrices arrays.
% Input: RotmArray - 3x3xN rotation matrices (global←segment)
% Output: alpha - sagittal elevation angles (XY plane) in DEGREE
%         beta  - frontal elevation angle (XZ plane) in DEGREE
% Usage:
%   [alpha_thigh, beta_thigh] = isc_rotm2elev(R_thigh)
%   [alpha_shank, beta_shank] = isc_rotm2elev(R_shank)
%   [alpha_foot, beta_foot] = isc_rotm2elev(R_foot)
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
beta  = rad2deg(atan2(segmentAxis(:,3), segmentAxis(:,2))); % frontal (YZ in ISB lab frame)
end
