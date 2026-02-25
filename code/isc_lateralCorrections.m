function [Qout] = isc_lateralCorrections(Q,side,opts)
%isc_lateralCorrections function takes input joint
%angles/moments/velocities and corrects them to Left and Right sides.
%   NOTICE: in the transformation matrices, we replace ankle inversion and
%   rotation orders. This DOES NOT mean we do so in our definitions of the
%   vectors: q, qdot, tau. This is why in lateral corrections we we DO NOT 
%   replace them.
%
arguments
    Q (:,12) double
    side (1,1) char {mustBeMember(side,{'L','R'})}
    opts.offsetAnkleDorsiflex double = pi/2
end

[pelvicTilt, pelvicObliquity, pelvicRotation] = deal(Q(:,1),Q(:,2),Q(:,3));
[hipFlexion, hipAdduction, hipRotation] = deal(Q(:,4),Q(:,5),Q(:,6));
[kneeFlexion, kneeAdduction, kneeRotation] = deal(Q(:,7),Q(:,8),Q(:,9));
[ankleDorsiflexion, ankleInversion, ankleRotation] = deal(Q(:,10),Q(:,11),Q(:,12));
if strcmp(side,'L')
    % Here there is no ankle replacement of rotation with inversion
    Qout = [-pelvicTilt, pelvicObliquity, pelvicRotation, ...
        hipFlexion, -hipAdduction, -hipRotation, ...
        -kneeFlexion, -kneeAdduction, -kneeRotation, ...
        opts.offsetAnkleDorsiflex + ankleDorsiflexion, ankleInversion, -ankleRotation]; 
elseif  strcmp(side,'R')
    % Here there is no ankle replacement of rotation with inversion
    Qout = [-pelvicTilt, -pelvicObliquity, pelvicRotation, ...
        hipFlexion, hipAdduction, hipRotation, ...
        -kneeFlexion, kneeAdduction, kneeRotation, ...
        opts.offsetAnkleDorsiflex + ankleDorsiflexion, -ankleInversion, ankleRotation]; 
    
end
end