function [elevationsMoment,Malpha] = isc_jointTorque2ElevationMoments(Q,Tau,side)
%isc_jointTorque2ElevationMoments computes elevation moments of p,t,s,f in
%three axes using the full 12x12 jacobian and its pseudo inverse.
% Inputs:
%   Q - pelvis, hip, knee, ankle arrays as they come from MOCAP in DEGREE.
%       [Ntimesteps,Ncols,Mtrials] = size(Q); 
%       Ncols should be 12; -> 3 per joint.
%   Tau - hip, knee, ankle arrays as they come from MOCAP.
%         [Ntimesteps,Ncols,Mtrials] = size(Q); 
%         Ncols should be 9; -> 3 per (hip, knee, ankle) joint
% Outputs:
%   elevationsMoment - [Malpha; Mbeta; Mgamma] in same units as the input. 
%       pelvis, thigh, shank, foot respectively. 
%       [Mtrials,Ntimesteps,12] = size(elevationsDot)
%   Malpha  - elevation moments only for the thigh, shank, foot
%       [Mtrials,Ntimesteps,3] = size(elevationsDot)
arguments
    Q (:,12,:)
    Tau (:,9,:)
    side (1,1) char {mustBeMember(side,{'L','R'})}
end
[Ntimesteps,Ncols,Mtrials] = size(Tau); % Ncols should be 12; -> 3 per joint.

% Init the output array with nan values (Mtrials x Ntimesteps x 12).
elevationsMoment = nan(Mtrials,Ntimesteps,12);
Malpha = nan(Mtrials,Ntimesteps,3);

for m = 1:Mtrials
    % Extract joint angles
    pelvis = squeeze(Q(:,1:3,m));
    hip = squeeze(Q(:,4:6,m));
    knee = squeeze(Q(:,7:9,m));
    ankle = squeeze(Q(:,10:12,m)); % Nx3 joint angles matrices

    % Extract joint torques
    tau_hip = squeeze(Tau(:,1:3,m));
    tau_knee = squeeze(Tau(:,4:6,m));
    tau_ankle = squeeze(Tau(:,7:9,m)); % Nx3 joint angles matrices
    
    % Make corrections to joint angles and velocities from MOCAP outputs.
    QL = isc_lateralCorrections([pelvis,hip,knee,ankle], side);
    QR = isc_lateralCorrections([pelvis,hip,knee,ankle], side);

    % Crucial: do not add an ankle offsets to torques. 
    % To fit with corrections function tau_hip is inserted twice only to be
    % neglected later on..
    QLtau = isc_lateralCorrections([zeros(size(tau_hip)),tau_hip,tau_knee,tau_ankle], side, "offsetAnkleDorsiflex", 0);
    QRtau = isc_lateralCorrections([zeros(size(tau_hip)),tau_hip,tau_knee,tau_ankle], side, "offsetAnkleDorsiflex", 0);

    switch side
        case "L"
            for t = 1:Ntimesteps
                q = QL(t,:);
                tau = QLtau(t,:); % 12x1 vector (with zeros in pelvis..)

                Jalpha = isc_elevationJacobian(deg2rad(q));
                Jinv = pinv(Jalpha);
                Malpha(m,t,1:3) = Jinv(:,1:3).' * tau.';
            end
        case "R"
            for t = 1:Ntimesteps
                q = QR(t,:);
                tau = QRtau(t,:); % 12x1 vector (no pelvis..)

                Jalpha = isc_elevationJacobian(deg2rad(q));
                Jinv = pinv(Jalpha);
                Malpha(m,t,1:3) = Jinv(:,1:3).' * tau.';
            end
    end
end
end