# ISC3d
A toolbox to streamline the computation and analysis of Intersegmental Coordination in human gait studies. The core algorithm computes 3D elevation angles based solely on 3D anatomical joint angles. Intersegmental Coordination metrics can then be easily computed. Reference frames adhered to ISB recommendations.

Designed in the Neurorobotics and Bionic Limbs Lab (eNaBLe) in the Mechanical Engineering faculty, Technion - IIT.

This work was funded by the Israel Science Foundation grant 2937/24.
## Elevation Angles
This example demonstrates how to compute elevation angles from joint angles, and then compute intersegmental coordination metrics.
```matlab
% Compute the orientations of the thigh, shank, and foot segments, given the pelvis, hip, knee, ankle joint angles and the side (left or right).
[R_thigh, R_shank, R_foot] = isc_joint2rotm(pelvis_deg, hip_deg, knee_deg, ankle_deg, side)

% Compute the elevation angles (alpha, beta) for each segment, given the rotation matrices array.
[alpha_thigh, beta_thigh] = isc_rotm2elev(R_thigh)
[alpha_shank, beta_shank] = isc_rotm2elev(R_shank)
[alpha_foot, beta_foot] = isc_rotm2elev(R_foot)

% Quantify the intersegmental coordination metrics
planarityIndex = isc_quantify(alpha_thigh,alpha_shank,alpha_foot,type="PI")
u3t = isc_quantify(alpha_thigh,alpha_shank,alpha_foot,type="u3t")
PVPC2 = isc_quantify(alpha_thigh,alpha_shank,alpha_foot,type="PVPC2")
```

## The Elevation Jacobian
It is important to input the joint angles in their correct signs for left and right. 
In addition, mind the ankle replacement of rotation with inversion as in ISB standards.
Assuming the angles are given in VICON format, the usage is as follows:


```matlab
% Define the joint coordinate vector for single vector of the angles (q) (1 x 12 angles).
q_left =  [ -pelvicTilt,                pelvicObliquity,      pelvicRotation, ...
             hipFlexion,               -hipAdduction,        -hipRotation, ...
            -kneeFlexion,              -kneeAdduction,       -kneeRotation, ...
             ankleDorsiflexion + pi/2,  ankleInversion,      -ankleRotation]; 
q_right = [ -pelvicTilt,               -pelvicObliquity,      pelvicRotation, ...
             hipFlexion,                hipAdduction,         hipRotation, ...
            -kneeFlexion,               kneeAdduction,        kneeRotation, ...
             ankleDorsiflexion + pi/2, -ankleInversion,       ankleRotation]; 

% Alternatively, use the lateral_corrections function for matrix of angles (Q) (Nsamples x 12 angles)
Q_left = isc_lateralCorrections([pelvis,hip,knee,ankle], 'L');
Q_right = isc_lateralCorrections([pelvis,hip,knee,ankle], 'R');

% Compute the elevation Jacobian - (3,12) matrix, at a given configuration (q), defined as the mapping from pelvis, hip, knee, ankle angular velocities to the angular velocities of the elevation angles in the sagittal plane (alpha_dot).
J_alpha_left = isc_elevationJacobian(q_left)
J_alpha_right = isc_elevationJacobian(q_right)

% Compute the full elevation Jacobian - (12,12) matrix, at a given configuration (q), defined as the mapping from pelvis, hip, knee, ankle angular velocities to the angular velocities of the pelvis, thigh, shank, and foot segments in the elevation space (alpha_dot,beta_dot,gamma_dot).
J_full_left = isc_fullJacobian(q_left)
J_full_right = isc_fullJacobian(q_right)
```

## Elevation Space Moments (ESM)
This example demonstrates how to compute the elevation space moments (ESM) given the joint torques and angles.
```matlab
% Compute ESM for an array of angles (Q) and torques (Tau).
% Note that "Malpha" is the moment in the sagittal plane, and is the one relevant for traditional intersegmental coordination analysis. The "elevationMoments3D" (Malpha; Mbeta; Mgamma) are in the frontal and transverse planes, respectively, and may be used to extend the analysis of intersegmental coordination to 3D.
[elevationMoments_left,Malpha_left] = isc_jointTorque2ElevationMoments(Q_left,Tau_left,'L');
[elevationMoments_right,Malpha_right] = isc_jointTorque2ElevationMoments(Q_right,Tau_right,'R');

% size(Q) = [Ntimesteps, 12, Mtrials] ;
% size(Tau) = [Ntimesteps, 9, Mtrials] ;
% size(Malpha_left) = [Mtrials,Ntimesteps,3]
% size(elevationMoments_left) = [Mtrials,Ntimesteps,12]
```

## Citation
If you use this toolbox in your research, please cite:
**Extending the Law of Intersegmental Coordination: Implications for Powered Prosthetic Controls** *Elad Siman Tov & Nili E. Krausz (2026)* [Read in arXiv](https://arxiv.org/pdf/2602.02181.pdf).

```bibtex
@misc{tov2026extendinglawintersegmentalcoordination,
      title        = {Extending the Law of Intersegmental Coordination: Implications for Powered Prosthetic Controls},
      author       = {Elad Siman Tov and Nili E. Krausz},
      year         = {2026},
      eprint       = {2602.02181},
      archivePrefix= {arXiv},
      primaryClass = {cs.RO},
      url          = {https://arxiv.org/abs/2602.02181}
}
```