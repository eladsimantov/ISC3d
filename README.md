# ISC3d
A toolbox to streamline the computation and analysis of Intersegmental Coordination in human gait studies. The core algorithm computes 3D elevation angles based solely on 3D anatomical joint angles. Intersegmental Coordination metrics can then be easily computed. Reference frames adhered to ISB recommendations.

Designed in the Neurorobotics and Bionic Limbs Lab (eNaBLe) in the Mechanical Engineering faculty, Technion - IIT.

This work was funded by the Israel Science Foundation grant 2937/24.
## Elevation Angles
To compute the rotation matrix per segment, use the following command

```matlab
[R_thigh, R_shank, R_foot] = isc_joint2rotm(pelvis, hip, knee, ankle, side)
```

## The Elevation Jacobian
It is important to input the joint angles in their correct signs for left and right. 
In addition, mind the ankle replacement of rotation with inversion as in ISB standards.
Assuming the angles are given in VICON format, the usage is as follows:


```matlab
% Define the joint coordinate vector
q_left = [-pelvicTilt, pelvicObliquity, pelvicRotation, hipFlexion, -hipAdduction, -hipRotation, -kneeFlexion, -kneeAdduction, -kneeRotation, offsetAnkleDorsiflex + ankleDorsiflexion, -ankleRotation, -ankleInversion]; 
q_right = [-pelvicTilt, -pelvicObliquity, pelvicRotation, hipFlexion, hipAdduction, hipRotation, -kneeFlexion, kneeAdduction, kneeRotation, offsetAnkleDorsiflex + ankleDorsiflexion, ankleRotation, -ankleInversion];

% Compute the elevation Jacobians - A (3,12) matrix, at a given configuration (q).
J_alpha_left = isc_elevationJacobian(q_left)
J_alpha_right = isc_elevationJacobian(q_right)
```


## Citation
If you use this toolbox in your research, please cite:

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