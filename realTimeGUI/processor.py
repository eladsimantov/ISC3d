import numpy as np

class GaitProcessor:
    """
    Real-time processor for gait cycle kinematics using sliding window PCA.
    Computes intersegmental coordination (ISC) planar constraints exactly.
    """
    def __init__(self, window_size):
        self.window_size = window_size
        self.buffer = []
        
        # Initialize running statistics for Welford's algorithm
        self.mu = np.zeros(3)
        self.C = np.zeros((3, 3))
        self.n = 0

        # Initialize previous eigenvectors for sign consistency to avoid sudden flips of PCs
        self.prev_evecs = np.eye(3) # Track previous vectors

    def update(self, thigh, shank, foot):
        """
        Updates the sliding window with a new sample and computes the exact PCA.
        
        Uses Welford's online algorithm to recursively add the newest sample 
        and remove the oldest sample from the running unnormalized covariance. 
        This update avoids recomputing the covariance matrix from scratch.
        
        Args:
            thigh (float): Thigh elevation angle.
            shank (float): Shank elevation angle.
            foot (float): Foot elevation angle.
            
        Returns:
            dict: Contains the current raw angles and the updated PCA results 
                  (mean, principal components, and variance explained).
        """
        x = np.array([thigh, shank, foot])
        self.buffer.append(x)
        
        # 1. Add newest sample to running stats
        self.n += 1
        dx = x - self.mu
        self.mu = self.mu + dx / self.n
        self.C = self.C + np.outer(dx, x - self.mu)
        
        # 2. Remove oldest sample if window capacity is exceeded
        if self.n > self.window_size:
            x_old = self.buffer.pop(0)
            self.n -= 1
            # Downdate the mean and covariance
            dx_old = x_old - self.mu
            self.mu = self.mu - dx_old / self.n
            self.C = self.C - np.outer(dx_old, x_old - self.mu)

        # Default fallback output before minimum samples are reached
        pca_data = {"mean": [0,0,0], "pc1": [0,0,0], "pc2": [0,0,0], "pc3": [0,0,0], "var": [0,0,0], "scores": [0,0,0]}
        
        # 3. Compute exact PCA from the running unnormalized covariance
        if self.n >= 3:
            # np.linalg.eigh is fast for symmetric matrices; returns ascending
            evals, evecs = np.linalg.eigh(self.C / (self.n - 1))
            
            # Sort in descending order of variance
            idx = np.argsort(evals)[::-1]
            evecs = evecs[:, idx]
            evals = evals[idx]
            
            var_exp = (evals / np.sum(evals)) * 100
            
            # Prevent sudden flips by aligning with previous vectors
            signs = np.sign(np.sum(self.prev_evecs * evecs, axis=0))
            signs[signs == 0] = 1 # Handle rare exact orthogonality 
            evecs = evecs * signs
            self.prev_evecs = evecs.copy()

            pca_data = {
                "mean": self.mu.tolist(),
                "pc1": evecs[:, 0].tolist(),
                "pc2": evecs[:, 1].tolist(),
                "pc3": evecs[:, 2].tolist(),
                "var": var_exp.tolist(),
                "scores": ((x - self.mu) @ evecs).tolist()
            }

        return {"angles": {"thigh": thigh, "shank": shank, "foot": foot}, "pca": pca_data}


class ExponentialGaitProcessor:
    """
    Real-time processor for gait cycle kinematics using exponential forgetting PCA
    and a calibration prior. Computes planar constraints without window-exit jumps.
    """
    def __init__(self, window_size, calib_data=None, lambda_prior=0.05):
        """
        Args:
            window_size (int): Converted internally to an exponential forgetting factor.
            calib_data (np.ndarray, optional): (N, 3) array of a baseline gait cycle.
            lambda_prior (float): Weight of the prior to prevent rank deficiency.
        """
        self.window_size = window_size
        self.lambda_prior = lambda_prior
        
        # Convert window size to equivalent exponential decay factor
        # e.g., window_size=150 -> alpha=0.9933
        self.alpha = 1.0 - (1.0 / window_size) if window_size > 0 else 0.98
        
        # Initialize previous eigenvectors for sign consistency to avoid sudden flips of PCs
        self.prev_evecs = np.eye(3) # Track previous vectors

        # 1. Setup calibration prior
        if calib_data is not None:
            self.mu_calib = np.mean(calib_data, axis=0)
            self.C_calib = np.cov(calib_data, rowvar=False)
        else:
            # Safe default prior if no data is provided
            self.mu_calib = np.zeros(3)
            self.C_calib = np.eye(3) 
            
        # 2. Initialize running states
        self.mu = self.mu_calib.copy()
        self.C = self.C_calib.copy()
        self.n = 0

    def update(self, thigh, shank, foot):
        """
        Updates the exponential PCA with a new sample.
        Returns the exact same dictionary structure as the sliding window processor.
        """
        x = np.array([thigh, shank, foot])
        self.n += 1
        
        # 1. Recursive mean update
        dx = x - self.mu
        self.mu = self.mu + (1 - self.alpha) * dx
        
        # 2. Recursive covariance update
        self.C = self.alpha * self.C + (1 - self.alpha) * np.outer(dx, dx)
        
        # 3. Regularize with the prior to anchor the plane
        C_eff = self.C + (self.lambda_prior * self.C_calib)
        
        pca_data = {"mean": [0,0,0], "pc1": [0,0,0], "pc2": [0,0,0], "pc3": [0,0,0], "var": [0,0,0], "scores": [0,0,0]}
        
        # 4. Exact PCA on the regularized covariance
        # We can compute immediately because the prior guarantees rank >= 3
        if self.n > 0:
            evals, evecs = np.linalg.eigh(C_eff)
            
            # Sort descending
            idx = np.argsort(evals)[::-1]
            evecs = evecs[:, idx]
            evals = evals[idx]
            
            var_sum = np.sum(evals)
            var_exp = (evals / var_sum) * 100 if var_sum > 0 else np.zeros(3)
            
            # Prevent sudden flips by aligning with previous vectors
            signs = np.sign(np.sum(self.prev_evecs * evecs, axis=0))
            signs[signs == 0] = 1 # Handle rare exact orthogonality 
            evecs = evecs * signs
            self.prev_evecs = evecs.copy()

            pca_data = {
                "mean": self.mu.tolist(),
                "pc1": evecs[:, 0].tolist(),
                "pc2": evecs[:, 1].tolist(),
                "pc3": evecs[:, 2].tolist(),
                "var": var_exp.tolist(),
                "scores": ((x - self.mu) @ evecs).tolist()
            }

        return {"angles": {"thigh": thigh, "shank": shank, "foot": foot}, "pca": pca_data}