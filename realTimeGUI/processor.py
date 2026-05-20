import numpy as np

class GaitProcessor:
    def __init__(self, window_size):
        self.window_size = window_size
        self.buffer = []

    def update(self, thigh, shank, foot):
        angles = [thigh, shank, foot]
        self.buffer.append(angles)
        
        if len(self.buffer) > self.window_size:
            self.buffer.pop(0)
            
        pca_data = {"mean": [0,0,0], "pc1": [0,0,0], "pc2": [0,0,0], "pc3": [0,0,0], "var": [0,0,0]}
        
        # Compute PCA when buffer is full
        if len(self.buffer) == self.window_size:
            mat = np.array(self.buffer)
            mean_vals = np.mean(mat, axis=0)
            mat_centered = mat - mean_vals
            cov = np.cov(mat_centered.T)
            
            if not np.any(np.isnan(cov)):
                evals, evecs = np.linalg.eigh(cov)
                # Sort descending 
                idx = np.argsort(evals)[::-1]
                evecs = evecs[:, idx]
                evals = evals[idx]
                
                # Calculate % Variance Explained
                var_exp = ((evals / np.sum(evals)) * 100).tolist()
                
                pca_data = {
                    "mean": mean_vals.tolist(),
                    "pc1": evecs[:, 0].tolist(),
                    "pc2": evecs[:, 1].tolist(),
                    "pc3": evecs[:, 2].tolist(),
                    "var": var_exp
                }

        return {"angles": {"thigh": thigh, "shank": shank, "foot": foot}, "pca": pca_data}