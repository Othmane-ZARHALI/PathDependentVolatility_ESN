"""Known-mean Gaussian controls fitted at the complete-replication level."""
import numpy as np


def gaussian_control_rows(normals,weights,basis):
    """Return w-1, w U'z, w((U'z)^2-1); each has target expectation zero.

    z is the full whitened target N(0,I) coordinate, w=f/q, and U has
    orthonormal columns. This formula is for random draws from q, with a plain
    mean, not deterministic allocation weights. No observed normalization.
    """
    z,w,U=np.asarray(normals,float),np.asarray(weights,float),np.asarray(basis,float)
    if (z.ndim!=2 or w.shape!=(len(z),) or U.ndim!=2 or U.shape[0]!=z.shape[1]
            or not 1<=U.shape[1]<=U.shape[0] or not np.isfinite(z).all()
            or not np.isfinite(w).all() or np.any(w<=0) or not np.isfinite(U).all()
            or not np.allclose(U.T@U,np.eye(U.shape[1]),rtol=0,atol=2e-12)):
        raise ValueError('Require finite full states, positive density weights and orthonormal basis')
    y=z@U
    return np.column_stack((w-1,w[:,None]*y,w[:,None]*(y*y-1)))


def fit_replication_controls(outcomes,controls,ridge=1.):
    """Fit slopes from independent complete pilot replications, no intercept use.

    Center pilot arrays to estimate covariance, standardize each control by
    sample SD, and solve ridge regression with fixed dimensionless ridge.
    The checking estimator subtracts uncentered zero-mean controls times these
    frozen slopes. Centering here never changes their exact checking mean.
    """
    Y,C=np.asarray(outcomes,float),np.asarray(controls,float)
    if (Y.ndim!=2 or C.ndim!=2 or len(Y)!=len(C) or len(Y)<4 or C.shape[1]<1
            or not np.isfinite(Y).all() or not np.isfinite(C).all()
            or not np.isfinite(ridge) or ridge<=0):
        raise ValueError('Require >=4 complete finite pilot replications and positive ridge')
    centered=C-C.mean(axis=0);target=Y-Y.mean(axis=0)
    sd=np.sqrt(np.sum(centered*centered,axis=0)/(len(C)-1));sd=np.where(sd>0,sd,1.)
    X=centered/sd
    gram=X.T@X/(len(C)-1);rhs=X.T@target/(len(C)-1)
    return np.linalg.solve(gram+ridge*np.eye(C.shape[1]),rhs)/sd[:,None]
