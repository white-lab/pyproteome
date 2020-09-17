import numpy as np

def vip(model):
    '''
    Calculate VIP scores for a PLSR model.

    Parameters
    ----------
    model

    Returns
    -------
    :class:`numpy.array`
    '''
    try:
        t = model.x_scores_
        w = model.x_weights_
        q = model.y_loadings_
    except:
        t = model.x_scores
        w = model.x_weights
        q = model.y_loadings
    
    p, h = w.shape
    vips = np.zeros((p,))
    s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
    total_s = np.sum(s)
    
    for i in range(p):
        weight = np.array([
            (w[i, j] / np.linalg.norm(w[:, j])) ** 2
            for j in range(h)
        ])
        vips[i] = np.sqrt(p*(s.T @ weight) / total_s)

    return vips