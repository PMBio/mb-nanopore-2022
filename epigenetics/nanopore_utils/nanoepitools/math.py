import numpy as np
from scipy.stats import rankdata
import scipy


def llr_to_p(llr, prior=0.5):
    """
    Convert log-likelihood ratios log(p(x|a)/p(x|~a)) to posterior
    probabilty p(a|x) given a prior p(a). For unbiased prediction,
    leave prior at 0.5
    """
    return 1 / (1 + np.exp(-llr) * (1 / prior - 1))


def p_to_llr(p, prior=0.5, eps=1e-16):
    """
    Converts the posterior probability p(a|x) into a log-likelihood ratio
    log(p(x|a)/p(x|~a)) given a prior p(a)
    """
    p = np.clip(p, eps, 1-eps) # Avoid numerical errors when p is 0 or 1
    return -np.log(prior * (1 - p) / (p * (1 - prior)))


def llr_to_uncertainty(llr, method="linear"):
    if method == "linear":
        p = llr_to_p(llr)
        return 0.5 - np.abs(0.5 - p)


def fdr_from_pvals(p_vals: np.ndarray) -> np.ndarray:
    """
    Computes FDR from p-values using the Benjamini-Hochberg method.
    :param p_vals: numpy array of p-values
    :return: numpy array of adjusted p-values
    """
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    
    return fdr


def bs_from_llrs(llrs: np.ndarray, thres: float = 1, min_reads: int = 1) -> float:
    """
    Computes methylation beta score from a list of log-likelihood ratios
    :param llrs: Log-likelihood ratio array
    :param thres: threshold for absolute llr - excluding all llrs with an absolute llr lower than this threshold
                  (default: 1.0)
    :param min_reads: return np.nan if length of llrs after threshold filtering is less than min_reads (default: 1)
    :return: methylation beta score
    """
    llrs_used = llrs[np.abs(llrs) > thres]
    if len(llrs_used) < min_reads:
        return np.nan
    return (llrs_used > 0).sum() / len(llrs_used)

def __ensure_numpy(x) -> np.ndarray:
    if not isinstance(x, np.ndarray):
        x = np.array(x)
    return x

def nangmean(x: np.ndarray) -> float:
    """ Computes geometric mean while ignoring NaNs """
    x = __ensure_numpy(x)
    x = x[~np.isnan(x)]
    return scipy.stats.gmean(x)

def maxabs(x: np.ndarray) -> float:
    x = __ensure_numpy(x)
    """ Returns the value with the maximum magnitude """
    return x[np.unravel_index(np.argmax(np.abs(x)), x.shape)]

def nanpearsonr(a: np.ndarray, b:np.ndarray):
    idx = ~np.isnan(a) & ~np.isnan(b)
    return scipy.stats.pearsonr(a[idx], b[idx])