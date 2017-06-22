import numpy as np
from scipy import stats




def summarize_data(data, method='mean', srange=(0,-1), axis=1):
    avail_methods = ['mean', 'median', 'max', 'min', 'sum']
    if method not in avail_methods:
        raise ValueError("Method must be one of "+", ".join(avail_methods))
    
    if method == 'median' and isinstance(data, np.ma.MaskedArray):
        return np.ma.median(data, axis=axis, keepdims=True)
    else:
        return getattr(data, method)(axis=axis)


def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh
#end is_outlier()


def compute_error(data, method, axis=0, ci=0.95):
    """Compute the error 
    
    Parameters:
        data: array-like
        method: (string) error method to use, one of {sem, std, ci} for Standard Error, Standard Deviation or Confidence Interval
        axis: (int) Axis of the array to operate over (0 for column-wise, 1 for row-wise)
        ci: (float) Only used if method = ci. Confidence interval to compute
        
    Returns:
        two tuple, upper and lower error. i.e. (0.01, 0.01)
    """
    std = np.std(data, axis=axis)
    n = data.shape[axis]
    sem = std / np.sqrt(n)
    
    if method == 'sem':
        return (sem, sem)
    elif method == 'std':
        return (std, std)
    else: #ci
        h = sem * stats.t._ppf((1 + float(ci)) / 2., n - 1)
        return (h, h)
#end compute_error()


def correct_invalid_data(data, method):
    """Corrects invalid data points in a numpy array
    
    Method controls the method of invalid value correction:
        ignore masks invalid values
        zero sets invalid values to zero
    
    Parameters:
        data: array of data
        method: one of {ignore, zero}
        
    Returns:
        data with invalid values corrected
    """
    if method == "ignore":
        data = np.ma.masked_array(data, np.isnan(data))
    elif method == "zero":
        data = data[np.isnan(data)] = 0
    return data
#end correct_out_of_bounds_data()


def format_poly_equation(poly, variable="x", precision=3):
    """Produce a human readable string representation of a polynomial 
       
    Parameters:
        poly: coefficients produced by numpy.polyfit
        variable: string variable to use in the polynomail string
        precision: precision of the resulting polynomial
        
    Returns:
        String representation of the polynomial equation
    """
    # joiner[first, negative] = str
    joiner = {
        (True, True): '-',
        (True, False): '',
        (False, True): ' - ',
        (False, False): ' + '
    }
    cf = ':0.'+str(int(precision))+'f'

    result = []
    for power, coeff in enumerate(reversed(poly)):
        j = joiner[not result, coeff < 0]
        coeff = abs(coeff)
        if coeff == 1 and power != 0:
            coeff = ''

        f = {0: '{}{'+cf+'}', 1: '{}{'+cf+'}{}'}.get(power, '{}{'+cf+'}{}^{}')
        result.append(f.format(j, coeff, variable, power))

    return ''.join(result) or '0'
#end format_poly_equation()





