import numpy as np
from scipy import stats




def summarize_data(data, method='mean', srange=None, axis=1):
    """
    Retruns a summary of the supplied data along the axis requested
    
    Method must be one of [mean, median, max, min, sum]. A special value of method is 
    the value "all", in which case this function will return a dict with keys representing
    summary methods with their corresponding values.
    
    Parameters:
        data: an n-d array of data to summarize
        method: the summary statistic to compute
        range: tuple with min/max index of the data subrange to summarize'
        axis: axis to summarize data along.
        
    
    """
    avail_methods = ['mean', 'median', 'max', 'min', 'sum']
    if method == 'all':
        results = {}
        for m in avail_methods:
            results[m] = summarize_data(data, m, srange, axis)
    else:
        if method not in avail_methods:
            raise ValueError("Method must be one of "+", ".join(avail_methods))
        
        if srange is not None:
            d = data[:,srange[0]:srange[1]]
        else:
            d = data
        
        if method == 'median' and isinstance(data, np.ma.MaskedArray):
            return np.ma.median(d, axis=axis)
        else:
            return getattr(d, method)(axis=axis)


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


def max_exclude_outliers(data, threshold=3.5):
    return data[(~is_outlier(data, threshold))].max()


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
#end correct_invalid_data()


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

def make_random_shuffle(df, columns, n=1, axis=0):
    """
    Given a pandas dataframe, returns a new (copied) dataframe 
    where the columns specificied are replaced with the origional
    data that has been randomly shuffled.
    
    Parameters:
        df: pandas dataframe to operate on
        columns: iterable of column indicies to operate on
        
    Returns:
        copy of df with `columns` shuffled
    """
    rand_df = df.copy()
    for c in columns:
        rand_df[:, c] = np.random.shuffle(rand_df.ix[:, c])
    return rand_df
#end make_random_shuffle()

def make_random_gaussian(df, columns):
    """
    Given a pandas dataframe, returns a new (copied) dataframe 
    where the columns specificied are replace with random
    gaussian data with the same distribution of the origional
    data (same mean, std dev, and count).
    
    Parameters:
        df: pandas dataframe to operate on
        columns: iterable of column indicies to operate on
        
    Returns:
        copy of df with `columns` randomized
    """
    rand_df = df.copy()
    for c in columns:
        mu = df.ix[:, c].mean()
        sigma = df.ix[:, c].std()
        count = df.ix[:, c].count()
        rand_df.ix[:, c] = np.random.normal(mu, sigma, count)
    return rand_df
#end make_random_gaussian()





