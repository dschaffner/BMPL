from matplotlib import pylab
import numpy


def get_corr(t, sig1, sig2, detrend=pylab.detrend_linear,
             mode='same', normalized=1, optimize=1):
    """
    lag,corr = get_corr(t,sig1,sig2,detrend=pylab.detrend_linear,
                             mode='same',normalized=1,):
    NORMALIZED uses (N-1)*sigma1*sigma2 to normalize 
        correlation function (standard deviation normalization)
    OPTIMIZE calls OPTLENGTH on both signals so that the correlation
        runs quickly.  NOTE: Do NOT run without this on raw probe 
        signals.  The number of points is absurdly un-optimized (odd,
        maybe close to prime, etc) and the traces are huge (1-2 Msamples).

    """
    # optimize lengths
    #sig1 = optlength(sig1)
    #sig2 = optlength(sig2)
    # detrend
    # pylab.plot(sig1)
    #sig1 = detrend(sig1)
    # print 'This one!!!'
    meansig1 = numpy.mean(sig1)
    sig1_sub = sig1-meansig1
    # pylab.plot(sig1)
    #sig2 = detrend(sig2)
    meansig2 = numpy.mean(sig2)
    sig2_sub = sig2-meansig2
    n = len(sig1)
    # correlate
    corr = numpy.correlate(sig1_sub, sig2_sub, mode=mode)
    #corr = numpy.correlate(sig1,sig2,mode=mode)
    # normalize
    if normalized:
        corr /= (t.size - 1)*sig1.std()*sig2.std()
    # if normalized: #added by David Schaffner 6-18-2015
    #    corr /= (1.0/(n-1))*sig1.std()*sig2.std()
    # calculate lag array
    dt = t[1]-t[0]
    tau = dt*(numpy.arange(corr.size) - corr.size/2)
    # integer division makes this agree with correlate
    # correlate leaves off the most positive lag value if the number of points is even

    return tau, corr


def get_corr_wmean(t, sig1, sig2, detrend=pylab.detrend_linear,
                   mode='same', normalized=1, optimize=1):
    """
    lag,corr = get_corr(t,sig1,sig2,detrend=pylab.detrend_linear,
                             mode='same',normalized=1,):
    NORMALIZED uses (N-1)*sigma1*sigma2 to normalize 
        correlation function (standard deviation normalization)
    OPTIMIZE calls OPTLENGTH on both signals so that the correlation
        runs quickly.  NOTE: Do NOT run without this on raw probe 
        signals.  The number of points is absurdly un-optimized (odd,
        maybe close to prime, etc) and the traces are huge (1-2 Msamples).

    """
    # optimize lengths
    #sig1 = optlength(sig1)
    #sig2 = optlength(sig2)
    # detrend///
    # pylab.plot(sig1)
    #sig1 = detrend(sig1)
    meansig1 = numpy.mean(sig1)
    sig1_sub = sig1-meansig1
    # pylab.plot(sig1)
    #sig2 = detrend(sig2)
    meansig2 = numpy.mean(sig2)
    sig2_sub = sig2-meansig2
    n = len(sig1)
    # correlate
    corr = numpy.correlate(sig1_sub, sig2_sub, mode=mode)
    # normalize
    if normalized:
        corr /= (t.size - 1)*sig1.std()*sig2.std()
    # if normalized: #added by David Schaffner 6-18-2015
    #    corr /= (1.0/(n-1))*sig1.std()*sig2.std()
    # calculate lag array
    dt = t[1]-t[0]
    tau = dt*(numpy.arange(corr.size) - corr.size/2)
    # integer division makes this agree with correlate
    # correlate leaves off the most positive lag value if the number of points is even

    return tau, corr, meansig1, meansig2
