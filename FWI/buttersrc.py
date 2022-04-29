import numpy as np
import scipy
import matplotlib.pyplot as plt
from devito import configuration
configuration['log-level'] = 'WARNING'
from scipy import signal
from devito.tools import Pickable

from scipy import interpolate
from cached_property import cached_property
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    plt = None

from devito.types import SparseTimeFunction

__all__ = ['PointSource', 'Receiver', 'Shot', 'WaveletSource',
           'RickerSource', 'GaborSource', 'DGaussSource', 'TimeAxis']


class TimeAxis(object):
    """
    Data object to store the TimeAxis. Exactly three of the four key arguments
    must be prescribed. Because of remainder values it is not possible to create
    a TimeAxis that exactly adhears to the inputs therefore start, stop, step
    and num values should be taken from the TimeAxis object rather than relying
    upon the input values.

    The four possible cases are:
    start is None: start = step*(1 - num) + stop
    step is None: step = (stop - start)/(num - 1)
    num is None: num = ceil((stop - start + step)/step);
                 because of remainder stop = step*(num - 1) + start
    stop is None: stop = step*(num - 1) + start

    Parameters
    ----------
    start : float, optional
        Start of time axis.
    step : float, optional
        Time interval.
    num : int, optional
        Number of values (Note: this is the number of intervals + 1).
        Stop value is reset to correct for remainder.
    stop : float, optional
        End time.
    """
    def __init__(self, start=None, step=None, num=None, stop=None):
        try:
            if start is None:
                start = step*(1 - num) + stop
            elif step is None:
                step = (stop - start)/(num - 1)
            elif num is None:
                num = int(np.ceil((stop - start + step)/step))
                stop = step*(num - 1) + start
            elif stop is None:
                stop = step*(num - 1) + start
            else:
                raise ValueError("Only three of start, step, num and stop may be set")
        except:
            raise ValueError("Three of args start, step, num and stop may be set")

        if not isinstance(num, int):
            raise TypeError("input argument must be of type int")

        self.start = start
        self.stop = stop
        self.step = step
        self.num = num

    def __str__(self):
        return "TimeAxis: start=%g, stop=%g, step=%g, num=%g" % \
               (self.start, self.stop, self.step, self.num)

    def _rebuild(self):
        return TimeAxis(start=self.start, stop=self.stop, num=self.num)

    @cached_property
    def time_values(self):
        return np.linspace(self.start, self.stop, self.num)


class PointSource(SparseTimeFunction):
    """Symbolic data object for a set of sparse point sources

    Parameters
    ----------
    name : str
        Name of the symbol representing this source.
    grid : Grid
        The computational domain.
    time_range : TimeAxis
        TimeAxis(start, step, num) object.
    npoint : int, optional
        Number of sparse points represented by this source.
    data : ndarray, optional
        Data values to initialise point data.
    coordinates : ndarray, optional
        Point coordinates for this source.
    space_order : int, optional
        Space discretization order.
    time_order : int, optional
        Time discretization order (defaults to 2).
    dtype : data-type, optional
        Data type of the buffered data.
    dimension : Dimension, optional
        Represents the number of points in this source.
    """

    @classmethod
    def __args_setup__(cls, *args, **kwargs):
        kwargs['nt'] = kwargs['time_range'].num

        # Either `npoint` or `coordinates` must be provided
        npoint = kwargs.get('npoint')
        if npoint is None:
            coordinates = kwargs.get('coordinates', kwargs.get('coordinates_data'))
            if coordinates is None:
                raise TypeError("Need either `npoint` or `coordinates`")
            kwargs['npoint'] = coordinates.shape[0]

        return args, kwargs

    def __init_finalize__(self, *args, **kwargs):
        time_range = kwargs.pop('time_range')
        data = kwargs.pop('data', None)

        kwargs.setdefault('time_order', 2)
        super(PointSource, self).__init_finalize__(*args, **kwargs)

        self._time_range = time_range._rebuild()

        # If provided, copy initial data into the allocated buffer
        if data is not None:
            self.data[:] = data

    @cached_property
    def time_values(self):
        return self._time_range.time_values

    @property
    def time_range(self):
        return self._time_range

    def resample(self, dt=None, num=None, rtol=1e-5, order=3):
        # Only one of dt or num may be set.
        if dt is None:
            assert num is not None
        else:
            assert num is None

        start, stop = self._time_range.start, self._time_range.stop
        dt0 = self._time_range.step

        if dt is None:
            new_time_range = TimeAxis(start=start, stop=stop, num=num)
            dt = new_time_range.step
        else:
            new_time_range = TimeAxis(start=start, stop=stop, step=dt)

        if np.isclose(dt, dt0):
            return self

        nsamples, ntraces = self.data.shape

        new_traces = np.zeros((new_time_range.num, ntraces))

        for i in range(ntraces):
            tck = interpolate.splrep(self._time_range.time_values,
                                     self.data[:, i], k=order)
            new_traces[:, i] = interpolate.splev(new_time_range.time_values, tck)

        # Return new object
        return PointSource(name=self.name, grid=self.grid, data=new_traces,
                           time_range=new_time_range, coordinates=self.coordinates.data)

    # Pickling support
    _pickle_kwargs = SparseTimeFunction._pickle_kwargs + ['time_range']
    _pickle_kwargs.remove('nt')  # `nt` is inferred from `time_range`


Receiver = PointSource
Shot = PointSource


class WaveletSource(PointSource):

    """
    Abstract base class for symbolic objects that encapsulate a set of
    sources with a pre-defined source signal wavelet.

    Parameters
    ----------
    name : str
        Name for the resulting symbol.
    grid : Grid
        The computational domain.
    f0 : float
        Peak frequency for Ricker wavelet in kHz.
    time_values : TimeAxis
        Discretized values of time in ms.
    a : float, optional
        Amplitude of the wavelet (defaults to 1).
    t0 : float, optional
        Firing time (defaults to 1 / f0)
    """

    @classmethod
    def __args_setup__(cls, *args, **kwargs):
        kwargs.setdefault('npoint', 1)

        return super(WaveletSource, cls).__args_setup__(*args, **kwargs)

    def __init_finalize__(self, *args, **kwargs):
        super(WaveletSource, self).__init_finalize__(*args, **kwargs)

        self.f0 = kwargs.get('f0')
        self.a = kwargs.get('a')
        self.t0 = kwargs.get('t0')
        self.fc = kwargs.get('fc')
        self.fs = kwargs.get('fs')
        for p in range(kwargs['npoint']):
            self.data[:, p] = self.wavelet

    @property
    def wavelet(self):
        """
        Return a wavelet with a peak frequency ``f0`` at time ``t0``.
        """
        raise NotImplementedError('Wavelet not defined')

    def show(self, idx=0, wavelet=None, wh=None):
        """
        Plot the wavelet of the specified source.

        Parameters
        ----------
        idx : int
            Index of the source point for which to plot wavelet.
        wavelet : ndarray or callable
            Prescribed wavelet instead of one from this symbol.
        """
        wavelet = wavelet or self.data[:, idx]
        plt.figure()
        plt.plot(self.time_values, wavelet)
        plt.xlabel('Time (ms)')
        plt.ylabel('Amplitude')
        plt.tick_params()
        plt.show()
    # Pickling support
    _pickle_kwargs = PointSource._pickle_kwargs + ['f0', 'a', 'f0']


class RickerSource(WaveletSource):

    """
    Symbolic object that encapsulate a set of sources with a
    pre-defined Ricker wavelet:

    http://subsurfwiki.org/wiki/Ricker_wavelet

    Parameters
    ----------
    name : str
        Name for the resulting symbol.
    grid : Grid
        The computational domain.
    f0 : float
        Peak frequency for Ricker wavelet in kHz.
    time : TimeAxis
        Discretized values of time in ms.

    Returns
    ----------
    A Ricker wavelet.
    """

    @property
    def wavelet(self):
        t0 = self.t0 or 1 / self.f0
        a = self.a or 1
        r = (np.pi * self.f0 * (self.time_values - t0))
        return a * (1-2.*r**2)*np.exp(-r**2)


class GaborSource(WaveletSource):

    """
    Symbolic object that encapsulate a set of sources with a
    pre-defined Gabor wavelet:

    https://en.wikipedia.org/wiki/Gabor_wavelet

    Parameters
    ----------
    name : str
        Name for the resulting symbol.
    grid : Grid
        defining the computational domain.
    f0 : float
        Peak frequency for Ricker wavelet in kHz.
    time : TimeAxis
        Discretized values of time in ms.

    Returns
    -------
    A Gabor wavelet.
    """

    @property
    def wavelet(self):
        agauss = 0.5 * self.f0
        tcut = self.t0 or 1.5 / agauss
        s = (self.time_values - tcut) * agauss
        a = self.a or 1
        return a * np.exp(-2*s**2) * np.cos(2 * np.pi * s)


class DGaussSource(WaveletSource):

    """
    Symbolic object that encapsulate a set of sources with a
    pre-defined 1st derivative wavelet of a Gaussian Source.

    Notes
    -----
    For visualizing the second or third order derivative
    of Gaussian wavelets, the convention is to use the
    negative of the normalized derivative. In the case
    of the second derivative, scaling by -1 produces a
    wavelet with its main lobe in the positive y direction.
    This scaling also makes the Gaussian wavelet resemble
    the Mexican hat, or Ricker, wavelet. The validity of
    the wavelet is not affected by the -1 scaling factor.

    Parameters
    ----------
    name : str
        Name for the resulting symbol.
    grid : Grid
        The computational domain.
    f0 : float
        Peak frequency for wavelet in kHz.
    time : TimeAxis
        Discretized values of time in ms.

    Returns
    -------
    The 1st order derivative of the Gaussian wavelet.
    """

    @property
    def wavelet(self):
        t0 = self.t0 or 1 / self.f0
        a = self.a or 1
        time = (self.time_values - t0)
        return -2 * a * time * np.exp(- a * time**2)

# Make source bandlimitied from 10 Hz - 600 Hz 
class ButterSource(WaveletSource):
    @property
    def wavelet(self):
        t0 = self.t0 
        a = self.a or 1
        fc = self.fc        # cut off frequency
        fs = self.fs        # sampling frequency
        f0 = self.f0
        
        # Ricker Wavelet
        r = (np.pi * f0 * (self.time_values - 1./f0))
        wave =  a * (1-2.*r**2)*np.exp(-r**2)
        
        # Butterworth Filter
        sos = signal.butter(5, self.fc, 'high', fs=fs, output='sos')
        filtered_src = signal.sosfilt(sos, wave)
        
#        b,c = signal.butter(5, self.fc, 'high', analog=True)
#        w, h = signal.freqs(b,c)

#        plt.semilogx(w, 20*np.log10(abs(h)))
#        plt.semilogx(w, h)
#        plt.title('Butterworth filter frequency response')
#        plt.xlabel('Frequency in kHz')
#        plt.ylabel('Amplitude in dB')
#        plt.grid(which='both', axis='both')
#        plt.axvline(fc, color='green')
#        plt.show()

        return filtered_src


# Make special geometry class for banlimited source
class AcquisitionGeometry(Pickable):
    def __init__(self, model, rec_positions, src_positions, t0, tn, **kwargs):
        """
        In practice would be __init__(segyfile) and all below parameters
        would come from a segy_read (at property call rather than at init)
        """
        src_positions = np.reshape(src_positions, (-1, model.dim))
        rec_positions = np.reshape(rec_positions, (-1, model.dim))
        self.rec_positions = rec_positions
        self._nrec = rec_positions.shape[0]
        self.src_positions = src_positions
        self._nsrc = src_positions.shape[0]
        self._src_type = kwargs.get('src_type')
        assert (self.src_type in sources or self.src_type is None)
        self._f0 = kwargs.get('f0')
        self._fs = kwargs.get('fs')
        self._fc = kwargs.get('fc')
        self._a = kwargs.get('a', None)
        self._t0w = kwargs.get('t0w', None)
        if self._src_type is not None and self._f0 is None:
            error("Peak frequency must be provided in KH" +
                  " for source of type %s" % self._src_type)
        self._grid = model.grid
        self._model = model
        self._dt = model.critical_dt
        self._t0 = t0
        self._tn = tn

    def resample(self, dt):
        self._dt = dt
        return self

    @property
    def time_axis(self):
        return TimeAxis(start=self.t0, stop=self.tn, step=self.dt)

    @property
    def src_type(self):
        return self._src_type

    @property
    def grid(self):
        return self._grid

    @property
    def model(self):
        warning("Model is kept for backward compatibility but should not be"
                "obtained from the geometry")
        return self._model

    @property
    def f0(self):
        return self._f0
    
    @property
    def fs(self):
        return self._fs

    @property
    def fc(self):
        return self._fc

    @property
    def tn(self):
        return self._tn

    @property
    def t0(self):
        return self._t0

    @property
    def dt(self):
        return self._dt

    @property
    def nt(self):
        return self.time_axis.num

    @property
    def nrec(self):
        return self._nrec

    @property
    def nsrc(self):
        return self._nsrc

    @property
    def dtype(self):
        return self.grid.dtype

    @property
    def rec(self):
        return self.new_rec()
    
    def new_rec(self, name='rec'):
        return Receiver(name=name, grid=self.grid,
                        time_range=self.time_axis, npoint=self.nrec,
                        coordinates=self.rec_positions)

    @property
    def adj_src(self):
        if self.src_type is None:
            warning("No source type defined, returning uninitiallized (zero) shot record")
            return self.new_rec()
        adj_src = sources[self.src_type](name='rec', grid=self.grid, f0=self.f0,
                                         time_range=self.time_axis, npoint=self.nrec,
                                         coordinates=self.rec_positions,
                                         t0=self._t0w, a=self._a)
        # Revert time axis to have a proper shot record and not compute on zeros
        for i in range(self.nrec):
            adj_src.data[:, i] = adj_src.wavelet[::-1]
        return adj_src

    @property
    def src(self):
        return self.new_src()

    def new_src(self, name='src', src_type='self'):
        if self.src_type is None or src_type is None:
            warning("No surce type defined, returning uninistiallized (zero) source")
            return PointSource(name=name, grid=self.grid,
                               time_range=self.time_axis, npoint=self.nsrc,
                               coordinates=self.src_positions)
        elif self.src_type == 'Butter':
            return ButterSource(name=name, grid=self.grid, 
                                time_range=self.time_axis, npont=self.nsrc,
                                coordinates=self.src_positions, f0=self.f0, fs=self.fs,fc=self.fc)
        else:
            return sources[self.src_type](name=name, grid=self.grid, f0=self.f0,
                                          time_range=self.time_axis, npoint=self.nsrc,
                                          coordinates=self.src_positions,
                                          t0=self._t0w, a=self._a)

    _pickle_args = ['grid', 'rec_positions', 'src_positions', 't0', 'tn']
    _pickle_kwargs = ['f0', 'src_type', 'fs', 'fc' ]

sources = {'Wavelet' : WaveletSource, 'Ricker' : RickerSource, 'Gabor' : GaborSource, 'Butter': ButterSource  }
