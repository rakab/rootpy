from __future__ import absolute_import

import math
import operator
import re
import itertools

from .. import ROOT
from .canvas import _PadBase
from .hist import _Hist, Hist, HistStack
from .graph import _Graph1DBase, Graph
from .profile import Profile
from . import base
from .. import log; log = log[__name__]
from ..context import preserve_current_canvas, do_nothing
from ..extern.six.moves import range

__all__ = [
    'draw',
    'get_limits',
    'get_band',
    'canvases_with',
    'find_all_primitives',
    'tick_length_pixels',
]


def draw(plottables, pad=None, same=False,
         xaxis=None, yaxis=None,
         xtitle=None, ytitle=None,
         xlimits=None, ylimits=None,
         xdivisions=None, ydivisions=None,
         logx=False, logy=False,
         **kwargs):
    """
    Draw a list of histograms, stacks, and/or graphs.

    Parameters
    ----------
    plottables : Hist, Graph, HistStack, or list of such objects
        List of objects to draw.

    pad : Pad or Canvas, optional (default=None)
        The pad to draw onto. If None then use the current global pad.

    same : bool, optional (default=False)
        If True then use 'SAME' draw option for all objects instead of
        all but the first. Use this option if you are drawing onto a pad
        that already holds drawn objects.

    xaxis : TAxis, optional (default=None)
        Use this x-axis or use the x-axis of the first plottable if None.

    yaxis : TAxis, optional (default=None)
        Use this y-axis or use the y-axis of the first plottable if None.

    xtitle : str, optional (default=None)
        Set the x-axis title.

    ytitle : str, optional (default=None)
        Set the y-axis title.

    xlimits : tuple, optional (default=None)
        Set the x-axis limits with a 2-tuple of (min, max)

    ylimits : tuple, optional (default=None)
        Set the y-axis limits with a 2-tuple of (min, max)

    xdivisions : int, optional (default=None)
        Set the number of divisions for the x-axis

    ydivisions : int, optional (default=None)
        Set the number of divisions for the y-axis

    logx : bool, optional (default=False)
        If True, then set the x-axis to log scale.

    logy : bool, optional (default=False)
        If True, then set the y-axis to log scale.

    kwargs : dict
        All extra arguments are passed to get_limits when determining the axis
        limits.

    Returns
    -------
    (xaxis, yaxis), (xmin, xmax, ymin, ymax) : tuple
        The axes and axes bounds.

    See Also
    --------
    get_limits

    """
    context = preserve_current_canvas if pad else do_nothing
    if not isinstance(plottables, (tuple, list)):
        plottables = [plottables]
    elif not plottables:
        raise ValueError("plottables is empty")
    with context():
        if pad is not None:
            pad.cd()
        # get the axes limits
        xmin, xmax, ymin, ymax = get_limits(plottables,
                                            logx=logx, logy=logy,
                                            **kwargs)
        if xlimits is not None:
            xmin, xmax = xlimits
        if ylimits is not None:
            ymin, ymax = ylimits
        if not same:
            obj = plottables.pop(0)
            if isinstance(obj, ROOT.THStack):
                obj.SetMinimum(ymin)
                obj.SetMaximum(ymax)
            obj.Draw()
            xaxis = obj.xaxis
            yaxis = obj.yaxis
        # draw the plottables
        for i, obj in enumerate(plottables):
            if i == 0 and isinstance(obj, ROOT.THStack):
                # use SetMin/Max for y-axis
                obj.SetMinimum(ymin)
                obj.SetMaximum(ymax)
                # ROOT: please fix this...
            obj.Draw('SAME')
        # set the axes limits and titles
        if xaxis is not None:
            xaxis.SetLimits(xmin, xmax)
            xaxis.SetRangeUser(xmin, xmax)
            if xtitle is not None:
                xaxis.SetTitle(xtitle)
            if xdivisions is not None:
                xaxis.SetNdivisions(xdivisions)
        if yaxis is not None:
            yaxis.SetLimits(ymin, ymax)
            yaxis.SetRangeUser(ymin, ymax)
            if ytitle is not None:
                yaxis.SetTitle(ytitle)
            if ydivisions is not None:
                yaxis.SetNdivisions(ydivisions)
        if pad is None:
            pad = ROOT.gPad
        pad.SetLogx(bool(logx))
        pad.SetLogy(bool(logy))
        # redraw axes on top
        # axes ticks sometimes get hidden by filled histograms
        pad.RedrawAxis()
    return (xaxis, yaxis), (xmin, xmax, ymin, ymax)


multiadd = lambda a, b: map(operator.add, a, b)
multisub = lambda a, b: map(operator.sub, a, b)


def _limits_helper(x1, x2, a, b, snap=False):
    """
    Given x1, x2, a, b, where:

        x1 - x0         x3 - x2
    a = ------- ,   b = -------
        x3 - x0         x3 - x0

    determine the points x0 and x3:

    x0         x1                x2       x3
    |----------|-----------------|--------|

    """
    if x2 < x1:
        raise ValueError("x2 < x1")
    if a + b >= 1:
        raise ValueError("a + b >= 1")
    if a < 0:
        raise ValueError("a < 0")
    if b < 0:
        raise ValueError("b < 0")
    if snap:
        if x1 >= 0:
            x1 = 0
            a = 0
        elif x2 <= 0:
            x2 = 0
            b = 0
        if x1 == x2 == 0:
            # garbage in garbage out
            return 0., 1.
    elif x1 == x2:
        # garbage in garbage out
        return x1 - 1., x1 + 1.
    if a == 0 and b == 0:
        return x1, x2
    elif a == 0:
        return x1, (x2 - b * x1) / (1 - b)
    elif b == 0:
        return (x1 - a * x2) / (1 - a), x2
    x0 = ((b / a) * x1 + x2 - (x2 - x1) / (1 - a - b)) / (1 + b / a)
    x3 = (x2 - x1) / (1 - a - b) + x0
    return x0, x3


def get_limits(plottables,
               xpadding=0,
               ypadding=0.1,
               xerror_in_padding=True,
               yerror_in_padding=True,
               snap=True,
               logx=False,
               logy=False,
               logx_crop_value=1E-5,
               logy_crop_value=1E-5,
               logx_base=10,
               logy_base=10):
    """
    Get the axes limits that should be used for a 1D histogram, graph, or stack
    of histograms.

    Parameters
    ----------

    plottables : Hist, Graph, HistStack, or list of such objects
        The object(s) for which visually pleasing plot boundaries are
        requested.

    xpadding : float or 2-tuple, optional (default=0)
        The horizontal padding as a fraction of the final plot width.

    ypadding : float or 2-tuple, optional (default=0.1)
        The vertical padding as a fraction of the final plot height.

    xerror_in_padding : bool, optional (default=True)
        If False then exclude the x error bars from the calculation of the plot
        width.

    yerror_in_padding : bool, optional (default=True)
        If False then exclude the y error bars from the calculation of the plot
        height.

    snap : bool, optional (default=True)
        Make the minimum or maximum of the vertical range the x-axis depending
        on if the plot maximum and minimum are above or below the x-axis. If
        the plot maximum is above the x-axis while the minimum is below the
        x-axis, then this option will have no effect.

    logx : bool, optional (default=False)
        If True, then the x-axis is log scale.

    logy : bool, optional (default=False)
        If True, then the y-axis is log scale.

    logx_crop_value : float, optional (default=1E-5)
        If an x-axis is using a logarithmic scale then crop all non-positive
        values with this value.

    logy_crop_value : float, optional (default=1E-5)
        If the y-axis is using a logarithmic scale then crop all non-positive
        values with this value.

    logx_base : float, optional (default=10)
        The base used for the logarithmic scale of the x-axis.

    logy_base : float, optional (default=10)
        The base used for the logarithmic scale of the y-axis.

    Returns
    -------

    xmin, xmax, ymin, ymax : tuple of plot boundaries
        The computed x and y-axis ranges.

    """
    try:
        import numpy as np
        use_numpy = True
    except ImportError:
        use_numpy = False

    if not isinstance(plottables, (list, tuple)):
        plottables = [plottables]

    xmin = float('+inf')
    xmax = float('-inf')
    ymin = float('+inf')
    ymax = float('-inf')

    for h in plottables:

        if isinstance(h, HistStack):
            h = h.sum

        if not isinstance(h, (_Hist, _Graph1DBase)):
            raise TypeError(
                "unable to determine plot axes ranges "
                "from object of type `{0}`".format(
                    type(h)))

        if use_numpy:
            y_array_min = y_array_max = np.array(list(h.y()))
            if yerror_in_padding:
                y_array_min = y_array_min - np.array(list(h.yerrl()))
                y_array_max = y_array_max + np.array(list(h.yerrh()))
            _ymin = y_array_min.min()
            _ymax = y_array_max.max()
        else:
            y_array_min = y_array_max = list(h.y())
            if yerror_in_padding:
                y_array_min = multisub(y_array_min, list(h.yerrl()))
                y_array_max = multiadd(y_array_max, list(h.yerrh()))
            _ymin = min(y_array_min)
            _ymax = max(y_array_max)

        if isinstance(h, _Graph1DBase):
            if use_numpy:
                x_array_min = x_array_max = np.array(list(h.x()))
                if xerror_in_padding:
                    x_array_min = x_array_min - np.array(list(h.xerrl()))
                    x_array_max = x_array_max + np.array(list(h.xerrh()))
                _xmin = x_array_min.min()
                _xmax = x_array_max.max()
            else:
                x_array_min = x_array_max = list(h.x())
                if xerror_in_padding:
                    x_array_min = multisub(x_array_min, list(h.xerrl()))
                    x_array_max = multiadd(x_array_max, list(h.xerrh()))
                _xmin = min(x_array_min)
                _xmax = max(x_array_max)
        else:
            _xmin = h.xedgesl(1)
            _xmax = h.xedgesh(h.nbins(0))

        if logy:
            _ymin = max(logy_crop_value, _ymin)
            _ymax = max(logy_crop_value, _ymax)
        if logx:
            _xmin = max(logx_crop_value, _xmin)
            _xmax = max(logx_crop_value, _xmax)

        if _xmin < xmin:
            xmin = _xmin
        if _xmax > xmax:
            xmax = _xmax
        if _ymin < ymin:
            ymin = _ymin
        if _ymax > ymax:
            ymax = _ymax

    if isinstance(xpadding, (list, tuple)):
        if len(xpadding) != 2:
            raise ValueError("xpadding must be of length 2")
        xpadding_left = xpadding[0]
        xpadding_right = xpadding[1]
    else:
        xpadding_left = xpadding_right = xpadding

    if isinstance(ypadding, (list, tuple)):
        if len(ypadding) != 2:
            raise ValueError("ypadding must be of length 2")
        ypadding_top = ypadding[0]
        ypadding_bottom = ypadding[1]
    else:
        ypadding_top = ypadding_bottom = ypadding

    if logx:
        x0, x3 = _limits_helper(
            math.log(xmin, logx_base), math.log(xmax, logx_base),
            xpadding_left, xpadding_right)
        xmin = logx_base ** x0
        xmax = logx_base ** x3
    else:
        xmin, xmax = _limits_helper(
            xmin, xmax, xpadding_left, xpadding_right)

    if logy:
        y0, y3 = _limits_helper(
            math.log(ymin, logy_base), math.log(ymax, logy_base),
            ypadding_bottom, ypadding_top, snap=False)
        ymin = logy_base ** y0
        ymax = logy_base ** y3
    else:
        ymin, ymax = _limits_helper(
            ymin, ymax, ypadding_bottom, ypadding_top, snap=snap)

    return xmin, xmax, ymin, ymax


def get_band(low_hist, high_hist, middle_hist=None):
    """
    Convert the low and high histograms into a TGraphAsymmErrors centered at
    the middle histogram if not None otherwise the middle between the low and
    high points, to be used to draw a (possibly asymmetric) error band.
    """
    npoints = low_hist.nbins(0)
    band = Graph(npoints)
    for i in range(npoints):
        center = low_hist.x(i + 1)
        width = low_hist.xwidth(i + 1)
        low, high = low_hist.y(i + 1), high_hist.y(i + 1)
        if middle_hist is not None:
            middle = middle_hist.y(i + 1)
        else:
            middle = (low + high) / 2.
        yerrh = max(high - middle, low - middle, 0)
        yerrl = abs(min(high - middle, low - middle, 0))
        band.SetPoint(i, center, middle)
        band.SetPointError(i, width / 2., width / 2.,
                           yerrl, yerrh)
    return band


def canvases_with(drawable):
    """
    Return a list of all canvases where `drawable` has been painted.

    Note: This function is inefficient because it inspects all objects on all
          canvases, recursively. Avoid calling it if you have a large number of
          canvases and primitives.
    """
    return [c for c in ROOT.gROOT.GetListOfCanvases()
            if drawable in find_all_primitives(c)]


def find_all_primitives(pad):
    """
    Recursively find all primities on a pad, even those hiding behind a
    GetListOfFunctions() of a primitive
    """
    result = []
    for primitive in pad.GetListOfPrimitives():
        result.append(primitive)
        if hasattr(primitive, "GetListOfFunctions"):
            result.extend(primitive.GetListOfFunctions())
        if hasattr(primitive, "GetHistogram"):
            p = primitive.GetHistogram()
            if p:
                result.append(p)
        if isinstance(primitive, ROOT.TPad):
            result.extend(find_all_primitives(primitive))
    return result


def tick_length_pixels(pad, xaxis, yaxis, xlength, ylength=None):
    """
    Set the axes tick lengths in pixels
    """
    if ylength is None:
        ylength = xlength
    xaxis.SetTickLength(xlength / float(pad.height_pixels))
    yaxis.SetTickLength(ylength / float(pad.width_pixels))


class Compare(object):
    """
    Compare arbitrary number of plotables to each other.
    Performs any kind of user defined mathematical operation on arbitrary
    groups of input plotables (Hist, Profile, Graph).
    Currently supports only 1D objects and mathematical functions, which can be
    calculated bin-wise, i.e. no integration or differentiation.
    Function for the error calculation can be manually indicated as well or if
    left empty, assuming variables are uncorrelated it will be automatically
    calculated from main user defined function.

    Parameters
    ----------

    plotables : list or dict, mandatory
              A list or dict of input plotable objects. All of them should be
              of the same type (1D Hist, Profile or Graph) and should have same
              bin number and size. If list is being passed, first item will be
              assigned to the key x0, all others - to x1 and x1/x0 will be used
              as a comparisson function. If the user defined function is being
              desired, pass plotables as python dict of lists, where keys have
              a following form 'x0', 'x1', 'x2',... Later on this keys can be
              used to define function for the comparison and/or for error
              calculation.

    func : string, optional (default=None)
         A string defining function to use for the comparison. If list is
         being passed as plotables, it should be left empty and default x1/x0
         will be used, otherwise it's mandatory to indicate some function. User
         is free to choose any kind of mathematical operation, as long as it
         can be performed bin-wise on input plotables. To refer to some
         particular object from the input, use keys defined in plotables dict
         as variable names (x0,x1,x2...).

    errfunc : string, optional (default=None)
            A string defining function to use for the error calculation. If
            left empty, assuming there is no correlation between input
            plotables, it will be automatically calculated from func.
            If custom function is being desired user is free to choose any kind
            of mathematical operation, as long as it can be performed bin-wise
            on input plotables. To refer to some particular object from the
            input, use keys defined in plotables dict as variable names
            (x0,x1,x2...). To access the error of some particular object use
            the letter 'd' in front of the corresponding key name
            (dx0,dx1,dx2...).

    keytoclone: string, optional (default=None)
              A string indicating to the object, from which the style of the
              resulting plotable will be copied. It should be a key name
              defined in plotables ('x0','x1','x2'...). If left empty key 'xN'
              with highest N from plotables (N=0,1,2...) will be used by
              default.

    Returns
    -------
    The object to Compare class, from which various attributes and list of the
    resulting plotable objects can be accessed.

    Example
    -------
        result =  Compare({'x1':[h10,h11],'x2':[h20,h21,h22]}, '(x1-x2)/x1')

    This will calculate:
        (h10-h20)/h10, (h10-h21)/h10, (h10-h22)/h10, (h11-h20)/h11,
        (h11-h21)/h11 and (h11-h22)/h11

    To acces the results involving h10 in the calculation, use:
        result[h10]

    To further filter the results and get a list of objects involving h10 and
    h21 together in the calculation, use:
        result[h10, h21]

    """

    def __init__(self, plotables, func=None, errfunc=None, keytoclone=None):
        try:
            import sympy
        except ImportError:
            raise ImportError((
                    "`simpy` module can not be imported. Compare class depends"
                    " on SymPy library, please install it and try again..."))
        else:
            self.sympy = sympy

        self.plotables = {}
        if not isinstance(plotables,dict):
            if not isinstance(plotables,list):
                raise TypeError("Please pass at least a list of 2 plotables.")
            else:
                log.info("Assigning first plotable in the list to x0, all others to x1")
                self.plotables['x0'] = plotables[0]
                for item in plotables[1:]:
                    self.plotables['x1'] = item
                if func is not None or errfunc is not None:
                    raise RuntimeError(("Keys weren't assigned to plotables,"
                        " but function is explicitly provided.\n"
                        "Please either pass plotables as dict with key"
                        " names of form" " xN (N=0,1,2,...) or leave\n"
                        "func empty, in which case it will be set to a x1/x0."))
                else:
                    log.info("Setting compare function x1/x0")
                    log.warning(("If you want to set your own compare function"
                        " pass plotables\nas dict with key names of form"
                        " xN (N=0,1,2,...)."))
                    self.func = 'x1/x0'
        else:
            for key in plotables:
                if not isinstance(plotables[key],list):
                    self.plotables[key] = [plotables[key]]
                else:
                    self.plotables[key] = plotables[key]

            if not isinstance(func,str):
                raise TypeError(
                        "Plotables were assigned to keys but func is not provided at all or is not a string")
            else:
                self.func = func

        #Plotable objects consistency check with each other
        flat_list = list(itertools.chain.from_iterable(list(self.plotables.values())))
        for i, item in enumerate(flat_list):
            if i == 0:
                self.type = item.__class__.__name__
                if self.type not in ['Hist','Profile','Graph']:
                    raise TypeError(
                            "Plotables should be either rootpy Hist, Profile or Graph objects.")
            elif item.__class__.__name__ != self.type:
                raise TypeError(
                        "All plotables should be of same type either rootpy Hist, Profile or Graph.")
            elif (self.type is not 'Graph') and (flat_list[0].compatible(item, True) is False):
                raise ValueError(
                        "Number of bins and bin edges must be the same for every {0}".format(self.type))
            elif self.type is 'Graph':
                check = False
                if item.num_points == flat_list[0].num_points:
                    check = True
                    for g1, g2 in zip(item,flat_list[0]):
                        if g1.x.value != g2.x.value:
                            check = False
                            break
                        if g1.x.error != g2.x.error:
                            check = False
                            break
                if check is False:
                    raise TypeError(
                            "Number of points, exlow and exhigh must be the same for every Graph")
            if base.dim(item) != 1:
                raise TypeError(
                        "Currently we support comparisons only for 1D plotable objects.")
        #End of object consistency checks

        self.symfunc = self.sympy.sympify(self.func)

        for item in self.plotables:
            if re.match(r"x\d+", item) is None:
                raise KeyError("Key {0} must have a form of xN (N=0,1,2,...).".format(item))

        for item in self.symfunc.free_symbols:
            if (str(item) in self.plotables.keys()) is False:
                raise KeyError("{0} is being used in func, but not found in plotable keys.".format(str(item)))

        if (errfunc is not None) and (not isinstance(errfunc,str)):
            raise ValueError("errfunc must be of type string.")
        else:
            self.errfunc = errfunc

        if self.errfunc is not None:
            self.errfunc = self.sympy.sympify(self.errfunc)
            for item in self.errfunc.free_symbols:
                item = str(item)
                if item[0] == 'd':
                    item = item[1:]
                elif item[0] != 'x':
                    raise KeyError(("{0} has an incorrect form. Value from"
                        " plotable should be indicated as xN, error - as dxN (N=0,1,2,...)".format(item)))
                if (item in self.plotables.keys()) is False:
                    raise KeyError("{0} is being used in errfunc, but not found in plotable keys.".format(item))

        else:
            #Calculate standard error from func, assuming variables are independent
            self.errfunc = self.sympy.sympify(0)
            for item in self.symfunc.free_symbols:
                d = (self.sympy.diff(self.symfunc,item)*self.sympy.var('d'+str(item)))**2
                self.errfunc += d
            self.errfunc = self.errfunc**0.5

        #List of output plotables
        self.result = []
        #Resulting plotable will have the same style as input plotable with the
        #biggest N in xN keys of input plotables
        self.symkeys = self.symfunc.free_symbols
        self.keys = [int(str(key)[1:]) for key in self.symkeys]
        self.keys.sort()

        if keytoclone is None:
            if len(self.keys) > 0:
                keytoclone = 'x'+str(self.keys[-1])
            else:
                keytoclone = self.plotables.keys()[-1]
                log.warning(
                    "func has evaluated to constant, using style from {0} for resulting plotables".format(keytoclone))
        elif not isinstance(keytoclone,str):
            raise TypeError(
                    "keytoclone must be of type string.")
        elif keytoclone not in self.plotables.keys():
            raise TypeError(
                    "{0} is not found in the keys of input plotables.".format(keytoclone))
        self.keytoclone = keytoclone

        #Perform calculation
        self.prepare_output()

    @staticmethod
    def product_dict(**kwargs):
        """
        Create Cartesian product from entries of dict, while retaining the
        original key structure
        """
        keys = kwargs.keys()
        vals = kwargs.values()
        for instance in itertools.product(*vals):
            yield dict(zip(keys, instance))

    def prepare_output(self):
        """
        Determine all combination of input plotables and fill the output list
        """
        self.combinations = list(self.product_dict(**self.plotables))
        for comb in self.combinations:
            name = self.func
            for k in comb.keys():
                name = re.sub(k,comb[k].name,name)
            log.info("Performing calculation of {0}...".format(name))
            result = comb[self.keytoclone].clone()
            result.name = name
            #Add every plotable object in the given combination as a parents
            #attribute of the given resulting object
            result.parents = comb.values()
            #Perform actual evaluation of func and errorfunc
            self.evaluate(comb,result)
            #Save result in the output list
            self.result.append(result)


    def evaluate(self, comb, result):
        """
        Perform func and errfunc evaluation for given plotable combination and
        fill corresponding output plotable object
        """
        if self.type is not 'Graph':
            for i, ibin in enumerate(result):
                subs_func = {}
                subs_errfunc = {}
                for key in self.symkeys:
                    subs_func[key] = comb[str(key)][i].value
                for key in self.errfunc.free_symbols:
                    if str(key)[0] != 'd':
                        subs_errfunc[key] = comb[str(key)][i].value
                    else:
                        subs_errfunc[key] = comb[str(key)[1:]][i].error
                ibin.value = self.symfunc.evalf(subs=subs_func)
                ibin.error = self.errfunc.evalf(subs=subs_errfunc)
        else:
            for i, point in enumerate(result):
                subs_func = {}
                subs_errfunc_hi = {}
                subs_errfunc_low = {}
                for key in self.symkeys:
                    subs_func[key] = comb[str(key)][i].y.value
                for key in self.errfunc.free_symbols:
                    if str(key)[0] != 'd':
                        subs_errfunc_low[key] = comb[str(key)][i].y.value
                        subs_errfunc_hi[key] = subs_errfunc_low[key]
                    else:
                        subs_errfunc_low[key] = comb[str(key)[1:]][i].y.error_low
                        subs_errfunc_hi[key] = comb[str(key)[1:]][i].y.error_hi
                point.y.value = self.symfunc.evalf(subs=subs_func)
                point.y.error_low = self.errfunc.evalf(subs=subs_errfunc_low)
                point.y.error_hi = self.errfunc.evalf(subs=subs_errfunc_hi)

    def __getitem__(self, item):
        """
        Get list of all the resulting plotables for which h1 input was used
        using calls like:

            result[h1]

        Multiple input objects cab be passed as a selection as well:

            result[h1,h2,h3]
        """
        if isinstance(item,tuple) or isinstance(item,list):
            return [obj for obj in self.result if set(item).issubset(set(obj.parents))]
        return [obj for obj in self.result if item in obj.parents]
