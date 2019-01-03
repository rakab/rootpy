import re
import sympy
import itertools

from .. import log; log = log[__name__]
from . import Hist, Graph

class Compare(object):
    """

    Compare arbitrary number of plottables

    """

    def __init__(self, plottables, func=None, errfunc=None):

        self.plottables = {}
        if not isinstance(plottables,dict):
            if not isinstance(plottables,list):
                raise TypeError("Please pass at least a list of 2 plottables.")
            else:
                log.info("Assigning first plottable in the list to x0, all others to x1")
                self.plottables['x0'] = plottables[0]
                for item in plottables[1:]:
                    self.plottables['x1'] = item
                if func is not None or errfunc is not None:
                    raise RuntimeError(("Keys weren't assigned to plottables,"
                        " but function is explicitly provided.\n"
                        "Please either pass plottables as dict with key"
                        " names of form" " xN (N=0,1,2,...) or leave\n"
                        "func empty, in which case it will be set to a x1/x0."))
                else:
                    log.info("Setting compare function x1/x0")
                    log.warning(("If you want to set your own compare function"
                        " pass plottables\nas dict with key names of form"
                        " xN (N=0,1,2,...)."))
                    self.func = 'x1/x0'
        else:
            self.plottables = plottables

            if not isinstance(func,str):
                raise TypeError("func must be of type string")
            else:
                self.func = func

        self.symfunc = sympy.sympify(self.func)

        for item in self.plottables:
            if re.match(r"x\d+", item) is None:
                raise KeyError("Key {0} must have a form of xN (N=0,1,2,...).".format(item))

        for item in self.symfunc.free_symbols:
            if (str(item) in self.plottables.keys()) is False:
                raise KeyError("{0} is being used in func, but not found in plottable keys.".format(str(item)))

        if (errfunc is not None) and (not isinstance(errfunc,str)):
            raise ValueError("errfunc must be of type string.")
        else:
            self.errfunc = errfunc

        if self.errfunc is not None:
            self.errfunc = sympy.sympify(self.errfunc)
            for item in self.errfunc.free_symbols:
                item = str(item)
                if item[0] == 'd':
                    item = item[1:]
                elif item[0] != 'x':
                    raise KeyError(("{0} has an incorrect form. Value from"
                        " plottable should be indicated as xN, error - as dxN (N=0,1,2,...)".format(item)))
                if (item in self.plottables.keys()) is False:
                    raise KeyError("{0} is being used in errfunc, but not found in plottable keys.".format(item))

        else:
            #Calculate standard error from func, assuming variables are independent
            self.errfunc = 0
            for item in self.symfunc.free_symbols:
                d = (sympy.diff(self.symfunc,item)*sympy.var('d'+str(item)))**2
                self.errfunc += d
            self.errfunc = self.errfunc**0.5

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
        Determine all combination of input plottables and fill the output dict
        """
        self.combinations = list(self.product_dict(**self.plottables))
        for comb in self.combinations:
            name = self.func
            for k in comb.keys():
                name = re.sub(k,str(comb[k]),name)
            log.info("Performing calculation of {0}...".format(name))
            self.evaluate_func()
            self.evaluate_errfunc()


    def evaluate_func(self):
        """
        Perform func evaluation and fill central values of output histos
        """
        pass

    def evaluate_errfunc(self):
        """
        Perform errfunc evaluation and fill errors of output histos
        """
        pass
