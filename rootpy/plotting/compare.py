import re
import sympy
import itertools

from .. import log; log = log[__name__]
from . import base, Hist, Profile, Graph

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

        self.symfunc = sympy.sympify(self.func)

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
            self.errfunc = sympy.sympify(self.errfunc)
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
            self.errfunc = sympy.sympify(0)
            for item in self.symfunc.free_symbols:
                d = (sympy.diff(self.symfunc,item)*sympy.var('d'+str(item)))**2
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
