#!/usr/bin/python
# vim:fileencoding=utf-8
# -*- coding: utf-8 -*-

# Copyright 2013 by Carsten Richter and Robert Mietrach
# Contact: carsten.richter@desy.de and
#          carsten.richter@physik.tu-freiberg.de
#


def wrap_for_fit(func, argdict, varlist, unpack=True):
    """
    Wraps func so that it only expects a tuple containing the values for
    the variables specified in varlist. The rest of the variables is
    supplied from argdict.

    func -- original function
    argdict -- complete dictionary of func's arguments with default values
    varlist -- variables that that returned function is still expecting


    returns a (wrapped_function, start_values) pair

    start_values -- list of values from argdict for the keys in varlist
    wrapped_function -- function which takes tuple of values corresponding
                        to the variables in varlist
    
    unpack : bool
        when True, unpacks the dictionary of parameters when passing it
        to the cost function. used by default 
    """

    def helper(t):
        try: t = (t.item(),)
        except: pass
        helper_kwargs = dict(zip(varlist, t))
        func_kwargs = argdict.copy()
        func_kwargs.update(helper_kwargs)
        if unpack:
            return func(**func_kwargs)
        else:
            return func(func_kwargs)
    return helper, [argdict[key] for key in varlist]
