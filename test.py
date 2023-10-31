"""
Tests for passing Python callables when constructing TFX classes.

This feature is not implemented by a PyROOT pythonization, but by a converter of
Cppyy that creates a C++ wrapper to invoke the Python callable.
"""

import unittest
import math

import ROOT


def pyf_tf1_identity(x, p):
    return x[0]


def pyf_tf1_params(x, p):
    return p[0] * x[0] + p[1]


class pyf_tf1_callable:
    def __call__(self, x, p):
        return p[0] * x[0] + p[1]


def pyf_tf1_gauss(x, p):
    return p[0] * 1.0 / math.sqrt(2.0 * math.pi * p[2]**2) * math.exp(-(x[0] - p[1])**2 / 2.0 / p[2]**2)


f = ROOT.TF1("tf1_identity", pyf_tf1_identity, 0.0, 1.0)


