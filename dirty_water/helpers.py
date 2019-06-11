#!/usr/bin/env python3

import math

def round_to_pipet(x):
    return '{:.2f}'.format(x)

def with_extra(x, extra=0.1):
    return math.ceil(x + x * extra)

class UserInputError (ValueError):
    pass
