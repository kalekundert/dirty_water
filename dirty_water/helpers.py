#!/usr/bin/env python3

def round_to_pipet(x):
    return '{:.2f}'.format(x)


class UserInputError (ValueError):
    pass
