#!/usr/bin/env python3

from dirty_water.helpers import round_to_pipet

def test_round_to_pipet():
    assert round_to_pipet(1) == '1.00'
    assert round_to_pipet(1.2) == '1.20'
    assert round_to_pipet(1.23) == '1.23'
    assert round_to_pipet(1.234) == '1.23'

    assert round_to_pipet(12) == '12.00'
    assert round_to_pipet(12.3) == '12.30'
    assert round_to_pipet(12.34) == '12.34'
    assert round_to_pipet(12.345) == '12.35'

    assert round_to_pipet(123) == '123.00'
    assert round_to_pipet(123.4) == '123.40'
    assert round_to_pipet(123.45) == '123.45'
    assert round_to_pipet(123.456) == '123.46'
