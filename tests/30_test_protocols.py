#!/usr/bin/env python3

import dirty_water

def test_protocol_str_steps():
    protocol = dirty_water.Protocol()
    protocol += 'hello'
    protocol += 'world'

    assert str(protocol) == """\
1. hello

2. world"""

def test_protocol_iterable_steps():
    protocol = dirty_water.Protocol()
    protocol += 'hello', 'world'

    assert str(protocol) == """\
1. hello

2. world"""


def test_pcr_default_args():
    pcr = dirty_water.Pcr()
    assert pcr.num_reactions == 1
    assert pcr.extra_master_mix == 0
    assert pcr.template_in_master_mix == True
    assert pcr.primers_in_master_mix == False
    assert pcr.annealing_temp == 60
    assert pcr.extension_time == 120
    assert pcr.num_cycles == 35
    assert pcr.steps[0] == '''\
Setup 1 PCR reaction and 1 negative control:

Reagent             Conc  Each Rxn  Master Mix
──────────────────────────────────────────────
Q5 master mix         2x  25.00 μL    50.00 μL
primer mix           10x   5.00 μL            
template DNA   100 pg/μL   1.00 μL     2.00 μL
water                     19.00 μL    38.00 μL
──────────────────────────────────────────────
                          50.00 μL    45.00 μL/rxn'''
    assert pcr.steps[1] == '''\
Run the following thermocycler protocol:

98°C → 98°C → 60°C → 72°C → 72°C → 12°C
0:30   0:10   0:20   2:00   2:00    ∞
      └──────────────────┘
               35x'''

def test_pcr_num_reactions():
    pcr = dirty_water.Pcr()
    pcr.num_reactions = 2
    assert pcr.num_reactions == 2
    assert pcr.steps[0] == '''\
Setup 2 PCR reactions and 1 negative control:

Reagent             Conc  Each Rxn  Master Mix
──────────────────────────────────────────────
Q5 master mix         2x  25.00 μL    75.00 μL
primer mix           10x   5.00 μL            
template DNA   100 pg/μL   1.00 μL     3.00 μL
water                     19.00 μL    57.00 μL
──────────────────────────────────────────────
                          50.00 μL    45.00 μL/rxn'''

    pcr.num_reactions = 3
    assert pcr.num_reactions == 3
    assert pcr.steps[0] == '''\
Setup 3 PCR reactions and 1 negative control:

Reagent             Conc  Each Rxn  Master Mix
──────────────────────────────────────────────
Q5 master mix         2x  25.00 μL   100.00 μL
primer mix           10x   5.00 μL            
template DNA   100 pg/μL   1.00 μL     4.00 μL
water                     19.00 μL    76.00 μL
──────────────────────────────────────────────
                          50.00 μL    45.00 μL/rxn'''

def test_pcr_extra_master_mix():
    pcr = dirty_water.Pcr()
    pcr.extra_master_mix = 10
    assert pcr.extra_master_mix == 10
    assert pcr.steps[0] == '''\
Setup 1 PCR reaction and 1 negative control:

Reagent             Conc  Each Rxn  Master Mix
──────────────────────────────────────────────
Q5 master mix         2x  25.00 μL    55.00 μL
primer mix           10x   5.00 μL            
template DNA   100 pg/μL   1.00 μL     2.20 μL
water                     19.00 μL    41.80 μL
──────────────────────────────────────────────
                          50.00 μL    45.00 μL/rxn'''

def test_pcr_master_mix_components():
    pcr = dirty_water.Pcr()
    pcr.template_in_master_mix = False
    pcr.primers_in_master_mix = False
    assert not pcr.template_in_master_mix
    assert not pcr.primers_in_master_mix
    assert pcr.steps[0] == '''\
Setup 1 PCR reaction and 1 negative control:

Reagent             Conc  Each Rxn  Master Mix
──────────────────────────────────────────────
Q5 master mix         2x  25.00 μL    50.00 μL
primer mix           10x   5.00 μL            
template DNA   100 pg/μL   1.00 μL            
water                     19.00 μL    38.00 μL
──────────────────────────────────────────────
                          50.00 μL    44.00 μL/rxn'''

    pcr.template_in_master_mix = True
    pcr.primers_in_master_mix = False
    assert pcr.template_in_master_mix
    assert not pcr.primers_in_master_mix
    assert pcr.steps[0] == '''\
Setup 1 PCR reaction and 1 negative control:

Reagent             Conc  Each Rxn  Master Mix
──────────────────────────────────────────────
Q5 master mix         2x  25.00 μL    50.00 μL
primer mix           10x   5.00 μL            
template DNA   100 pg/μL   1.00 μL     2.00 μL
water                     19.00 μL    38.00 μL
──────────────────────────────────────────────
                          50.00 μL    45.00 μL/rxn'''

    pcr.template_in_master_mix = False
    pcr.primers_in_master_mix = True
    assert not pcr.template_in_master_mix
    assert pcr.primers_in_master_mix
    assert pcr.steps[0] == '''\
Setup 1 PCR reaction and 1 negative control:

Reagent             Conc  Each Rxn  Master Mix
──────────────────────────────────────────────
Q5 master mix         2x  25.00 μL    50.00 μL
primer mix           10x   5.00 μL    10.00 μL
template DNA   100 pg/μL   1.00 μL            
water                     19.00 μL    38.00 μL
──────────────────────────────────────────────
                          50.00 μL    49.00 μL/rxn'''

    pcr.template_in_master_mix = True
    pcr.primers_in_master_mix = True
    assert pcr.template_in_master_mix
    assert pcr.primers_in_master_mix
    assert pcr.steps[0] == '''\
Setup 1 PCR reaction and 1 negative control:

Reagent             Conc  Each Rxn  Master Mix
──────────────────────────────────────────────
Q5 master mix         2x  25.00 μL    50.00 μL
primer mix           10x   5.00 μL    10.00 μL
template DNA   100 pg/μL   1.00 μL     2.00 μL
water                     19.00 μL    38.00 μL
──────────────────────────────────────────────
                          50.00 μL    50.00 μL/rxn'''

def test_pcr_annealing_temp():
    pcr = dirty_water.Pcr()
    pcr.annealing_temp = 62
    assert pcr.steps[1] == '''\
Run the following thermocycler protocol:

98°C → 98°C → 62°C → 72°C → 72°C → 12°C
0:30   0:10   0:20   2:00   2:00    ∞
      └──────────────────┘
               35x'''

def test_pcr_extension_time():
    pcr = dirty_water.Pcr()
    pcr.extension_time = 59
    assert pcr.steps[1] == '''\
Run the following thermocycler protocol:

98°C → 98°C → 60°C → 72°C → 72°C → 12°C
0:30   0:10   0:20   0:59   2:00    ∞
      └──────────────────┘
               35x'''

    pcr.extension_time = 60
    assert pcr.steps[1] == '''\
Run the following thermocycler protocol:

98°C → 98°C → 60°C → 72°C → 72°C → 12°C
0:30   0:10   0:20   1:00   2:00    ∞
      └──────────────────┘
               35x'''

    pcr.extension_time = 61
    assert pcr.steps[1] == '''\
Run the following thermocycler protocol:

98°C → 98°C → 60°C → 72°C → 72°C → 12°C
0:30   0:10   0:20   1:01   2:00    ∞
      └──────────────────┘
               35x'''

def test_pcr_num_cycles():
    pcr = dirty_water.Pcr()
    pcr.num_cycles = 25
    assert pcr.steps[1] == '''\
Run the following thermocycler protocol:

98°C → 98°C → 60°C → 72°C → 72°C → 12°C
0:30   0:10   0:20   2:00   2:00    ∞
      └──────────────────┘
               25x'''

