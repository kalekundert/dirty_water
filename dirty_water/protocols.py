#!/usr/bin/env python3

class Protocol:

    def __init__(self):
        self.steps = []

    def __iter__(self):
        yield from self.steps

    def __iadd__(self, step_or_steps):
        from nonstdlib import MagicFormatter

        if isinstance(step_or_steps, str):
            step_or_steps = (step_or_steps,)

        for step in step_or_steps:
            self.steps.append(str(step) | MagicFormatter(level=2))

        return self

    def __str__(self):
        from textwrap import indent
        formatted_steps = []

        for i, step in enumerate(self.steps, 1):
            number = "{}. ".format(i)
            padding = ' ' * len(number)
            formatted_steps.append(
                    indent(number + step, ' ' * len(number)).strip())

        return '\n\n'.join(formatted_steps)


class Pcr(Protocol):

    def __init__(self, num_reactions=1, ta=60, tx=120, nc=35):
        from .reaction import Reaction
        self.reaction = Reaction()
        self.primer_mix = Reaction()
        self.num_reactions = num_reactions
        self.make_primer_mix = False
        self.extra_master_mix = 0
        self.annealing_temp = ta
        self.extension_time = tx
        self.num_cycles = nc

        self.primer_mix['water'].std_volume = 38, 'μL'
        self.primer_mix.show_master_mix = False

        self.primer_mix['forward primer'].std_volume = 1, 'μL'
        self.primer_mix['forward primer'].std_stock_conc = 200, 'μM'

        self.primer_mix['reverse primer'].std_volume = 1, 'μL'
        self.primer_mix['reverse primer'].std_stock_conc = 200, 'μM'

        self.reaction['water'].std_volume = 19, 'μL'
        self.reaction['water'].master_mix = True

        self.reaction['primer mix'].std_volume = 5, 'μL'
        self.reaction['primer mix'].std_stock_conc = '10x'
        self.reaction['primer mix'].master_mix = False

        self.reaction['template DNA'].std_volume = 1, 'μL'
        self.reaction['template DNA'].std_stock_conc = 100, 'pg/μL'
        self.reaction['template DNA'].master_mix = True

        self.reaction['Q5 master mix'].std_volume = 25, 'μL'
        self.reaction['Q5 master mix'].std_stock_conc = '2x'
        self.reaction['Q5 master mix'].master_mix = True

    def __getitem__(self, key):
        return self.reaction[key]

    @property
    def steps(self):
        s = 's' if self.num_reactions != 1 else ''
        ta = self.annealing_temp
        tx = self.extension_time
        tx = '{}:{:02d}'.format(tx // 60, tx % 60)
        pad = ' ' * (7 - len(tx))
        nc = self.num_cycles

        from .reaction import Reaction

        primer_step = """\
Prepare each 10x primer mix:

{}""".format(self.primer_mix)

        setup_step = """\
Setup {self.num_reactions} PCR reaction{s} and 1 negative control:

{self.reaction}""".format(**locals())

        thermocycler_step = """\
Run the following thermocycler protocol:

98°C → 98°C → {ta}°C → 72°C → 72°C → 12°C
0:30   0:10   0:20   {tx}{pad}2:00    ∞
      └──────────────────┘
               {nc}x""".format(**locals())

        if self.make_primer_mix:
            return primer_step, setup_step, thermocycler_step
        else:
            return setup_step, thermocycler_step

    @property
    def num_reactions(self):
        return self.reaction.num_reactions - 1

    @num_reactions.setter
    def num_reactions(self, value):
        self.reaction.num_reactions = value + 1

    @property
    def extra_master_mix(self):
        return self.reaction.extra_master_mix

    @extra_master_mix.setter
    def extra_master_mix(self, value):
        self.reaction.extra_master_mix = value

    @property
    def template_in_master_mix(self):
        return self.reaction['template DNA'].master_mix

    @template_in_master_mix.setter
    def template_in_master_mix(self, value):
        self.reaction['template DNA'].master_mix = value

    @property
    def primers_in_master_mix(self):
        return self.reaction['primer mix'].master_mix

    @primers_in_master_mix.setter
    def primers_in_master_mix(self, value):
        self.reaction['primer mix'].master_mix = value

    @primers_in_master_mix.setter
    def additives_in_master_mix(self, value):
        if 'DMSO' in self.reaction:
            self.reaction['DMSO'].master_mix = value
        if 'Betaine' in self.reaction:
            self.reaction['Betaine'].master_mix = value

    @property
    def dmso(self):
        return 'DMSO' in self.reaction

    @dmso.setter
    def dmso(self, dmso):
        if dmso:
            percent = 2 if dmso is True else dmso
            volume = 50 * percent / 100
            
            self.reaction['DMSO'].std_volume = volume, 'μL'
            self.reaction['DMSO'].std_stock_conc = '100%'
            self.reaction['DMSO'].master_mix = True
            self.reaction['water'].std_volume -= volume
        else:
            del self.reaction['DMSO']

    @property
    def betaine(self):
        return 'Betaine' in self.reaction

    @betaine.setter
    def betaine(self, betaine):
        if betaine:
            percent = 2 if betaine is True else betaine
            volume = 50 * percent / 100
            
            self.reaction['Betaine'].std_volume = volume, 'μL'
            self.reaction['Betaine'].std_stock_conc = '5M'
            self.reaction['Betaine'].master_mix = True
            self.reaction['water'].std_volume -= volume
        else:
            del self.reaction['Betaine']

