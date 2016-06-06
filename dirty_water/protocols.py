#!/usr/bin/env python3

class Protocol:

    def __init__(self):
        self.steps = []

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


class Pcr:

    def __init__(self, num_reactions=1, ta=60, tx=120, nc=35):
        from .reaction import Reaction
        self.reaction = Reaction()
        self.num_reactions = num_reactions
        self.extra_master_mix = 0
        self.annealing_temp = ta
        self.extension_time = tx
        self.num_cycles = nc

        self.reaction['Q5 master mix'].std_volume = 25, 'μL'
        self.reaction['Q5 master mix'].std_stock_conc = '2x'
        self.reaction['Q5 master mix'].master_mix = True

        self.reaction['primer mix'].std_volume = 5, 'μL'
        self.reaction['primer mix'].std_stock_conc = '10x'
        self.reaction['primer mix'].master_mix = False

        self.reaction['template DNA'].std_volume = 1, 'μL'
        self.reaction['template DNA'].std_stock_conc = 100, 'pg/μL'
        self.reaction['template DNA'].master_mix = True

        self.reaction['water'].std_volume = 19
        self.reaction['water'].volume_unit = 'μL'
        self.reaction['water'].master_mix = True

    def __iter__(self):
        yield from self.steps

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

        setup_step = """\
Setup {self.num_reactions} PCR reaction{s} and 1 negative control:

{self.reaction}""".format(**locals())

        thermocycler_step = """\
Run the following thermocycler protocol:

98°C → 98°C → {ta}°C → 72°C → 72°C → 12°C
0:30   0:10   0:20   {tx}{pad}2:00    ∞
      └──────────────────┘
               {nc}x""".format(**locals())

        return [setup_step, thermocycler_step]

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

