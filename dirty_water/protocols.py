#!/usr/bin/env python3

from nonstdlib import plural

class Protocol:
    """
    Typical usage:

    >>> p = dirty_water.Protocol()
    >>> p += "First step..."
    >>> p += "Second step..."
    >>> print(p)
    """

    def __init__(self):
        self.steps = []
        self.notes = ""

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

    thermocycler_protocols = {
            # This needs to be in a config file somewhere...
            'q5': {
                'initial_denature_temp': 98,
                'initial_denature_time': 30,
                'denature_temp': 98,
                'denature_time': 10,
                'anneal_temp': 60,
                'anneal_time': 20,
                'extend_temp': 72,
                'extend_time': 120,
                'final_extend_temp': 72,
                'final_extend_time': 120,
                'num_cycles': 35,
                'hold': 4,
            },
            'ssoadv': {
                'initial_denature_temp': 95,
                'initial_denature_time': 30,
                'denature_temp': 95,
                'denature_time': 10,
                'anneal_temp': 60,
                'anneal_time': 15,
                'melt_curve_low_temp': 65,
                'melt_curve_high_temp': 95,
                'melt_curve_temp_step': 0.5,
                'melt_curve_time_step': 5,
                'num_cycles': 40,
                'two_step': True,
            },
    }
    thermocycler_aliases = {
            'ta': 'annealing_temp',
            'tx': 'extension_time',
            'nc': 'num_cycles',
    }
    polymerase_reagents = {
            'q5': {
                'name': 'Q5 master mix',
                'std_volume': (25, 'µL'),
                'std_stock_conc': '2x',
            },
            'ssoadv': {
                'name': 'SsoAdvanced mix',
                'std_volume': (25, 'µL'),
                'std_stock_conc': '2x',
            },
    }

    def __init__(self, polymerase='q5', num_reactions=1, **thermocycler_kwargs):
        from .reaction import Reaction
        self.reaction = Reaction()
        self.primer_mix = Reaction()
        self.thermocycler_protocol = self.thermocycler_protocols[polymerase]
        self.num_reactions = num_reactions
        self.make_primer_mix = False
        self.extra_master_mix = 0

        self.thermocycler_protocol.update({
                thermocycler_aliases.get(k, k): v
                for k, v in thermocycler_kwargs.items()
        })

        self.primer_mix['water'].std_volume = 36, 'μL'
        self.primer_mix.show_master_mix = False

        self.primer_mix['forward primer'].std_volume = 2, 'μL'
        self.primer_mix['forward primer'].std_stock_conc = 100, 'μM'

        self.primer_mix['reverse primer'].std_volume = 2, 'μL'
        self.primer_mix['reverse primer'].std_stock_conc = 100, 'μM'

        polymerase_reagent = self.polymerase_reagents[polymerase]

        self.reaction['water'].std_volume = 44 - polymerase_reagent['std_volume'][0], 'μL'
        self.reaction['water'].master_mix = True

        self.reaction['primer mix'].std_volume = 5, 'μL'
        self.reaction['primer mix'].std_stock_conc = '10x'
        self.reaction['primer mix'].master_mix = False

        self.reaction['template DNA'].std_volume = 1, 'μL'
        self.reaction['template DNA'].std_stock_conc = 100, 'pg/μL'
        self.reaction['template DNA'].master_mix = True

        self.reaction[polymerase_reagent['name']].std_volume = polymerase_reagent['std_volume']
        self.reaction[polymerase_reagent['name']].std_stock_conc = polymerase_reagent['std_stock_conc']
        self.reaction[polymerase_reagent['name']].master_mix = True

    def __getitem__(self, key):
        return self.reaction[key]

    @property
    def steps(self):
        primer_step = f"""\
Prepare each 10x primer mix:

{self.primer_mix}"""

        setup_step = f"""\
Setup {self.num_reactions} PCR {plural(self.num_reactions):reaction/s} and 1 negative control:

{self.reaction}"""

        def time(x):
            if x < 60:
                return f'{x}s'
            elif x % 60:
                return f'{x//60}m{x%60:02}'
            else:
                return f'{x//60} min'

        def has_step(protocol, step, params=['temp', 'time']):
            return all((f'{step}_{param}' in protocol) for param in params)

        p = self.thermocycler_protocol
        three_step = not p.get('two_step', False) and has_step(p, 'extend')

        thermocycler_steps = [
                f"- {p['initial_denature_temp']}°C for {time(p['initial_denature_time'])}",
                f"- Repeat {p['num_cycles']}x:",
                f"  - {p['denature_temp']}°C for {time(p['denature_time'])}",
                f"  - {p['anneal_temp']}°C for {time(p['anneal_time'])}",
        ]
        if three_step:
            thermocycler_steps += [
                f"  - {p['extend_temp']}°C for {time(p['extend_time'])}",
            ]

        if has_step(p, 'final_extend'):
            thermocycler_steps += [
                f"- {p['final_extend_temp']}°C for {time(p['final_extend_time'])}",
            ]

        if has_step(p, 'melt_curve', 'low_temp high_temp temp_step time_step'.split()):
            thermocycler_steps += [
                f"- {p['melt_curve_low_temp']}-{p['melt_curve_high_temp']}°C in {time(p['melt_curve_time_step'])} steps of {p['melt_curve_temp_step']}°C",
            ]

        if 'hold' in p:
            thermocycler_steps += [
                f"- {p['hold']}°C hold",
            ]

        br = '\n'
        thermocycler_step = f"""\
Run the following thermocycler protocol:

{br.join(thermocycler_steps)}"""

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
    def annealing_temp(self):
        return self.thermocycler_protocol['anneal_temp']

    @annealing_temp.setter
    def annealing_temp(self, value):
        self.thermocycler_protocol['anneal_temp'] = value

    @property
    def extension_time(self):
        return self.thermocycler_protocol['extend_time']

    @extension_time.setter
    def extension_time(self, value):
        self.thermocycler_protocol['extend_time'] = value

    @property
    def num_cycles(self):
        return self.thermocycler_protocol['extend_time']

    @num_cycles.setter
    def num_cycles(self, value):
        self.thermocycler_protocol['num_cycles'] = value

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

