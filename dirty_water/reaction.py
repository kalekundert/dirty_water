#!/usr/bin/env python3

import re, numbers, collections
from .helpers import round_to_pipet, UserInputError

class Reaction:

    def __init__(self, std_reagents=None):
        self.reagents = collections.OrderedDict()
        self._num_reactions = 1
        self._extra_master_mix = 10
        self.show_each_rxn = True
        self.show_master_mix = True
        self.show_totals = True

        if std_reagents:
            self.load_std_reagents(std_reagents)

    def __repr__(self):
        return 'Reaction()'

    def __getitem__(self, name):
        if name == 'master mix':
            return MasterMix(self)

        if name not in self.reagents:
            if name.lower() == 'water':
                self.reagents[name] = Water(self, name)
            else:
                self.reagents[name] = Reagent(self, name)

        return self.reagents[name]

    def __iter__(self):
        yield from self.reagents.values()

    def __str__(self):
        return self.reagent_table

    @property
    def num_reactions(self):
        return self._num_reactions

    @num_reactions.setter
    def num_reactions(self, num_reactions):
        self._num_reactions = int(num_reactions)

    @property
    def extra_master_mix(self):
        return self._extra_master_mix

    @extra_master_mix.setter
    def extra_master_mix(self, percent):
        self._extra_master_mix = float(percent)

    @property
    def scale(self):
        return self.num_reactions * (1 + self.extra_master_mix / 100)

    @property
    def volume(self):
        return sum(x.std_volume for x in self)

    @volume.setter
    def volume(self, volume):
        ratio = volume / self.volume
        for reagent in self:
            reagent.std_volume *= ratio

    @property
    def volume_unit(self):
        return next(iter(self)).volume_unit

    @property
    def volume_str(self):
        return '{} {}'.format(round_to_pipet(self.volume), self.volume_unit)

    @property
    def reagent_table(self):

        def cols(*cols):
            """
            Eliminate columns the user doesn't want to see.
            """
            cols = list(cols)
            if not self.show_master_mix:
                del cols[3]
            if not self.show_each_rxn:
                del cols[2]
            return cols

        # Figure out how big the table should be.

        column_titles = cols(
                "Reagent",
                "Conc",
                "Each Rxn",
                "Master Mix",
        )
        column_footers = cols(
                '',
                '',
                self.volume_str,
                self['master mix'].volume_str,
        )
        column_getters = cols(
                lambda x: x.name,
                lambda x: x.stock_conc_str,
                lambda x: x.volume_str,
                lambda x: x.scaled_volume_str,
        )
        column_widths = [
                max(map(len,
                    [title, footer] + [getter(x) for x in self]))
                for title, footer,getter in 
                    zip(column_titles, column_footers, column_getters)
        ]
        column_alignments = '<>>>'
        row_template = '  '.join(
                '{{:{}{}}}'.format(column_alignments[i], column_widths[i])
                for i in range(len(column_titles))
        )

        # Print the table
        
        rule = '─' * (sum(column_widths) + 2 * len(column_widths) - 2)
        rows = [
            row_template.format(*column_titles),
            rule,
        ] + [
            row_template.format(
                *[getter(reagent) for getter in column_getters])
            for reagent in self
        ]
        if self.show_totals and (self.show_each_rxn or self.show_master_mix):
            rows += [
                rule,
                row_template.format(*column_footers),
            ]
        return '\n'.join(rows) + '/rxn'

    def load_std_reagents(self, std_reagents):
        lines = std_reagents.strip().split('\n')

        # Examine the second line of the table to determine where each column 
        # starts and stops.

        column_slices = [
                slice(x.start(), x.end())
                for x in re.finditer('[-=]+', lines[1])
        ]
        if len(column_slices) != 4:
            raise UserInputError("Expected to find 4 columns, delineated by '=' or '-' in the second line.")

        def split_columns(line):
            return tuple(line[x].strip() for x in column_slices)

        # Parse standard concentrations and volumes from each line in the 
        # table.

        def parse_amount(x, unit_required=True):
            try:
                amount, unit = x.split()
                return float(amount), unit
            except ValueError:
                if unit_required: raise
                else: return x or None

        for line in lines[2:]:
            reagent, stock_conc, volume, master_mix = split_columns(line)
            self[reagent].std_stock_conc = parse_amount(stock_conc, False)
            self[reagent].std_volume = parse_amount(volume)
            self[reagent].master_mix = (master_mix == 'yes')


class Reagent:

    def __init__(self, reaction, name):
        self.reaction = reaction
        self.name = name
        self._std_volume = None
        self._std_stock_conc = None
        self.volume_unit = None
        self.conc_unit = None
        self.master_mix = False
        self._stock_conc = None
        self._conc = None

    def __repr__(self):
        return 'Reagent({0.name})'.format(self)

    def __str__(self):
        return '{0.volume_str} {0.name}'.format(self)

    @property
    def std_stock_conc(self):
        return self._std_stock_conc

    @std_stock_conc.setter
    def std_stock_conc(self, conc):
        if isinstance(conc, tuple):
            self._std_stock_conc, self.conc_unit = conc
        else:
            self._std_stock_conc = conc

    @property
    def std_conc(self):
        return self.std_stock_conc * self.std_volume / self.reaction.volume

    @property
    def std_volume(self):
        return self._std_volume

    @std_volume.setter
    def std_volume(self, volume):
        if isinstance(volume, tuple):
            self._std_volume, self.volume_unit = volume
        else:
            self._std_volume = float(volume)

    @property
    def stock_conc(self):
        return self._stock_conc or self.std_stock_conc

    @stock_conc.setter
    def stock_conc(self, stock_conc):
        if not isinstance(self.std_stock_conc, numbers.Real):
            raise UserInputError("The 'std_stock_conc' for '{}' isn't numeric, so 'stock_conc' can't be changed.".format(self.name))
        self._stock_conc = float(stock_conc)

    @property
    def stock_conc_str(self):
        if self.stock_conc and self.conc_unit:
            return '{0.stock_conc:.0f} {0.conc_unit}'.format(self)
        elif self.stock_conc:
            return '{0.stock_conc}'.format(self)
        else:
            return ''

    @property
    def conc(self):
        return self._conc or self.std_conc

    @conc.setter
    def conc(self, conc):
        if not isinstance(self.std_stock_conc, numbers.Real):
            raise UserInputError("The 'std_stock_conc' for '{}' isn't numeric, so 'conc' can't be changed.".format(self.name))
        self._conc = float(conc)

    @property
    def volume(self):
        # It isn't possible to calculate std_conc for reagents for which 
        # std_stock_conc is a string (i.e. '10x') or left undefined (i.e.  
        # water).  So don't calculate a new volume if std_volume will do.
        if self._conc is None:
            return self.std_volume
        else:
            return self.reaction.volume * self.conc / self.stock_conc

    @property
    def volume_str(self):
        return '{} {}'.format(round_to_pipet(self.volume), self.volume_unit)

    @property
    def scaled_volume_str(self):
        if not self.master_mix:
            return ''
        else:
            return '{} {}'.format(
                    round_to_pipet(self.volume * self.reaction.scale),
                    self.volume_unit)


class Water (Reagent):

    @property
    def volume(self):
        volume = self.reaction.volume

        for reagent in self.reaction:
            if reagent is not self:
                volume -= reagent.volume

        if volume < 0:
            from warnings import warn
            warn("Reaction volume exceeds {}".format(self.reaction.volume_str))

        return volume


class MasterMix (Reagent):

    def __init__(self, reaction):
        super().__init__(reaction, 'master mix')
        self.master_mix = True

    def __iter__(self):
        yield from (x for x in self.reaction if x.master_mix)

    @property
    def volume(self):
        return sum(x.volume for x in self)

    @property
    def volume_unit(self):
        return self.reaction.volume_unit

    @volume_unit.setter
    def volume_unit(self, unit):
        pass



