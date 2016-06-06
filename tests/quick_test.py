#!/usr/bin/env python3

import dirty_water

rxn = dirty_water.Reaction()
rxn.num_reactions = 10

rxn['Cas9'].std_volume = 1, 'μL'
rxn['Cas9'].std_stock_conc = 20, 'μM'
rxn['Cas9'].master_mix = True
rxn['Cas9'].conc = 4

rxn['buffer'].std_volume = 1, 'μL'
rxn['buffer'].std_stock_conc = '10x'
rxn['buffer'].master_mix = True

rxn['sgRNA'].std_volume = 1, 'μL'
rxn['sgRNA'].std_stock_conc = 32, 'μM'

rxn['trypsin'].std_volume = 1, 'μL'
rxn['trypsin'].std_stock_conc = '10x'

rxn['water'].std_volume = 6, 'μL'
rxn['water'].master_mix = True

print(rxn)
