'''

'''
import sciris as sc

class HumanMeta(sc.prettyobj):
    ''' Defines all keys used for human agents '''

    def __init__(self):

        # Set the properties of a person
        self.agent = [
            'uid',              # Int
            'age',              # Float
            'sex',              # Float
            'symp_prob',        # Float
            'severe_prob',      # Float
            'death_prob',       # Float
            'rel_trans',        # Float
            'rel_sus',          # Float
            'viral_load',       # Float
            'rescaled_vl',      # Float
            'n_infections',     # Int
            'n_breakthroughs',  # Int
        ]

        # Set the states that a person can be in: these are all booleans per person -- used in people.py
        self.states = [
            'susceptible',
            'naive',
            'exposed',
            'infectious',
            'symptomatic',
            'severe',
            'tested',
            'diagnosed',
            'recovered',
            'known_dead',
            'dead',
            'quarantined',
            'vaccinated',
        ]

        # Variant states -- these are ints
        self.variant_states = [
            'exposed_variant',
            'infectious_variant',
            'recovered_variant',
        ]

        # Variant states -- these are ints, by variant
        self.by_variant_states = [
            'exposed_by_variant',
            'infectious_by_variant',
        ]

        # Immune states, by variant
        self.imm_states = [
            'sus_imm',  # Float, by variant
            'symp_imm', # Float, by variant
            'sev_imm',  # Float, by variant
        ]

        # Neutralizing antibody states
        self.nab_states = [
            'peak_nab',    # Float, peak neutralization titre relative to convalescent plasma
            'nab',         # Float, current neutralization titre relative to convalescent plasma
            't_nab_event', # Int, time since nab-conferring event
        ]

        # Additional vaccination states
        self.vacc_states = [
            'doses',          # Number of doses given per person
            'vaccine_source', # index of vaccine that individual received
        ]

        # Set the dates various events took place: these are floats per person -- used in people.py
        self.dates = [f'date_{state}' for state in self.states] # Convert each state into a date
        self.dates.append('date_pos_test') # Store the date when a person tested which will come back positive
        self.dates.append('date_end_quarantine') # Store the date when a person comes out of quarantine
        # Duration of different states: these are floats per person -- used in people.py
        self.durs = [
            'dur_exp2inf',
            'dur_inf2sym',
            'dur_sym2sev',
            'dur_sev2crit',
            'dur_disease',
        ]

        # Timing of control points for viral load
        self.vl_points = [
            'x_p_inf',
            'y_p_inf',
            'x_p1',
            'y_p1',
            'x_p2',
            'y_p2',
            'x_p3',
            'y_p3',
        ]

        self.all_states = self.person + self.states + self.variant_states + self.by_variant_states + self.imm_states + self.nab_states + self.vacc_states + self.dates + self.durs + self.vl_points

        # Validate
        self.state_types = ['person', 'states', 'variant_states', 'by_variant_states', 'imm_states',
                            'nab_states', 'vacc_states', 'dates', 'durs', 'all_states']
        for state_type in self.state_types:
            states = getattr(self, state_type)
            n_states        = len(states)
            n_unique_states = len(set(states))
            if n_states != n_unique_states: # pragma: no cover
                errormsg = f'In {state_type}, only {n_unique_states} of {n_states} state names are unique'
                raise ValueError(errormsg)

        return