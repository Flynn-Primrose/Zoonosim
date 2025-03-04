'''
Defines the properties and states of individually modeled agents
'''
import sciris as sc

__all__ = ['IndividualMeta']

class individualMeta(sc.prettyobj):
    ''' For storing all the keys relating to individually modelled agents '''

    def __init__(self):

        # Set the properties of a person
        self.person = [
            'uid',              # Int
            'role',             # string? (i.e. owner, worker, inspector, visitor, none)
            'age',              # Float
            'sex',              # Float
            'symp_prob',        # Float
            'crit_prob',        # Float
            'death_prob',       # Float
            'abs_symp_prob',        # Float
            'abs_crit_prob',        # Float
            'abs_death_prob',       # Float
            'rel_trans',        # Float
            'rel_sus',          # Float
            'viral_load',       # Float
            'n_infections',     # Int
            'n_breakthroughs',  # Int
            'is_coinfected'
        ]

        # Set the states that a person can be in: these are all booleans per person -- used in people.py
        self.states = [
            'susceptible',
            'exposed',
            'infectious',
            'symptomatic',
            'critical',
            'tested',
            'diagnosed',
            'isolated',
            'recovered',
            'known_dead',
            'dead',
            'known_contact',
            'quarantined',
            'vaccinated',
        ]

        #Each field would be initialized as an matrix NxK where N is the number of pathogens in the simulation, K is the number of agents in the simulation  
        self.pathogen_states = [
            'p_susceptible',    #susceptible with a specific pathogen
            'p_naive',          #naive for a specific pathogen
            'p_exposed',        #exposed by a specific pathogen
            'p_infectious',     #infectious with a specific pathogen
            'p_symptomatic',    #symptomatic with a specific pathogen
            'p_critical',       #in critical condition caused by specific pathogen
            'p_tested',         #tested for a specific pathogen
            'p_diagnosed',      #diagnosed for a specific pathogen
            'p_recovered',      #recovered after infectious with a specific pathogen
            'p_dead',           #dead due to a specific pathogen
        ]
        
        #Information on variants of pathogens (variant states).
        self.pathogen_variants =[
            'p_exposed_variant',    #exposed to variant: matrix NxK, int
            'p_infectious_variant', #infectious with variant: matrix NxK, int
            'p_recovered_variant',  #recovered from variant: matrix NxK, int

            'p_exposed_by_variant',     #Records of exposed to a variant: matrix NxNkxK
            'p_infectious_by_variant',  #Records of infectiou with variant: matrix NxNkxK
]

         
        # Immune states, by pathogen by variant 
        self.imm_states = [
            'sus_imm',  # Float, by pathogen, by variant, by people,(Matrix NxNkxp, where N is the number of pathogens, and Nk the number of variants of N, whrere p is the population size)
            'symp_imm', # Float, by pathogen, by variant, by people,(Matrix NxNkxp, where N is the number of pathogens, and Nk the number of variants of N, whrere p is the population size)
        ]

        # Neutralizing antibody states, 1D array of length: number of pathogens in simulation
        self.nab_states = [
            'peak_nab',    # Float, peak neutralization titre relative to convalescent plasma
            'nab',         # Float, current neutralization titre relative to convalescent plasma
            't_nab_event', # Int, time since nab-conferring event
        ]

        self.imm_levels = [
            'imm_level',
            't_peak_imm',
            't_min_imm', 
            'curr_min',
            'decay_rate',
            'growth_rate',
            'imm_min',
            'imm_peak',
            'growth_start_date',
            'decay_start_date'
            ]

        # Additional vaccination states
        self.vacc_states = [
            'doses',          # Number of doses given per person
            'vaccine_source', # index of vaccine that individual received
            'vaccine_path', # index of pathogen vaccine is for
        ]

        # Set the dates various events took place: these are floats per person -- used in people.py
        self.dates = [f'date_{state}' for state in self.states] # Convert each state into a date
        self.dates.append('date_pos_test') # Store the date when a person tested which will come back positive 
        
        self.dates.append('date_end_quarantine') # Store the date when a person comes out of quarantine
        self.dates.append('date_end_isolation') # Store the date when a person comes out of isolation
        
        self.pathogen_dates = [f'date_{state}' for state in self.pathogen_states] # Convert each state into a date, arrays of NxP where N is num of pathogen and P is num of people
        self.pathogen_dates.append('p_date_pos_test')
        #self.dates.append('date_pos_test') # Store the date when a person tested which will come back positive
        #self.dates.append('date_end_quarantine') # Store the date when a person comes out of quarantine

        # Duration of different states: these are floats per person -- used in people.py
       
       #each field is a 1D array of length: number of pathogens in simulation
        self.durs = [
            'dur_exp2inf',
            'dur_inf2sym',
            'dur_sym2crit',
            'dur_disease',
        ]

        # Timing of control points for viral load; each field is a 1D array of length: number of pathogens in simulation
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

        self.population_sampling = [
            'provides_sample_prob', #float 
            'IgG_level' #float
        ]

        
        self.all_states = self.person + self.states + self.pathogen_states + self.pathogen_variants+ self.imm_states + self.nab_states + self.imm_levels+ self.vacc_states + self.dates+ self.pathogen_dates + self.durs + self.vl_points + self.population_sampling

        self.result_stocks = {
            'susceptible': 'Number susceptible',
            'exposed':     'Number exposed',
            'infectious':  'Number infectious',
            'symptomatic': 'Number symptomatic',
            'critical':    'Number of critical cases',
            'recovered':   'Number recovered',
            'dead':        'Number dead',
            'diagnosed':   'Number of confirmed cases',
            'isolated':    'Number in isolation',
            'known_dead':  'Number of confirmed deaths',
            'quarantined': 'Number in quarantine',
            'vaccinated':  'Number of people vaccinated',
        }

        self.result_stocks_by_variant = {
            'exposed_by_variant':    'Number exposed by variant',
            'infectious_by_variant': 'Number infectious by variant',
        }

        # The types of result that are counted as flows -- used in sim.py; value is the label suffix
        self.result_flows = {
            'infections':   'infections',
            'reinfections': 'reinfections',
            'infectious':   'infectious',
            'symptomatic':  'symptomatic cases',
            'critical':     'critical cases',
            'recoveries':   'recoveries',
            'deaths':       'deaths',
            'tests':        'tests',
            'diagnoses':    'diagnoses',
            'isolated':     'isolated people',
            'known_deaths': 'known deaths',
            'quarantined':  'quarantines started',
            'doses':        'vaccine doses',
            'vaccinated':   'vaccinated people'
        }

        self.result_flows_by_variant = {
            'infections_by_variant':  'infections by variant',
            'symptomatic_by_variant': 'symptomatic by variant',
            'infectious_by_variant':  'infectious by variant',
            'diagnoses_by_variant':   'diagnoses by variant'
        }

        # Define new and cumulative flows
        self.new_result_flows = [f'new_{key}' for key in self.result_flows.keys()]
        self.cum_result_flows = [f'cum_{key}' for key in self.result_flows.keys()]
        self.new_result_flows_by_variant = [f'new_{key}' for key in self.result_flows_by_variant.keys()]
        self.cum_result_flows_by_variant = [f'cum_{key}' for key in self.result_flows_by_variant.keys()]

        # Validate 
        self.state_types = ['person', 'states','pathogen_states', 'pathogen_variants', 'imm_states',
                            'nab_states', 'vacc_states', 'dates', 'pathogen_dates', 'durs', 'all_states', 'imm_levels', 'population_sampling']
        for state_type in self.state_types:
            states = getattr(self, state_type)
            n_states        = len(states)
            n_unique_states = len(set(states))
            if n_states != n_unique_states: # pragma: no cover
                errormsg = f'In {state_type}, only {n_unique_states} of {n_states} state names are unique'
                raise ValueError(errormsg)

        return
    
