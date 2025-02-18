import sciris as sc

__all__ = ['FlocskMeta']

class FlocksMeta(sc.prettyobj):
    ''' For storing all the keys relating to a person and people '''

    def __init__(self):

        # Set the properties of a flock of birds
        self.flock = [
            'uid',              # Int
            'role',             # string? (i.e. broiler, layer, breeder)
            'age',              # Float (not sure if we need this)
            'head',             # Int (headcount of flock)
            'susceptible',      # Int (number of susceptible birds)
            'exposed',          # Int (number of exposed birds)
            'infectious',       # Int (number of infectious birds)
        ]

        # Set the states that a flock can be in
        self.states = [
            'infectious',
            'symptomatic',
            'culled'
        ]

        #Each field would be initialized as an matrix NxK where N is the number of pathogens in the simulation, K is the number of flocks in the simulation  
        self.pathogen_states = [
            'p_infectious',     #infectious with a specific pathogen
        ]
        
        #Information on variants of pathogens (variant states).
        self.pathogen_variants =[
            'p_infectious_variant', #infectious with variant: matrix NxK, int
        ]

        # Set the dates various events took place: these are floats per flock
        self.dates = [f'date_{state}' for state in self.states] # Convert each state into a date
        
        self.pathogen_dates = [f'date_{state}' for state in self.pathogen_states] # Convert each state into a date, arrays of NxP where N is num of pathogen and P is num of people
  
        self.all_states = self.flock + self.states + self.pathogen_states + self.pathogen_variants + self.dates+ self.pathogen_dates

        self.result_stocks = {
            'infectious':  'Number infectious',
            'symptomatic': 'Number symptomatic',
            'culled':      'Number culled',
        }

        self.result_stocks_by_variant = {
            'infectious_by_variant': 'Number infectious by variant',
            'symptomatic_by_variant': 'Number symptomatic by variant',
            'culled_by_variant': 'Number culled by variant',
        }

        # The types of result that are counted as flows -- used in sim.py; value is the label suffix
        self.result_flows = {
            'infectious':   'infectious',
            'symptomatic':  'symptomatic',
            'culled':       'culled',
        }

        self.result_flows_by_variant = {
            'infectious_by_variant':  'infectious by variant',
            'symptomatic_by_variant': 'symptomatic by variant',
            'culled_by_variant':      'culled by variant',
        }

        # Define new and cumulative flows
        self.new_result_flows = [f'new_{key}' for key in self.result_flows.keys()]
        self.cum_result_flows = [f'cum_{key}' for key in self.result_flows.keys()]
        self.new_result_flows_by_variant = [f'new_{key}' for key in self.result_flows_by_variant.keys()]
        self.cum_result_flows_by_variant = [f'cum_{key}' for key in self.result_flows_by_variant.keys()]

        # Validate 
        self.state_types = ['flock', 'states','pathogen_states', 'pathogen_variants', 'dates', 'pathogen_dates', 'all_states']
        for state_type in self.state_types:
            states = getattr(self, state_type)
            n_states        = len(states)
            n_unique_states = len(set(states))
            if n_states != n_unique_states: # pragma: no cover
                errormsg = f'In {state_type}, only {n_unique_states} of {n_states} state names are unique'
                raise ValueError(errormsg)

        return