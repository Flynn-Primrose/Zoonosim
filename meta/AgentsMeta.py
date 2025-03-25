'''
Keys that are common amongst all agent types.
'''

import sciris as sc

class AgentsMeta(sc.prettyobj):
    '''
    For storing all keys that are common across all agent types
    '''

    def __init__(self):

        self.agent = [
            'uid', #Int
            'agent_type', #string? or int?
        ]

        self.states = [ # all boolean
            'susceptible', #?
            'exposed', #?
            'infectious', #?

        ]

        self.variant_states =[
            'exposed_variant', # int
            'infectious_variant' # int
        ]

        self.by_variant_states = [
            'exposed_by_variant',
            'infectious_by_variant',
        ]

        self.dates = [f'date_{state}' for state in self.states]

        self.all_states = self.agent + self.states +self.variant_states + self.by_variant_states + self.dates

        # Validate
        self.state_types = ['agent', 'states', 'variant_states', 'by_variant_states', 'dates', 'all_states']
        for state_type in self.state_types:
            states = getattr(self, state_type)
            n_states        = len(states)
            n_unique_states = len(set(states))
            if n_states != n_unique_states: # pragma: no cover
                errormsg = f'In {state_type}, only {n_unique_states} of {n_states} state names are unique'
                raise ValueError(errormsg)

        return