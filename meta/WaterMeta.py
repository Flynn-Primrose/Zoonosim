'''

'''

import sciris as sc

class WaterMeta(sc.prettyobj):
    ''' Defines all keys that are used by barns '''

    def __init__(self):
        
        self.agent = [
            'uid', # int
            'temperature'
        ]

        self.states = [

        ]

        self.variant_states = [

        ]

        self.by_variant_states = [

        ]

                # Validate
        self.state_types = ['egg', 'states', 'variant_states', 'by_variant_states']
        for state_type in self.state_types:
            states = getattr(self, state_type)
            n_states        = len(states)
            n_unique_states = len(set(states))
            if n_states != n_unique_states: # pragma: no cover
                errormsg = f'In {state_type}, only {n_unique_states} of {n_states} state names are unique'
                raise ValueError(errormsg)

        return