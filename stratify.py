import numpy as np
import operator

def get_indices_to_track(sim, indices):
    '''
        Returns the list of indices of the people to track results for 
        (satisfying the stratification constraints set by passing a constraint dict to the Sim object, which is parsed by stratify.set_stratified_indices)
    '''
    if not sim.pars['enable_stratifications']:  
        return indices
    if indices is None:
        return sim.stratification_indices
    return np.intersect1d(indices, sim.stratification_indices)

def set_stratified_indices(indices, constraint, sim):
    """
        Function to that parses constraint and chooses which indices to track results for results in result_keys_to_stratify, 
        based on if a person's attributes satisfy the constraint condition 

        Parameters: 
        sim(Sim) = A reference to the simulation class.
        indices(array of int) = indices of the list of people to track a certain results prior to filtering which people have characteristics which satisfy the stratification constraints.

         
        Example to stratify some results: 
             
           constraint = {
                'age': [('>', 25), ('<', 50)],
                'has_watch': [('==', True)],
                'income': [('>=', 10000)]
            }

            sim = zoonosim.Sim(pars, enable_stratifications = True, stratification_pars = constraint)

    """  
   
    # Check for valid constraints
    operators = {'>': operator.gt, '<': operator.lt, '>=': operator.ge, '<=': operator.le, '==': operator.eq, '=': operator.eq,'!=': operator.ne}
     
    mask = np.copy(indices)
    
    for attr, const_list in constraint.items():
        for (operator_str, value) in const_list:

            if operator_str not in operators:
                raise ValueError(f"Invalid operator: {operator_str}")
            if not attr in sim.people.keys():
                raise ValueError(f"attribute {attr} does not exist in people object") 

            indices_with_cond = (operators[operator_str](sim.people[attr], value))

            for ind, i in enumerate(indices_with_cond): 
                if i == True:
                    if not (operators[operator_str](sim.people[attr][ind], value)): 
                        raise Exception()

            if isinstance(indices_with_cond, np.ndarray):
                nonzero_indices = indices_with_cond.nonzero()[0]
                mask = np.intersect1d(mask, nonzero_indices)
            else: 
                mask = np.intersect1d(mask, [])
    
     
    # Apply mask to results and return indices satisfying condition
    return indices[mask] 