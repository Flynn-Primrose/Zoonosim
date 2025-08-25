'''
Additional analysis functions that are not part of the core Zoonosim workflow,
but which are useful for particular investigations.
'''

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from . import interventions as zni
from . import defaults as znd
from . import plotting as znplt
from . import options as zno


__all__ = ['Analyzer', 'snapshot', 'biography']


class Analyzer(sc.prettyobj):
    '''
    Base class for analyzers. Based on the Intervention class. Analyzers are used
    to provide more detailed information about a simulation than is available by
    default -- for example, pulling states out of sim.agents on a particular timestep
    before it gets updated in the next timestep.

    To retrieve a particular analyzer from a sim, use sim.get_analyzer().

    Args:
        label (str): a label for the Analyzer (used for ease of identification)
    '''

    def __init__(self, label=None):
        if label is None:
            label = self.__class__.__name__ # Use the class name if no label is supplied
        self.label = label # e.g. "Record ages"
        self.initialized = False
        self.finalized = False
        return


    def __call__(self, *args, **kwargs):
        # Makes Analyzer(sim) equivalent to Analyzer.apply(sim)
        if not self.initialized:
            errormsg = f'Analyzer (label={self.label}, {type(self)}) has not been initialized'
            raise RuntimeError(errormsg)
        return self.apply(*args, **kwargs)


    def initialize(self, sim=None):
        '''
        Initialize the analyzer, e.g. convert date strings to integers.
        '''
        self.initialized = True
        self.finalized = False
        return


    def finalize(self, sim=None):
        '''
        Finalize analyzer

        This method is run once as part of `sim.finalize()` enabling the analyzer to perform any
        final operations after the simulation is complete (e.g. rescaling)
        '''
        if self.finalized:
            raise RuntimeError('Analyzer already finalized')  # Raise an error because finalizing multiple times has a high probability of producing incorrect results e.g. applying rescale factors twice
        self.finalized = True
        return


    def apply(self, sim):
        '''
        Apply analyzer at each time point. The analyzer has full access to the
        sim object, and typically stores data/results in itself. This is the core
        method which each analyzer object needs to implement.

        Args:
            sim: the Sim instance
        '''
        raise NotImplementedError


    def shrink(self, in_place=False):
        '''
        Remove any excess stored data from the intervention; for use with sim.shrink().

        Args:
            in_place (bool): whether to shrink the intervention (else shrink a copy)
        '''
        if in_place:
            return self
        else:
            return sc.dcp(self)


    def to_json(self):
        '''
        Return JSON-compatible representation

        Custom classes can't be directly represented in JSON. This method is a
        one-way export to produce a JSON-compatible representation of the
        intervention. This method will attempt to JSONify each attribute of the
        intervention, skipping any that fail.

        Returns:
            JSON-serializable representation
        '''
        # Set the name
        json = {}
        json['analyzer_name'] = self.label if hasattr(self, 'label') else None
        json['analyzer_class'] = self.__class__.__name__

        # Loop over the attributes and try to process
        attrs = self.__dict__.keys()
        for attr in attrs:
            try:
                data = getattr(self, attr)
                try:
                    attjson = sc.jsonify(data)
                    json[attr] = attjson
                except Exception as E:
                    json[attr] = f'Could not jsonify "{attr}" ({type(data)}): "{str(E)}"'
            except Exception as E2:
                json[attr] = f'Could not jsonify "{attr}": "{str(E2)}"'
        return json


def validate_recorded_dates(sim, requested_dates, recorded_dates, die=True):
    '''
    Helper method to ensure that dates recorded by an analyzer match the ones
    requested.
    '''
    requested_dates = sorted(list(requested_dates))
    recorded_dates = sorted(list(recorded_dates))
    if recorded_dates != requested_dates: # pragma: no cover
        errormsg = f'The dates {requested_dates} were requested but only {recorded_dates} were recorded: please check the dates fall between {sim.date(sim["start_day"])} and {sim.date(sim["start_day"])} and the sim was actually run'
        if die:
            raise RuntimeError(errormsg)
        else:
            print(errormsg)
    return

class snapshot(Analyzer):
    '''
    Analyzer that takes a "snapshot" of the sim.agents array at specified points
    in time, and saves them to itself. To retrieve them, you can either access
    the dictionary directly, or use the get() method.

    Args:
        days   (list): list of ints/strings/date objects, the days on which to take the snapshot
        args   (list): additional day(s)
        die    (bool): whether or not to raise an exception if a date is not found (default true)
        kwargs (dict): passed to Analyzer()


    **Example**::

        sim = zn.Sim(analyzers=cv.snapshot('2020-04-04', '2020-04-14'))
        sim.run()
        snapshot = sim['analyzers'][0]
        agents = snapshot.snapshots[0]            # Option 1
        agents = snapshot.snapshots['2020-04-04'] # Option 2
        agents = snapshot.get('2020-04-14')       # Option 3
        agents = snapshot.get(34)                 # Option 4
        agents = snapshot.get()                   # Option 5
    '''

    def __init__(self, days, *args, die=True, **kwargs):
        super().__init__(**kwargs) # Initialize the Analyzer object
        days = sc.promotetolist(days) # Combine multiple days
        days.extend(args) # Include additional arguments, if present
        self.days      = days # Converted to integer representations
        self.die       = die  # Whether or not to raise an exception
        self.dates     = None # String representations
        self.start_day = None # Store the start date of the simulation
        self.snapshots = sc.odict() # Store the actual snapshots
        return


    def initialize(self, sim):
        self.start_day = sim['start_day'] # Store the simulation start day
        self.days, self.dates = zni.process_days(sim, self.days, return_dates=True) # Ensure days are in the right format
        max_snapshot_day = self.days[-1]
        max_sim_day = sim.day(sim['end_day'])
        if max_snapshot_day > max_sim_day: # pragma: no cover
            errormsg = f'Cannot create snapshot for {self.dates[-1]} (day {max_snapshot_day}) because the simulation ends on {self.end_day} (day {max_sim_day})'
            raise ValueError(errormsg)
        self.initialized = True
        return


    def apply(self, sim):
        for ind in zni.find_day(self.days, sim.t):
            date = self.dates[ind]
            self.snapshots[date] = sc.dcp(sim.agents) # Take snapshot!


    def finalize(self, sim):
        super().finalize()
        validate_recorded_dates(sim, requested_dates=self.dates, recorded_dates=self.snapshots.keys(), die=self.die)
        return


    def get(self, key=None):
        ''' Retrieve a snapshot from the given key (int, str, or date) '''
        if key is None:
            key = self.days[0]
        day  = sc.day(key, start_date=self.start_day)
        date = sc.date(day, start_date=self.start_day, as_date=False)
        if date in self.snapshots:
            snapshot = self.snapshots[date]
        else: # pragma: no cover
            dates = ', '.join(list(self.snapshots.keys()))
            errormsg = f'Could not find snapshot date {date} (day {day}): choices are {dates}'
            raise sc.KeyNotFoundError(errormsg)
        return snapshot
    
class biography(Analyzer):
    '''
    Analyzer that records the state of a given agent at given points of time. To retrieve them, you can either access
    the dictionary directly, or use the get() method.

    Args:
        uid    (Int): the uid of the agent to be tracked (optional, if agent_type is specified)
        agent_type (str): the type of agent to be tracked (optional, if uid is specified)
        days   (list): list of ints/strings/date objects, the days on which to take the snapshot
        args   (list): additional day(s)
        die    (bool): whether or not to raise an exception if a date is not found (default true)
        kwargs (dict): passed to Analyzer()
    '''
    def __init__(self, uid, agent_type, days, *args, die=True, **kwargs):
        super().__init__(**kwargs) # Initialize the Analyzer object
        days = sc.promotetolist(days) # Combine multiple days
        days.extend(args) # Include additional arguments, if present
        self.uid       = uid
        self.agent_type= agent_type
        self.days      = days # Converted to integer representations
        self.die       = die  # Whether or not to raise an exception
        self.dates     = None # String representations
        self.start_day = None # Store the start date of the simulation
        self.bio = sc.odict() # Store the actual snapshots
        self.df = None # Dataframe of results
        return
    
    def initialize(self, sim):
        self.start_day = sim['start_day'] # Store the simulation start day
        self.days, self.dates = zni.process_days(sim, self.days, return_dates=True) # Ensure days are in the right format
        max_bio_day = self.days[-1]
        max_sim_day = sim.day(sim['end_day'])
        if max_bio_day > max_sim_day: # pragma: no cover
            errormsg = f'Cannot create biography for {self.dates[-1]} (day {max_bio_day}) because the simulation ends on {self.end_day} (day {max_sim_day})'
            raise ValueError(errormsg)
        
        if self.agent_type is None and self.uid is None:
            errormsg = f'At least one of agent_type or uid arguments must be supplied.'
            raise ValueError(errormsg)
        elif self.uid is None:
            self.uid = np.random.choice(sim.agents[self.agent_type]['uid'])
        elif self.agent_type is None:
            self.agent_type = sim.agents.agent_type[sim.agents.uid == self.uid]
        else:
            if self.agent_type != sim.agents.agent_type[sim.agents.uid == self.uid]:
                errormsg = f'The specified agent type is inconsistent with the specified uid.'
                raise ValueError(errormsg)
            
        self.initialized = True
        return
    
    def apply(self, sim):
        for ind in zni.find_day(self.days, sim.t):
            date = self.dates[ind]
            bio_record = {}
            agent_ind = np.where(sim.agents[self.agent_type]['uid'] == self.uid)[0]
            for state in sim.agents[self.agent_type]['meta'].all_recordable_states:
                bio_record[state] = sim.agents[self.agent_type][state][agent_ind][0] # Get the state for the agent
            self.bio[date] = bio_record 
        return
    
    def finalize(self, sim):
        super().finalize()
        validate_recorded_dates(sim, requested_dates=self.dates, recorded_dates=self.bio.keys(), die=self.die)
        return
    
    def get(self, key=None):
        ''' Retrieve a bio record from the given key (int, str, or date) '''
        if key is None:
            key = self.days[0]
        day  = sc.day(key, start_date=self.start_day)
        date = sc.date(day, start_date=self.start_day, as_date=False)
        if date in self.bio:
            bio_record = self.bio[date]
        else: # pragma: no cover
            dates = ', '.join(list(self.bio.keys()))
            errormsg = f'Could not find bio record date {date} (day {day}): choices are {dates}'
            raise sc.KeyNotFoundError(errormsg)
        return bio_record
    
    def to_df(self):
        '''Create dataframe totals for each day'''
        df = pd.DataFrame()
        for date, k in self.bio.items():
            df_ = pd.DataFrame(k, index=[0])  # Convert the record to a DataFrame
            df_['date'] = date
            df = pd.concat((df, df_))
        cols = list(df.columns.values)
        cols = [cols[-1]] + [cols[-2]] + cols[:-2]
        self.df = df[cols]
        return self.df
    
    def plot(self, props_to_plot, do_show=None, fig_args=None, axis_args=None, plot_args=None,
             dateformat=None, **kwargs):
        '''
        Plot the results.

        Args:
            props_to_plot (list): list of properties to plot. Note that these must be valid properties for the agent type.
            do_show   (bool): whether to show the plot
            fig_args  (dict): passed to pl.figure()
            axis_args (dict): passed to pl.subplots_adjust()
            plot_args (dict): passed to pl.plot()
            dateformat (str): the format to use for the x-axes (only used for time series)
            kwargs    (dict): passed to ``zn.options.with_style()``
        '''
        if self.df is None:
            self.to_df()

        fig_args  = sc.mergedicts(dict(figsize=(18,11)), fig_args)
        axis_args = sc.mergedicts(dict(left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.25, hspace=0.4), axis_args)
        plot_args = sc.mergedicts(dict(lw=2, alpha=0.5, marker='o'), plot_args)

        with zno.options.with_style(**kwargs):
            fig, axs = plt.subplots(nrows=1, ncols=1, **fig_args)
            plt.subplots_adjust(**axis_args)
            colors = sc.vectocolor(len(props_to_plot))
            axs.set_title('Biography of agent uid={}'.format(self.uid))
            axs.set_xlabel('Date')
            axs.set_ylabel('Value')
            for i, prop in enumerate(props_to_plot):
                if prop not in self.df.columns:
                    errormsg = f'Property "{prop}" not found in the biography data. Available properties: {list(self.df.columns)}'
                    raise ValueError(errormsg)
                axs.plot(self.df['date'], self.df[prop], **plot_args, color=colors[i % len(colors)], label=prop)
                if dateformat is not None:
                    znplt.format_date_axis(axs, dateformat=dateformat)
                axs.legend(loc='best')

        return znplt.handle_show_return(fig=fig, do_show=do_show)