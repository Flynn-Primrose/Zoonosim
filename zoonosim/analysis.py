'''
Additional analysis functions that are not part of the core Zoonosim workflow,
but which are useful for particular investigations.
'''

import os
import sciris as sc
import numpy as np
import pylab as pl
import pandas as pd
import matplotlib.pyplot as plt
from . import interventions as zni
from . import misc as znm
from . import plotting as znplt
from . import options as zno
from . import run as znr


__all__ = ['Analyzer', 'snapshot', 'biography', 'Fit' ,'Calibration' ]


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

class Fit(Analyzer):
    '''
    A class for calculating the fit between the model and the data. Note the
    following terminology is used here:

        - fit: nonspecific term for how well the model matches the data
        - difference: the absolute numerical differences between the model and the data (one time series per result)
        - goodness-of-fit: the result of passing the difference through a statistical function, such as mean squared error
        - loss: the goodness-of-fit for each result multiplied by user-specified weights (one time series per result)
        - mismatches: the sum of all the losses (a single scalar value per time series)
        - mismatch: the sum of the mismatches -- this is the value to be minimized during calibration

    Args:
        sim (Sim): the sim object
        weights (dict): the relative weight to place on each result (by default: 10 for deaths, 5 for diagnoses, 1 for everything else)
        keys (list): the keys to use in the calculation
        custom (dict): a custom dictionary of additional data to fit; format is e.g. {'my_output':{'data':[1,2,3], 'sim':[1,2,4], 'weights':2.0}}
        compute (bool): whether to compute the mismatch immediately
        verbose (bool): detail to print
        die (bool): whether to raise an exception if no data are supplied
        label (str): the label for the analyzer
        kwargs (dict): passed to zn.compute_gof() -- see this function for more detail on goodness-of-fit calculation options

    **Example**::

        sim = zn.Sim(datafile='my-data-file.csv')
        sim.run()
        fit = sim.compute_fit()
        fit.plot()
    '''

    def __init__(self, sim, weights=None, keys=None, custom=None, compute=True, verbose=False, die=True, label=None, **kwargs):
        super().__init__(label=label) # Initialize the Analyzer object

        # Handle inputs
        self.weights    = weights
        self.custom     = sc.mergedicts(custom)
        self.verbose    = verbose
        self.weights    = sc.mergedicts({'cum_deaths':10, 'cum_diagnoses':5}, weights)
        self.keys       = keys
        self.gof_kwargs = kwargs
        self.die        = die

        # Copy data
        if sim.data is None: # pragma: no cover
            errormsg = 'Model fit cannot be calculated until data are loaded'
            if self.die:
                raise RuntimeError(errormsg)
            else:
                znm.warn(errormsg)
                sim.data = pd.DataFrame() # Use an empty dataframe
        self.data = sim.data

        # Copy sim results
        if not sim.results_ready: # pragma: no cover
            errormsg = 'Model fit cannot be calculated until results are run'
            if self.die: raise RuntimeError(errormsg)
            else:        znm.warn(errormsg)
        self.sim_results = sc.objdict()
        for key in sim.result_keys() + ['t', 'date']:
            self.sim_results[key] = sim.results[key]
        self.sim_npts = sim.npts # Number of time points in the sim

        # Copy other things
        self.sim_dates = sim.datevec.tolist()

        # These are populated during initialization
        self.inds         = sc.objdict() # To store matching indices between the data and the simulation
        self.inds.sim     = sc.objdict() # For storing matching indices in the sim
        self.inds.data    = sc.objdict() # For storing matching indices in the data
        self.date_matches = sc.objdict() # For storing matching dates, largely for plotting
        self.pair         = sc.objdict() # For storing perfectly paired points between the data and the sim
        self.diffs        = sc.objdict() # Differences between pairs
        self.gofs         = sc.objdict() # Goodness-of-fit for differences
        self.losses       = sc.objdict() # Weighted goodness-of-fit
        self.mismatches   = sc.objdict() # Final mismatch values
        self.mismatch     = None # The final value

        if compute:
            self.compute()

        return


    def compute(self):
        ''' Perform all required computations '''
        self.reconcile_inputs() # Find matching values
        self.compute_diffs() # Perform calculations
        self.compute_gofs()
        self.compute_losses()
        self.compute_mismatch()
        return self.mismatch


    def reconcile_inputs(self):
        ''' Find matching keys and indices between the model and the data '''

        data_cols = self.data.columns
        if self.keys is None:
            sim_keys = [k for k in self.sim_results.keys() if k.startswith('cum_')] # Default sim keys, only keep cumulative keys if no keys are supplied
            intersection = list(set(sim_keys).intersection(data_cols)) # Find keys in both the sim and data
            self.keys = [key for key in sim_keys if key in intersection] # Maintain key order
            if not len(self.keys): # pragma: no cover
                errormsg = f'No matches found between simulation result keys:\n{sc.strjoin(sim_keys)}\n\nand data columns:\n{sc.strjoin(data_cols)}'
                if self.die: raise sc.KeyNotFoundError(errormsg)
                else:        znm.warn(errormsg)
        mismatches = [key for key in self.keys if key not in data_cols]
        if len(mismatches): # pragma: no cover
            mismatchstr = ', '.join(mismatches)
            errormsg = f'The following requested key(s) were not found in the data: {mismatchstr}'
            if self.die: raise sc.KeyNotFoundError(errormsg)
            else:        znm.warn(errormsg)

        for key in self.keys: # For keys present in both the results and in the data
            self.inds.sim[key]  = []
            self.inds.data[key] = []
            self.date_matches[key] = []
            count = -1
            for d, datum in self.data[key].items():
                count += 1
                if np.isfinite(datum):
                    if d in self.sim_dates:
                        self.date_matches[key].append(d)
                        self.inds.sim[key].append(self.sim_dates.index(d))
                        self.inds.data[key].append(count)
            self.inds.sim[key]  = np.array(self.inds.sim[key])
            self.inds.data[key] = np.array(self.inds.data[key])

        # Convert into paired points
        matches = 0 # Count how many data points match
        for key in self.keys:
            self.pair[key] = sc.objdict()
            sim_inds = self.inds.sim[key]
            data_inds = self.inds.data[key]
            n_inds = len(sim_inds)
            self.pair[key].sim  = np.zeros(n_inds)
            self.pair[key].data = np.zeros(n_inds)
            for i in range(n_inds):
                matches += 1
                self.pair[key].sim[i]  = self.sim_results[key].values[sim_inds[i]]
                self.pair[key].data[i] = self.data[key].values[data_inds[i]]

        # Process custom inputs
        self.custom_keys = list(self.custom.keys())
        for key in self.custom.keys():
            matches += 1 # If any of these exist, count it as  amatch

            # Initialize and do error checking
            custom = self.custom[key]
            c_keys = list(custom.keys())
            if 'sim' not in c_keys or 'data' not in c_keys:
                errormsg = f'Custom input must have "sim" and "data" keys, not {c_keys}'
                raise sc.KeyNotFoundError(errormsg)
            c_data = custom['data']
            c_sim  = custom['sim']
            try:
                assert len(c_data) == len(c_sim)
            except: # pragma: no cover
                errormsg = f'Custom data and sim must be arrays, and be of the same length: data = {c_data}, sim = {c_sim} could not be processed'
                raise ValueError(errormsg)
            if key in self.pair: # pragma: no cover
                errormsg = f'You cannot use a custom key "{key}" that matches one of the existing keys: {self.pair.keys()}'
                raise ValueError(errormsg)

            # If all tests pass, simply copy the data
            self.pair[key] = sc.objdict()
            self.pair[key].sim  = c_sim
            self.pair[key].data = c_data

            # Process weight, if available
            wt = custom.get('weight', 1.0) # Attempt to retrieve key 'weight', or use the default if not provided
            wt = custom.get('weights', wt) # ...but also try "weights"
            self.weights[key] = wt # Set the weight

        if matches == 0:
            errormsg = 'No paired data points were found between the supplied data and the simulation; please check the dates for each'
            if self.die: raise ValueError(errormsg)
            else:        znm.warn(errormsg)

        return


    def compute_diffs(self, absolute=False):
        ''' Find the differences between the sim and the data '''
        for key in self.pair.keys():
            self.diffs[key] = self.pair[key].sim - self.pair[key].data
            if absolute:
                self.diffs[key] = np.abs(self.diffs[key])
        return


    def compute_gofs(self, **kwargs):
        ''' Compute the goodness-of-fit '''
        kwargs = sc.mergedicts(self.gof_kwargs, kwargs)
        for key in self.pair.keys():
            actual    = sc.dcp(self.pair[key].data)
            predicted = sc.dcp(self.pair[key].sim)
            self.gofs[key] = znm.compute_gof(actual, predicted, **kwargs)
        return


    def compute_losses(self):
        ''' Compute the weighted goodness-of-fit '''
        for key in self.gofs.keys():
            if key in self.weights:
                weight = self.weights[key]
                if sc.isiterable(weight): # It's an array
                    len_wt = len(weight)
                    len_sim = self.sim_npts
                    len_match = len(self.gofs[key])
                    if len_wt == len_match: # If the weight already is the right length, do nothing
                        pass
                    elif len_wt == len_sim: # Most typical case: it's the length of the simulation, must trim
                        weight = weight[self.inds.sim[key]] # Trim to matching indices
                    else: # pragma: no cover
                        errormsg = f'Could not map weight array of length {len_wt} onto simulation of length {len_sim} or data-model matches of length {len_match}'
                        raise ValueError(errormsg)
            else:
                weight = 1.0
            self.losses[key] = self.gofs[key]*weight
        return


    def compute_mismatch(self, use_median=False):
        ''' Compute the final mismatch '''
        for key in self.losses.keys():
            if use_median:
                self.mismatches[key] = np.median(self.losses[key])
            else:
                self.mismatches[key] = np.sum(self.losses[key])
        self.mismatch = self.mismatches[:].sum()
        return self.mismatch


    def plot(self, keys=None, width=0.8, fig_args=None, axis_args=None, plot_args=None,
             date_args=None, do_show=None, fig=None, **kwargs):
        '''
        Plot the fit of the model to the data. For each result, plot the data
        and the model; the difference; and the loss (weighted difference). Also
        plots the loss as a function of time.

        Args:
            keys      (list):  which keys to plot (default, all)
            width     (float): bar width
            fig_args  (dict):  passed to ``pl.figure()``
            axis_args (dict):  passed to ``pl.subplots_adjust()``
            plot_args (dict):  passed to ``pl.plot()``
            date_args (dict):  passed to ``zn.plotting.reset_ticks()`` (handle date format, rotation, etc.)
            do_show   (bool):  whether to show the plot
            fig       (fig):   if supplied, use this figure to plot in
            kwargs    (dict):  passed to ``zn.options.with_style()``

        Returns:
            Figure object
        '''

        fig_args  = sc.mergedicts(dict(figsize=(18,11)), fig_args)
        axis_args = sc.mergedicts(dict(left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.3, hspace=0.3), axis_args)
        plot_args = sc.mergedicts(dict(lw=2, alpha=0.5, marker='o'), plot_args)
        date_args = sc.mergedicts(sc.objdict(as_dates=True, dateformat=None, rotation=None, start=None, end=None), date_args)

        if keys is None:
            keys = self.keys + self.custom_keys
        n_keys = len(keys)

        loss_ax = None
        colors = sc.gridcolors(n_keys)
        n_rows = 4

        # Plot
        with zno.options.with_style(**kwargs):
            if fig is None:
                fig = pl.figure(**fig_args)
            pl.subplots_adjust(**axis_args)
            main_ax1 = pl.subplot(n_rows, 2, 1)
            main_ax2 = pl.subplot(n_rows, 2, 2)
            bottom = sc.objdict() # Keep track of the bottoms for plotting cumulative
            bottom.daily = np.zeros(self.sim_npts)
            bottom.cumul = np.zeros(self.sim_npts)
            for k,key in enumerate(keys):
                if key in self.keys: # It's a time series, plot with days and dates
                    days      = self.inds.sim[key] # The "days" axis (or not, for custom keys)
                    daylabel  = 'Date'
                else: #It's custom, we don't know what it is
                    days      = np.arange(len(self.losses[key])) # Just use indices
                    daylabel  = 'Index'

                # Cumulative totals can't mix daily and non-daily inputs, so skip custom keys
                if key in self.keys:
                    for i,ax in enumerate([main_ax1, main_ax2]):

                        if i == 0:
                            data = self.losses[key]
                            ylabel = 'Daily mismatch'
                            title = 'Daily total mismatch'
                        else:
                            data = np.cumsum(self.losses[key])
                            ylabel = 'Cumulative mismatch'
                            title = f'Cumulative mismatch: {self.mismatch:0.3f}'

                        dates = self.sim_results['date'][days] # Show these with dates, rather than days, as a reference point
                        ax.bar(dates, data, width=width, bottom=bottom[i][self.inds.sim[key]], color=colors[k], label=f'{key}')

                        if i == 0:
                            bottom.daily[self.inds.sim[key]] += self.losses[key]
                        else:
                            bottom.cumul = np.cumsum(bottom.daily)

                        if k == len(self.keys)-1:
                            ax.set_xlabel('Date')
                            ax.set_ylabel(ylabel)
                            ax.set_title(title)
                            znplt.reset_ticks(ax=ax, date_args=date_args, start_day=self.sim_results['date'][0])
                            ax.legend()

                ts_ax = pl.subplot(n_rows, n_keys, k+1*n_keys+1)
                ts_ax.plot(days, self.pair[key].data, c='k', label='Data', **plot_args)
                ts_ax.plot(days, self.pair[key].sim, c=colors[k], label='Simulation', **plot_args)
                ts_ax.set_title(key)
                if k == 0:
                    ts_ax.set_ylabel('Time series (counts)')
                    ts_ax.legend()

                diff_ax = pl.subplot(n_rows, n_keys, k+2*n_keys+1)
                diff_ax.bar(days, self.diffs[key], width=width, color=colors[k], label='Difference')
                diff_ax.axhline(0, c='k')
                if k == 0:
                    diff_ax.set_ylabel('Differences (counts)')
                    diff_ax.legend()

                loss_ax = pl.subplot(n_rows, n_keys, k+3*n_keys+1, sharey=loss_ax)
                loss_ax.bar(days, self.losses[key], width=width, color=colors[k], label='Losses')
                loss_ax.set_xlabel(daylabel)
                loss_ax.set_title(f'Total loss: {self.losses[key].sum():0.3f}')
                if k == 0:
                    loss_ax.set_ylabel('Losses')
                    loss_ax.legend()

                if daylabel == 'Date':
                    for ax in [ts_ax, diff_ax, loss_ax]:
                        znplt.reset_ticks(ax=ax, date_args=date_args, start_day=self.sim_results['date'][0])

        return znplt.handle_show_return(fig=fig, do_show=do_show)

# Helper functions for calibration

def import_optuna():
    ''' A helper function to import Optuna, which is an optional dependency '''
    try:
        import optuna as op # Import here since it's slow
    except ModuleNotFoundError as E: # pragma: no cover
        errormsg = f'Optuna import failed ({str(E)}), please install first (pip install optuna)'
        raise ModuleNotFoundError(errormsg)
    return op

def compare_pars(dict_new, dict_orig):
    """
    Compare two dictionaries recursively.
    Returns:
        common  -> keys (and subkeys) present in both dict_new and dict_orig
        missing -> keys (and subkeys) present in dict_new but not in dict_orig
    """
    common = {}
    missing = {}
    for key, val in dict_new.items():
        if key in dict_orig:
            if isinstance(val, dict) and isinstance(dict_orig[key], dict):
                # Recurse into subdicts
                sub_common, sub_missing = compare_pars(val, dict_orig[key])
                if sub_common:
                    common[key] = sub_common
                if sub_missing:
                    missing[key] = sub_missing
            else:
                # Key exists in both, keep dict_new's value
                common[key] = val
        else:
            # Key missing from dict_orig
            missing[key] = val
    return common, missing

def is_empty(d):
    ''' Check if a dictionary (possibly nested) is empty '''
    if not isinstance(d, dict):
        return False
    return all(is_empty(v) for v in d.values()) if d else True

def pars_sampler(trial, calib_pars, par_samplers):
    '''
    A helper function to sample parameters from Optuna trials
    '''
    sampled_pars = {}
    for key, value in calib_pars.items():
        if isinstance(value, list) and len(value) == 3:
            best, low, high = value
            if key in par_samplers: # If a custom sampler is used, get it now
                try:
                    sampler_fn = getattr(trial, par_samplers[key])
                except Exception as E:
                    errormsg = 'The requested sampler function is not found: ensure it is a valid attribute of an Optuna Trial object'
                    raise AttributeError(errormsg) from E
            else:
                sampler_fn = trial.suggest_uniform
            sampled_pars[key] = sampler_fn(key, low, high) # Sample from values within this range
        elif isinstance(value, list):
            errormsg = f'Parameter "{key}" must be a list of [best, low, high]'
            raise ValueError(errormsg)
        elif isinstance(value, dict):
            sampled_pars[key] = pars_sampler(trial, value, par_samplers.get(key, {})) # Recurse into sub-dictionaries
        else:
            errormsg = f'Parameter "{key}" must be a list of [best, low, high] or a dictionary for nested parameters'
            raise ValueError(errormsg)
    return sampled_pars
            
def pars_parser(calib_pars):
    '''
    A helper function to parse calibration parameters into initial, and bounds.
    '''
    initial_pars = {}
    par_bounds = {}
    for key, value in calib_pars.items():
        if isinstance(value, list) and len(value) == 3:
            best, low, high = value
            initial_pars[key] = best
            par_bounds[key] = [low, high]
        elif isinstance(value, dict):
            sub_best, sub_bounds = pars_parser(value) # Recurse into sub-dictionaries
            initial_pars[key] = sub_best
            par_bounds[key] = sub_bounds
        else:
            errormsg = f'Parameter "{key}" must be a list of [best, low, high] or a dictionary for nested parameters'
            raise ValueError(errormsg)
    return initial_pars, par_bounds

class Calibration(Analyzer):
    '''
    A class to handle calibration of Covasim simulations. Uses the Optuna hyperparameter
    optimization library (optuna.org), which must be installed separately (via
    pip install optuna).

    Note: running a calibration does not guarantee a good fit! You must ensure that
    you run for a sufficient number of iterations, have enough free parameters, and
    that the parameters have wide enough bounds. Please see the tutorial on calibration
    for more information.

    Args:
        sim          (Sim)  : the simulation to calibrate
        calib_pars   (dict) : a dictionary of the parameters to calibrate of the format dict(key1=[best, low, high])
        fit_args     (dict) : a dictionary of options that are passed to sim.compute_fit() to calculate the goodness-of-fit
        par_samplers (dict) : an optional mapping from parameters to the Optuna sampler to use for choosing new points for each; by default, suggest_uniform
        custom_fn    (func) : a custom function for modifying the simulation; receives the sim and calib_pars as inputs, should return the modified sim
        n_trials     (int)  : the number of trials per worker
        n_workers    (int)  : the number of parallel workers (default: maximum
        total_trials (int)  : if n_trials is not supplied, calculate by dividing this number by n_workers)
        name         (str)  : the name of the database (default: 'zoonosim_calibration')
        db_name      (str)  : the name of the database file (default: 'zoonosim_calibration.db')
        keep_db      (bool) : whether to keep the database after calibration (default: false)
        storage      (str)  : the location of the database (default: sqlite)
        label        (str)  : a label for this calibration object
        die          (bool) : whether to stop if an exception is encountered (default: false)
        verbose      (bool) : whether to print details of the calibration
        kwargs       (dict) : passed to zn.Calibration()

    Returns:
        A Calibration object

    **Example**::

        sim = zn.Sim(datafile='data.csv')
        calib_pars = dict(beta=[0.015, 0.010, 0.020])
        calib = zn.Calibration(sim, calib_pars, total_trials=100)
        calib.calibrate()
        calib.plot()

 
    '''

    def __init__(self, sim, calib_pars=None, fit_args=None, custom_fn=None, par_samplers=None,
                 n_trials=None, n_workers=None, total_trials=None, name=None, db_name=None,
                 keep_db=None, storage=None, label=None, die=False, verbose=True):
        super().__init__(label=label) # Initialize the Analyzer object

        import multiprocessing as mp # Import here since it's also slow

        # Handle run arguments
        if n_trials  is None: n_trials  = 20
        if n_workers is None: n_workers = mp.cpu_count()
        if name      is None: name      = 'zoonosim_calibration'
        if db_name   is None: db_name   = f'{name}.db'
        if keep_db   is None: keep_db   = False
        if storage   is None: storage   = f'sqlite:///{db_name}'
        if total_trials is not None: n_trials = total_trials/n_workers
        self.run_args   = sc.objdict(n_trials=int(n_trials), n_workers=int(n_workers), name=name, db_name=db_name, keep_db=keep_db, storage=storage)

        self.run_args.setdefault('parallelizer', 'concurrent')

        # if self.run_args['parallelizer'] != 'concurrent' and verbose >= 1:
        #     print(f"Warning: Parallelizer is explicitly set to '{self.run_args['parallelizer']}', this may cause issues on Windows machines. Consider setting parallelizer to 'concurrent'")

        # Handle other inputs
        self.sim          = sim
        self.calib_pars   = calib_pars
        self.fit_args     = sc.mergedicts(fit_args)
        self.par_samplers = sc.mergedicts(par_samplers)
        self.custom_fn    = custom_fn
        self.die          = die
        self.verbose      = verbose
        self.calibrated   = False

        # Handle if the sim has already been run
        if self.sim.complete:
            warnmsg = 'Sim has already been run; re-initializing, but in future, use a sim that has not been run'
            znm.warn(warnmsg)
            self.sim = self.sim.copy()
            self.sim.initialize()

        return




    def run_sim(self, calib_pars, label=None, return_sim=False):
        ''' Create and run a simulation '''
        sim = self.sim.copy()
        if label: sim.label = label
        valid_pars, invalid_pars = compare_pars(calib_pars, sim.pars)
        sim.update_pars(valid_pars, recursive=True)
        if self.custom_fn:
            sim = self.custom_fn(sim, calib_pars)
        else:
            if not is_empty(invalid_pars):
                errormsg = f'The following parameters are not part of the sim, nor is a custom function specified to use them: {invalid_pars}'
                raise ValueError(errormsg)
        try:
            sim.run()
            sim.compute_fit(**self.fit_args)
            if return_sim:
                return sim
            else:
                return sim.fit.mismatch
        except Exception as E:
            if self.die:
                raise E
            else:
                warnmsg = f'Encountered error running sim!\nParameters:\n{valid_pars}\nTraceback:\n{sc.traceback()}'
                znm.warn(warnmsg)
                output = None if return_sim else np.inf
                return output


    def run_trial(self, trial):
        ''' Define the objective for Optuna '''
        try:
            pars = pars_sampler(trial, self.calib_pars, self.par_samplers)
            mismatch = self.run_sim(pars)
        except Exception as E:
            errormsg = f'Error during trial sampling or simulation run: {str(E)}'
            raise RuntimeError(errormsg) from E
        return mismatch


    def worker(self):
        ''' Run a single worker '''
        op = import_optuna()
        if self.verbose:
            op.logging.set_verbosity(op.logging.DEBUG)
        else:
            op.logging.set_verbosity(op.logging.ERROR)
        try:
            study = op.load_study(storage=self.run_args.storage, study_name=self.run_args.name)
            output = study.optimize(self.run_trial, n_trials=self.run_args.n_trials)
        except Exception as E:
            errormsg = f'Optuna study could not be loaded or optimized; please ensure the database path is correct: {self.run_args.storage}'
            raise RuntimeError(errormsg) from E
        return output


    def run_workers(self):
        ''' Run multiple workers in parallel '''
        if self.run_args.n_workers > 1: # Normal use case: run in parallel
            try:
                output = sc.parallelize(self.worker, iterarg=self.run_args.n_workers, parallelizer=self.run_args.parallelizer)
            except Exception as E:
                if isinstance(E, RuntimeError):
                    if 'freeze_support' in E.args[0]: # For this error, add additional information
                        errormsg = '''
                                Uh oh! It appears you are trying to run with multiprocessing on Windows outside
                                of the __main__ block; please see https://docs.python.org/3/library/multiprocessing.html
                                for more information. The correct syntax to use is e.g.
                            
                                    import zoonosim as zn
                                    sim = zn.Sim(data_file='data.csv')
                                    calib = zn.Calibration(sim, calib_pars, n_workers=4, total_trials=100)
                            
                                    if __name__ == '__main__':
                                        calib.calibrate()
                            
                                Alternatively, to run without multiprocessing, set n_workers = 1.
                                '''
                    raise RuntimeError(errormsg) from E
                else: # For all other runtime errors, raise the original exception
                    raise E
        else: # Special case: just run one
            output = [self.worker()]
        return output


    def remove_db(self):
        '''
        Remove the database file if keep_db is false and the path exists.
        '''
        if os.path.exists(self.run_args.db_name):
            os.remove(self.run_args.db_name)
            if self.verbose:
                print(f'Removed existing calibration {self.run_args.db_name}')
        return


    def make_study(self):
        ''' Make a study, deleting one if it already exists '''
        op = import_optuna()
        if not self.run_args.keep_db:
            self.remove_db()
        output = op.create_study(storage=self.run_args.storage, study_name=self.run_args.name)
        return output


    def calibrate(self, calib_pars=None, verbose=True, **kwargs):
        '''
        Actually perform calibration.

        Args:
            calib_pars (dict): if supplied, overwrite stored calib_pars
            verbose (bool): whether to print output from each trial
            kwargs (dict): if supplied, overwrite stored run_args (n_trials, n_workers, etc.)
        '''
        op = import_optuna()

        # Load and validate calibration parameters
        if calib_pars is not None:
            self.calib_pars = calib_pars
        if self.calib_pars is None:
            errormsg = 'You must supply calibration parameters either when creating the calibration object or when calling calibrate().'
            raise ValueError(errormsg)
        self.run_args.update(kwargs) # Update optuna settings


        # Run the optimization
        t0 = sc.tic()
        self.make_study()
        self.run_workers()
        self.study = op.load_study(storage=self.run_args.storage, study_name=self.run_args.name)
        #self.best_pars = sc.objdict(self.study.best_params)
        self.elapsed = sc.toc(t0, output=True)

        # Compare the results
        initial_pars, par_bounds = pars_parser(self.calib_pars)
        self.initial_pars  = sc.objdict(initial_pars)
        self.par_bounds    = sc.objdict(par_bounds)
        #self.initial_pars = sc.objdict({k:v[0] for k,v in self.calib_pars.items()})
        #self.par_bounds   = sc.objdict({k:np.array([v[1], v[2]]) for k,v in self.calib_pars.items()})
        self.before = self.run_sim(calib_pars=self.initial_pars, label='Before calibration', return_sim=True)
        #self.after  = self.run_sim(calib_pars=self.best_pars,    label='After calibration',  return_sim=True)
        self.parse_study()

        # Tidy up
        self.calibrated = True
        if not self.run_args.keep_db:
            self.remove_db()
        if verbose:
            self.summarize()

        return self


    def summarize(self):
        ''' Print out results from the calibration '''
        if self.calibrated:
            print(f'Calibration for {self.run_args.n_workers*self.run_args.n_trials} total trials completed in {self.elapsed:0.1f} s.')
            before = self.before.fit.mismatch
            after = self.after.fit.mismatch
            print('\nInitial parameter values:')
            print(self.initial_pars)
            print('\nBest parameter values:')
            print(self.best_pars)
            print(f'\nMismatch before calibration: {before:n}')
            print(f'Mismatch after calibration:  {after:n}')
            print(f'Percent improvement:         {((before-after)/before)*100:0.1f}%')
            return before, after
        else:
            print('Calibration not yet run; please run calib.calibrate()')
            return


    def parse_study(self):
        '''Parse the study into a data frame -- called automatically '''
        best = self.best_pars

        print('Making results structure...')
        results = []
        n_trials = len(self.study.trials)
        failed_trials = []
        for trial in self.study.trials:
            data = {'index':trial.number, 'mismatch': trial.value}
            for key,val in trial.params.items():
                data[key] = val
            if data['mismatch'] is None:
                failed_trials.append(data['index'])
            else:
                results.append(data)
        print(f'Processed {n_trials} trials; {len(failed_trials)} failed')

        keys = ['index', 'mismatch'] + list(best.keys())
        data = sc.objdict().make(keys=keys, vals=[])
        for i,r in enumerate(results):
            for key in keys:
                if key not in r:
                    warnmsg = f'Key {key} is missing from trial {i}, replacing with default'
                    znm.warn(warnmsg)
                    r[key] = best[key]
                data[key].append(r[key])
        self.data = data
        self.df = pd.DataFrame.from_dict(data)

        return


    def to_json(self, filename=None):
        '''
        Convert the data to JSON.
        '''
        order = np.argsort(self.df['mismatch'])
        json = []
        for o in order:
            row = self.df.iloc[o,:].to_dict()
            rowdict = dict(index=row.pop('index'), mismatch=row.pop('mismatch'), pars={})
            for key,val in row.items():
                rowdict['pars'][key] = val
            json.append(rowdict)
        if filename:
            sc.savejson(filename, json, indent=2)
        else:
            return json


    def plot_sims(self, **kwargs):
        '''
        Plot sims, before and after calibration.
        '''
        msim = znr.MultiSim([self.before, self.after])
        fig = msim.plot(**kwargs)
        return znplt.handle_show_return(fig=fig)


    def plot_trend(self, best_thresh=2):
        '''
        Plot the trend in best mismatch over time.
        '''
        mismatch = sc.dcp(self.df['mismatch'].values)
        best_mismatch = np.zeros(len(mismatch))
        for i in range(len(mismatch)):
            best_mismatch[i] = mismatch[:i+1].min()
        smoothed_mismatch = sc.smooth(mismatch)
        fig = pl.figure(figsize=(16,12), dpi=120)

        ax1 = pl.subplot(2,1,1)
        pl.plot(mismatch, alpha=0.2, label='Original')
        pl.plot(smoothed_mismatch, lw=3, label='Smoothed')
        pl.plot(best_mismatch, lw=3, label='Best')

        ax2 = pl.subplot(2,1,2)
        max_mismatch = mismatch.min()*best_thresh
        inds = sc.findinds(mismatch<=max_mismatch)
        pl.plot(best_mismatch, lw=3, label='Best')
        pl.scatter(inds, mismatch[inds], c=mismatch[inds], label='Usable indices')
        for ax in [ax1, ax2]:
            pl.sca(ax)
            pl.grid(True)
            pl.legend()
            sc.setylim()
            sc.setxlim()
            pl.xlabel('Trial number')
            pl.ylabel('Mismatch')
        return znplt.handle_show_return(fig=fig)


    def plot_all(self): # pragma: no cover
        '''
        Plot every point in the calibration. Warning, very slow for more than a few hundred trials.
        '''
        g = pairplotpars(self.data, color_column='mismatch', bounds=self.par_bounds)
        return g


    def plot_best(self, best_thresh=2): # pragma: no cover
        ''' Plot only the points with lowest mismatch. New in version 3.1.1. '''
        max_mismatch = self.df['mismatch'].min()*best_thresh
        inds = sc.findinds(self.df['mismatch'].values <= max_mismatch)
        g = pairplotpars(self.data, inds=inds, color_column='mismatch', bounds=self.par_bounds)
        return g


    def plot_stride(self, npts=200): # pragma: no cover
        '''
        Plot a fixed number of points in order across the results.
        '''
        npts = min(len(self.df), npts)
        inds = np.linspace(0, len(self.df)-1, npts).round()
        g = pairplotpars(self.data, inds=inds, color_column='mismatch', bounds=self.par_bounds)
        return g

def pairplotpars(data, inds=None, color_column=None, bounds=None, cmap='parula', bins=None, edgecolor='w', facecolor='#F8A493', figsize=(20,16)): # pragma: no cover
    ''' Plot scatterplots, histograms, and kernel densities for calibration results '''
    try:
        import seaborn as sns # Optional import
    except ModuleNotFoundError as E:
        errormsg = 'Calibration plotting requires Seaborn; please install with "pip install seaborn"'
        raise ModuleNotFoundError(errormsg) from E

    data = sc.odict(sc.dcp(data))

    # Create the dataframe
    df = pd.DataFrame.from_dict(data)
    if inds is not None:
        df = df.iloc[inds,:].copy()

    # Choose the colors
    if color_column:
        colors = sc.vectocolor(df[color_column].values, cmap=cmap)
    else:
        colors = [facecolor for i in range(len(df))]
    df['color_column'] = [sc.rgb2hex(rgba[:-1]) for rgba in colors]

    # Make the plot
    grid = sns.PairGrid(df)
    grid = grid.map_lower(pl.scatter, **{'facecolors':df['color_column']})
    grid = grid.map_diag(pl.hist, bins=bins, edgecolor=edgecolor, facecolor=facecolor)
    grid = grid.map_upper(sns.kdeplot)
    grid.fig.set_size_inches(figsize)
    grid.fig.tight_layout()

    # Set bounds
    if bounds:
        for ax in grid.axes.flatten():
            xlabel = ax.get_xlabel()
            ylabel = ax.get_ylabel()
            if xlabel in bounds:
                ax.set_xlim(bounds[xlabel])
            if ylabel in bounds:
                ax.set_ylim(bounds[ylabel])

    return grid    