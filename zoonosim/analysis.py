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


def import_optuna():
    ''' A helper function to import Optuna, which is an optional dependency '''
    try:
        import optuna as op # Import here since it's slow
    except ModuleNotFoundError as E: # pragma: no cover
        errormsg = f'Optuna import failed ({str(E)}), please install first (pip install optuna)'
        raise ModuleNotFoundError(errormsg)
    return op
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
        name         (str)  : the name of the database (default: 'covasim_calibration')
        db_name      (str)  : the name of the database file (default: 'covasim_calibration.db')
        keep_db      (bool) : whether to keep the database after calibration (default: false)
        storage      (str)  : the location of the database (default: sqlite)
        label        (str)  : a label for this calibration object
        die          (bool) : whether to stop if an exception is encountered (default: false)
        verbose      (bool) : whether to print details of the calibration
        kwargs       (dict) : passed to cv.Calibration()

    Returns:
        A Calibration object

    **Example**::

        sim = cv.Sim(datafile='data.csv')
        calib_pars = dict(beta=[0.015, 0.010, 0.020])
        calib = cv.Calibration(sim, calib_pars, total_trials=100)
        calib.calibrate()
        calib.plot()

    New in version 3.0.3.
    '''

    def __init__(self, sim, calib_pars=None, fit_args=None, custom_fn=None, par_samplers=None,
                 n_trials=None, n_workers=None, total_trials=None, name=None, db_name=None,
                 keep_db=None, storage=None, label=None, die=False, verbose=True):
        super().__init__(label=label) # Initialize the Analyzer object

        import multiprocessing as mp # Import here since it's also slow

        # Handle run arguments
        if n_trials  is None: n_trials  = 20
        if n_workers is None: n_workers = mp.cpu_count()
        if name      is None: name      = 'covasim_calibration'
        if db_name   is None: db_name   = f'{name}.db'
        if keep_db   is None: keep_db   = False
        if storage   is None: storage   = f'sqlite:///{db_name}'
        if total_trials is not None: n_trials = total_trials/n_workers
        self.run_args   = sc.objdict(n_trials=int(n_trials), n_workers=int(n_workers), name=name, db_name=db_name, keep_db=keep_db, storage=storage)

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
        valid_pars = {k:v for k,v in calib_pars.items() if k in sim.pars}
        sim.update_pars(valid_pars)
        if self.custom_fn:
            sim = self.custom_fn(sim, calib_pars)
        else:
            if len(valid_pars) != len(calib_pars):
                extra = set(calib_pars.keys()) - set(valid_pars.keys())
                errormsg = f'The following parameters are not part of the sim, nor is a custom function specified to use them: {sc.strjoin(extra)}'
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
        pars = {}
        for key, (best,low,high) in self.calib_pars.items():
            if key in self.par_samplers: # If a custom sampler is used, get it now
                try:
                    sampler_fn = getattr(trial, self.par_samplers[key])
                except Exception as E:
                    errormsg = 'The requested sampler function is not found: ensure it is a valid attribute of an Optuna Trial object'
                    raise AttributeError(errormsg) from E
            else:
                sampler_fn = trial.suggest_uniform
            pars[key] = sampler_fn(key, low, high) # Sample from values within this range
        mismatch = self.run_sim(pars)
        return mismatch


    def worker(self):
        ''' Run a single worker '''
        op = import_optuna()
        if self.verbose:
            op.logging.set_verbosity(op.logging.DEBUG)
        else:
            op.logging.set_verbosity(op.logging.ERROR)
        study = op.load_study(storage=self.run_args.storage, study_name=self.run_args.name)
        output = study.optimize(self.run_trial, n_trials=self.run_args.n_trials)
        return output


    def run_workers(self):
        ''' Run multiple workers in parallel '''
        if self.run_args.n_workers > 1: # Normal use case: run in parallel
            output = sc.parallelize(self.worker, iterarg=self.run_args.n_workers)
        else: # Special case: just run one
            output = [self.worker()]
        return output


    def remove_db(self):
        '''
        Remove the database file if keep_db is false and the path exists.

        New in version 3.1.0.
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
        self.best_pars = sc.objdict(self.study.best_params)
        self.elapsed = sc.toc(t0, output=True)

        # Compare the results
        self.initial_pars = sc.objdict({k:v[0] for k,v in self.calib_pars.items()})
        self.par_bounds   = sc.objdict({k:np.array([v[1], v[2]]) for k,v in self.calib_pars.items()})
        self.before = self.run_sim(calib_pars=self.initial_pars, label='Before calibration', return_sim=True)
        self.after  = self.run_sim(calib_pars=self.best_pars,    label='After calibration',  return_sim=True)
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

        New in version 3.1.1.
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

        New in version 3.1.1: renamed from plot() to plot_sims().
        '''
        msim = znr.MultiSim([self.before, self.after])
        fig = msim.plot(**kwargs)
        return znplt.handle_show_return(fig=fig)


    def plot_trend(self, best_thresh=2):
        '''
        Plot the trend in best mismatch over time.

        New in version 3.1.1.
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

        New in version 3.1.1.
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

        New in version 3.1.1.
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