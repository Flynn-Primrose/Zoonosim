'''
ParsObj inherits from FlexPretty and adds an attribute self.pars, which is a dictionary of parameters.
It also defines methods for manipulating the parameters.

BaseSim inherits from ParsObj and defines methods for manipulating the simulation object that are not
directly related to simulating the epidemic. This includes methods for saving and loading the simulation,
exporting results, and getting interventions and analyzers.

Also defines the function set_metadata(), which sets the metadata for the simulation object.
'''

import sciris as sc
import numpy as np
import pandas as pd
import datetime as dt

from . import FlexPretty
from . import Result
from .. import version as znv
from .. import misc as znm
from .. import utils as znu

__all__ = ['BaseSim']

class ParsObj(FlexPretty):
    '''
    A class based around performing operations on a self.pars dict.
    '''

    def __init__(self, pars):
        self.update_pars(pars, create=True)
        return


    def __getitem__(self, key):
        ''' Allow sim['par_name'] instead of sim.pars['par_name'] '''
        try:
            return self.pars[key]
        except:
            all_keys = '\n'.join(list(self.pars.keys()))
            errormsg = f'Key "{key}" not found; available keys:\n{all_keys}'
            raise sc.KeyNotFoundError(errormsg)


    def __setitem__(self, key, value):
        ''' Ditto '''
        if key in self.pars:
            self.pars[key] = value
        else:
            all_keys = '\n'.join(list(self.pars.keys()))
            errormsg = f'Key "{key}" not found; available keys:\n{all_keys}'
            raise sc.KeyNotFoundError(errormsg)
        return


    def update_pars(self, pars=None, create=False):
        '''
        Update internal dict with new pars.

        Args:
            pars (dict): the parameters to update (if None, do nothing)
            create (bool): if create is False, then raise a KeyNotFoundError if the key does not already exist
        '''
        if pars is not None:
            if not isinstance(pars, dict):
                raise TypeError(f'The pars object must be a dict; you supplied a {type(pars)}')
            if not hasattr(self, 'pars'):
                self.pars = pars
            if not create:
                available_keys = list(self.pars.keys())
                mismatches = [key for key in pars.keys() if key not in available_keys]
                if len(mismatches):
                    errormsg = f'Key(s) {mismatches} not found; available keys are {available_keys}'
                    raise sc.KeyNotFoundError(errormsg)
            self.pars.update(pars)
        return
    
def set_metadata(obj, **kwargs):
    ''' Set standard metadata for an object '''
    obj.created = kwargs.get('created', sc.now())
    obj.version = kwargs.get('version', znv.__version__)
    obj.git_info = kwargs.get('git_info', znm.git_info())
    return

class BaseSim(ParsObj):
    '''
    The BaseSim class stores various methods useful for the Sim that are not directly
    related to simulating the epidemic. It is not used outside of the Sim object,
    so the separation of methods into the BaseSim and Sim classes is purely to keep
    each one of manageable size.
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs) # Initialize and set the parameters as attributes
        return


    def _disp(self):
        '''
        Print a verbose display of the sim object. Used by repr(). See sim.disp()
        for the user version. Equivalent to sc.prettyobj().
        '''
        return sc.prepr(self)


    def _brief(self):
        '''
        Return a one-line description of a sim -- used internally and by repr();
        see sim.brief() for the user version.
        '''
        # Try to get a detailed description of the sim...
        try:
            if self.results_ready:

                results = f'Print Results here' # Placeholder until results are finalized. TODO: Update this
            else:
                results = 'not run'

            # Set label string
            labelstr = f'"{self.label}"' if self.label else '<no label>'

            start = sc.date(self['start_day'], as_date=False)
            if self['end_day']:
                end = sc.date(self['end_day'], as_date=False)
            else:
                end = sc.date(self['n_days'], start_date=start)

            #TODO: Include populations af all species?
            string   = f'Sim({labelstr}; {start} to {end}; epi: {results})'

        # ...but if anything goes wrong, return the default with a warning
        except Exception as E: # pragma: no cover
            string = sc.objectid(self)
            string += f'Warning, sim appears to be malformed; use sim.disp() for details:\n{str(E)}'

        return string


    def update_pars(self, pars=None, create=False, **kwargs):
        ''' Ensure that metaparameters get used properly before being updated '''

        # Merge everything together
        pars = sc.mergedicts(pars, kwargs)
        if pars:

            # Handle other special parameters
            # if pars.get('some_parameter'):
            #     self.some_method(pars['some_parameter'])

            # Call update_pars() for ParsObj
            super().update_pars(pars=pars, create=create)

        return


    def set_metadata(self, simfile):
        ''' Set the metadata for the simulation -- creation time and filename '''
        set_metadata(self)
        if simfile is None:
            self.simfile = 'zoonosim.sim'
        return


    def set_seed(self, seed=-1):
        '''
        Set the seed for the random number stream from the stored or supplied value

        Args:
            seed (None or int): if no argument, use current seed; if None, randomize; otherwise, use and store supplied seed
        '''
        # Unless no seed is supplied, reset it
        if seed != -1:
            self['rand_seed'] = seed
        znu.set_seed(self['rand_seed'])
        return

    @property
    def n(self):
        ''' Count the number of people -- if it fails, assume none '''
        # Im not sure how this needs to change to reflect multiple species.
        # TODO: Revisit this
        return -1 # Placeholder until this is updated
        # try: # By default, the length of the people dict
        #     return len(self.people)
        # except:  # pragma: no cover # If it's None or missing
        #     return 0

    @property
    def scaled_pop_size(self):
        ''' Get the total population size, i.e. the number of agents times the scale factor -- if it fails, assume none '''
        try:
            return self['pop_size']*self['pop_scale']
        except:  # pragma: no cover # If it's None or missing
            return 0

    @property
    def npts(self):
        ''' Count the number of time points '''
        try:
            return int(self['n_days'] + 1)
        except: # pragma: no cover
            return 0

    @property
    def tvec(self):
        ''' Create a time vector '''
        try:
            return np.arange(self.npts)
        except: # pragma: no cover
            return np.array([])

    @property
    def datevec(self):
        '''
        Create a vector of dates

        Returns:
            Array of `datetime` instances containing the date associated with each
            simulation time step

        '''
        try:
            return self['start_day'] + self.tvec * dt.timedelta(days=1)
        except: # pragma: no cover
            return np.array([])


    def day(self, day, *args):
        '''
        Convert a string, date/datetime object, or int to a day (int).

        Args:
            day (str, date, int, or list): convert any of these objects to a day relative to the simulation's start day

        Returns:
            days (int or str): the day(s) in simulation time

        **Example**::

            sim.day('2020-04-05') # Returns 35
        '''
        return sc.day(day, *args, start_date=self['start_day'])


    def date(self, ind, *args, dateformat=None, as_date=False):
        '''
        Convert one or more integer days of simulation time to a date/list of dates --
        by default returns a string, or returns a datetime Date object if as_date is True.
        See also cv.date(), which provides a partly overlapping set of date conversion
        features.

        Args:
            ind (int, list, or array): the index day(s) in simulation time (NB: strings and date objects are accepted, and will be passed unchanged)
            args (list): additional day(s)
            dateformat (str): the format to return the date in
            as_date (bool): whether to return as a datetime date instead of a string

        Returns:
            dates (str, Date, or list): the date(s) corresponding to the simulation day(s)

        **Examples**::

            sim = cv.Sim()
            sim.date(34) # Returns '2020-04-04'
            sim.date([34, 54]) # Returns ['2020-04-04', '2020-04-24']
            sim.date([34, '2020-04-24']) # Returns ['2020-04-04', '2020-04-24']
            sim.date(34, 54, as_date=True) # Returns [datetime.date(2020, 4, 4), datetime.date(2020, 4, 24)]
        '''

        # Handle inputs
        if not isinstance(ind, list): # If it's a number, string, or dateobj, convert it to a list
            ind = sc.promotetolist(ind)
        ind.extend(args)
        if dateformat is None:
            dateformat = '%Y-%m-%d'

        # Do the conversion
        dates = []
        for raw in ind:
            if sc.isnumber(raw):
                date_obj = sc.date(self['start_day'], as_date=True) + dt.timedelta(days=int(raw))
            else:
                date_obj = sc.date(raw, as_date=True)
            if as_date:
                dates.append(date_obj)
            else:
                dates.append(date_obj.strftime(dateformat))

        # Return a string rather than a list if only one provided
        if len(ind)==1:
            dates = dates[0]

        return dates


    def result_keys(self, which='main'):
        '''
        Get the actual results objects, not other things stored in sim.results.

        If which is 'main', return only the main results keys. If 'variant', return
        only variant keys. If 'all', return all keys.

        '''
        keys = []
        choices = ['main', 'variant', 'all']
        if which in ['main', 'all']:
            keys += [key for key,res in self.results.items() if isinstance(res, Result)]
        if which in ['variant', 'all'] and 'variant' in self.results:
            keys += [key for key,res in self.results['variant'].items() if isinstance(res, Result)]
        if which not in choices: # pragma: no cover
            errormsg = f'Choice "which" not available; choices are: {sc.strjoin(choices)}'
            raise ValueError(errormsg)
        return keys


    def copy(self):
        ''' Returns a deep copy of the sim '''
        return sc.dcp(self)


    def export_results(self, for_json=True, filename=None, indent=2, *args, **kwargs):
        '''
        Convert results to dict -- see also to_json().

        The results written to Excel must have a regular table shape, whereas
        for the JSON output, arbitrary data shapes are supported.

        Args:
            for_json (bool): if False, only data associated with Result objects will be included in the converted output
            filename (str): filename to save to; if None, do not save
            indent (int): indent (int): if writing to file, how many indents to use per nested level
            args (list): passed to savejson()
            kwargs (dict): passed to savejson()

        Returns:
            resdict (dict): dictionary representation of the results

        '''

        if not self.results_ready: # pragma: no cover
            errormsg = 'Please run the sim before exporting the results'
            raise RuntimeError(errormsg)

        resdict = {}
        resdict['t'] = self.results['t'] # Assume that there is a key for time

        if for_json:
            resdict['timeseries_keys'] = self.result_keys()
        for key,res in self.results.items():
            if isinstance(res, Result):
                resdict[key] = res.values
                if res.low is not None:
                    resdict[key+'_low'] = res.low
                if res.high is not None:
                    resdict[key+'_high'] = res.high
            elif for_json:
                if key == 'date':
                    resdict[key] = [str(d) for d in res] # Convert dates to strings
                else:
                    resdict[key] = res
        if filename is not None:
            sc.savejson(filename=filename, obj=resdict, indent=indent, *args, **kwargs)
        return resdict


    def export_pars(self, filename=None, indent=2, *args, **kwargs):
        '''
        Return parameters for JSON export -- see also to_json().

        This method is required so that interventions can specify
        their JSON-friendly representation.

        Args:
            filename (str): filename to save to; if None, do not save
            indent (int): indent (int): if writing to file, how many indents to use per nested level
            args (list): passed to savejson()
            kwargs (dict): passed to savejson()

        Returns:
            pardict (dict): a dictionary containing all the parameter values
        '''
        pardict = {}
        for key in self.pars.keys():
            if key == 'interventions':
                pardict[key] = [intervention.to_json() for intervention in self.pars[key]]
            elif key == 'start_day':
                pardict[key] = str(self.pars[key])
            else:
                pardict[key] = self.pars[key]
        if filename is not None:
            sc.savejson(filename=filename, obj=pardict, indent=indent, *args, **kwargs)
        return pardict


    def to_json(self, filename=None, keys=None, tostring=False, indent=2, verbose=False, *args, **kwargs):
        '''
        Export results and parameters as JSON.

        Args:
            filename (str): if None, return string; else, write to file
            keys (str or list): attributes to write to json (default: results, parameters, and summary)
            tostring (bool): if not writing to file, whether to write to string (alternative is sanitized dictionary)
            indent (int): if writing to file, how many indents to use per nested level
            verbose (bool): detail to print
            args (list): passed to savejson()
            kwargs (dict): passed to savejson()

        Returns:
            A unicode string containing a JSON representation of the results,
            or writes the JSON file to disk

        **Examples**::

            json = sim.to_json()
            sim.to_json('results.json')
            sim.to_json('summary.json', keys='summary')
        '''

        # Handle keys
        if keys is None:
            keys = ['results', 'pars', 'summary']
        keys = sc.promotetolist(keys)

        # Convert to JSON-compatible format
        d = {}
        for key in keys:
            if key == 'results':
                resdict = self.export_results(for_json=True)
                d['results'] = resdict
            elif key in ['pars', 'parameters']:
                pardict = self.export_pars()
                d['parameters'] = pardict
            elif key == 'summary':
                d['summary'] = dict(sc.dcp(self.summary))
            else: # pragma: no cover
                try:
                    d[key] = sc.sanitizejson(getattr(self, key))
                except Exception as E:
                    errormsg = f'Could not convert "{key}" to JSON: {str(E)}; continuing...'
                    print(errormsg)

        if filename is None:
            output = sc.jsonify(d, tostring=tostring, indent=indent, verbose=verbose, *args, **kwargs)
        else:
            output = sc.savejson(filename=filename, obj=d, indent=indent, *args, **kwargs)

        return output


    def to_df(self, date_index=False):
        '''
        Export results to a pandas dataframe

        Args:
            date_index  (bool): if True, use the date as the index
        '''
        resdict = self.export_results(for_json=False)
        df = pd.DataFrame.from_dict(resdict)
        df['date'] = self.datevec
        new_columns = ['t','date'] + df.columns[1:-1].tolist() # Get column order
        df = df.reindex(columns=new_columns) # Reorder so 't' and 'date' are first
        if date_index:
            df = df.set_index('date')
        return df


    def to_excel(self, filename=None, skip_pars=None):
        '''
        Export parameters and results as Excel format

        Args:
            filename  (str): if None, return string; else, write to file
            skip_pars (list): if provided, a custom list parameters to exclude

        Returns:
            An sc.Spreadsheet with an Excel file, or writes the file to disk
        '''
        if skip_pars is None:
            skip_pars = ['variant_map', 'vaccine_map'] # These include non-string keys so fail at sc.flattendict()

        # Export results
        result_df = self.to_df(date_index=True)

        # Export parameters
        pars = {k:v for k,v in self.pars.items() if k not in skip_pars}
        par_df = pd.DataFrame.from_dict(sc.flattendict(pars, sep='_'), orient='index', columns=['Value'])
        par_df.index.name = 'Parameter'

        # Convert to spreadsheet
        spreadsheet = sc.Spreadsheet()
        spreadsheet.freshbytes()
        with pd.ExcelWriter(spreadsheet.bytes, engine='xlsxwriter') as writer:
            result_df.to_excel(writer, sheet_name='Results')
            par_df.to_excel(writer, sheet_name='Parameters')
        spreadsheet.load()

        if filename is None:
            output = spreadsheet
        else:
            output = spreadsheet.save(filename)

        return output


    def shrink(self, skip_attrs=None, in_place=True):
        '''
        "Shrinks" the simulation by removing the people and other memory-intensive
        attributes (e.g., some interventions and analyzers), and returns a copy of
        the "shrunken" simulation. Used to reduce the memory required for RAM or
        for saved files.

        Args:
            skip_attrs (list): a list of attributes to skip (remove) in order to perform the shrinking; default "people"
            in_palce (bool): whether to perform the shrinking in place (default), or return a shrunken copy instead

        Returns:
            shrunken (Sim): a Sim object with the listed attributes removed
        '''
        # from . import interventions as zni # To avoid circular imports
        # from . import analysis as zna

        # By default, skip people (~90% of memory), the popdict (which is usually empty anyway), and _orig_pars (which is just a backup)
        if skip_attrs is None:
            skip_attrs = ['popdict', 'people', '_orig_pars'] # TODO: Modify these to reflect the actual attributes to skip. Ideally this should be set in defaults.py
        # Create the new object, and copy original dict, skipping the skipped attributes
        if in_place:
            shrunken = self
            for attr in skip_attrs:
                setattr(self, attr, None)
        else:
            shrunken = object.__new__(self.__class__)
            shrunken.__dict__ = {k:(v if k not in skip_attrs else None) for k,v in self.__dict__.items()}

        # Shrink interventions and analyzers, with a lot of checking along the way
        # for key in ['interventions', 'analyzers']:
        #     ias = self.pars[key] # List of interventions or analyzers
        #     shrunken_ias = [ia.shrink(in_place=in_place) for ia in ias if isinstance(ia, (cvi.Intervention, cva.Analyzer))]
        #     self.pars[key] = shrunken_ias # Actually shrink, and re-store

        # Don't return if in place
        if in_place:
            return
        else:
            return shrunken


    def save(self, filename=None, keep_people=None, skip_attrs=None, **kwargs):
        '''
        Save to disk as a gzipped pickle.

        Args:
            filename (str or None): the name or path of the file to save to; if None, uses stored
            kwargs: passed to sc.makefilepath()

        Returns:
            filename (str): the validated absolute path to the saved file

        **Example**::

            sim.save() # Saves to a .sim file
        '''

        # Set keep_people based on whether or not we're in the middle of a run
        if keep_people is None:
            if self.initialized and not self.results_ready:
                keep_people = True
            else:
                keep_people = False

        # Handle the filename
        if filename is None:
            filename = self.simfile
        filename = sc.makefilepath(filename=filename, **kwargs)
        self.filename = filename # Store the actual saved filename

        # Handle the shrinkage and save
        if skip_attrs or not keep_people:
            obj = self.shrink(skip_attrs=skip_attrs, in_place=False)
        else:
            obj = self
        znm.save(filename=filename, obj=obj)

        return filename


    @staticmethod
    def load(filename, *args, **kwargs):
        '''
        Load from disk from a gzipped pickle.

        Args:
            filename (str): the name or path of the file to load from
            kwargs: passed to cv.load()

        Returns:
            sim (Sim): the loaded simulation object

        **Example**::

            sim = cv.Sim.load('my-simulation.sim')
        '''
        sim = znm.load(filename, *args, **kwargs)
        if not isinstance(sim, BaseSim): # pragma: no cover
            errormsg = f'Cannot load object of {type(sim)} as a Sim object'
            raise TypeError(errormsg)
        return sim


    def _get_ia(self, which, label=None, partial=False, as_list=False, as_inds=False, die=True, first=False):
        ''' Helper method for get_interventions() and get_analyzers(); see get_interventions() docstring '''

        # Handle inputs
        if which not in ['interventions', 'analyzers']: # pragma: no cover
            errormsg = f'This method is only defined for interventions and analyzers, not "{which}"'
            raise ValueError(errormsg)

        ia_list = sc.tolist(self.pars[which]) # List of interventions or analyzers
        n_ia = len(ia_list) # Number of interventions/analyzers

        if label == 'summary': # Print a summary of the interventions
            df = pd.DataFrame(columns=['ind', 'label', 'type'])
            for ind,ia_obj in enumerate(ia_list):
                df = df.append(dict(ind=ind, label=str(ia_obj.label), type=type(ia_obj)), ignore_index=True)
            print(f'Summary of {which}:')
            print(df)
            return

        else: # Standard usage case
            position = 0 if first else -1 # Choose either the first or last element
            if label is None: # Get all interventions if no label is supplied, e.g. sim.get_interventions()
                label = np.arange(n_ia)
            if isinstance(label, np.ndarray): # Allow arrays to be provided
                label = label.tolist()
            labels = sc.promotetolist(label)

            # Calculate the matches
            matches = []
            match_inds = []
            for label in labels:
                if sc.isnumber(label):
                    matches.append(ia_list[label]) # This will raise an exception if an invalid index is given
                    label = n_ia + label if label<0 else label # Convert to a positive number
                    match_inds.append(label)
                elif sc.isstring(label) or isinstance(label, type):
                    for ind,ia_obj in enumerate(ia_list):
                        if sc.isstring(label) and ia_obj.label == label or (partial and (label in str(ia_obj.label))):
                            matches.append(ia_obj)
                            match_inds.append(ind)
                        elif isinstance(label, type) and isinstance(ia_obj, label):
                            matches.append(ia_obj)
                            match_inds.append(ind)
                else: # pragma: no cover
                    errormsg = f'Could not interpret label type "{type(label)}": should be str, int, list, or {which} class'
                    raise TypeError(errormsg)

            # Parse the output options
            if as_inds:
                output = match_inds
            elif as_list: # Used by get_interventions()
                output = matches
            else:
                if len(matches) == 0: # pragma: no cover
                    if die:
                        errormsg = f'No {which} matching "{label}" were found'
                        raise ValueError(errormsg)
                    else:
                        output = None
                else:
                    output = matches[position] # Return either the first or last match (usually), used by get_intervention()

            return output


    def get_interventions(self, label=None, partial=False, as_inds=False):
        '''
        Find the matching intervention(s) by label, index, or type. If None, return
        all interventions. If the label provided is "summary", then print a summary
        of the interventions (index, label, type).

        Args:
            label (str, int, Intervention, list): the label, index, or type of intervention to get; if a list, iterate over one of those types
            partial (bool): if true, return partial matches (e.g. 'beta' will match all beta interventions)
            as_inds (bool): if true, return matching indices instead of the actual interventions

        **Examples**::

            tp = cv.test_prob(symp_prob=0.1)
            cb1 = cv.change_beta(days=5, changes=0.3, label='NPI')
            cb2 = cv.change_beta(days=10, changes=0.3, label='Masks')
            sim = cv.Sim(interventions=[tp, cb1, cb2])
            cb1, cb2 = sim.get_interventions(cv.change_beta)
            tp, cb2 = sim.get_interventions([0,2])
            ind = sim.get_interventions(cv.change_beta, as_inds=True) # Returns [1,2]
            sim.get_interventions('summary') # Prints a summary
        '''
        return self._get_ia('interventions', label=label, partial=partial, as_inds=as_inds, as_list=True)


    def get_intervention(self, label=None, partial=False, first=False, die=True):
        '''
        Like get_interventions(), find the matching intervention(s) by label,
        index, or type. If more than one intervention matches, return the last
        by default. If no label is provided, return the last intervention in the list.

        Args:
            label (str, int, Intervention, list): the label, index, or type of intervention to get; if a list, iterate over one of those types
            partial (bool): if true, return partial matches (e.g. 'beta' will match all beta interventions)
            first (bool): if true, return first matching intervention (otherwise, return last)
            die (bool): whether to raise an exception if no intervention is found

        **Examples**::

            tp = cv.test_prob(symp_prob=0.1)
            cb = cv.change_beta(days=5, changes=0.3, label='NPI')
            sim = cv.Sim(interventions=[tp, cb])
            cb = sim.get_intervention('NPI')
            cb = sim.get_intervention('NP', partial=True)
            cb = sim.get_intervention(cv.change_beta)
            cb = sim.get_intervention(1)
            cb = sim.get_intervention()
            tp = sim.get_intervention(first=True)
        '''
        return self._get_ia('interventions', label=label, partial=partial, first=first, die=die, as_inds=False, as_list=False)


    def get_analyzers(self, label=None, partial=False, as_inds=False):
        '''
        Same as get_interventions(), but for analyzers.
        '''
        return self._get_ia('analyzers', label=label, partial=partial, as_list=True, as_inds=as_inds)


    def get_analyzer(self, label=None, partial=False, first=False, die=True):
        '''
        Same as get_intervention(), but for analyzers.
        '''
        return self._get_ia('analyzers', label=label, partial=partial, first=first, die=die, as_inds=False, as_list=False)