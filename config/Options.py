#%% Define the options class
'''
Define options for Zoonosim, mostly plotting and Numba options. All options should
be set using set() or directly, e.g.::

    zn.options(font_size=18)

To reset default options, use::

    zn.options('default')

Note: "options" is used to refer to the choices available (e.g., DPI), while "settings"
is used to refer to the choices made (e.g., DPI=150).
'''

import os
import sciris as sc

__all__ = ['options']

class Options(sc.objdict):
    '''
    Set options for Zoonosim.

    Use ``zn.options.set('defaults')`` to reset all values to default, or ``zn.options.set(dpi='default')``
    to reset one parameter to default. See ``zn.options.help(detailed=True)`` for
    more information.

    Options can also be saved and loaded using ``zn.options.save()`` and ``zn.options.load()``.
    See ``zn.options.context()`` and ``zn.options.with_style()`` to set options
    temporarily.

    Common options are (see also ``zn.options.help(detailed=True)``):

        - verbose:        default verbosity for simulations to use
        - warnings:       how to handle warnings (e.g. print, raise as errors, ignore)

    **Examples**::
        zn.options('defaults') # Reset to default options
    '''

    def __init__(self):
        super().__init__()
        optdesc, options = self.get_orig_options() # Get the options
        self.update(options) # Update this object with them
        self.setattribute('optdesc', optdesc) # Set the description as an attribute, not a dict entry
        self.setattribute('orig_options', sc.dcp(options)) # Copy the default options
        return


    def __call__(self, *args, **kwargs):
        '''Allow ``zn.options(dpi=150)`` instead of ``zn.options.set(dpi=150)`` '''
        return self.set(*args, **kwargs)


    def to_dict(self):
        ''' Pull out only the settings from the options object '''
        return {k:v for k,v in self.items()}


    def __repr__(self):
        ''' Brief representation '''
        output = sc.objectid(self)
        output += 'Zoonosim options (see also zn.options.disp()):\n'
        output += sc.pp(self.to_dict(), output=True)
        return output


    def __enter__(self):
        ''' Allow to be used in a with block '''
        return self


    def __exit__(self, *args, **kwargs):
        ''' Allow to be used in a with block '''
        try:
            reset = {}
            for k,v in self.on_entry.items():
                if self[k] != v: # Only reset settings that have changed
                    reset[k] = v
            self.set(**reset)
            self.delattribute('on_entry')
        except AttributeError as E:
            errormsg = 'Please use zn.options.context() if using a with block'
            raise AttributeError(errormsg) from E
        return


    def disp(self):
        ''' Detailed representation '''
        output = 'Zoonosim options (see also zn.options.help()):\n'
        keylen = 14 # Maximum key length  -- "numba_parallel"
        for k,v in self.items():
            keystr = sc.colorize(f'  {k:>{keylen}s}: ', fg='cyan', output=True)
            reprstr = sc.pp(v, output=True)
            reprstr = sc.indent(n=keylen+4, text=reprstr, width=None)
            output += f'{keystr}{reprstr}'
        print(output)
        return


    @staticmethod
    def get_orig_options():
        '''
        Set the default options for Zoonosim -- not to be called by the user, use
        ``zn.options.set('defaults')`` instead.
        '''

        # Options acts like a class, but is actually an objdict for simplicity
        optdesc = sc.objdict() # Help for the options
        options = sc.objdict() # The options

        optdesc.verbose = 'Set default level of verbosity (i.e. logging detail): e.g., 0.1 is an update every 10 simulated days'
        options.verbose = float(os.getenv('ZOONOSIM_VERBOSE', 0.1))

        optdesc.warnings = 'How warnings are handled: options are "warn" (default), "print", and "error"'
        options.warnings = str(os.getenv('ZOONOSIM_WARNINGS', 'warn'))

        optdesc.sep = 'Set thousands seperator for text output'
        options.sep = str(os.getenv('ZOONOSIM_SEP', ','))

        optdesc.precision = 'Set arithmetic precision for Numba -- 32-bit by default for efficiency'
        options.precision = int(os.getenv('ZOONOSIM_PRECISION', 32))

        optdesc.numba_parallel = 'Set Numba multithreading -- none, safe, full; full multithreading is ~20% faster, but results become nondeterministic'
        options.numba_parallel = str(os.getenv('ZOONOSIM_NUMBA_PARALLEL', 'none'))

        optdesc.numba_cache = 'Set Numba caching -- saves on compilation time; disabling is not recommended'
        options.numba_cache = bool(int(os.getenv('ZOONOSIM_NUMBA_CACHE', 1)))

        return optdesc, options


    def set(self, key=None, value=None, **kwargs):
        '''
        Actually change the style. See ``zn.options.help()`` for more information.

        Args:
            key    (str):    the parameter to modify, or 'defaults' to reset everything to default values
            value  (varies): the value to specify; use None or 'default' to reset to default
            kwargs (dict):   if supplied, set multiple key-value pairs

        **Example**::

            zn.options.set(dpi=50) # Equivalent to zn.options(dpi=50)
        '''

        reload_required = False

        # Reset to defaults
        if key in ['default', 'defaults']:
            kwargs = self.orig_options # Reset everything to default

        # Handle other keys
        elif key is not None:
            kwargs = sc.mergedicts(kwargs, {key:value})


        # Reset options
        for key,value in kwargs.items():

            if key not in self:
                keylist = self.orig_options.keys()
                keys = '\n'.join(keylist)
                errormsg = f'Option "{key}" not recognized; options are "defaults" or:\n{keys}\n\nSee help(zn.options.set) for more information.'
                raise sc.KeyNotFoundError(errormsg)
            else:
                if value in [None, 'default']:
                    value = self.orig_options[key]
                self[key] = value
                numba_keys = ['precision', 'numba_parallel', 'numba_cache'] # Specify which keys require a reload
                if key in numba_keys:
                    reload_required = True

        if reload_required:
            reload_numba()

        return


    def context(self, **kwargs):
        '''
        Alias to set() for non-plotting options, for use in a "with" block.

        Note: for plotting options, use ``zn.options.with_style()``, which is linked
        to Matplotlib's context manager. If you set plotting options with this,
        they won't have any effect.

        **Examples**::

            # Silence all output
            with zn.options.context(verbose=0):
                zn.Sim().run()

            # Convert warnings to errors
            with zn.options.context(warnings='error'):
                zn.Sim(location='not a location').initialize()

            # Use with_style(), not context(), for plotting options
            with zn.options.with_style(dpi=50):
                zn.Sim().run().plot()

        New in version 3.1.2.
        '''

        # Store current settings
        on_entry = {k:self[k] for k in kwargs.keys()}
        self.setattribute('on_entry', on_entry)

        # Make changes
        self.set(**kwargs)
        return self


    def get_default(self, key):
        ''' Helper function to get the original default options '''
        return self.orig_options[key]


    def changed(self, key):
        ''' Check if current setting has been changed from default '''
        if key in self.orig_options:
            return self[key] != self.orig_options[key]
        else:
            return None


    def help(self, detailed=False, output=False):
        '''
        Print information about options.

        Args:
            detailed (bool): whether to print out full help
            output (bool): whether to return a list of the options

        **Example**::

            zn.options.help(detailed=True)
        '''

        # If not detailed, just print the docstring for zn.options
        if not detailed:
            print(self.__doc__)
            return

        n = 15 # Size of indent
        optdict = sc.objdict()
        for key in self.orig_options.keys():
            entry = sc.objdict()
            entry.key = key
            entry.current = sc.indent(n=n, width=None, text=sc.pp(self[key], output=True)).rstrip()
            entry.default = sc.indent(n=n, width=None, text=sc.pp(self.orig_options[key], output=True)).rstrip()
            if not key.startswith('rc'):
                entry.variable = f'ZOONOSIM_{key.upper()}' # NB, hard-coded above!
            else:
                entry.variable = 'No environment variable'
            entry.desc = sc.indent(n=n, text=self.optdesc[key])
            optdict[key] = entry

        # Convert to a dataframe for nice printing
        print('Zoonosim global options ("Environment" = name of corresponding environment variable):')
        for k, key, entry in optdict.enumitems():
            sc.heading(f'{k}. {key}', spaces=0, spacesafter=0)
            changestr = '' if entry.current == entry.default else ' (modified)'
            print(f'          Key: {key}')
            print(f'      Current: {entry.current}{changestr}')
            print(f'      Default: {entry.default}')
            print(f'  Environment: {entry.variable}')
            print(f'  Description: {entry.desc}')

        sc.heading('Methods:', spacesafter=0)
        print('''
    zn.options(key=value) -- set key to value
    zn.options[key] -- get or set key
    zn.options.set() -- set option(s)
    zn.options.get_default() -- get default setting(s)
    zn.options.load() -- load settings from file
    zn.options.save() -- save settings to file
    zn.options.to_dict() -- convert to dictionary
    zn.options.style() -- create style context for plotting
''')

        if output:
            return optdict
        else:
            return


    def load(self, filename, verbose=True, **kwargs):
        '''
        Load current settings from a JSON file.

        Args:
            filename (str): file to load
            kwargs (dict): passed to ``sc.loadjson()``
        '''
        json = sc.loadjson(filename=filename, **kwargs)
        current = self.to_dict()
        new = {k:v for k,v in json.items() if v != current[k]} # Don't reset keys that haven't changed
        self.set(**new)
        if verbose: print(f'Settings loaded from {filename}')
        return


    def save(self, filename, verbose=True, **kwargs):
        '''
        Save current settings as a JSON file.

        Args:
            filename (str): file to save to
            kwargs (dict): passed to ``sc.savejson()``
        '''
        json = self.to_dict()
        output = sc.savejson(filename=filename, obj=json, **kwargs)
        if verbose: print(f'Settings saved to {filename}')
        return output




def reload_numba():
    '''
    Apply changes to Numba functions -- reloading modules is necessary for
    changes to propagate. Not necessary to call directly if zn.options.set() is used.

    **Example**::

        import zoonosim as zn
        zn.options.set(precision=64)
        sim = zn.Sim()
        sim.run()
        assert sim.people.rel_trans.dtype == np.float64
    '''
    print('Reloading Zoonosim so changes take effect...')
    import importlib
    import Zoonosim as zn
    importlib.reload(zn.config.DataTypes)
    importlib.reload(zn.utils)
    importlib.reload(zn)
    print("Reload complete. Note: for some options to take effect, you may also need to delete Zoonosim's __pycache__ folder.")
    return



# Create the options on module load, and load the fonts
options = Options()