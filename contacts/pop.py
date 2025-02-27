"""
This module provides the main class for interacting with SynthPops, the Pop class.
"""

import numpy as np
import sciris as sc

from ..config import Defaults as znd
from ..utils import sample_ops as znso

from .. import base as znb

__all__ = ['Pop']

class Pop(sc.prettyobj):

    def __init__(self,
                 n_people=None,
                 n_poultry=None,
                 max_farm_contacts=None,
                 farm_pars=None,
                 max_water_contacts=None,
                 water_pars=None,
                 max_equipment_contacts=None,
                 equipment_pars=None,
                 rand_seed=None,
                 smooth_ages=False,
                 window_length=7,
                 do_make=True
                 ):
        '''
        Make a full population network including both people, poultry, and
        contacts.

        Args:
            n_people (int)                          : The number of human agents to be created
            n_poultry (int)                         : The number of poultry flocks to be created
            max_contacts (dict)                     : A dictionary for maximum number of contacts per layer..
            rand_seed (int)                         : Start point random sequence is generated from.
            smooth_ages (bool)                      : If True, use smoothed out age distribution.
            window_length (int)                     : length of window over which to average or smooth out age distribution
            do_make (bool)                          : whether to make the population

        Returns:
            network (dict): A dictionary of the full population with ages, connections, and other attributes.
        '''

        # General parameters
        if n_people is None:
            #log.warning(f"Pop size n not given, generating a population with a default size of {defaults.default_pop_size} people.")
            n_people = znd.default_human_pop
        if n_poultry is None:
            n_poultry = znd.default_poultry_pop

        # Assign all the variables
        self.farm_pars          = sc.objdict()
        self.water_pars         = sc.objdict()
        self.equipment_pars     = sc.objdict()

        self.n_people           = int(n_people)
        self.n_poultry          = int(n_poultry)
        self.max_contacts       = sc.mergedicts({'F': max_farm_contacts}, 
                                                {'W': max_water_contacts},
                                                {'E': max_equipment_contacts}) 
        self.rand_seed          = rand_seed

        # Age distribution parameters
        self.smooth_ages                                 = smooth_ages
        self.window_length                               = window_length

        # Farm parameters
        # If any parameters are supplied as a dict to override defaults, merge them in now
        self.farm_pars = sc.objdict(sc.mergedicts(self.farm_pars, farm_pars))

        # Water parameters
        # If any parameters are supplied as a dict to override defaults, merge them in now
        self.water_pars = sc.objdict(sc.mergedicts(self.water_pars, water_pars))

        # Equipment parameters
        # If any parameters are supplied as a dict to override defaults, merge them in now
        self.equipment_pars = sc.objdict(sc.mergedicts(self.equipment_pars, equipment_pars))

        # what are the layers generated?
        self.layers = ['F', 'W', 'E']
        self.layer_mappings = dict(F = 'Farms', W = 'Water', E = 'Equipment')

        # Handle the seed
        if self.rand_seed is not None:
            znso.set_seed(self.rand_seed) # TODO: Decide where this method should be stored

        # Heavy lift: make the contacts and their connections
        population = self.generate()
        self.popdict = population

        # Add summaries post hoc  --- TBD: summaries during generation
        self.compute_information()  # compute full information
        self.compute_summary()  # then compute condensed summary

        return

    def generate(self):
        """
        Actually generate the network.

        Returns:
            network (dict): A dictionary of the full population with ages, connections, and other attributes.
        """

        # TODO: unpack variables -- to be refactored to pass parameters directly

        # General parameters

        n_people                        = self.n_people
        n_poultry                       = self.n_poultry
        max_contacts                    = self.max_contacts


        # Age distribution parameters
        smooth_ages                     = self.smooth_ages
        window_length                   = self.window_length

        # Household parameters
        household_method                = self.household_method


        # School parameters
        with_school_types               = self.school_pars.with_school_types
        school_mixing_type              = self.school_pars.school_mixing_type
        average_class_size              = self.school_pars.average_class_size
        inter_grade_mixing              = self.school_pars.inter_grade_mixing
        average_student_teacher_ratio   = self.school_pars.average_student_teacher_ratio
        average_teacher_teacher_degree  = self.school_pars.average_teacher_teacher_degree
        teacher_age_min                 = self.school_pars.teacher_age_min
        teacher_age_max                 = self.school_pars.teacher_age_max
        with_non_teaching_staff         = self.school_pars.with_non_teaching_staff
        average_student_all_staff_ratio = self.school_pars.average_student_all_staff_ratio
        average_additional_staff_degree = self.school_pars.average_additional_staff_degree
        staff_age_min                   = self.school_pars.staff_age_min
        staff_age_max                   = self.school_pars.staff_age_max

        # Load and store the expected age distribution of the population
        age_bracket_dist = spdata.read_age_bracket_distr(**loc_pars)  # age distribution defined by bins or age brackets
        expected_age_dist = spdata.get_smoothed_single_year_age_distr(**loc_pars, window_length=self.window_length)
        self.expected_age_dist = expected_age_dist
        expected_age_dist_values = [expected_age_dist[a] for a in expected_age_dist]
        self.expected_age_dist_values = expected_age_dist_values

        # Load and store the age brackets
        age_brackets = spdata.get_census_age_brackets(**loc_pars)
        self.age_brackets = age_brackets
        # mapping
        age_by_brackets = spb.get_age_by_brackets(age_brackets)
        self.age_by_brackets = age_by_brackets

        # Load the contact matrix
        contact_matrices = spdata.get_contact_matrices(datadir, sheet_name=sheet_name)
        # Store expected contact matrices
        self.contact_matrices = contact_matrices

        # Load age brackets, and mapping dictionary that matches contact matrices
        contact_matrix_shape = contact_matrices[list(contact_matrices.keys())[0]].shape
        contact_matrix_row = contact_matrix_shape[0]

        cm_age_brackets = spdata.get_census_age_brackets(**loc_pars, nbrackets=contact_matrix_row)
        self.cm_age_brackets = cm_age_brackets
        cm_age_by_brackets = spb.get_age_by_brackets(cm_age_brackets)
        self.cm_age_by_brackets = cm_age_by_brackets

        # Generate an age count for the population --- this will get passed around to methods generating the different layers where people live: long term care facilities, households, agricultural living quarters, other group living arrangements
        age_count = znhh.generate_age_count_multinomial(n, expected_age_dist_values)

        # Ages left to assign to a residence
        ages_left_to_assign = sc.dcp(age_count)

        # Generate LTCFs and remove some people from the age count of people left to place in a resident by age
        n_nonltcf, ltcf_adjusted_age_dist, ltcf_adjusted_age_dist_values, ages_left_to_assign, facilities = spltcf.generate_ltcfs(n, with_facilities, loc_pars, expected_age_dist, ages_left_to_assign)

        # Generate households
        household_size_dist = spdata.get_household_size_distr(**loc_pars)
        hh_sizes = znhh.generate_household_size_count_from_fixed_pop_size(n_nonltcf, household_size_dist)
        hha_brackets = spdata.get_head_age_brackets(**loc_pars)
        hha_by_size = spdata.get_head_age_by_size_distr(**loc_pars)

        if household_method == 'fixed_ages':

            homes_dic, homes = znhh.generate_all_households_fixed_ages(n_nonltcf, hh_sizes, hha_by_size, hha_brackets, cm_age_brackets, cm_age_by_brackets, contact_matrices, ages_left_to_assign)

        else:
            homes_dic, homes = znhh.generate_all_households_infer_ages(n, n_nonltcf, hh_sizes, hha_by_size, hha_brackets, cm_age_brackets, cm_age_by_brackets, contact_matrices, ltcf_adjusted_age_dist, ages_left_to_assign)

        # Handle homes and facilities
        homes_by_uids, age_by_uid = znhh.assign_uids_by_homes(homes)  # include facilities to assign ids
        age_by_uid_arr = np.array([age_by_uid[i] for i in range(self.n)], dtype=int)
        self.age_by_uid = age_by_uid_arr

        facilities_by_uid_lists = homes_by_uids[0:len(facilities)]

        # Generate school sizes
        school_sizes_dist_by_brackets = spdata.get_school_size_distr_by_brackets(**loc_pars)  # without school type
        school_size_brackets = spdata.get_school_size_brackets(**loc_pars)  # for right now the size distribution for all school types will use the same brackets or bins

        # Figure out who's going to school as a student with enrollment rates (gets called inside sp.get_uids_in_school)
        uids_in_school, uids_in_school_by_age, ages_in_school_count = znsch.get_uids_in_school(datadir, n_nonltcf, location, state_location, country_location, age_by_uid, homes_by_uids, use_default=use_default)  # this will call in school enrollment rates

        if with_school_types:
            school_size_distr_by_type = spdata.get_school_size_distr_by_type(**loc_pars)

            school_type_age_ranges = spdata.get_school_type_age_ranges(**loc_pars)

            school_types_distr_by_age = znsch.get_school_types_distr_by_age(school_type_age_ranges)
            school_type_by_age = znsch.get_school_types_by_age_single(school_types_distr_by_age)

            student_age_lists, student_uid_lists, school_types = znsch.send_students_to_school_with_school_types(school_size_distr_by_type,
                                                                                                                 school_size_brackets,
                                                                                                                 uids_in_school,
                                                                                                                 uids_in_school_by_age,
                                                                                                                 ages_in_school_count,
                                                                                                                 school_types_distr_by_age,
                                                                                                                 school_type_age_ranges)

        else:
            # Get school sizes
            school_sizes = znsch.generate_school_sizes(school_sizes_dist_by_brackets, school_size_brackets, uids_in_school)

            # Assign students to school using contact matrix method - generic schools
            student_age_lists, student_uid_lists, school_types = znsch.send_students_to_school(school_sizes,
                                                                                               uids_in_school,
                                                                                               uids_in_school_by_age,
                                                                                               ages_in_school_count,
                                                                                               cm_age_brackets,
                                                                                               cm_age_by_brackets,
                                                                                               contact_matrices)

            school_type_by_age = None

        # Get employment rates
        employment_rates = spdata.get_employment_rates(**loc_pars)

        # Find people who can be workers (removing everyone who is currently a student)
        uids_by_age = spb.get_ids_by_age(age_by_uid)  # Make a dictionary listing out uids of people by their age
        potential_worker_uids, potential_worker_uids_by_age, potential_worker_ages_left_count = znwp.get_uids_potential_workers(student_uid_lists,
                                                                                                                               employment_rates,
                                                                                                                               age_by_uid)
        workers_by_age_to_assign_count = znwp.get_workers_by_age_to_assign(employment_rates, potential_worker_ages_left_count, uids_by_age)

        # Assign teachers and update school lists
        teacher_age_lists, teacher_uid_lists, potential_worker_uids, potential_worker_uids_by_age, workers_by_age_to_assign_count = znsch.assign_teachers_to_schools(student_age_lists,
                                                                                                                                                                     student_uid_lists,
                                                                                                                                                                     employment_rates,
                                                                                                                                                                     workers_by_age_to_assign_count,
                                                                                                                                                                     potential_worker_uids,
                                                                                                                                                                     potential_worker_uids_by_age,
                                                                                                                                                                     potential_worker_ages_left_count,
                                                                                                                                                                     average_student_teacher_ratio=average_student_teacher_ratio,
                                                                                                                                                                     teacher_age_min=teacher_age_min,
                                                                                                                                                                     teacher_age_max=teacher_age_max)
        # Assign non teaching staff and update who's available to work at other places
        non_teaching_staff_uid_lists, potential_worker_uids, potential_worker_uids_by_age, workers_by_age_to_assign_count = znsch.assign_additional_staff_to_schools(student_uid_lists,
                                                                                                                                                                     teacher_uid_lists,
                                                                                                                                                                     workers_by_age_to_assign_count,
                                                                                                                                                                     potential_worker_uids,
                                                                                                                                                                     potential_worker_uids_by_age,
                                                                                                                                                                     potential_worker_ages_left_count,
                                                                                                                                                                     average_student_teacher_ratio=average_student_teacher_ratio,
                                                                                                                                                                     average_student_all_staff_ratio=average_student_all_staff_ratio,
                                                                                                                                                                     staff_age_min=staff_age_min,
                                                                                                                                                                     staff_age_max=staff_age_max,
                                                                                                                                                                     with_non_teaching_staff=with_non_teaching_staff)


        # Generate non-school workplace sizes needed to send everyone to work
        workplace_size_brackets = spdata.get_workplace_size_brackets(**loc_pars)
        workplace_size_distr_by_brackets = spdata.get_workplace_size_distr_by_brackets(**loc_pars)
        workplace_sizes = znwp.generate_workplace_sizes(workplace_size_distr_by_brackets, workplace_size_brackets, workers_by_age_to_assign_count)

        # Assign all workers who are not staff at schools to workplaces
        workplace_age_lists, workplace_uid_lists, potential_worker_uids, potential_worker_uids_by_age, workers_by_age_to_assign_count = spw.assign_rest_of_workers(workplace_sizes,
                                                                                                                                                                   potential_worker_uids,
                                                                                                                                                                   potential_worker_uids_by_age,
                                                                                                                                                                   workers_by_age_to_assign_count,
                                                                                                                                                                   age_by_uid,
                                                                                                                                                                   cm_age_brackets,
                                                                                                                                                                   cm_age_by_brackets,
                                                                                                                                                                   contact_matrices)


        population = spcnx.make_contacts(self,
                                         age_by_uid=age_by_uid,
                                         homes_by_uids=homes_by_uids,
                                         students_by_uid_lists=student_uid_lists,
                                         teachers_by_uid_lists=teacher_uid_lists,
                                         non_teaching_staff_uid_lists=non_teaching_staff_uid_lists,
                                         workplace_by_uid_lists=workplace_uid_lists,
                                         use_two_group_reduction=use_two_group_reduction,
                                         with_school_types=with_school_types,
                                         school_mixing_type=school_mixing_type,
                                         average_class_size=average_class_size,
                                         inter_grade_mixing=inter_grade_mixing,
                                         average_student_teacher_ratio=average_student_teacher_ratio,
                                         average_teacher_teacher_degree=average_teacher_teacher_degree,
                                         average_student_all_staff_ratio=average_student_all_staff_ratio,
                                         average_additional_staff_degree=average_additional_staff_degree,
                                         school_type_by_age=school_type_by_age,
                                         max_contacts=max_contacts)

        # Change types
        for key, person in population.items():
            for layerkey in population[key]['contacts'].keys():
                population[key]['contacts'][layerkey] = list(population[key]['contacts'][layerkey])

        school_mixing_types = [self.schools_in_groups[ns]['school_mixing_type'] for ns in range(len(self.schools_in_groups))]

        # temporarily store some information
        self.homes_by_uids = homes_by_uids
        self.workplace_uid_lists = workplace_uid_lists
        self.student_uid_lists = student_uid_lists
        self.teacher_uid_lists = teacher_uid_lists
        self.non_teaching_staff_uid_lists = non_teaching_staff_uid_lists
        self.school_types = school_types
        self.school_mixing_types = school_mixing_types

        self.set_layer_classes()
        self.clean_up_layer_info()

        return population

    def set_layer_classes(self):
        """Add layer classes."""
        self.initialize_households_list()
        self.populate_households(self.homes_by_uids, self.age_by_uid)
        self.initialize_workplaces_list()
        self.populate_workplaces(self.workplace_uid_lists)

        self.initialize_schools_list()
        self.populate_schools(self.student_uid_lists, self.teacher_uid_lists,
                              self.non_teaching_staff_uid_lists, self.age_by_uid,
                              self.school_types, self.school_mixing_types)

        self.populate_all_classrooms(self.schools_in_groups)
        return

    def clean_up_layer_info(self):
        """
        Clean up temporary data from the pop object after storing them in specific layer classes.
        """
        for key in ['workplace_uid_lists', 'student_uid_lists', 'teacher_uid_lists',
                    'non_teaching_staff_uid_lists', 'school_types',
                    'school_mixing_types', 'schools_in_groups',
                    ]:
            self.pop_item(key)
        return

    def pop_item(self, key):
        """Pop key from self."""
        self.__dict__.pop(key, None)  # pop checks if the key exists as an attribute and removes it in that case. Returns a default value of None if the key does not exist

    def to_dict(self):
        """
        Export to a dictionary -- official way to get the popdict.

        **Example**::

            popdict = pop.to_dict()
        """
        return sc.dcp(self.popdict)

    def to_json(self, filename, indent=2, **kwargs):
        """
        Export to a JSON file.

        **Example**::

            pop.to_json('my-pop.json')
        """
        return sc.savejson(filename, self.popdict, indent=indent, **kwargs)

    def save(self, filename, **kwargs):
        """
        Save population to an binary, gzipped object file.

        **Example**::

            pop.save('my-pop.pop')
        """
        return sc.saveobj(filename, self, **kwargs)

    @staticmethod
    def load(filename, *args, **kwargs):
        """
        Load from disk from a gzipped pickle.

        Args:
            filename (str): the name or path of the file to load from
            kwargs: passed to sc.loadobj()

        **Example**::

            pop = sp.Pop.load('my-pop.pop')
        """
        pop = sc.loadobj(filename, *args, **kwargs)
        if not isinstance(pop, Pop):
            errormsg = f'Cannot load object of {type(pop)} as a Pop object'
            raise TypeError(errormsg)
        return pop

    def initialize_households_list(self):
        """Initialize a new households list."""
        self.households = []
        return

    def initialize_empty_households(self, n_households=None):
        """
        Create a list of empty households.

        Args:
            n_households (int) : the number of households to initialize
        """
        znhh.initialize_empty_households(self, n_households)
        return

    def populate_households(self, households, age_by_uid):
        """
        Populate all of the households. Store each household at the index corresponding to it's hhid.

        Args:
            households (list) : list of lists where each sublist represents a household and contains the ids of the household members
            age_by_uid (dict) : dictionary mapping each person's id to their age
        """
        znhh.populate_households(self, households, age_by_uid)
        return

    def get_household(self, hhid):
        """
        Return household with id: hhid.

        Args:
            hhid (int) : household id number

        Returns:
            sp.Household: A populated household.
        """
        return znhh.get_household(self, hhid)

    def add_household(self, household):
        """
        Add a household to the list of households.

        Args:
            household (sp.Household): household with at minimum the hhid, member_uids, member_ages, reference_uid, and reference_age.
        """
        znhh.add_household(self, household)
        return

    def initialize_workplaces_list(self):
        """Initialize a new workplaces list."""
        self.workplaces = []
        return

    def initialize_empty_workplaces(self, n_workplaces=None):
        """
        Create a list of empty workplaces.

        Args:
            n_households (int) : the number of workplaces to initialize
        """
        znhh.initialize_empty_workplaces(self, n_workplaces)
        return

    def populate_workplaces(self, workplaces):
        """
        Populate all of the workplaces. Store each workplace at the index corresponding to it's wpid.

        Args:
            workplaces (list) : list of lists where each sublist represents a workplace and contains the ids of the workplace members
            age_by_uid (dict) : dictionary mapping each person's id to their age
        """
        znwp.populate_workplaces(self, workplaces)
        return

    def get_workplace(self, wpid):
        """
        Return workplace with id: wpid.

        Args:
            wpid (int) : workplace id number

        Returns:
            sp.Workplace: A populated workplace.
        """
        return znwp.get_workplace(self, wpid)

    def add_workplace(self, workplace):
        """
        Add a workplace to the list of workplaces.

        Args:
            workplace (sp.Workplace): workplace with at minimum the wpid, member_uids, member_ages, reference_uid, and reference_age.
        """
        znwp.add_workplace(self, workplace)
        return

    def initialize_schools_list(self):
        """Initialize a new schools list."""
        self.schools = []
        return

    def initialize_empty_schools(self, n_schools=None):
        """
        Create a list of empty schools.

        Args:
            n_schools (int) : the number of schools to initialize
        """
        znsch.initialize_empty_schools(self, n_schools)
        return

    def populate_schools(self, student_lists, teacher_lists, non_teaching_staff_lists, age_by_uid, school_types=None, school_mixing_types=None):
        """
        Populate all of the schools. Store each school at the index corresponding to it's scid.

        Args:
            student_lists (list)            : list of lists where each sublist represents a school and contains the ids of the students
            teacher_lists (list)            : list of lists where each sublist represents a school and contains the ids of the teachers
            non_teaching_staff_lists (list) : list of lists where each sublist represents a school and contains the ids of the non teaching staff
            age_by_uid (dict)               : dictionary mapping each person's id to their age
            school_types (list)             : list of the school types
            school_mixing_types (list)      : list of the school mixing types
        """
        znsch.populate_schools(self, student_lists, teacher_lists, non_teaching_staff_lists, age_by_uid, school_types, school_mixing_types)
        return

    def get_school(self, scid):
        """
        Return school with id: scid.

        Args:
            scid (int) : school id number

        Returns:
            sp.School: A populated school.
        """
        return znsch.get_school(self, scid)

    def add_school(self, school):
        """
        Add a school to the list of schools.

        Args:
            school (sp.School): school
        """
        znsch.add_school(self, school)
        return

    def populate_all_classrooms(self, schools_in_groups):
        """
        Populate all of the classrooms in schools for each school that has
        school_mixing_type equal to 'age_and_class_clustered'. Each classroom
        will be indexed at id clid.

        Args:
            schools_in_groups (dict) : a dictionary representing each school in terms of student_groups and teacher_groups corresponding to classrooms
        """
        for ns in range(self.n_schools):
            znsch.initialize_empty_classrooms(self.schools[ns], len(schools_in_groups[ns]['student_groups']))
            znsch.populate_classrooms(self.schools[ns], schools_in_groups[ns]['student_groups'], schools_in_groups[ns]['teacher_groups'], self.age_by_uid)
        return

    def get_classroom(self, scid, clid):
        """
        Return classroom with id: clid from school with id: scid.

        Args:
            scid (int) : school id number
            clid (int) : classroom id number

        Returns:
            sp.Classroom : A populated classroom.
        """
        return znsch.get_classroom(self, scid, clid)

    def compute_information(self):
        """Computing an advanced description of the population."""
        self.information = sc.objdict()
        self.information.age_count = self.count_pop_ages()
        self.information.layer_degrees = dict()
        self.information.layer_stats = dict()
        self.information.layer_degree_description = dict()

        for layer in self.layers:
            self.information.layer_degrees[layer] = spcnx.count_layer_degree(self, layer=layer)
            self.information.layer_stats[layer] = self.information.layer_degrees[layer].describe()[['age', 'degree']]
            self.information.layer_degree_description[layer] = self.information.layer_degrees[layer].groupby('age')['degree'].describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])  # default percentiles to include

        self.information.household_sizes = self.get_household_sizes()
        self.information.household_size_count = self.count_household_sizes()

        self.information.household_heads = self.get_household_heads()
        self.information.household_head_ages = self.get_household_head_ages()
        self.information.household_head_age_count = self.count_household_head_ages()
        self.information.household_head_ages_by_size_count = self.get_household_head_ages_by_size()

        self.information.ltcf_sizes = self.get_ltcf_sizes()
        self.information.ltcf_size_count = self.count_ltcf_sizes()

        self.information.enrollment_by_age = self.count_enrollment_by_age()
        self.information.enrollment_by_school_type = self.count_enrollment_by_school_type()

        self.information.employment_by_age = self.count_employment_by_age()
        self.information.workplace_sizes = self.get_workplace_sizes()
        self.information.workplace_size_count = self.count_workplace_sizes()

        return

    def compute_summary(self):
        """Compute summaries and add to pop post generation."""
        self.summary = sc.objdict()
        self.summary.mean_age = spb.calculate_mean_from_count(self.information.age_count)
        self.summary.std_age = spb.calculate_std_from_count(self.information.age_count)

        self.summary.layers = dict()
        for layer in self.layers:
            self.summary.layers[layer] = dict()

        percentiles = [5, 95]

        self.summary.layers['H']['mean'] = np.mean(list(self.information.household_sizes.values()))
        self.summary.layers['H']['std'] = np.std(list(self.information.household_sizes.values()))
        for p in percentiles:
            self.summary.layers['H'][p] = np.percentile(list(self.information.household_sizes.values()), q=p)

        sizes = []
        for s in self.information.enrollment_by_school_type.keys():
            sizes.extend(self.information.enrollment_by_school_type[s])
        self.summary.layers['S']['mean'] = np.mean(sizes)
        self.summary.layers['S']['std'] = np.std(sizes)
        for p in percentiles:
            self.summary.layers['S'][p] = np.percentile(sizes, q=p)

        self.summary.layers['W']['mean'] = np.mean(list(self.information.workplace_sizes.values()))
        self.summary.layers['W']['std'] = np.std(list(self.information.workplace_sizes.values()))
        for p in percentiles:
            self.summary.layers['W'][p] = np.percentile(list(self.information.workplace_sizes.values()), q=p)

    def summarize(self, return_msg=False):
        """Print and optionally return a brief summary string of the pop."""
        msg = ""
        msg += f"This networked population is created to resemble the population of {self.location + ',' if self.location is not None else ''} {self.state_location + ',' if self.state_location is not None else ''} {self.country_location if self.country_location is not None else ''}.\n"
        msg += f"The number of people is {self.n:.0f}.\n"
        msg += f"The mean age is {self.summary.mean_age:.2f} +/- {self.summary.std_age:.2f} years old.\n"
        msg += "\n"

        for layer in self.layers:
            s = self.information.layer_stats[layer]

            msg += f"Layer {layer}: {self.layer_mappings[layer]}\n"
            msg += f"   Number of people: {len(self.information.layer_degrees[layer]):.0f}\n"
            msg += f"   Number of edges: {self.n * s.loc[s.index == 'mean']['degree'][0] * 2:.0f} ({s.loc[s.index == 'mean']['degree'][0]:.1f} ± {s.loc[s.index == 'std']['degree'][0]:.1f} per person)\n"
            msg += f"   Age (years): {s.loc[s.index == 'mean']['age'][0]:.1f} ({s.loc[s.index == 'min']['age'][0]:.0f}-{s.loc[s.index == 'max']['age'][0]:.0f})\n"

            if layer in ['H', 'S', 'W']:
                msg += f"   {self.layer_mappings[layer].title()} size: {self.summary.layers[layer]['mean']:.1f} ± {self.summary.layers[layer]['std']:.1f} people (range is {self.summary.layers[layer][5]:.1f}-{self.summary.layers[layer][95]:.1f}).\n"

            msg += "\n"

        msg += f"The rand_seed used to generate this population is {self.rand_seed}."

        print(msg)
        if return_msg:
            return msg
        else:
            return

    def count_pop_ages(self):
        """
        Create an age count of the generated population post generation.

        Returns:
            dict: Dictionary of the age count of the generated population.
        """
        return spb.count_ages(self.popdict)

    # convert to work on array
    def get_household_sizes(self):
        """
        Create household sizes in the generated population post generation.

        Returns:
            dict: Dictionary of household size by household id (hhid).
        """
        return znhh.get_household_sizes(self.popdict)

    # convert to work on array
    def count_household_sizes(self):
        """
        Count of household sizes in the generated population.

        Returns:
            dict: Dictionary of the count of household sizes.
        """
        return spb.count_values(self.information.household_sizes)

    # convert to work on array
    def get_household_heads(self):
        """Get the ids of the head of households in the generated population post generation."""
        return znhh.get_household_heads(self.popdict)

    def get_household_head_ages(self):
        """Get the age of the head of each household in the generated population post generation."""
        return {hhid: self.popdict[head_id]['age'] for hhid, head_id in self.information.household_heads.items()}

    def count_household_head_ages(self, bins=None):
        """
        Count of household head ages in the generated population.

        Args:
            bins (array) : If supplied, use this to create a binned count of the household head ages. Otherwise, count discrete household head ages.

        Returns:
            dict: Dictionary of the count of household head ages.
        """
        if bins is None:
            return spb.count_values(self.information.household_head_ages)
        else:
            head_ages = list(self.information.household_head_ages.values())
            hist, bins = np.histogram(head_ages, bins=bins, density=0)
            return {i: hist[i] for i in range(len(hist))}

    def get_household_head_ages_by_size(self):
        """
        Get the count of households by size and the age of the head of the
        household, assuming the minimal household members id is the id of the
        head of the household.

        Returns:
            np.ndarray: An array with row as household size and columns as
            household head age brackets.
        """
        return znhh.get_household_head_ages_by_size(self)

    def count_enrollment_by_age(self):
        """
        Create enrollment count by age for students in the generated population post generation.

        Returns:
            dict: Dictionary of the count of enrolled students by age in the generated population.
        """
        return znsch.count_enrollment_by_age(self.popdict)

    @property
    def enrollment_rates_by_age(self):
        """
        Enrollment rates by age for students in the generated population.

        Returns:
            dict: Dictionary of the enrollment rates by age for students in the generated population.
        """
        return {k: self.information.enrollment_by_age[k]/self.information.age_count[k] if self.information.age_count[k] > 0 else 0 for k in range(defaults.settings.max_age)}

    def count_enrollment_by_school_type(self, *args, **kwargs):
        """
        Create enrollment sizes by school types in the generated population post generation.

        Returns:
            list: List of generated enrollment sizes by school type.
        """
        enrollment_by_school_type = znsch.count_enrollment_by_school_type(self.popdict, *args, **kwargs)
        return enrollment_by_school_type

    def count_employment_by_age(self):
        """
        Create employment count by age for workers in the generated population post generation.

        Returns:
            dict: Dictionary of the count of employed workers by age in the generated population.
        """
        return znwp.count_employment_by_age(self.popdict)

    @property
    def employment_rates_by_age(self):
        """
        Employment rates by age for workers in the generated population.

        Returns:
            dict: Dictionary of the employment rates by age for workers in the generated population.
        """
        return {k: self.information.employment_by_age[k]/self.information.age_count[k] if self.information.age_count[k] > 0 else 0 for k in range(defaults.settings.max_age)}

    # convert to work on array
    def get_workplace_sizes(self):
        """
        Create workplace sizes in the generated population post generation.

        Returns:
            dict: Dictionary of workplace size by workplace id (wpid).
        """
        return znwp.get_workplace_sizes(self.popdict)

    # convert to work on array
    def count_workplace_sizes(self):
        """
        Count of workplace sizes in the generated population.

        Returns:
            dict:Dictionary of the count of workplace sizes.
        """
        return spb.count_values(self.information.workplace_sizes)

    def get_contact_counts_by_layer(self, layer='S', **kwargs):
        """
        Get the number of contacts by layer.

        Returns:
            dict: Dictionary of the count of contacts in the layer for the
            different people types in the layer. See
            sp.contact_networks.get_contact_counts_by_layer() for method details.
        """
        return spcnx.get_contact_counts_by_layer(self.popdict, layer, **kwargs)

    def to_people(self):
        ''' Convert to the alternative People representation of a population '''
        ppl = spp.make_people(popdict=self.popdict, rand_seed=self.rand_seed)  # Create the corresponding population
        return ppl
