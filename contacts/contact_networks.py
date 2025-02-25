"""
This module generates contact networks.
"""

import sciris as sc
import numpy as np
import pandas as pd
import networkx as nx

from .. import utils as znu

def make_all_contacts(pop, structs, pars):
    """
    From microstructure objects (dictionary mapping ID to age, lists of lists in different settings, etc.), create a dictionary of individuals.
    Each key is the ID of an individual which maps to a dictionary for that individual with attributes such as their age, household ID (hhid),
    school ID (scid), workplace ID (wpid), workplace industry code (wpindcode) if available, and contacts in different layers.

    Args:
    structs is a Sciris objdict, which contains the following structural information:


    pars is a Sciris objdict, which specifies the following parameters:



    Returns:
        A popdict of people with attributes.

    Notes:
        Methods to trim large groups of contacts down to better approximate a sense of close contacts (such as classroom sizes or
        smaller work groups are available via sp.trim_contacts() or sp.create_reduced_contacts_with_group_types(): see these methods for more details).

        If with_school_types==False, completely random schools will be generated with respect to the average_class_size,
        but other parameters such as average_additional_staff_degree will not be used.
    """
    pop.popdict = init_popdict_skele(structs, sexes=np.random.randint(2, size=len(structs.age_by_uid)))


def init_popdict_skele(structs):
    # Return population dictionary skeleton.
    popdict = {}
    layer_keys = [] 

    for uid in structs.age_by_uid:
        popdict[uid] = {}
        popdict[uid]['age'] = int(structs.age_by_uid[uid])

        popdict[uid]['contacts'] = {}

        for k in layer_keys:
            popdict[uid]['contacts'][k] = set()
    return popdict



def get_contact_counts_by_layer(popdict, layer='S', with_layer_ids=False):
    """
    Method to count the number of contacts for individuals in the population
    based on their role in a layer and the role of their contacts. For example,
    in schools this method can distinguish the number of contacts between
    students, teachers, and non teaching staff in the population, as well as
    return the number of contacts between all individuals present in a school.
    In a population with a school layer and roles defined as students, teachers,
    and non teaching staff, this method will return the number of contacts or
    edges for sc_students, sc_teachers, and sc_staff to sc_student, sc_teacher,
    sc_staff, all_staff, all. all_staff is the combination of sc_teacher and
    sc_staff, and all is all kinds of people in schools.

    Args:
        popdict (dict)        : popdict of a Pop object, Dictionary keys are the IDs of individuals in the population and the values are a dictionary
        layer (str)           : name of the physical contact layer: H for households, S for schools, W for workplaces, C for community, etc.
        with_layer_ids (bool) : If True, return additional dictionary on contacts by layer group id

    Returns:
        If with_layer_ids is False: A dictionary with keys = people_types
        (default to ['sc_student', 'sc_teacher', 'sc_staff']) and each value is
        a dictionary which stores the list of counts for each type of contact:
        default to ['sc_student', 'sc_teacher', 'sc_staff', 'all_staff', 'all']
        for example: contact_counter['sc_teacher']['sc_teacher'] store the
        counts of each teacher's contacts or edges to other teachers. If
        with_layer_ids is True: additionally return a dictionary with keys =
        layer_id (for example: scid, wpid...), and value is list of contact
        contacts.

    """
    layer = layer.upper()
    # layer keys are used to identify the people in that layer
    layer_keys = {"F": "fid",
                  "R": "rid",
                  "A": "aid"}

    # for all layers, 'all' contact_types will store counts for all contacts but
    # based on each different layer, there can be more contact_types.
    if layer == 'F':
        people_types = ['sc_student', 'sc_teacher', 'sc_staff']
        contact_types = people_types + ['all_staff', 'all']
    elif layer == "LTCF":
        people_types = ['ltcf_res', 'ltcf_staff']
        contact_types = people_types + ['all']
    elif layer in ["W", "H"]:
        people_types = [layer_keys[layer]]
        contact_types = ['all']
    else:
        raise NotImplementedError(f"layer {layer} not supported.")

    # initialize the contact counter between each people type and contact type as empty list
    contact_counter = {k: dict(zip(contact_types, ([] for _ in contact_types))) for k in
                       dict.fromkeys(people_types)}
    # index_switcher is a case-switch selector for the person selected by its type
    index_switcher = {k: contact_counter[k] for k in people_types}

    # also store all contacts count per layer id in contacts_counter_by_id
    contacts_counter_by_id = dict()
    for uid, person in popdict.items():
        if person[layer_keys[layer]] is not None:
            # count_switcher is a case-switch selector for contact counts by type
            count_switcher = {
                'sc_student': len([c for c in person["contacts"]["S"] if popdict[c]['sc_student']]),
                'sc_teacher': len([c for c in person["contacts"]["S"] if popdict[c]['sc_teacher']]),
                'sc_staff': len([c for c in person["contacts"]["S"] if popdict[c]['sc_staff']]),
                'ltcf_res': len([c for c in person["contacts"]["LTCF"] if popdict[c]['ltcf_res']]),
                'ltcf_staff': len([c for c in person["contacts"]["LTCF"] if popdict[c]['ltcf_staff']]),
                'all_staff': len([c for c in person["contacts"]["S"] if popdict[c]['sc_teacher']]) + len([c for c in person["contacts"]["S"] if popdict[c]['sc_staff']]),
                'all': len([c for c in person["contacts"][layer]])
            }

            contacts_counter_by_id.setdefault(person[layer_keys[layer]], [])
            for k1 in people_types:
                # if this person does not belong to a particular key, we don't need to store the counts under this key
                if person.get(k1) is not None:
                    # store sc_teacher, sc_student, sc_staff, all_staff and all below
                    if layer == "S":
                        for k2 in people_types:
                            index_switcher.get(k1)[k2].append(count_switcher.get(k2))
                        index_switcher.get(k1)["all_staff"].append(
                            count_switcher.get('sc_teacher') + count_switcher.get('sc_staff'))
                    # for other types, only all contacts are stored
                    index_switcher.get(k1)["all"].append(count_switcher.get('all'))
            if with_layer_ids:
                contacts_counter_by_id[person[layer_keys[layer]]].append(count_switcher.get('all'))
    if with_layer_ids:
        return contact_counter, contacts_counter_by_id
    else:
        return contact_counter


def filter_people(pop, ages=None, uids=None):
    """
    Helper function to filter people based on their uid and age.

    Args:
        pop (sp.Pop)         : population
        ages (list or array) : ages of people to include
        uids (list or array) : ids of people to include

    Returns:
        array: An array of the ids of people to include for further analysis.
    """
    output = np.arange(pop.n)
    if uids is not None:  # catch instance where the only uids supplied is the first one, 0
        output = np.intersect1d(output, uids)

    if ages is not None:  # catch instance where the only ages supplied is age 0
        output = np.intersect1d(output, sc.findinds(np.isin(pop.age_by_uid, ages)))
    return output


def count_layer_degree(pop, layer='H', ages=None, uids=None, uids_included=None):
    """
    Create a dataframe from the population of people in the layer, including
    their uid, age, degree, and the ages of contacts in the layer.

    Args:
        pop (sp.Pop)                 : population
        layer (str)                  : name of the physical contact layer: F for farms, R for shared resources, A for Abbatoirs.
        ages (list or array)         : ages of people to include
        uids (list or array)         : ids of people to include
        uids_included (list or None) : pre-calculated mask of people to include

    Returns:
        pandas.DataFrame: A pandas DataFrame of people in the layer including uid, age,
        degree, and the ages of contacts in the layer.
    """
    if uids_included is None:
        uids_included = filter_people(pop, ages=ages, uids=uids)

    layerid_mapping = {'F': 'fid', 'R': 'rid', 'A': 'aid'}

    degree_dicts = []

    for i in uids_included:
        a = pop.age_by_uid[i]

        if pop.popdict[i][layerid_mapping[layer]] is not None:
            nc = len(pop.popdict[i]['contacts'][layer])
            ca = [pop.age_by_uid[j] for j in pop.popdict[i]['contacts'][layer]]
            degree_dicts.append({'uid': i, 'age': a, 'degree': nc, 'contact_ages': ca})

    degree_df = pd.DataFrame(degree_dicts)

    return degree_df


def compute_layer_degree_description(pop, layer='H', ages=None, uids=None, uids_included=None, degree_df=None, percentiles=None):
    """
    Compute a description of the statistics for the degree distribution by age
    for a layer in the population contact network. See
    pandas.Dataframe.describe() for more details on all of the statistics
    included by default.

    Args:
        pop (sp.Pop)         : population
        layer (str)  : name of the physial contact layer: H for households, S for schools, W for workplaces, C for community or other
        ages (list or array) : ages of people to include
        uids (list or array) : ids of people to include
        uids_included (list or None): pre-calculated mask of people to include
        degree_df (dataframe) : pandas dataframe of people in the layer and their uid, age, degree, and ages of their contacts in the layer
        percentiles (list) : list of the percentiles to include as statistics

    Returns:
        pandas.DataFrame: A pandas DataFrame of the statistics for the layer
        degree distribution by age.
    """
    if degree_df is None:
        degree_df = count_layer_degree(pop, layer, ages, uids, uids_included)

    if percentiles is None:
        percentiles = [0.05, 0.25, 0.5, 0.75, 0.95]

    d = degree_df.groupby('age')['degree'].describe(percentiles=percentiles)
    return d


def random_graph_model(uids, average_degree, seed=None):
    """
    Generate edges for a group of individuals given their ids from an Erdos-Renyi
    random graph model given the expected average degree.

    Args:
        uids (list, np.ndarray) : a list or array of the ids of people in the graph
        average_degree (float)  : the average degree in the generated graph

    Returns:
        nx.Graph : Fast implementation of the Erdos-Renyi random graph model.
    """
    N = len(uids)
    if N == 0:
        raise ValueError(f"Expected uids to a non-empty list or array. Instead, the length of uids is {len(uids)}.")

    if average_degree >= N:
        #log.debug(f"Desired average degree is greater than or equal to the number of nodes. This method does not support multi-edges; returning a fully connected graph.")
        G = nx.complete_graph(N)

    else:
        p = average_degree / N
        G = nx.fast_gnp_random_graph(N, p, seed=seed)

    return G


def get_expected_density(average_degree, n_nodes):
    """
    Calculate the expected density of an undirected graph with no self-loops
    given graph properties. The expected density of an undirected graph with
    no self-loops is defined as the number of edges as a fraction of the
    number of maximal edges possible.

    Reference: Newman, M. E. J. (2010). Networks: An Introduction (pp 134-135).
    Oxford University Press.

    Args:
        average_degree (float) : average expected degree
        n_nodes (int) : number of nodes in the graph

    Returns:
        float: The expected graph density.
    """
    E = n_nodes * average_degree / 2
    Emax = n_nodes * (n_nodes - 1) / 2
    density = min(E / Emax, 1)  # capture when the average density is greater than the number of nodes - 1
    return density