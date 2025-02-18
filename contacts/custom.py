import numpy as np
import sciris as sc
import random

def group_uids_by_income_age_brackets(age_bracs, inc_bracs, age_by_uid, fam_income_by_uid):
    # Group UIDs by income and age brackets. 
    ppl_income_age = {}
    for i in range(len(inc_bracs)):
        ppl_income_age[i] = {i:[] for i in range(len(age_bracs))} # Each cell contains a list of people we'll append to. 

    for uid in age_by_uid:
        person = uid

        # Get brackets. 
        income = fam_income_by_uid[uid]
        age = age_by_uid[uid]

        # See which bracket the income calls into. 
        b_i = None
        for i, b_inc in enumerate(inc_bracs):
            if income < b_inc:
                b_i = i
                break
        # If still not assigned, assign to last bracket. 
        if b_i is None:
            b_i = len(inc_bracs) - 1
        
        # same for age.
        b_a = None
        for i, b_age in enumerate(age_bracs):
            if age < b_age:
                b_a = i
                break
        if b_a is None:
            b_a = len(age_bracs) - 1

        ppl_income_age[b_i][b_a].append(person)
    
    # Now, shuffle each list. 
    for b_i in ppl_income_age:
        for b_a in ppl_income_age[b_i]:
            random.shuffle(ppl_income_age[b_i][b_a])

    return ppl_income_age


def allocate_household_incomes(pars, age_by_uid, homes_by_uids, homes):
    ret = dict()
    A_I_joint = pars.age_income_dist
    p_age = np.sum(A_I_joint, axis=1)
    all_incs = pars.ai_inc_bracs    

    for i, home_ages in enumerate(homes): # TODO: I think homes contains the age of each person per home.
        # TODO: I think that home is the key. 

        # Assume the age of the reference individual is simply the highest age in the household.
        hh_holder_age = max(home_ages)

        # Calculate P(I|A); a vector of probabilities of each income, given the age.
        cur_joint = A_I_joint[hh_holder_age,:]
        p_inc_given_age = cur_joint / p_age[hh_holder_age]

        # Sample an income from this distribution.
        hh_income = np.random.choice(all_incs, p=p_inc_given_age)

        for uid in homes_by_uids[i]:
            ret[uid] = hh_income
    
    return ret 

def update_pars(default_pars, pars):
    for key, item in pars.items():
        default_pars[key] = item
    return default_pars

def make_pars():
    # Make default pars for the behaviour model. 
    ret = sc.objdict(
            n=None,
            max_contacts=None,
            ltcf_pars=None,
            school_pars=None,
            with_industry_code=False,
            with_facilities=False,
            use_default=False,
            with_school_types=False,
            school_mixing_type='random',
            average_class_size=20,
            inter_grade_mixing=0.1,
            average_student_teacher_ratio=20,
            average_teacher_teacher_degree=3,
            teacher_age_min=25,
            teacher_age_max=75,
            with_non_teaching_staff=False,
            average_student_all_staff_ratio=15,
            average_additional_staff_degree=20,
            staff_age_min=20,
            staff_age_max=75,
            com_contacts=20,
            com_dispersion=None,
            rand_seed=None,
            country_location=None,
            state_location=None,
            location=None,
            sheet_name=None,
            household_method='infer_ages',
            smooth_ages=False,
            window_length=7,
            do_make=True,
            as_region = False,
            base_uid=0,)
    return ret