# Questions

## About poultry farms

- What roles do people play on the farm?
  - Owner
  - Worker
  - Vet?
  - Inspector?
- How many people work on a typical poultry farm?
- How closely does each role work with the flocks?
  - How much time on a typical day would each role spend in close proximity with the flock?
- What type of poultry flocks are there?
  - broiler
  - laying
  - incubator
  - What are the differences in how these flocks are handled?
    - I assume layers have their eggs collected regularly.
    - Removal of fecal matter?
    - Feeding?
    - Inspections?
    - Differences in housing conditions?
      - Are incubators more/less isolated from the environment than broilers?
      - Backyard flocks vs Commercial?
      - Free-run, Free-range?
      - To quote Beaudoin "Researchers commonly distinguish four broad types of production systems for poultry
production ... . Two are low-biosecurity systems: backyard or village pro-
duction, and minimally biosecure open houses with or without netting. In the former,
poultry roam freely in domestic settings; in the latter, poultry may be allowed out-
side of the housing, as with free-grazing duck flocks. Two other systems have greater
biosecurity: contract farms (in which poultry are in closed facilities with basic physical
barriers) and biosecure closed systems (with all aspects of production and transport
conducted internally)"
- Do farms have multiple types of flocks?
  - If so, do they come into contact with each other?
- How many incubators does the typical farm work with?
  - Is each poultry flock replaced entirely? Or is there mixing of old and new stock?
  - Do farms consistently work with the same incubators?
- How frequently are flocks sent to the abattoir?
  - Is the whole flock sent at once?
  - Is the same abattoir always used?
- What happens if a farmer suspects their flock is infected?
- How frequently do outside parties visit the farm?
  - Inspectors?
  - Vets?
  - Buyers/Sellers?
  - Nosey scientists trying to learn about poultry farming?
- What equipment is used and is it shared between farms?

## About H5N1

- For each species (currently human and poultry flocks),
and pathogen (currently H5N1), we need estimates of the following:
  - Time from exposure to infectious
  - Time from infectious to symptomatic
  - Probability of being symptomatic
- For poultry:
  - How does the disease progress in a poultry flock?
    - i.e. what is the best way to model the progression considering:
      - Time from first exposure to infectious
        - Ward suggests this is ~4 days in humans
      - Duration of infectious period
      - Probability of an infection going undetected (time dependent?)
  - What are the first signs of infection a farmer might notice?
- For humans:
  - Time from symptomatic to critical
    - Lai suggests this is ~4 days
    - Hui suggests ~4 days from onset of symptoms to hospitalization
  - Probability of critical symptoms
    - Lai suggests 90% of cases will eventually be hospitalized
  - Time from asymptomatic to recovered
  - Time from mild to recovered
  - Time from Critical to recovered
    - Lai suggests ~5 days
  - Time from critical to death
    - Lai suggests ~5 days
- Are there any vaccines available for each species type?
  - If so, what is there efficacy?
- Other types of immunity (e.g. from previous infection)
- Immunity progression
  - What is the peak immunity achieved
  - When is peak immunity achieved
  - What is the final level of immunity achieved?
  - When is the final level of immunity achieved?
- Cross-species infection
  - probability of poultry infecting human
  - probability of human infecting poultry

## About Quebec

- How many poultry farms are there in Quebec?
  - Martine suggested 700+ Broiler flocks, 200 laying flocks, and 99 incubator flocks.
- How many abattoirs are in Quebec?

## General questions

- What do we include in the model and what can be excluded?
- Whats the best way to model immunity (e.g. should we use the nab framework?)
- Relative strength of contact types? A key aspect of the model is that contacts between agents come in a variety of types each of which modifies the probability of transmission. For example in pathosim contacts are broken down into 'work', 'home', 'school', and 'community'. Home contacts might have a weight of 2.0 where work contacts have a weight of 1.0, this would indicate that an agent is twice as likely to be infected at home as at work, all else being equal. This reflects the inherent differences in the types of contact that happen in these locations.
  - In our case the contact types will be different but we are still interested in their relative importance. Our contact types could be:
    - different contact type for different activities related to poultry (e.g. 'feeding', 'collecting eggs', 'removing fecal matter'...). We will want to have a rough idea of how likely transmission is during each of these activities relative to each other.

## Sources

- Global epidemiology of avian influenza A H5N1 virus infection in humans, 1997–2015: a systematic review of individual case data
    Lai, Shengjie et al. The Lancet Infectious Diseases, Volume 16, Issue 7, e108 - e118
- Estimates of epidemiological parameters for H5N1 influenza in humans: a rapid review
    Jack Ward, Joshua W. Lambert, Timothy W. Russell, James M. Azam, Adam J. Kucharski, Sebastian Funk, Billy J. Quilty, Oswaldo Gressani, Niel Hens, W. John Edmunds
    medRxiv 2024.12.11.24318702; doi: [link](https://doi.org/10.1101/2024.12.11.24318702)
- Hui DS. Review of clinical symptoms and spectrum in humans with influenza A/H5N1 infection. Respirology. 2008 Mar;13 Suppl 1:S10-3. doi: 10.1111/j.1440-1843.2008.01247.x. PMID: 18366521.
- Beaudoin, Amanda, and Alan G. Isaac. "Direct and indirect transmission of avian influenza: results from a calibrated agent-based model." Journal of Economic Interaction and Coordination 18.1 (2023): 191-212.
