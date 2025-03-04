# Read Me (meta)

Unlike Covasim and Pathosim, Zoonosim includes several different agent types. Some agent's (i.e. human agents) are modeled explicitly at the individual level.
Others (e.g. poultry flocks) are modeled as a single agent which represents a collection of many individual constituents (A flock is a single agent that contains many birds).
Additionally, it is sometimes useful to treat certain inanimate objects as agents so there contact history can be easily recorded. All of these agent types have different
properties, potential states, and values that we may wish to track. This module defines a 'meta' object for each agent type which dictates their internal states, properties, and tracked values.
