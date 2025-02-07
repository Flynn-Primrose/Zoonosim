# Base Module

Base classes for Zoonosim. These classes handle a lot of the boilerplate of the
derived classes (e.g. loading, saving, key lookups, etc.), so those classes
can be focused on the disease-specific functionality.

This class system is inherited from Covasim and is designed to separate the mundane tasks from disease specific functionality.
I'm unsure at the time of writing (2025-02-07) how much of this will have to change to accommodate the addition of multiple species. I'm hoping not much.
I have moved the sub-module into it's own subdirectory and split it into separate files for each class because I find that more intuitive.

The major changes include:

- BasePeople class is now BaseRoster, This class will be used as the base for agent rosters of all species.
