# Run-time instructions for example files

This directory contains examples of both static and dynamic usage of KETCHUP. There are python files, YAML option files to help specify a problem, and solver option files, as well as the data in the `data` dir used in these examples.

The following instructions lists commands to run from the command window to run KETCHUP example files.

## Basic examples

For the static main example, run
```console
python KETCHUP_example.py
```
which is set up for k-ecoli74 data analysis.

For the dynamic KETCHUP examples, run
```console
python KETCHUP_dynamic.py
```
which is set up for the FDH data analysis.

These use the default values with a ketchup_user_options function at the beginning that defines the problem, including where files are located.

## Command line options

KETCHUP has three command line options:
 * --seedvalue (-s)
   * seedvalue sets the random seed used in the random number generator. It can be an integer or 'time' to use the system clock to set the seed.
 * --solver-options (-so)
   * solver-options is defined based on what solver Pyomo is using.
 * --program-options (-po)
   * program-options is a YAML file that contains user-defined KETCHUP options that override the ketchup_user_options function in the python code.

All of these options can be combined. The --seedvalue option has the highest priority of setting the random seed value.

## Command line arguments examples

To include user-defined details, run
```console
python KETCHUP_example.py -po options_example.yml
```
If seed value and/or solver option files need to be overwritten (e.g., set new seedvalue to 1 and solver option file to ipopt_hsl.opt) from the YAML options file, run
```console
python KETCHUP_example.py -po options_example.yml -s 1 -so ipopt_hsl.opt
```

As noted earlier seedvalue can also be 'time' to use clock time: run
```console
python KETCHUP_example.py -s time
```

For a dynamic example using a program option file, run
```console
python KETCHUP_dynamic.py -po options_dynamic_example_FDH_flat_diff.yml
```
to analyze the FDH data, but this time using a flat TSV form of the data.
Although the provided YAML file otherwise contains all the same information as in the ketchup_user_options function, you could modify the YAML file to remove all but the item(s) that override parts of the options set in the ketchup_user_options function. 
In that case, not every item needs to be in the program option file, as it just changes things from their default values or those set by the ketchup_user_options function.

## Defining all options in YAML file

Alternatively, if you don't wish to define any user options within the python file directly, you can use the example KETCHUP file that only reads from the YAML program options file. Here, everything needs to be included therein or on the command line. The above data analyses can be performed by

```console
python KETCHUP_yaml_only.py -po options_example.yml -s time -so ipopt.opt
```
for the static k-ecoli74 example and 
```console
python KETCHUP_yaml_only.py -po options_dynamic_example_FDH.yml
```
for the dynamic FDH example, and
```console
python KETCHUP_yaml_only.py -po options_dynamic_example_FDH_flat.yml
```
using the flat TSV experimental data for the FDH example.

This python program can also be a good way to check that everything is correctly defined for your problem within the YAML file, since the only items otherwise set are the KETCHUP default values.
