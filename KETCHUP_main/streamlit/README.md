#KETCHUP GUI via Streamlit
Running KETCHUP as a browser interface instructions

#Quickstart
KETCHUP GUI is accessible via streamlit where a user-interface will pop-up via browser upon command. To start the GUI, open command prompt in the streamlit folder and type: 
```console
streamlit run main.pt
```
The GUI will take in user-defined YAML files and populate the file input boxes. The minimum information required to run KETCHUP is the model, mechanism, and data file. However, user can toggle the "Additional options" checkbox to select more available options (for more details see
[main KETCHUP Quickstart guide](../../doc/intallation.rst). 
Any changes made after an preset YAML file is uploaded will be active (i.e., the final settings/files selected upon pressing "Start parameterization" will be used).

