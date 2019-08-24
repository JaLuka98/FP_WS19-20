errors.py is a very simple python script.
It will help you to calculate and format the formulas for gaussian propagation of uncertainty.
Simply type your equation in line 22 in python and run the script.
The output is already formatted in latex.

Note_1: There is currently no support for values without an uncertainty in your equation. Python considers every variable in your equation to have an uncertainty.
Note_2: You can customize the names of the variables in line 20. Simply name them differently in 'var(...)' and they will appear with that name in the output.

Credits go to https://toolbox.pep-dortmund.org/files/archive/2013/AutomatisierenVonFehlerrechnung.html
