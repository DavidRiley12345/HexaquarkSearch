1 - Select your variables you would like to plot weighted in tree_filter.py
      also apply initial cuts at this stage
      The resulting tree is output to cut.root

2 - Populate run gensweight.C with the selected variables and then run it to create a a new tree with the variables plus the sweights in in weighted.root

3 - Then use /allvars/ to plot variables with signal and bkg
    due to some weirdness with the trees, only include 'sw_allvar_makesig_M();' the first time you run the program on a freshlyt created tree, as it wont let you overwrite it, comment it out otherwise

