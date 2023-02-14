1 - Select your variables you would like to plot weighted in tree_filter.py
      also apply initial cuts at this stage
      The resulting tree is output to cut.root

2 - Populate run gensweight.C with the selected variables and then run it to create a a new tree with the variables plus the sweights in in weighted.root

3 - Then use /allvars/ to plot variables with signal and bkg
    due to some weirdness with the trees, only include 'sw_allvar_makesig_M();' the first time you run the program on a freshlyt created tree, as it wont let you overwrite it, comment it out otherwise
    
    also from that point forward use the tree created in the allvars folder, as this contains the sigma_c invariant mass
    
And then use sigc_bkg_fit to fit just the sweighted background of the for later use


  

After that you can use sigmac_sw to create a fit and sweights for the sigma_c invariant mass plot
  This generates its own set of sweights and is not at all based on the sweights from the lambdac plots made above


