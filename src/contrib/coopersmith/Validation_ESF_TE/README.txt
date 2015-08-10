Validation directions:

To start, open up $GENIE/config/EventGeneratorListAssembler.xml

Change the following options:

  <param_set name="Default"> 
     <param type="int" name="NGenerators">   13                                 </param>
     <param type="alg" name="Generator-0">   genie::EventGenerator/QEL-CC       </param>
     <param type="alg" name="Generator-1">   genie::EventGenerator/QEL-NC       </param>
     <param type="alg" name="Generator-2">   genie::EventGenerator/RES-CC       </param>
     <param type="alg" name="Generator-3">   genie::EventGenerator/RES-NC       </param>
     <param type="alg" name="Generator-4">   genie::EventGenerator/DIS-CC       </param>
     <param type="alg" name="Generator-5">   genie::EventGenerator/DIS-NC       </param>
     <param type="alg" name="Generator-6">   genie::EventGenerator/COH-CC       </param>
     <param type="alg" name="Generator-7">   genie::EventGenerator/COH-NC       </param>
     <param type="alg" name="Generator-8">   genie::EventGenerator/DIS-CC-CHARM </param>
     <param type="alg" name="Generator-9">   genie::EventGenerator/QEL-CC-CHARM </param>
     <param type="alg" name="Generator-10">  genie::EventGenerator/NUE-EL       </param>
     <param type="alg" name="Generator-11">  genie::EventGenerator/IMD          </param>
     <param type="alg" name="Generator-12">  genie::EventGenerator/IMD-ANH      </param>
 </param_set>

To:

  <param_set name="Default"> 
     <param type="int" name="NGenerators">   1                                  </param>
     <param type="alg" name="Generator-0">   genie::EventGenerator/QEL-CC       </param>
 </param_set>

This is not strictly necessary, but will speed up the spline generation step a LOT.
(Or if you know another way to make gspl2root do only CCQE events, then do that.)



Now:

Configure GENIE to use Bodek-Ritchie and no Transverse Enhancement (default)
Run the following:

./make_splines.sh DEF

after the file DEF_QE_splines.xml is created (shouldn't take more than a minute or so)
Run the following:

./generate_CCQE_events.py DEF

You will probably want to run the second of these on nohup.





Now configure GENIE to use Bodek-Ritchie and Transverse Enhancement by changing
the following line in $GENIE/config/UserPhysicsOptions.xml :

change
<param type="bool" name="UseElFFTransverseEnhancement">  false                              </param>
to
<param type="bool" name="UseElFFTransverseEnhancement">  true                               </param>


Now run the following:

./make_splines.sh TE

after the file TE_QE_splines.xml is created (shouldn't take more than a minute or so)
Run the following:

./generate_CCQE_events.py TE

You will probably want to run the second of these on nohup.



Now configure GENIE to use Effective Spectral Function and no Transverse Enhancement by changing
the following lines in $GENIE/config/UserPhysicsOptions.xml :

change
<param type="alg"  name="NuclearModel">                 genie::FGMBodekRitchie/Default     </param>
to
<param type="alg"  name="NuclearModel">                 genie::EffectiveSF/Default </param>

change
<param type="bool" name="UseElFFTransverseEnhancement">  true                              </param>
back to
<param type="bool" name="UseElFFTransverseEnhancement">  false                               </param>


Now run the following:

./make_splines.sh EFF

after the file EFF_QE_splines.xml is created (shouldn't take more than couple minutes or so)
Run the following:

./generate_CCQE_events.py EFF

You will probably want to run the second of these on nohup.


Finally, run


root -b -q Validation.cpp+'()'

Validation plots will appear in the plots/ folder
