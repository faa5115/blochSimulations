# Bloch Simulations
  A fundamentals project showing how to break down slice-select excitation in Magnetic Resonance Imaging
## slice selection
This file path relies on the following files in the blochSimulations directory:
  1. getSincEnvelope_fa.m:  creates a time-varying b1 field for a given gyromagnetic ratio, specified duration, number of timepoints within that duration, 
                            time-badwidth-product, and flip angle. the output for each timepoint is in tesla. 
  2. func_sliceSelection.m:  carryies out the slice selection on a given set of linear spatially sorted spins along the slice's direction.  It relies on
                            the following input:
                            B1e (from getSincEnvelope_fa.m)  the  (1 x number of RF timepoints)  array which holds the B1 strength for each timepoint. 
                            sliceDir:  a (1 x number of slice points) array detaliing the slice location (in mm) for each spin).  
                            RFBW:  the excitation bandwidth of the B1e in Hz. 
                            sliceThickness:  the thickness we aim to excite. 
                            RFDur (same as for the getSincEnvelope_fa function):  RF's duration.
                            M0:  (number of slice points x 3) array that detail the initizl magnetization  for each spin (each spatial point has Mx, My, and Mz). If
                            we were starting from thermal equilibrium, each slice point would have its magnetization be [0 0 Mequilibirium]' where Meq is the     equilibrium       value. 
                            gamma:  gyromagnetic ratio (%rad/s/T.)
                            sliceOffCenter (mm):  the off isocenter excitation.
                            RFPhase0:  the Phase of the RF pulse.  
                            Meq:  similar to M0, but detailing the equilibrium magnetization for each spin.  
                            T1/T2: (seconds) recovery/relaxation constants.
                            dFArr:  (1 x number of slice points ) array detailing the off resonant frequencies along the slice select direction in Hz. 
                           
This function first will add a linear phase ramp to the B1e array: 
  B1ea = B1e;
  for n = 1 : size(B1e,2)
      B1ea(1,n) = B1e(1,n) * (cos(RFPhase0 + gamma*Gslice*sliceOffCenter*(tRF(1,n))) + 1i*(sin(RFPhase0 + gamma*Gslice*sliceOffCenter*(tRF(1,n)))));
  end

  B1e = B1ea ;
  
  So B1ea will be the RF pulse having the set excitaiton phase and specifying which location along the slice select axis to excite.  
  Below is an at isocenter example and another example exciting 4mm away from isocenter:
  
  ![GitHub Logo](/images/0mmOffIsoCenterB1Plot.jpeg)
  ![alt text](<./sliceSelectionDemo/0mmOffIsoCenterB1Plot.jpeg>)
  ![alt text](https://github.com/faa5115/blochSimulations/sliceSelectionDemo/neg4mmOffIsoCenterB1Plot.jpeg)
