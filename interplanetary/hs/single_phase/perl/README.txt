Files seq*.pl and input.dat are in this folder for GTOC3 purposes.

They take data from input.dat and produce output in out/Ast1_Ast2_Ast3. 
The file output.dat stored in the last mentioned folder, contains in sequence the times and the final mass for each
one of the phases, plus the final value of the objective function.

- seqi.pl is the main automated script used for multiple ampl calls. It considers possible unfeasibilities.
- seqifast.pl does not perform a full optimization, but also considers unfeasibilities and gives (supposedly) a quite good approximation for the result of
	seqi.pl
- seq.pl is the original version, it ignores unfeasibilities, although they are noted in output.dat.
- seqguess.pl ignores them also, but can ALSO be called with a sequence of numbers, identifying the possible guess used.
  It can be called without numbers, but in that case it will look in the output.dat file of a previous run to reproduce the selected ones.
- refine.pl interpolates a previous solution. It modifies Solution.out
  file and main.mod to duplicate (2*N-1) the number of nodes. It changes the guess option in main.mod for guessrefine.inc. 
- guessrefine.inc allows to use the previous stored solution (without interpolation)
  in a call from the outside of ampl.
- smoothref.inc performs a mean average on the solution stored. Combined with guessrefine.inc
  may lead to better results reducing the number of peaks.   

SEQI.PL CALL:

After a call of SEQI.PL the next step will be to get the solution stored. For that, seqguess.pl can be called without arguments, but the input.dat file must be modified. 
For that purpose check the output.dat file. Check the different phases and the corresponding guesses. 
-If the guess selected has been executed just once, no modification has to be made. 
-If it has been executed three times, the period corresponding to that phase in the input.dat must be increased by 50 days. 
-If it has been executed just two times:
     If in the first execution you see POSSIBLE MISMATCHING WITH NEXT PHASE
            Decrease -100 in the initial date, -50 in the period
     Else
	    Decrease -50 in the period
		
            
SEQGUESS.PL CALL:
   seqguess.pl N1 N2 N3 N4
   N1, N2... are the numbers representing the guesses for the corresponding phases
   0 -"crazyguess.inc"
   1 -"guesslambert.inc"
   2 -"guesslambert2.inc"
   3 -"guesstangential.inc"
   4 -"guesstangential2.inc"
   5 -"guesstangential2m.inc"
   6.-"guesstangential3.inc"

   Example:
   
   seqguess.pl 5 5 4 1

   launches AMPL choosing the guesses (E-A1 t2m) (A1-A2 t2m) (A2-A3 t2) (A3-A4
   t3)
   
   It has been verified to worked properly with one of the sequences simulated
   by Dario, producing the same numbers!
   


