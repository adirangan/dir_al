# dir_al
directory of c functions related to al simulation

%%%%%%%%;
% Compilation. ;
%%%%%%%%;
make al ;
%%%%%%%%;
% Note that there are several older libraries referenced in the makefile. ;
% For example, the current makefile references libglut, libGL and libGLU. ;
%%%%%%%%;

%%%%%%%%;
% Input-Parameters. ;
%%%%%%%%;
I believe that the parameters within:
al_wilson6b_reliability.in
recapitulate those within the paper:
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002622

For example, you can compare the global-variables:
CS__
SPARSE__
SPARSE_OGLI__
to Eqs. 9 and 10 in the paper.

The input-parameters within
al_wilson6b_reliability.in
are set up so that:
(i) there is no visualization, and
(ii) the system is stimulated with multiple periodic input-pulses, and
(iii) multiple forms of data-collection are enabled (see the 'suite' data-structure).

The input-parameters within:
al_wilson6b_view_0.in 
are quite similar, with the differences:
(i) the visualization is enabled, and
(ii) the stimulus is far simpler (with pulses controlled by the viewer), and
(iii) most forms of data-collection are disabled (with the exception of the 'power' data-structure). 

%%%%%%%%;
% Importing connectivity: ;
%%%%%%%%;
To incorprate a particular connectivity, you can set up an ascii file with the following format:
%%%%;
Nra->ntypes;
Nra->lengthra[0],Nra->lengthra[1],...,Nra->lengthra[Nra->ntypes-1];
LINK= ntA,nrA,NUMBER_OF_CONNECTED_NEURONS,nt1,nr1,nt2,nr2,nt3,nr3,...,ntN,nrN;
LINK= ntB,nrB,NUMBER_OF_CONNECTED_NEURONS,nt1,nr1,nt2,nr2,nt3,nr3,...,ntN,nrN;
... 
END;
%%%%;

