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
To incorporate a particular connectivity, you can set up an ascii file with the following format:
%%%%;
Nra->ntypes;
Nra->lengthra[0],Nra->lengthra[1],...,Nra->lengthra[Nra->ntypes-1];
LINK= ntA,nrA,NUMBER_OF_CONNECTED_NEURONS,nt1,nr1,nt2,nr2,nt3,nr3,...,ntN,nrN;
LINK= ntB,nrB,NUMBER_OF_CONNECTED_NEURONS,nt1,nr1,nt2,nr2,nt3,nr3,...,ntN,nrN;
... 
END;
%%%%;

%%%%%%%%;
% Programming the time-course of the input (i.e., input pulses). ;
%%%%%%%%;
The function
system_monitor
found within
al_datahandling.c
checks for the flag
SUITE_BOTHER.
If this flag is set (e.g., SUITE_BOTHER= 7), then various input pulses are programmed.
Unfortunately, due to the large number of hyperparameters involved in structuring the input,
most of system_monitor is currently hardcoded.
As an example, see:
case 7
which runs multiple odors for multiple trials over various concentrations.
Some of these parameters can be easily changed by modifying the input file.
For example, the number of trials can be modified by changing:
SUITE_NINSTANCES
and so forth.
Other parameters are not yet pulled from the input file.
For example, the onset and length of each pulse are determined by
pulse_waitfor
pulse_howlong
which are currently hardcoded (i.e., 1024 and 512 respectively).

%%%%%%%%;
% Strobe data-structure. ;
%%%%%%%%;
The strobe data-structure (i.e., st) records system variables (of type double) across time.
The total time of the record is fixed by st->length and st->update_timestep.
There are two ways in which strobed variables are updated.
the function
strobeupdate_old
takes a snapshot of the variable every st->update_timestep.
This is reasonable for some variables (e.g., the time itself), but not others (e.g., varname_registry_spike_flag).
The function
strobeupdate
takes a local average of the variable (over the time-window associated with st->update_timestep).
This is used for most variables.
This record (stored in st->data) is accessed periodically,
and successive periods will be overwritten (starting from the beginning)
if the simulation time exceeds the total time of the record.
If requested (by st->cycle_bother) the recording periods will be accumulated and stored (st->cycledata).
If requested (by st->lpower_bother), then the (log) of the power-spectrum is accumulated (st->lpowerdata).

%%%%%%%%;
% Power data-structure. ;
%%%%%%%%;
The power data-structure makes use of the strobe data-structure.
The power data-structure stores strobe data-structures associated with multiple variables.
These are read from the global neuron array (linked to by p->Nra)
and stored in the array of array of strobes (linked to by p->strarara).

%%%%%%%%;
% Outputs for the Power and Strobe data-structures. ;
%%%%%%%%;
The function
powerdump
dumps each of the strobes in p->strarara to a file
(with name of the form powerstra_*_*_*).
The dump is determined by dump_type.
If dump_type==0, then a matlab readable (ascii) file is created.
If dump_type==1, then fig, jpg and pnm files are created.
The pnm files are typically by WritePNMfile_color found in d_llists.c,
while the fig and jpg files are created by ra2jpg found in d_llists.c.
The fig file is xfig readable (ascii), and requires fig2dev to convert into a jpg.



