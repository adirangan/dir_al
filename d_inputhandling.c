/* /\* Here are the input functions *\/ */
/* void readinput(); */
/* void updateglobals(char *); */
/* void dumpoutput(); */

void readinput()
{
  /* This reads the piped input file, or standard input.
     the variable names should not be longer than 64 characters */
  int verbose=0;
  char *vname= (char *) tmalloc(sizeof(char)*128),*equals= (char *) tmalloc(sizeof(char)*128);
  do{
    scanf("%[^=]",vname);scanf("%s",equals);updateglobals(vname);scanf("%c",equals);
    if (verbose){ printf("At this point variable name is (%s)\n",vname);}}
  while (strcmp(vname,"END")!=0);
  tfree(vname);tfree(equals);
}

void updateglobals(char *vname)
{
  /* given a variable name, scan the appropriate format into the appropriate global variable */
  int verbose=0;
  if (strcmp(vname,"ON_MY_COMPUTER")==0){ scanf("%d",&ON_MY_COMPUTER); if (verbose){ printf("%s read to be %d\n",vname,ON_MY_COMPUTER);}}
  if (strcmp(vname,"GLOBAL_STRING")==0){ scanf("%s",GLOBAL_STRING); if (verbose){ printf("%s read to be %s\n",vname,GLOBAL_STRING);}}
  if (strcmp(vname,"GLOBAL_RECORD_NUMBER")==0){ scanf("%d",&GLOBAL_RECORD_NUMBER); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_RECORD_NUMBER);}}
  if (strcmp(vname,"GLOBAL_SPIKEINPUT_RSEED")==0){ scanf("%d",&GLOBAL_SPIKEINPUT_RSEED); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_SPIKEINPUT_RSEED);}}
  if (strcmp(vname,"GLOBAL_CLEANUP")==0){ scanf("%d",&GLOBAL_CLEANUP); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_CLEANUP);}}
  if (strcmp(vname,"NSLICES")==0){ scanf("%d",&NSLICES); if (verbose){ printf("%s read to be %d\n",vname,NSLICES);}}
  if (strcmp(vname,"NPHASES")==0){ scanf("%d",&NPHASES); if (verbose){ printf("%s read to be %d\n",vname,NPHASES);}}
  if (strcmp(vname,"SUITE_BOTHER")==0){ scanf("%d",&SUITE_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,SUITE_BOTHER);}}
  if (strcmp(vname,"SUITE_NSECONDS")==0){ scanf("%d",&SUITE_NSECONDS); if (verbose){ printf("%s read to be %d\n",vname,SUITE_NSECONDS);}}
  if (strcmp(vname,"SUITE_DUMPEVERY")==0){ scanf("%d",&SUITE_DUMPEVERY); if (verbose){ printf("%s read to be %d\n",vname,SUITE_DUMPEVERY);}}
  if (strcmp(vname,"SUITE_BITBYBIT_REMOVE")==0){ scanf("%d",&SUITE_BITBYBIT_REMOVE); if (verbose){ printf("%s read to be %d\n",vname,SUITE_BITBYBIT_REMOVE);}}
  if (strcmp(vname,"SUITE_SINDEXMAX")==0){ scanf("%d",&SUITE_SINDEXMAX); if (verbose){ printf("%s read to be %d\n",vname,SUITE_SINDEXMAX);}}
  if (strcmp(vname,"SUITE_DINDEXMAX")==0){ scanf("%d",&SUITE_DINDEXMAX); if (verbose){ printf("%s read to be %d\n",vname,SUITE_DINDEXMAX);}}
  if (strcmp(vname,"SUITE_TINDEXMAX")==0){ scanf("%d",&SUITE_TINDEXMAX); if (verbose){ printf("%s read to be %d\n",vname,SUITE_TINDEXMAX);}}
  if (strcmp(vname,"LGN_BOTHER")==0){ scanf("%d",&LGN_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,LGN_BOTHER);}}
  if (strcmp(vname,"LGN_DUMP")==0){ scanf("%d",&LGN_DUMP); if (verbose){ printf("%s read to be %d\n",vname,LGN_DUMP);}}
  if (strcmp(vname,"LGN_TYPE_FLAG")==0){ scanf("%d",&LGN_TYPE_FLAG); if (verbose){ printf("%s read to be %d\n",vname,LGN_TYPE_FLAG);}}
  if (strcmp(vname,"LGN_DUMB")==0){ scanf("%d",&LGN_DUMB); if (verbose){ printf("%s read to be %d\n",vname,LGN_DUMB);}}
  if (strcmp(vname,"GLOBAL_ODOR")==0){ scanf("%d",&GLOBAL_ODOR); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_ODOR);}}
  if (strcmp(vname,"GLOBAL_ODOR_BASE")==0){ scanf("%d",&GLOBAL_ODOR_BASE); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_ODOR_BASE);}}
  if (strcmp(vname,"GLOBAL_ODOR_BACON")==0){ scanf("%d",&GLOBAL_ODOR_BACON); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_ODOR_BACON);}}
  if (strcmp(vname,"GLOBAL_FPS")==0){ scanf("%d",&GLOBAL_FPS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_FPS);}}
  if (strcmp(vname,"OBSCURED")==0){ scanf("%d",&OBSCURED); if (verbose){ printf("%s read to be %d\n",vname,OBSCURED);}}
  if (strcmp(vname,"GRATING_VS_LMI")==0){ scanf("%d",&GRATING_VS_LMI); if (verbose){ printf("%s read to be %d\n",vname,GRATING_VS_LMI);}}
  if (strcmp(vname,"CYCLE_LENGTH")==0){ scanf("%lf",&CYCLE_LENGTH); if (verbose){ printf("%s read to be %lf\n",vname,CYCLE_LENGTH);}}
  if (strcmp(vname,"GRATING_PULSE")==0){ scanf("%lf",&GRATING_PULSE); if (verbose){ printf("%s read to be %lf\n",vname,GRATING_PULSE);}}
  if (strcmp(vname,"GRATING_DRIFT")==0){ scanf("%lf",&GRATING_DRIFT); if (verbose){ printf("%s read to be %lf\n",vname,GRATING_DRIFT);}}
  if (strcmp(vname,"STIMULUS_ONSET_TIME")==0){ scanf("%lf",&STIMULUS_ONSET_TIME); if (verbose){ printf("%s read to be %lf\n",vname,STIMULUS_ONSET_TIME);}}
  if (strcmp(vname,"SQUARE_DURATION_TIME")==0){ scanf("%lf",&SQUARE_DURATION_TIME); if (verbose){ printf("%s read to be %lf\n",vname,SQUARE_DURATION_TIME);}}
  if (strcmp(vname,"SQUARE_DRAG_TIME")==0){ scanf("%lf",&SQUARE_DRAG_TIME); if (verbose){ printf("%s read to be %lf\n",vname,SQUARE_DRAG_TIME);}}
  if (strcmp(vname,"LINE_DELAY_TIME")==0){ scanf("%lf",&LINE_DELAY_TIME); if (verbose){ printf("%s read to be %lf\n",vname,LINE_DELAY_TIME);}}
  if (strcmp(vname,"LGN_BACKRATE")==0){ scanf("%lf",&LGN_BACKRATE); if (verbose){ printf("%s read to be %lf\n",vname,LGN_BACKRATE);}}
  if (strcmp(vname,"OTHERLAYER_BACKRATE")==0){ scanf("%lf",&OTHERLAYER_BACKRATE); if (verbose){ printf("%s read to be %lf\n",vname,OTHERLAYER_BACKRATE);}}
  if (strcmp(vname,"OTHERLAYER_INPUTRATE")==0){ scanf("%lf",&OTHERLAYER_INPUTRATE); if (verbose){ printf("%s read to be %lf\n",vname,OTHERLAYER_INPUTRATE);}}
  if (strcmp(vname,"INPUT_CONTRAST")==0){ scanf("%lf",&INPUT_CONTRAST); if (verbose){ printf("%s read to be %lf\n",vname,INPUT_CONTRAST);}}
  if (strcmp(vname,"INPUT_SPACEK")==0){ scanf("%lf",&INPUT_SPACEK); if (verbose){ printf("%s read to be %lf\n",vname,INPUT_SPACEK);}}
  if (strcmp(vname,"INPUT_SPACEANGLE")==0){ scanf("%lf",&INPUT_SPACEANGLE); if (verbose){ printf("%s read to be %lf\n",vname,INPUT_SPACEANGLE);}}
  if (strcmp(vname,"INPUT_SPACEANGLE_BACON")==0){ scanf("%lf",&INPUT_SPACEANGLE_BACON); if (verbose){ printf("%s read to be %lf\n",vname,INPUT_SPACEANGLE_BACON);}}
  if (strcmp(vname,"INPUT_SPACEPHASE")==0){ scanf("%lf",&INPUT_SPACEPHASE); if (verbose){ printf("%s read to be %lf\n",vname,INPUT_SPACEPHASE);}}
  if (strcmp(vname,"CORTEX_BOTHER")==0){ scanf("%d",&CORTEX_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,CORTEX_BOTHER);}}
  if (strcmp(vname,"AUTAPSES_OFF")==0){ scanf("%d",&AUTAPSES_OFF); if (verbose){ printf("%s read to be %d\n",vname,AUTAPSES_OFF);}}
  if (strcmp(vname,"LR_BOTHER")==0){ scanf("%d",&LR_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,LR_BOTHER);}}
  if (strcmp(vname,"EIFvsIF")==0){ scanf("%d",&EIFvsIF); if (verbose){ printf("%s read to be %d\n",vname,EIFvsIF);}}
  if (strcmp(vname,"ORDER")==0){ scanf("%d",&ORDER); if (verbose){ printf("%s read to be %d\n",vname,ORDER);}}
  if (strcmp(vname,"LGN_DETERMINISM")==0){ scanf("%d",&LGN_DETERMINISM); if (verbose){ printf("%s read to be %d\n",vname,LGN_DETERMINISM);}}
  if (strcmp(vname,"HANSHELLEY_FLAG")==0){ scanf("%d",&HANSHELLEY_FLAG); if (verbose){ printf("%s read to be %d\n",vname,HANSHELLEY_FLAG);}}
  if (strcmp(vname,"STARTING_HOMOGENIZATION")==0){ scanf("%lf",&STARTING_HOMOGENIZATION); if (verbose){ printf("%s read to be %lf\n",vname,STARTING_HOMOGENIZATION);}}
  if (strcmp(vname,"SPIKETOL")==0){ scanf("%d",&SPIKETOL); if (verbose){ printf("%s read to be %d\n",vname,SPIKETOL);}}
  if (strcmp(vname,"GLOBAL_SPACE_SMOOTHER")==0){ scanf("%d",&GLOBAL_SPACE_SMOOTHER); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_SPACE_SMOOTHER);}}
  if (strcmp(vname,"ARBOR_DIA")==0){ scanf("%d",&ARBOR_DIA); if (verbose){ printf("%s read to be %d\n",vname,ARBOR_DIA);}}
  if (strcmp(vname,"BIG_SYSTEM_FLAG")==0){ scanf("%d",&BIG_SYSTEM_FLAG); if (verbose){ printf("%s read to be %d\n",vname,BIG_SYSTEM_FLAG);}}
  if (strcmp(vname,"LR_ANGLE_DIA")==0){ scanf("%lf",&LR_ANGLE_DIA); if (verbose){ printf("%s read to be %lf\n",vname,LR_ANGLE_DIA);}}
  if (strcmp(vname,"LR_DIST_DIA_INPIES")==0){ scanf("%lf",&LR_DIST_DIA_INPIES); if (verbose){ printf("%s read to be %lf\n",vname,LR_DIST_DIA_INPIES);}}
  if (strcmp(vname,"LR_DIR_DIA")==0){ scanf("%lf",&LR_DIR_DIA); if (verbose){ printf("%s read to be %lf\n",vname,LR_DIR_DIA);}}
  if (strcmp(vname,"LR_RATE")==0){ scanf("%lf",&LR_RATE); if (verbose){ printf("%s read to be %lf\n",vname,LR_RATE);}}
  if (strcmp(vname,"LR_TO_AMPA")==0){ scanf("%lf",&LR_TO_AMPA); if (verbose){ printf("%s read to be %lf\n",vname,LR_TO_AMPA);}}
  if (strcmp(vname,"AXONAL_DELAY_PER_PIE")==0){ scanf("%lf",&AXONAL_DELAY_PER_PIE); if (verbose){ printf("%s read to be %lf\n",vname,AXONAL_DELAY_PER_PIE);}}
  if (strcmp(vname,"P_EC")==0){ scanf("%lf",&P_EC); if (verbose){ printf("%s read to be %lf\n",vname,P_EC);}}
  if (strcmp(vname,"P_ES")==0){ scanf("%lf",&P_ES); if (verbose){ printf("%s read to be %lf\n",vname,P_ES);}}
  if (strcmp(vname,"P_IC")==0){ scanf("%lf",&P_IC); if (verbose){ printf("%s read to be %lf\n",vname,P_IC);}}
  if (strcmp(vname,"P_IS")==0){ scanf("%lf",&P_IS); if (verbose){ printf("%s read to be %lf\n",vname,P_IS);}}
  if (strcmp(vname,"P_AMPA")==0){ scanf("%lf",&P_AMPA); if (verbose){ printf("%s read to be %lf\n",vname,P_AMPA);}}
  if (strcmp(vname,"P_NMDA")==0){ scanf("%lf",&P_NMDA); if (verbose){ printf("%s read to be %lf\n",vname,P_NMDA);}}
  if (strcmp(vname,"P_GABA")==0){ scanf("%lf",&P_GABA); if (verbose){ printf("%s read to be %lf\n",vname,P_GABA);}}
  if (strcmp(vname,"GABOR_DRIFT")==0){ scanf("%lf",&GABOR_DRIFT); if (verbose){ printf("%s read to be %lf\n",vname,GABOR_DRIFT);}}
  if (strcmp(vname,"LGNANGLE_DRIFT")==0){ scanf("%lf",&LGNANGLE_DRIFT); if (verbose){ printf("%s read to be %lf\n",vname,LGNANGLE_DRIFT);}}
  if (strcmp(vname,"LGN_STRENGTH")==0){ scanf("%lf",&LGN_STRENGTH); if (verbose){ printf("%s read to be %lf\n",vname,LGN_STRENGTH);}}
  if (strcmp(vname,"OTHERLAYER_STRENGTH")==0){ scanf("%lf",&OTHERLAYER_STRENGTH); if (verbose){ printf("%s read to be %lf\n",vname,OTHERLAYER_STRENGTH);}}
  if (strcmp(vname,"CS_ESES_A")==0){ scanf("%lf",&CS_ESES_A); if (verbose){ printf("%s read to be %lf\n",vname,CS_ESES_A);}}
  if (strcmp(vname,"CS_ESIS_A")==0){ scanf("%lf",&CS_ESIS_A); if (verbose){ printf("%s read to be %lf\n",vname,CS_ESIS_A);}}
  if (strcmp(vname,"CS_ESEC_A")==0){ scanf("%lf",&CS_ESEC_A); if (verbose){ printf("%s read to be %lf\n",vname,CS_ESEC_A);}}
  if (strcmp(vname,"CS_ESIC_A")==0){ scanf("%lf",&CS_ESIC_A); if (verbose){ printf("%s read to be %lf\n",vname,CS_ESIC_A);}}
  if (strcmp(vname,"CS_ECES_A")==0){ scanf("%lf",&CS_ECES_A); if (verbose){ printf("%s read to be %lf\n",vname,CS_ECES_A);}}
  if (strcmp(vname,"CS_ECIS_A")==0){ scanf("%lf",&CS_ECIS_A); if (verbose){ printf("%s read to be %lf\n",vname,CS_ECIS_A);}}
  if (strcmp(vname,"CS_ECEC_A")==0){ scanf("%lf",&CS_ECEC_A); if (verbose){ printf("%s read to be %lf\n",vname,CS_ECEC_A);}}
  if (strcmp(vname,"CS_ECIC_A")==0){ scanf("%lf",&CS_ECIC_A); if (verbose){ printf("%s read to be %lf\n",vname,CS_ECIC_A);}}
  if (strcmp(vname,"CS_ESES_N")==0){ scanf("%lf",&CS_ESES_N); if (verbose){ printf("%s read to be %lf\n",vname,CS_ESES_N);}}
  if (strcmp(vname,"CS_ESIS_N")==0){ scanf("%lf",&CS_ESIS_N); if (verbose){ printf("%s read to be %lf\n",vname,CS_ESIS_N);}}
  if (strcmp(vname,"CS_ESEC_N")==0){ scanf("%lf",&CS_ESEC_N); if (verbose){ printf("%s read to be %lf\n",vname,CS_ESEC_N);}}
  if (strcmp(vname,"CS_ESIC_N")==0){ scanf("%lf",&CS_ESIC_N); if (verbose){ printf("%s read to be %lf\n",vname,CS_ESIC_N);}}
  if (strcmp(vname,"CS_ECES_N")==0){ scanf("%lf",&CS_ECES_N); if (verbose){ printf("%s read to be %lf\n",vname,CS_ECES_N);}}
  if (strcmp(vname,"CS_ECIS_N")==0){ scanf("%lf",&CS_ECIS_N); if (verbose){ printf("%s read to be %lf\n",vname,CS_ECIS_N);}}
  if (strcmp(vname,"CS_ECEC_N")==0){ scanf("%lf",&CS_ECEC_N); if (verbose){ printf("%s read to be %lf\n",vname,CS_ECEC_N);}}
  if (strcmp(vname,"CS_ECIC_N")==0){ scanf("%lf",&CS_ECIC_N); if (verbose){ printf("%s read to be %lf\n",vname,CS_ECIC_N);}}
  if (strcmp(vname,"CS_ISES_G")==0){ scanf("%lf",&CS_ISES_G); if (verbose){ printf("%s read to be %lf\n",vname,CS_ISES_G);}}
  if (strcmp(vname,"CS_ISIS_G")==0){ scanf("%lf",&CS_ISIS_G); if (verbose){ printf("%s read to be %lf\n",vname,CS_ISIS_G);}}
  if (strcmp(vname,"CS_ISEC_G")==0){ scanf("%lf",&CS_ISEC_G); if (verbose){ printf("%s read to be %lf\n",vname,CS_ISEC_G);}}
  if (strcmp(vname,"CS_ISIC_G")==0){ scanf("%lf",&CS_ISIC_G); if (verbose){ printf("%s read to be %lf\n",vname,CS_ISIC_G);}}
  if (strcmp(vname,"CS_ICES_G")==0){ scanf("%lf",&CS_ICES_G); if (verbose){ printf("%s read to be %lf\n",vname,CS_ICES_G);}}
  if (strcmp(vname,"CS_ICIS_G")==0){ scanf("%lf",&CS_ICIS_G); if (verbose){ printf("%s read to be %lf\n",vname,CS_ICIS_G);}}
  if (strcmp(vname,"CS_ICEC_G")==0){ scanf("%lf",&CS_ICEC_G); if (verbose){ printf("%s read to be %lf\n",vname,CS_ICEC_G);}}
  if (strcmp(vname,"CS_ICIC_G")==0){ scanf("%lf",&CS_ICIC_G); if (verbose){ printf("%s read to be %lf\n",vname,CS_ICIC_G);}}
  if (strcmp(vname,"CS_ESLR")==0){ scanf("%lf",&CS_ESLR); if (verbose){ printf("%s read to be %lf\n",vname,CS_ESLR);}}
  if (strcmp(vname,"CS_ISLR")==0){ scanf("%lf",&CS_ISLR); if (verbose){ printf("%s read to be %lf\n",vname,CS_ISLR);}}
  if (strcmp(vname,"CS_ECLR")==0){ scanf("%lf",&CS_ECLR); if (verbose){ printf("%s read to be %lf\n",vname,CS_ECLR);}}
  if (strcmp(vname,"CS_ICLR")==0){ scanf("%lf",&CS_ICLR); if (verbose){ printf("%s read to be %lf\n",vname,CS_ICLR);}}
  if (strcmp(vname,"CS_LRES")==0){ scanf("%lf",&CS_LRES); if (verbose){ printf("%s read to be %lf\n",vname,CS_LRES);}}
  if (strcmp(vname,"CS_LRIS")==0){ scanf("%lf",&CS_LRIS); if (verbose){ printf("%s read to be %lf\n",vname,CS_LRIS);}}
  if (strcmp(vname,"CS_LREC")==0){ scanf("%lf",&CS_LREC); if (verbose){ printf("%s read to be %lf\n",vname,CS_LREC);}}
  if (strcmp(vname,"CS_LRIC")==0){ scanf("%lf",&CS_LRIC); if (verbose){ printf("%s read to be %lf\n",vname,CS_LRIC);}}
  if (strcmp(vname,"SPARSE_ESES")==0){ scanf("%d",&SPARSE_ESES); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ESES);}}
  if (strcmp(vname,"SPARSE_ESIS")==0){ scanf("%d",&SPARSE_ESIS); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ESIS);}}
  if (strcmp(vname,"SPARSE_ESEC")==0){ scanf("%d",&SPARSE_ESEC); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ESEC);}}
  if (strcmp(vname,"SPARSE_ESIC")==0){ scanf("%d",&SPARSE_ESIC); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ESIC);}}
  if (strcmp(vname,"SPARSE_ECES")==0){ scanf("%d",&SPARSE_ECES); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ECES);}}
  if (strcmp(vname,"SPARSE_ECIS")==0){ scanf("%d",&SPARSE_ECIS); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ECIS);}}
  if (strcmp(vname,"SPARSE_ECEC")==0){ scanf("%d",&SPARSE_ECEC); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ECEC);}}
  if (strcmp(vname,"SPARSE_ECIC")==0){ scanf("%d",&SPARSE_ECIC); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ECIC);}}
  if (strcmp(vname,"SPARSE_ISES")==0){ scanf("%d",&SPARSE_ISES); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ISES);}}
  if (strcmp(vname,"SPARSE_ISIS")==0){ scanf("%d",&SPARSE_ISIS); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ISIS);}}
  if (strcmp(vname,"SPARSE_ISEC")==0){ scanf("%d",&SPARSE_ISEC); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ISEC);}}
  if (strcmp(vname,"SPARSE_ISIC")==0){ scanf("%d",&SPARSE_ISIC); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ISIC);}}
  if (strcmp(vname,"SPARSE_ICES")==0){ scanf("%d",&SPARSE_ICES); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ICES);}}
  if (strcmp(vname,"SPARSE_ICIS")==0){ scanf("%d",&SPARSE_ICIS); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ICIS);}}
  if (strcmp(vname,"SPARSE_ICEC")==0){ scanf("%d",&SPARSE_ICEC); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ICEC);}}
  if (strcmp(vname,"SPARSE_ICIC")==0){ scanf("%d",&SPARSE_ICIC); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ICIC);}}
  if (strcmp(vname,"LR_PCONNECT")==0){ scanf("%lf",&LR_PCONNECT); if (verbose){ printf("%s read to be %lf\n",vname,LR_PCONNECT);}}
  if (strcmp(vname,"TAU_REF")==0){ scanf("%lf",&TAU_REF); if (verbose){ printf("%s read to be %lf\n",vname,TAU_REF);}}
  if (strcmp(vname,"TAU_AMPA")==0){ scanf("%lf",&TAU_AMPA); if (verbose){ printf("%s read to be %lf\n",vname,TAU_AMPA);}}
  if (strcmp(vname,"TAU_NMDA")==0){ scanf("%lf",&TAU_NMDA); if (verbose){ printf("%s read to be %lf\n",vname,TAU_NMDA);}}
  if (strcmp(vname,"TAU_GABA")==0){ scanf("%lf",&TAU_GABA); if (verbose){ printf("%s read to be %lf\n",vname,TAU_GABA);}}
  if (strcmp(vname,"DEPRESSION")==0){ scanf("%d",&DEPRESSION); if (verbose){ printf("%s read to be %d\n",vname,DEPRESSION);}}
  if (strcmp(vname,"DEPRESS_pA")==0){ scanf("%lf",&DEPRESS_pA); if (verbose){ printf("%s read to be %lf\n",vname,DEPRESS_pA);}}
  if (strcmp(vname,"DEPRESS_pN")==0){ scanf("%lf",&DEPRESS_pN); if (verbose){ printf("%s read to be %lf\n",vname,DEPRESS_pN);}}
  if (strcmp(vname,"DEPRESS_pG")==0){ scanf("%lf",&DEPRESS_pG); if (verbose){ printf("%s read to be %lf\n",vname,DEPRESS_pG);}}
  if (strcmp(vname,"TAU_DEPRESS_AMPA")==0){ scanf("%lf",&TAU_DEPRESS_AMPA); if (verbose){ printf("%s read to be %lf\n",vname,TAU_DEPRESS_AMPA);}}
  if (strcmp(vname,"TAU_DEPRESS_NMDA")==0){ scanf("%lf",&TAU_DEPRESS_NMDA); if (verbose){ printf("%s read to be %lf\n",vname,TAU_DEPRESS_NMDA);}}
  if (strcmp(vname,"TAU_DEPRESS_GABA")==0){ scanf("%lf",&TAU_DEPRESS_GABA); if (verbose){ printf("%s read to be %lf\n",vname,TAU_DEPRESS_GABA);}}
  if (strcmp(vname,"TAU_DEPRESS_LR")==0){ scanf("%lf",&TAU_DEPRESS_LR); if (verbose){ printf("%s read to be %lf\n",vname,TAU_DEPRESS_LR);}}
  if (strcmp(vname,"VOLTAGE_REST")==0){ scanf("%lf",&VOLTAGE_REST); if (verbose){ printf("%s read to be %lf\n",vname,VOLTAGE_REST);}}
  if (strcmp(vname,"VOLTAGE_RESET")==0){ scanf("%lf",&VOLTAGE_RESET); if (verbose){ printf("%s read to be %lf\n",vname,VOLTAGE_RESET);}}
  if (strcmp(vname,"VOLTAGE_THRESHOLD")==0){ scanf("%lf",&VOLTAGE_THRESHOLD); if (verbose){ printf("%s read to be %lf\n",vname,VOLTAGE_THRESHOLD);}}
  if (strcmp(vname,"VOLTAGE_EX")==0){ scanf("%lf",&VOLTAGE_EX); if (verbose){ printf("%s read to be %lf\n",vname,VOLTAGE_EX);}}
  if (strcmp(vname,"VOLTAGE_IN")==0){ scanf("%lf",&VOLTAGE_IN); if (verbose){ printf("%s read to be %lf\n",vname,VOLTAGE_IN);}}
  if (strcmp(vname,"VOLTAGE_THRESHOLD_EIF")==0){ scanf("%lf",&VOLTAGE_THRESHOLD_EIF); if (verbose){ printf("%s read to be %lf\n",vname,VOLTAGE_THRESHOLD_EIF);}}
  if (strcmp(vname,"VOLTAGE_TAKEOFF")==0){ scanf("%lf",&VOLTAGE_TAKEOFF); if (verbose){ printf("%s read to be %lf\n",vname,VOLTAGE_TAKEOFF);}}
  if (strcmp(vname,"VOLTAGE_DELTAT")==0){ scanf("%lf",&VOLTAGE_DELTAT); if (verbose){ printf("%s read to be %lf\n",vname,VOLTAGE_DELTAT);}}
  if (strcmp(vname,"CONDUCTANCE_LK")==0){ scanf("%lf",&CONDUCTANCE_LK); if (verbose){ printf("%s read to be %lf\n",vname,CONDUCTANCE_LK);}}
  if (strcmp(vname,"GLOBAL_TI")==0){ scanf("%lf",&GLOBAL_TI); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_TI);}}
  if (strcmp(vname,"GLOBAL_TF")==0){ scanf("%lf",&GLOBAL_TF); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_TF);}}
  if (strcmp(vname,"GLOBAL_DTmax")==0){ scanf("%lf",&GLOBAL_DTmax); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_DTmax);}}
  if (strcmp(vname,"GLOBAL_dtadapt")==0){ scanf("%d",&GLOBAL_dtadapt); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_dtadapt);}}
  if (strcmp(vname,"GLOBAL_verbose")==0){ scanf("%d",&GLOBAL_verbose); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_verbose);}}
  if (strcmp(vname,"FIDDLE_ROW")==0){ scanf("%d",&FIDDLE_ROW); if (verbose){ printf("%s read to be %d\n",vname,FIDDLE_ROW);}}
  if (strcmp(vname,"FIDDLE_COL")==0){ scanf("%d",&FIDDLE_COL); if (verbose){ printf("%s read to be %d\n",vname,FIDDLE_COL);}}
  if (strcmp(vname,"FIDDLE_PARAMETER")==0){ scanf("%d",&FIDDLE_PARAMETER); if (verbose){ printf("%s read to be %d\n",vname,FIDDLE_PARAMETER);}}
  if (strcmp(vname,"DRAW_FLAG")==0){ scanf("%d",&DRAW_FLAG); if (verbose){ printf("%s read to be %d\n",vname,DRAW_FLAG);}}
  if (strcmp(vname,"DRAW_FLAG2")==0){ scanf("%d",&DRAW_FLAG2); if (verbose){ printf("%s read to be %d\n",vname,DRAW_FLAG2);}}
  if (strcmp(vname,"STEPS_PER_DRAW")==0){ scanf("%d",&STEPS_PER_DRAW); if (verbose){ printf("%s read to be %d\n",vname,STEPS_PER_DRAW);}}
  if (strcmp(vname,"STD_VIEW")==0){ scanf("%lf",&STD_VIEW); if (verbose){ printf("%s read to be %lf\n",vname,STD_VIEW);}}
  if (strcmp(vname,"PTREE_VIEW")==0){ scanf("%lf",&PTREE_VIEW); if (verbose){ printf("%s read to be %lf\n",vname,PTREE_VIEW);}}
  if (strcmp(vname,"GRAYSCALE")==0){ scanf("%d",&GRAYSCALE); if (verbose){ printf("%s read to be %d\n",vname,GRAYSCALE);}}
  if (strcmp(vname,"SUPERGLOBAL_DRAW_FLAG")==0){ scanf("%d",&SUPERGLOBAL_DRAW_FLAG); if (verbose){ printf("%s read to be %d\n",vname,SUPERGLOBAL_DRAW_FLAG);}}
  if (strcmp(vname,"MOVIE_BOTHER")==0){ scanf("%d",&MOVIE_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,MOVIE_BOTHER);}}
  if (strcmp(vname,"MOVIE_TIME_REFRESH")==0){ scanf("%lf",&MOVIE_TIME_REFRESH); if (verbose){ printf("%s read to be %lf\n",vname,MOVIE_TIME_REFRESH);}}
  if (strcmp(vname,"MOVIE_START_RECORDING")==0){ scanf("%lf",&MOVIE_START_RECORDING); if (verbose){ printf("%s read to be %lf\n",vname,MOVIE_START_RECORDING);}}
  if (strcmp(vname,"PNM_BOTHER")==0){ scanf("%d",&PNM_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,PNM_BOTHER);}}
  if (strcmp(vname,"PNM_REMOVE")==0){ scanf("%d",&PNM_REMOVE); if (verbose){ printf("%s read to be %d\n",vname,PNM_REMOVE);}}
  if (strcmp(vname,"OUTPUT_DUMP_EVERY")==0){ scanf("%lf",&OUTPUT_DUMP_EVERY); if (verbose){ printf("%s read to be %lf\n",vname,OUTPUT_DUMP_EVERY);}}
  if (strcmp(vname,"RTC_BOTHER")==0){ scanf("%d",&RTC_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,RTC_BOTHER);}}
  if (strcmp(vname,"GLOBAL_RTC_LENGTH")==0){ scanf("%d",&GLOBAL_RTC_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_RTC_LENGTH);}}
  if (strcmp(vname,"GLOBAL_RTC_FRAMELENGTH")==0){ scanf("%d",&GLOBAL_RTC_FRAMELENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_RTC_FRAMELENGTH);}}
  if (strcmp(vname,"GLOBAL_RTC_NANGLES")==0){ scanf("%d",&GLOBAL_RTC_NANGLES); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_RTC_NANGLES);}}
  if (strcmp(vname,"GLOBAL_RTC_NPHASES")==0){ scanf("%d",&GLOBAL_RTC_NPHASES); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_RTC_NPHASES);}}
  if (strcmp(vname,"STROBETRACE_BOTHER")==0){ scanf("%d",&STROBETRACE_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,STROBETRACE_BOTHER);}}
  if (strcmp(vname,"GLOBAL_STROBETRACE_LENGTH")==0){ scanf("%lf",&GLOBAL_STROBETRACE_LENGTH); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_STROBETRACE_LENGTH);}}
  if (strcmp(vname,"GLOBAL_STROBETRACE_NANGLES")==0){ scanf("%d",&GLOBAL_STROBETRACE_NANGLES); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_STROBETRACE_NANGLES);}}
  if (strcmp(vname,"GLOBAL_STROBETRACE_AVALANCHE_UPDATE_EVERY")==0){ scanf("%lf",&GLOBAL_STROBETRACE_AVALANCHE_UPDATE_EVERY); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_STROBETRACE_AVALANCHE_UPDATE_EVERY);}}
  if (strcmp(vname,"GLOBAL_STROBETRACE_AVALANCHE_LOGDMIN")==0){ scanf("%d",&GLOBAL_STROBETRACE_AVALANCHE_LOGDMIN); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_STROBETRACE_AVALANCHE_LOGDMIN);}}
  if (strcmp(vname,"GLOBAL_STROBETRACE_AVALANCHE_LOGDMAX")==0){ scanf("%d",&GLOBAL_STROBETRACE_AVALANCHE_LOGDMAX); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_STROBETRACE_AVALANCHE_LOGDMAX);}}
  if (strcmp(vname,"TUNINGCURVE_BOTHER")==0){ scanf("%d",&TUNINGCURVE_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,TUNINGCURVE_BOTHER);}}
  if (strcmp(vname,"GLOBAL_TUNINGCURVE_NANGLES")==0){ scanf("%d",&GLOBAL_TUNINGCURVE_NANGLES); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_TUNINGCURVE_NANGLES);}}
  if (strcmp(vname,"GLOBAL_TUNINGCURVE_NRADIUS")==0){ scanf("%d",&GLOBAL_TUNINGCURVE_NRADIUS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_TUNINGCURVE_NRADIUS);}}
  if (strcmp(vname,"LMITRI_BOTHER")==0){ scanf("%d",&LMITRI_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,LMITRI_BOTHER);}}
  if (strcmp(vname,"GLOBAL_LMITRI_T0")==0){ scanf("%d",&GLOBAL_LMITRI_T0); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_LMITRI_T0);}}
  if (strcmp(vname,"GLOBAL_LMITRI_TIMELENGTH")==0){ scanf("%d",&GLOBAL_LMITRI_TIMELENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_LMITRI_TIMELENGTH);}}
  if (strcmp(vname,"GLOBAL_LMITRI_ROW_MAX")==0){ scanf("%lf",&GLOBAL_LMITRI_ROW_MAX); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_LMITRI_ROW_MAX);}}
  if (strcmp(vname,"GLOBAL_LMITRI_ROW_MIN")==0){ scanf("%lf",&GLOBAL_LMITRI_ROW_MIN); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_LMITRI_ROW_MIN);}}
  if (strcmp(vname,"PTREE_BOTHER")==0){ scanf("%d",&PTREE_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,PTREE_BOTHER);}}
  if (strcmp(vname,"GLOBAL_PTREE_BITBYBIT")==0){ scanf("%d",&GLOBAL_PTREE_BITBYBIT); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_BITBYBIT);}}
  if (strcmp(vname,"GLOBAL_PTREE_ZZZ")==0){ scanf("%d",&GLOBAL_PTREE_ZZZ); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_ZZZ);}}
  if (strcmp(vname,"GLOBAL_PTREE_NLEGS")==0){ scanf("%d",&GLOBAL_PTREE_NLEGS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_NLEGS);}}
  if (strcmp(vname,"GLOBAL_PTREE_LEGTIME")==0){ scanf("%d",&GLOBAL_PTREE_LEGTIME); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_LEGTIME);}}
  if (strcmp(vname,"GLOBAL_PTREE_NREGIONS")==0){ scanf("%d",&GLOBAL_PTREE_NREGIONS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_NREGIONS);}}
  if (strcmp(vname,"P_REWIRE")==0){ scanf("%lf",&P_REWIRE); if (verbose){ printf("%s read to be %lf\n",vname,P_REWIRE);}}
  if (strcmp(vname,"SPARSE_ROW_DIA")==0){ scanf("%d",&SPARSE_ROW_DIA); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_ROW_DIA);}}
  if (strcmp(vname,"SPARSE_COL_DIA")==0){ scanf("%d",&SPARSE_COL_DIA); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_COL_DIA);}}
  if (strcmp(vname,"SPARSE_T2S_ORDERING")==0){ scanf("%d",&SPARSE_T2S_ORDERING); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_T2S_ORDERING);}}
  if (strcmp(vname,"SPARSE_CONNECTION_TYPE")==0){ scanf("%d",&SPARSE_CONNECTION_TYPE); if (verbose){ printf("%s read to be %d\n",vname,SPARSE_CONNECTION_TYPE);}}
  if (strcmp(vname,"SPARSE_CS_SCALE")==0){ scanf("%lf",&SPARSE_CS_SCALE); if (verbose){ printf("%s read to be %lf\n",vname,SPARSE_CS_SCALE);}}
  if (strcmp(vname,"GLOBAL_PTREE_EVENT_THRESHOLD")==0){ scanf("%d",&GLOBAL_PTREE_EVENT_THRESHOLD); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_EVENT_THRESHOLD);}}
  if (strcmp(vname,"GLOBAL_PTREE_EVENT_WITHIN")==0){ scanf("%lf",&GLOBAL_PTREE_EVENT_WITHIN); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_PTREE_EVENT_WITHIN);}}
  if (strcmp(vname,"GLOBAL_PTREE_REGION_TYPE")==0){ scanf("%d",&GLOBAL_PTREE_REGION_TYPE); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_REGION_TYPE);}}
  if (strcmp(vname,"CLOSET_BOTHER")==0){ scanf("%d",&CLOSET_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,CLOSET_BOTHER);}}
  if (strcmp(vname,"YGGDRASIL_BOTHER")==0){ scanf("%d",&YGGDRASIL_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,YGGDRASIL_BOTHER);}}
  if (strcmp(vname,"GLOBAL_YGGDRASIL_PPNREGIONS")==0){ scanf("%d",&GLOBAL_YGGDRASIL_PPNREGIONS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_YGGDRASIL_PPNREGIONS);}}
  if (strcmp(vname,"GLOBAL_YGGDRASIL_PPNLEGS")==0){ scanf("%d",&GLOBAL_YGGDRASIL_PPNLEGS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_YGGDRASIL_PPNLEGS);}}
  if (strcmp(vname,"GLOBAL_YGGDRASIL_PPLEGTIME")==0){ scanf("%d",&GLOBAL_YGGDRASIL_PPLEGTIME); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_YGGDRASIL_PPLEGTIME);}}
  if (strcmp(vname,"GLOBAL_YGGDRASIL_WEIGHT_MINIMUM")==0){ scanf("%d",&GLOBAL_YGGDRASIL_WEIGHT_MINIMUM); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_YGGDRASIL_WEIGHT_MINIMUM);}}
  if (strcmp(vname,"BONSAI_BOTHER")==0){ scanf("%d",&BONSAI_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,BONSAI_BOTHER);}}
  if (strcmp(vname,"HYDRA_BOTHER")==0){ scanf("%d",&HYDRA_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,HYDRA_BOTHER);}}
  if (strcmp(vname,"GLOBAL_HYDRA_SWITCH_EVERY_TWO")==0){ scanf("%d",&GLOBAL_HYDRA_SWITCH_EVERY_TWO); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_HYDRA_SWITCH_EVERY_TWO);}}
  if (strcmp(vname,"GLOBAL_HYDRA_JUSTONTIME")==0){ scanf("%d",&GLOBAL_HYDRA_JUSTONTIME); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_HYDRA_JUSTONTIME);}}
  if (strcmp(vname,"GLOBAL_HYDRA_STAYONTIME")==0){ scanf("%d",&GLOBAL_HYDRA_STAYONTIME); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_HYDRA_STAYONTIME);}}
  if (strcmp(vname,"GLOBAL_HYDRA_DUMP_EVERY")==0){ scanf("%d",&GLOBAL_HYDRA_DUMP_EVERY); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_HYDRA_DUMP_EVERY);}}
  if (strcmp(vname,"LYAPUNOV_BOTHER")==0){ scanf("%d",&LYAPUNOV_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,LYAPUNOV_BOTHER);}}
  if (strcmp(vname,"GLOBAL_LYAPUNOV_UPDATE_EVERY")==0){ scanf("%lf",&GLOBAL_LYAPUNOV_UPDATE_EVERY); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_LYAPUNOV_UPDATE_EVERY);}}
  if (strcmp(vname,"GLOBAL_LYAPUNOV_JIGGLE")==0){ scanf("%lf",&GLOBAL_LYAPUNOV_JIGGLE); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_LYAPUNOV_JIGGLE);}}
  if (strcmp(vname,"POWER_BOTHER")==0){ scanf("%d",&POWER_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,POWER_BOTHER);}}
  if (strcmp(vname,"GLOBAL_POWER_LENGTH")==0){ scanf("%d",&GLOBAL_POWER_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_POWER_LENGTH);}}
  if (strcmp(vname,"GLOBAL_POWER_HOWMANY")==0){ scanf("%d",&GLOBAL_POWER_HOWMANY); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_POWER_HOWMANY);}}
  if (strcmp(vname,"TAOF_BOTHER")==0){ scanf("%d",&TAOF_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,TAOF_BOTHER);}}
  if (strcmp(vname,"GLOBAL_TAOF_LENGTH")==0){ scanf("%d",&GLOBAL_TAOF_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_TAOF_LENGTH);}}
  if (strcmp(vname,"GLOBAL_TAOF_STEP_EVERY")==0){ scanf("%d",&GLOBAL_TAOF_STEP_EVERY); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_TAOF_STEP_EVERY);}}
  if (strcmp(vname,"SEIDCORR_BOTHER")==0){ scanf("%d",&SEIDCORR_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,SEIDCORR_BOTHER);}}
  if (strcmp(vname,"GLOBAL_SEIDCORR_SPACE_BIN_SIZE")==0){ scanf("%d",&GLOBAL_SEIDCORR_SPACE_BIN_SIZE); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_SEIDCORR_SPACE_BIN_SIZE);}}
  if (strcmp(vname,"GLOBAL_SEIDCORR_TIME_BIN_SIZE")==0){ scanf("%d",&GLOBAL_SEIDCORR_TIME_BIN_SIZE); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_SEIDCORR_TIME_BIN_SIZE);}}
  if (strcmp(vname,"GLOBAL_SEIDCORR_LENGTH")==0){ scanf("%d",&GLOBAL_SEIDCORR_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_SEIDCORR_LENGTH);}}
  if (strcmp(vname,"GLOBAL_SEIDCORR_TIME_START")==0){ scanf("%lf",&GLOBAL_SEIDCORR_TIME_START); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_SEIDCORR_TIME_START);}}
  if (strcmp(vname,"END")==0){ /* do nothing */ if (verbose){ printf("end of input reached\n");}}
/*   if (strcmp(vname,"yy")==0){ scanf("%zz",&yy); if (verbose){ printf("%s read to be %zz\n",vname,yy);}} */
}

void dumpoutput(char *filename)
{
  /* prints an output file */
  char text[256];
  char vname[128];
  FILE *fp=NULL;
  sprintf(text,"./%s%s",filename,GLOBAL_STRING);
  if (filename==NULL || (fp = fopen(text, "w")) == NULL){ printf("dumpoutput to stdout\n"); fp = stdout;}
  sprintf(vname,"ON_MY_COMPUTER"); fprintf(fp,"%s= %d\n",vname,ON_MY_COMPUTER);
  sprintf(vname,"GLOBAL_STRING"); fprintf(fp,"%s= %s\n",vname,GLOBAL_STRING);
  sprintf(vname,"GLOBAL_RECORD_NUMBER"); fprintf(fp,"%s= %d\n",vname,GLOBAL_RECORD_NUMBER);
  sprintf(vname,"GLOBAL_SPIKEINPUT_RSEED"); fprintf(fp,"%s= %d\n",vname,GLOBAL_SPIKEINPUT_RSEED);
  sprintf(vname,"GLOBAL_CLEANUP"); fprintf(fp,"%s= %d\n",vname,GLOBAL_CLEANUP);
  sprintf(vname,"NSLICES"); fprintf(fp,"%s= %d\n",vname,NSLICES);
  sprintf(vname,"NPHASES"); fprintf(fp,"%s= %d\n",vname,NPHASES);
  sprintf(vname,"SUITE_BOTHER"); fprintf(fp,"%s= %d\n",vname,SUITE_BOTHER);
  sprintf(vname,"SUITE_NSECONDS"); fprintf(fp,"%s= %d\n",vname,SUITE_NSECONDS);
  sprintf(vname,"SUITE_DUMPEVERY"); fprintf(fp,"%s= %d\n",vname,SUITE_DUMPEVERY);
  sprintf(vname,"SUITE_BITBYBIT_REMOVE"); fprintf(fp,"%s= %d\n",vname,SUITE_BITBYBIT_REMOVE);
  sprintf(vname,"SUITE_SINDEXMAX"); fprintf(fp,"%s= %d\n",vname,SUITE_SINDEXMAX);
  sprintf(vname,"SUITE_DINDEXMAX"); fprintf(fp,"%s= %d\n",vname,SUITE_DINDEXMAX);
  sprintf(vname,"SUITE_TINDEXMAX"); fprintf(fp,"%s= %d\n",vname,SUITE_TINDEXMAX);
  sprintf(vname,"LGN_BOTHER"); fprintf(fp,"%s= %d\n",vname,LGN_BOTHER);
  sprintf(vname,"LGN_DUMP"); fprintf(fp,"%s= %d\n",vname,LGN_DUMP);
  sprintf(vname,"LGN_TYPE_FLAG"); fprintf(fp,"%s= %d\n",vname,LGN_TYPE_FLAG);
  sprintf(vname,"LGN_DUMB"); fprintf(fp,"%s= %d\n",vname,LGN_DUMB);
  sprintf(vname,"GLOBAL_ODOR"); fprintf(fp,"%s= %d\n",vname,GLOBAL_ODOR);
  sprintf(vname,"GLOBAL_ODOR_BASE"); fprintf(fp,"%s= %d\n",vname,GLOBAL_ODOR_BASE);
  sprintf(vname,"GLOBAL_ODOR_BACON"); fprintf(fp,"%s= %d\n",vname,GLOBAL_ODOR_BACON);
  sprintf(vname,"GLOBAL_FPS"); fprintf(fp,"%s= %d\n",vname,GLOBAL_FPS);
  sprintf(vname,"OBSCURED"); fprintf(fp,"%s= %d\n",vname,OBSCURED);
  sprintf(vname,"GRATING_VS_LMI"); fprintf(fp,"%s= %d\n",vname,GRATING_VS_LMI);
  sprintf(vname,"CYCLE_LENGTH"); fprintf(fp,"%s= %0.16lf\n",vname,CYCLE_LENGTH);
  sprintf(vname,"GRATING_PULSE"); fprintf(fp,"%s= %0.16lf\n",vname,GRATING_PULSE);
  sprintf(vname,"GRATING_DRIFT"); fprintf(fp,"%s= %.16lf\n",vname,GRATING_DRIFT);
  sprintf(vname,"STIMULUS_ONSET_TIME"); fprintf(fp,"%s= %.16lf\n",vname,STIMULUS_ONSET_TIME);
  sprintf(vname,"SQUARE_DURATION_TIME"); fprintf(fp,"%s= %.16lf\n",vname,SQUARE_DURATION_TIME);
  sprintf(vname,"SQUARE_DRAG_TIME"); fprintf(fp,"%s= %.16lf\n",vname,SQUARE_DRAG_TIME);
  sprintf(vname,"LINE_DELAY_TIME"); fprintf(fp,"%s= %.16lf\n",vname,LINE_DELAY_TIME);
  sprintf(vname,"LGN_BACKRATE"); fprintf(fp,"%s= %0.16lf\n",vname,LGN_BACKRATE);
  sprintf(vname,"OTHERLAYER_BACKRATE"); fprintf(fp,"%s= %0.16lf\n",vname,OTHERLAYER_BACKRATE);
  sprintf(vname,"OTHERLAYER_INPUTRATE"); fprintf(fp,"%s= %0.16lf\n",vname,OTHERLAYER_INPUTRATE);
  sprintf(vname,"INPUT_CONTRAST"); fprintf(fp,"%s= %0.16lf\n",vname,INPUT_CONTRAST);
  sprintf(vname,"INPUT_SPACEK"); fprintf(fp,"%s= %0.16lf\n",vname,INPUT_SPACEK);
  sprintf(vname,"INPUT_SPACEANGLE"); fprintf(fp,"%s= %0.16lf\n",vname,INPUT_SPACEANGLE);
  sprintf(vname,"INPUT_SPACEANGLE_BACON"); fprintf(fp,"%s= %0.16lf\n",vname,INPUT_SPACEANGLE_BACON);
  sprintf(vname,"INPUT_SPACEPHASE"); fprintf(fp,"%s= %0.16lf\n",vname,INPUT_SPACEPHASE);
  sprintf(vname,"CORTEX_BOTHER"); fprintf(fp,"%s= %d\n",vname,CORTEX_BOTHER);
  sprintf(vname,"AUTAPSES_OFF"); fprintf(fp,"%s= %d\n",vname,AUTAPSES_OFF);
  sprintf(vname,"LR_BOTHER"); fprintf(fp,"%s= %d\n",vname,LR_BOTHER);
  sprintf(vname,"EIFvsIF"); fprintf(fp,"%s= %d\n",vname,EIFvsIF);
  sprintf(vname,"ORDER"); fprintf(fp,"%s= %d\n",vname,ORDER);
  sprintf(vname,"LGN_DETERMINISM"); fprintf(fp,"%s= %d\n",vname,LGN_DETERMINISM);
  sprintf(vname,"HANSHELLEY_FLAG"); fprintf(fp,"%s= %d\n",vname,HANSHELLEY_FLAG);
  sprintf(vname,"STARTING_HOMOGENIZATION"); fprintf(fp,"%s= %0.16lf\n",vname,STARTING_HOMOGENIZATION);
  sprintf(vname,"SPIKETOL"); fprintf(fp,"%s= %d\n",vname,SPIKETOL);
  sprintf(vname,"GLOBAL_SPACE_SMOOTHER"); fprintf(fp,"%s= %d\n",vname,GLOBAL_SPACE_SMOOTHER);
  sprintf(vname,"ARBOR_DIA"); fprintf(fp,"%s= %d\n",vname,ARBOR_DIA);
  sprintf(vname,"BIG_SYSTEM_FLAG"); fprintf(fp,"%s= %d\n",vname,BIG_SYSTEM_FLAG);
  sprintf(vname,"LR_ANGLE_DIA"); fprintf(fp,"%s= %0.16lf\n",vname,LR_ANGLE_DIA);
  sprintf(vname,"LR_DIST_DIA_INPIES"); fprintf(fp,"%s= %0.16lf\n",vname,LR_DIST_DIA_INPIES);
  sprintf(vname,"LR_DIR_DIA"); fprintf(fp,"%s= %0.16lf\n",vname,LR_DIR_DIA);
  sprintf(vname,"LR_RATE"); fprintf(fp,"%s= %0.16lf\n",vname,LR_RATE);
  sprintf(vname,"LR_TO_AMPA"); fprintf(fp,"%s= %0.16lf\n",vname,LR_TO_AMPA);
  sprintf(vname,"AXONAL_DELAY_PER_PIE"); fprintf(fp,"%s= %0.16lf\n",vname,AXONAL_DELAY_PER_PIE);
  sprintf(vname,"P_EC"); fprintf(fp,"%s= %0.16lf\n",vname,P_EC);
  sprintf(vname,"P_ES"); fprintf(fp,"%s= %0.16lf\n",vname,P_ES);
  sprintf(vname,"P_IC"); fprintf(fp,"%s= %0.16lf\n",vname,P_IC);
  sprintf(vname,"P_IS"); fprintf(fp,"%s= %0.16lf\n",vname,P_IS);
  sprintf(vname,"P_AMPA"); fprintf(fp,"%s= %0.16lf\n",vname,P_AMPA);
  sprintf(vname,"P_NMDA"); fprintf(fp,"%s= %0.16lf\n",vname,P_NMDA);
  sprintf(vname,"P_GABA"); fprintf(fp,"%s= %0.16lf\n",vname,P_GABA);
  sprintf(vname,"GABOR_DRIFT"); fprintf(fp,"%s= %0.16lf\n",vname,GABOR_DRIFT);
  sprintf(vname,"LGNANGLE_DRIFT"); fprintf(fp,"%s= %0.16lf\n",vname,LGNANGLE_DRIFT);
  sprintf(vname,"LGN_STRENGTH"); fprintf(fp,"%s= %0.16lf\n",vname,LGN_STRENGTH);
  sprintf(vname,"OTHERLAYER_STRENGTH"); fprintf(fp,"%s= %0.16lf\n",vname,OTHERLAYER_STRENGTH);
  sprintf(vname,"CS_ESES_A"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ESES_A);
  sprintf(vname,"CS_ESIS_A"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ESIS_A);
  sprintf(vname,"CS_ESEC_A"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ESEC_A);
  sprintf(vname,"CS_ESIC_A"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ESIC_A);
  sprintf(vname,"CS_ECES_A"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ECES_A);
  sprintf(vname,"CS_ECIS_A"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ECIS_A);
  sprintf(vname,"CS_ECEC_A"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ECEC_A);
  sprintf(vname,"CS_ECIC_A"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ECIC_A);
  sprintf(vname,"CS_ESES_N"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ESES_N);
  sprintf(vname,"CS_ESIS_N"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ESIS_N);
  sprintf(vname,"CS_ESEC_N"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ESEC_N);
  sprintf(vname,"CS_ESIC_N"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ESIC_N);
  sprintf(vname,"CS_ECES_N"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ECES_N);
  sprintf(vname,"CS_ECIS_N"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ECIS_N);
  sprintf(vname,"CS_ECEC_N"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ECEC_N);
  sprintf(vname,"CS_ECIC_N"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ECIC_N);
  sprintf(vname,"CS_ISES_G"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ISES_G);
  sprintf(vname,"CS_ISIS_G"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ISIS_G);
  sprintf(vname,"CS_ISEC_G"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ISEC_G);
  sprintf(vname,"CS_ISIC_G"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ISIC_G);
  sprintf(vname,"CS_ICES_G"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ICES_G);
  sprintf(vname,"CS_ICIS_G"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ICIS_G);
  sprintf(vname,"CS_ICEC_G"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ICEC_G);
  sprintf(vname,"CS_ICIC_G"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ICIC_G);
  sprintf(vname,"CS_ESLR"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ESLR);
  sprintf(vname,"CS_ISLR"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ISLR);
  sprintf(vname,"CS_ECLR"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ECLR);
  sprintf(vname,"CS_ICLR"); fprintf(fp,"%s= %0.16lf\n",vname,CS_ICLR);
  sprintf(vname,"CS_LRES"); fprintf(fp,"%s= %0.16lf\n",vname,CS_LRES);
  sprintf(vname,"CS_LRIS"); fprintf(fp,"%s= %0.16lf\n",vname,CS_LRIS);
  sprintf(vname,"CS_LREC"); fprintf(fp,"%s= %0.16lf\n",vname,CS_LREC);
  sprintf(vname,"CS_LRIC"); fprintf(fp,"%s= %0.16lf\n",vname,CS_LRIC);
  sprintf(vname,"SPARSE_ESES"); fprintf(fp,"%s= %d\n",vname,SPARSE_ESES);
  sprintf(vname,"SPARSE_ESIS"); fprintf(fp,"%s= %d\n",vname,SPARSE_ESIS);
  sprintf(vname,"SPARSE_ESEC"); fprintf(fp,"%s= %d\n",vname,SPARSE_ESEC);
  sprintf(vname,"SPARSE_ESIC"); fprintf(fp,"%s= %d\n",vname,SPARSE_ESIC);
  sprintf(vname,"SPARSE_ECES"); fprintf(fp,"%s= %d\n",vname,SPARSE_ECES);
  sprintf(vname,"SPARSE_ECIS"); fprintf(fp,"%s= %d\n",vname,SPARSE_ECIS);
  sprintf(vname,"SPARSE_ECEC"); fprintf(fp,"%s= %d\n",vname,SPARSE_ECEC);
  sprintf(vname,"SPARSE_ECIC"); fprintf(fp,"%s= %d\n",vname,SPARSE_ECIC);
  sprintf(vname,"SPARSE_ISES"); fprintf(fp,"%s= %d\n",vname,SPARSE_ISES);
  sprintf(vname,"SPARSE_ISIS"); fprintf(fp,"%s= %d\n",vname,SPARSE_ISIS);
  sprintf(vname,"SPARSE_ISEC"); fprintf(fp,"%s= %d\n",vname,SPARSE_ISEC);
  sprintf(vname,"SPARSE_ISIC"); fprintf(fp,"%s= %d\n",vname,SPARSE_ISIC);
  sprintf(vname,"SPARSE_ICES"); fprintf(fp,"%s= %d\n",vname,SPARSE_ICES);
  sprintf(vname,"SPARSE_ICIS"); fprintf(fp,"%s= %d\n",vname,SPARSE_ICIS);
  sprintf(vname,"SPARSE_ICEC"); fprintf(fp,"%s= %d\n",vname,SPARSE_ICEC);
  sprintf(vname,"SPARSE_ICIC"); fprintf(fp,"%s= %d\n",vname,SPARSE_ICIC);
  sprintf(vname,"LR_PCONNECT"); fprintf(fp,"%s= %0.16lf\n",vname,LR_PCONNECT);
  sprintf(vname,"TAU_REF"); fprintf(fp,"%s= %0.16lf\n",vname,TAU_REF);
  sprintf(vname,"TAU_AMPA"); fprintf(fp,"%s= %0.16lf\n",vname,TAU_AMPA);
  sprintf(vname,"TAU_NMDA"); fprintf(fp,"%s= %0.16lf\n",vname,TAU_NMDA);
  sprintf(vname,"TAU_GABA"); fprintf(fp,"%s= %0.16lf\n",vname,TAU_GABA);
  sprintf(vname,"DEPRESSION"); fprintf(fp,"%s= %d\n",vname,DEPRESSION);
  sprintf(vname,"DEPRESS_pA"); fprintf(fp,"%s= %0.16lf\n",vname,DEPRESS_pA);
  sprintf(vname,"DEPRESS_pN"); fprintf(fp,"%s= %0.16lf\n",vname,DEPRESS_pN);
  sprintf(vname,"DEPRESS_pG"); fprintf(fp,"%s= %0.16lf\n",vname,DEPRESS_pG);
  sprintf(vname,"TAU_DEPRESS_AMPA"); fprintf(fp,"%s= %0.16lf\n",vname,TAU_DEPRESS_AMPA);
  sprintf(vname,"TAU_DEPRESS_NMDA"); fprintf(fp,"%s= %0.16lf\n",vname,TAU_DEPRESS_NMDA);
  sprintf(vname,"TAU_DEPRESS_GABA"); fprintf(fp,"%s= %0.16lf\n",vname,TAU_DEPRESS_GABA);
  sprintf(vname,"TAU_DEPRESS_LR"); fprintf(fp,"%s= %0.16lf\n",vname,TAU_DEPRESS_LR);
  sprintf(vname,"VOLTAGE_REST"); fprintf(fp,"%s= %0.16lf\n",vname,VOLTAGE_REST);
  sprintf(vname,"VOLTAGE_RESET"); fprintf(fp,"%s= %0.16lf\n",vname,VOLTAGE_RESET);
  sprintf(vname,"VOLTAGE_THRESHOLD"); fprintf(fp,"%s= %0.16lf\n",vname,VOLTAGE_THRESHOLD);
  sprintf(vname,"VOLTAGE_EX"); fprintf(fp,"%s= %0.16lf\n",vname,VOLTAGE_EX);
  sprintf(vname,"VOLTAGE_IN"); fprintf(fp,"%s= %0.16lf\n",vname,VOLTAGE_IN);
  sprintf(vname,"VOLTAGE_THRESHOLD_EIF"); fprintf(fp,"%s= %0.16lf\n",vname,VOLTAGE_THRESHOLD_EIF);
  sprintf(vname,"VOLTAGE_TAKEOFF"); fprintf(fp,"%s= %0.16lf\n",vname,VOLTAGE_TAKEOFF);
  sprintf(vname,"VOLTAGE_DELTAT"); fprintf(fp,"%s= %0.16lf\n",vname,VOLTAGE_DELTAT);
  sprintf(vname,"CONDUCTANCE_LK"); fprintf(fp,"%s= %0.16lf\n",vname,CONDUCTANCE_LK);
  sprintf(vname,"GLOBAL_TI"); fprintf(fp,"%s= %0.16lf\n",vname,GLOBAL_TI);
  sprintf(vname,"GLOBAL_TF"); fprintf(fp,"%s= %0.16lf\n",vname,GLOBAL_TF);
  sprintf(vname,"GLOBAL_DTmax"); fprintf(fp,"%s= %0.16lf\n",vname,GLOBAL_DTmax);
  sprintf(vname,"GLOBAL_dtadapt"); fprintf(fp,"%s= %d\n",vname,GLOBAL_dtadapt);
  sprintf(vname,"GLOBAL_verbose"); fprintf(fp,"%s= %d\n",vname,GLOBAL_verbose);
  sprintf(vname,"FIDDLE_ROW"); fprintf(fp,"%s= %d\n",vname,FIDDLE_ROW);
  sprintf(vname,"FIDDLE_COL"); fprintf(fp,"%s= %d\n",vname,FIDDLE_COL);
  sprintf(vname,"FIDDLE_PARAMETER"); fprintf(fp,"%s= %d\n",vname,FIDDLE_PARAMETER);
  sprintf(vname,"DRAW_FLAG"); fprintf(fp,"%s= %d\n",vname,DRAW_FLAG);
  sprintf(vname,"DRAW_FLAG2"); fprintf(fp,"%s= %d\n",vname,DRAW_FLAG2);
  sprintf(vname,"STEPS_PER_DRAW"); fprintf(fp,"%s= %d\n",vname,STEPS_PER_DRAW);
  sprintf(vname,"STD_VIEW"); fprintf(fp,"%s= %0.16lf\n",vname,STD_VIEW);
  sprintf(vname,"PTREE_VIEW"); fprintf(fp,"%s= %0.16lf\n",vname,PTREE_VIEW);
  sprintf(vname,"GRAYSCALE"); fprintf(fp,"%s= %d\n",vname,GRAYSCALE);
  sprintf(vname,"SUPERGLOBAL_DRAW_FLAG"); fprintf(fp,"%s= %d\n",vname,SUPERGLOBAL_DRAW_FLAG);
  sprintf(vname,"MOVIE_BOTHER"); fprintf(fp,"%s= %d\n",vname,MOVIE_BOTHER);
  sprintf(vname,"MOVIE_TIME_REFRESH"); fprintf(fp,"%s= %0.16lf\n",vname,MOVIE_TIME_REFRESH);
  sprintf(vname,"MOVIE_START_RECORDING"); fprintf(fp,"%s= %0.16lf\n",vname,MOVIE_START_RECORDING);
  sprintf(vname,"PNM_BOTHER"); fprintf(fp,"%s= %d\n",vname,PNM_BOTHER);
  sprintf(vname,"PNM_REMOVE"); fprintf(fp,"%s= %d\n",vname,PNM_REMOVE);
  sprintf(vname,"OUTPUT_DUMP_EVERY"); fprintf(fp,"%s= %0.16lf\n",vname,OUTPUT_DUMP_EVERY);
  sprintf(vname,"RTC_BOTHER"); fprintf(fp,"%s= %d\n",vname,RTC_BOTHER);
  sprintf(vname,"GLOBAL_RTC_LENGTH"); fprintf(fp,"%s= %d\n",vname,GLOBAL_RTC_LENGTH);
  sprintf(vname,"GLOBAL_RTC_FRAMELENGTH"); fprintf(fp,"%s= %d\n",vname,GLOBAL_RTC_FRAMELENGTH);
  sprintf(vname,"GLOBAL_RTC_NANGLES"); fprintf(fp,"%s= %d\n",vname,GLOBAL_RTC_NANGLES);
  sprintf(vname,"GLOBAL_RTC_NPHASES"); fprintf(fp,"%s= %d\n",vname,GLOBAL_RTC_NPHASES);
  sprintf(vname,"STROBETRACE_BOTHER"); fprintf(fp,"%s= %d\n",vname,STROBETRACE_BOTHER);
  sprintf(vname,"GLOBAL_STROBETRACE_LENGTH"); fprintf(fp,"%s= %0.16lf\n",vname,GLOBAL_STROBETRACE_LENGTH);
  sprintf(vname,"GLOBAL_STROBETRACE_NANGLES"); fprintf(fp,"%s= %d\n",vname,GLOBAL_STROBETRACE_NANGLES);
  sprintf(vname,"GLOBAL_STROBETRACE_AVALANCHE_UPDATE_EVERY"); fprintf(fp,"%s= %0.16lf\n",vname,GLOBAL_STROBETRACE_AVALANCHE_UPDATE_EVERY);
  sprintf(vname,"GLOBAL_STROBETRACE_AVALANCHE_LOGDMIN"); fprintf(fp,"%s= %d\n",vname,GLOBAL_STROBETRACE_AVALANCHE_LOGDMIN);
  sprintf(vname,"GLOBAL_STROBETRACE_AVALANCHE_LOGDMAX"); fprintf(fp,"%s= %d\n",vname,GLOBAL_STROBETRACE_AVALANCHE_LOGDMAX);
  sprintf(vname,"TUNINGCURVE_BOTHER"); fprintf(fp,"%s= %d\n",vname,TUNINGCURVE_BOTHER);
  sprintf(vname,"GLOBAL_TUNINGCURVE_NANGLES"); fprintf(fp,"%s= %d\n",vname,GLOBAL_TUNINGCURVE_NANGLES);
  sprintf(vname,"GLOBAL_TUNINGCURVE_NRADIUS"); fprintf(fp,"%s= %d\n",vname,GLOBAL_TUNINGCURVE_NRADIUS);
  sprintf(vname,"LMITRI_BOTHER"); fprintf(fp,"%s= %d\n",vname,LMITRI_BOTHER);
  sprintf(vname,"GLOBAL_LMITRI_T0"); fprintf(fp,"%s= %d\n",vname,GLOBAL_LMITRI_T0);
  sprintf(vname,"GLOBAL_LMITRI_TIMELENGTH"); fprintf(fp,"%s= %d\n",vname,GLOBAL_LMITRI_TIMELENGTH);
  sprintf(vname,"GLOBAL_LMITRI_ROW_MAX"); fprintf(fp,"%s= %0.16lf\n",vname,GLOBAL_LMITRI_ROW_MAX);
  sprintf(vname,"GLOBAL_LMITRI_ROW_MIN"); fprintf(fp,"%s= %0.16lf\n",vname,GLOBAL_LMITRI_ROW_MIN);
  sprintf(vname,"PTREE_BOTHER"); fprintf(fp,"%s= %d\n",vname,PTREE_BOTHER);
  sprintf(vname,"GLOBAL_PTREE_BITBYBIT"); fprintf(fp,"%s= %d\n",vname,GLOBAL_PTREE_BITBYBIT);
  sprintf(vname,"GLOBAL_PTREE_ZZZ"); fprintf(fp,"%s= %d\n",vname,GLOBAL_PTREE_ZZZ);
  sprintf(vname,"GLOBAL_PTREE_NLEGS"); fprintf(fp,"%s= %d\n",vname,GLOBAL_PTREE_NLEGS);
  sprintf(vname,"GLOBAL_PTREE_LEGTIME"); fprintf(fp,"%s= %d\n",vname,GLOBAL_PTREE_LEGTIME);
  sprintf(vname,"GLOBAL_PTREE_NREGIONS"); fprintf(fp,"%s= %d\n",vname,GLOBAL_PTREE_NREGIONS);
  sprintf(vname,"P_REWIRE"); fprintf(fp,"%s= %0.16lf\n",vname,P_REWIRE);
  sprintf(vname,"SPARSE_ROW_DIA"); fprintf(fp,"%s= %d\n",vname,SPARSE_ROW_DIA);
  sprintf(vname,"SPARSE_COL_DIA"); fprintf(fp,"%s= %d\n",vname,SPARSE_COL_DIA);
  sprintf(vname,"SPARSE_T2S_ORDERING"); fprintf(fp,"%s= %d\n",vname,SPARSE_T2S_ORDERING);
  sprintf(vname,"SPARSE_CONNECTION_TYPE"); fprintf(fp,"%s= %d\n",vname,SPARSE_CONNECTION_TYPE);
  sprintf(vname,"SPARSE_CS_SCALE"); fprintf(fp,"%s= %0.16lf\n",vname,SPARSE_CS_SCALE);
  sprintf(vname,"GLOBAL_PTREE_EVENT_THRESHOLD"); fprintf(fp,"%s= %d\n",vname,GLOBAL_PTREE_EVENT_THRESHOLD);
  sprintf(vname,"GLOBAL_PTREE_EVENT_WITHIN"); fprintf(fp,"%s= %0.16lf\n",vname,GLOBAL_PTREE_EVENT_WITHIN);
  sprintf(vname,"GLOBAL_PTREE_REGION_TYPE"); fprintf(fp,"%s= %d\n",vname,GLOBAL_PTREE_REGION_TYPE);
  sprintf(vname,"CLOSET_BOTHER"); fprintf(fp,"%s= %d\n",vname,CLOSET_BOTHER);
  sprintf(vname,"YGGDRASIL_BOTHER"); fprintf(fp,"%s= %d\n",vname,YGGDRASIL_BOTHER);
  sprintf(vname,"GLOBAL_YGGDRASIL_PPNREGIONS"); fprintf(fp,"%s= %d\n",vname,GLOBAL_YGGDRASIL_PPNREGIONS);
  sprintf(vname,"GLOBAL_YGGDRASIL_PPNLEGS"); fprintf(fp,"%s= %d\n",vname,GLOBAL_YGGDRASIL_PPNLEGS);
  sprintf(vname,"GLOBAL_YGGDRASIL_PPLEGTIME"); fprintf(fp,"%s= %d\n",vname,GLOBAL_YGGDRASIL_PPLEGTIME);
  sprintf(vname,"GLOBAL_YGGDRASIL_WEIGHT_MINIMUM"); fprintf(fp,"%s= %d\n",vname,GLOBAL_YGGDRASIL_WEIGHT_MINIMUM);
  sprintf(vname,"BONSAI_BOTHER"); fprintf(fp,"%s= %d\n",vname,BONSAI_BOTHER);
  sprintf(vname,"HYDRA_BOTHER"); fprintf(fp,"%s= %d\n",vname,HYDRA_BOTHER);
  sprintf(vname,"GLOBAL_HYDRA_SWITCH_EVERY_TWO"); fprintf(fp,"%s= %d\n",vname,GLOBAL_HYDRA_SWITCH_EVERY_TWO);
  sprintf(vname,"GLOBAL_HYDRA_JUSTONTIME"); fprintf(fp,"%s= %d\n",vname,GLOBAL_HYDRA_JUSTONTIME);
  sprintf(vname,"GLOBAL_HYDRA_STAYONTIME"); fprintf(fp,"%s= %d\n",vname,GLOBAL_HYDRA_STAYONTIME);
  sprintf(vname,"GLOBAL_HYDRA_DUMP_EVERY"); fprintf(fp,"%s= %d\n",vname,GLOBAL_HYDRA_DUMP_EVERY);
  sprintf(vname,"LYAPUNOV_BOTHER"); fprintf(fp,"%s= %d\n",vname,LYAPUNOV_BOTHER);
  sprintf(vname,"GLOBAL_LYAPUNOV_UPDATE_EVERY"); fprintf(fp,"%s= %0.16f\n",vname,GLOBAL_LYAPUNOV_UPDATE_EVERY);
  sprintf(vname,"GLOBAL_LYAPUNOV_JIGGLE"); fprintf(fp,"%s= %0.16f\n",vname,GLOBAL_LYAPUNOV_JIGGLE);
  sprintf(vname,"POWER_BOTHER"); fprintf(fp,"%s= %d\n",vname,POWER_BOTHER);
  sprintf(vname,"GLOBAL_POWER_LENGTH"); fprintf(fp,"%s= %d\n",vname,GLOBAL_POWER_LENGTH);
  sprintf(vname,"GLOBAL_POWER_HOWMANY"); fprintf(fp,"%s= %d\n",vname,GLOBAL_POWER_HOWMANY);
  sprintf(vname,"TAOF_BOTHER"); fprintf(fp,"%s= %d\n",vname,TAOF_BOTHER);
  sprintf(vname,"GLOBAL_TAOF_LENGTH"); fprintf(fp,"%s= %d\n",vname,GLOBAL_TAOF_LENGTH);
  sprintf(vname,"GLOBAL_TAOF_STEP_EVERY"); fprintf(fp,"%s= %d\n",vname,GLOBAL_TAOF_STEP_EVERY);
  sprintf(vname,"SEIDCORR_BOTHER"); fprintf(fp,"%s= %d\n",vname,SEIDCORR_BOTHER);
  sprintf(vname,"GLOBAL_SEIDCORR_SPACE_BIN_SIZE"); fprintf(fp,"%s= %d\n",vname,GLOBAL_SEIDCORR_SPACE_BIN_SIZE);
  sprintf(vname,"GLOBAL_SEIDCORR_TIME_BIN_SIZE"); fprintf(fp,"%s= %d\n",vname,GLOBAL_SEIDCORR_TIME_BIN_SIZE);
  sprintf(vname,"GLOBAL_SEIDCORR_LENGTH"); fprintf(fp,"%s= %d\n",vname,GLOBAL_SEIDCORR_LENGTH);
  sprintf(vname,"GLOBAL_SEIDCORR_TIME_START"); fprintf(fp,"%s= %0.16f\n",vname,GLOBAL_SEIDCORR_TIME_START);
  sprintf(vname,"END"); fprintf(fp,"END= 0");
  if (fp!=stdout){ fclose(fp);}
/*   sprintf(vname,"yy"); fprintf(fp,"%s= %zz\n",vname,yy); */
}
