/* /\* Here are the input functions *\/ */
/* void readinput(); */
/* void updateglobals(char *); */
/* void dumpoutput(); */

void readinput()
{
  /* This reads the piped input file, or standard input.
     the variable names should not be longer than 128 characters */
  int verbose=0;
  char vname[128],equals[128],space[128],semicolon[128];
  do{
    scanf("%[^=]",vname);scanf("%s",equals);scanf("%c",space);updateglobals(vname);scanf("%c",semicolon);
    if (verbose){ printf("At this point variable name is (%s), equals is (%s), semicolon is (%s)\n",vname,equals,semicolon);} 
    scanf("%c",semicolon);}
  while (strcmp(vname,"END")!=0);
}

void updateglobals(char *vname)
{
  /* given a variable name, scan the appropriate format into the appropriate global variable */
  int verbose=GLOBAL_verbose;
  char comma_vs_semicolon[1],**strra=NULL;
  int length=0,nv=0;
  if (strcmp(vname,"GLOBAL_verbose")==0){ scanf("%d",&GLOBAL_verbose); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_verbose);}}
  if (strcmp(vname,"ON_MY_COMPUTER")==0){ scanf("%d",&ON_MY_COMPUTER); if (verbose){ printf("%s read to be %d\n",vname,ON_MY_COMPUTER);}}
  if (strcmp(vname,"GLOBAL_STRING")==0){ scanf("%[^,;]",GLOBAL_STRING); if (verbose){ printf("%s read to be %s\n",vname,GLOBAL_STRING);}}
  if (strcmp(vname,"GLOBAL_RECORD_NUMBER")==0){ scanf("%d",&GLOBAL_RECORD_NUMBER); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_RECORD_NUMBER);}}
  if (strcmp(vname,"GLOBAL_SPIKEINPUT_RSEED")==0){ scanf("%d",&GLOBAL_SPIKEINPUT_RSEED); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_SPIKEINPUT_RSEED);}}
  if (strcmp(vname,"RUN_DONE")==0){ scanf("%d",&RUN_DONE); if (verbose){ printf("%s read to be %d\n",vname,RUN_DONE);}}
  if (strcmp(vname,"GLOBAL_NTYPES")==0){ scanf("%d",&GLOBAL_NTYPES); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_NTYPES);}}
  if (strcmp(vname,"GLOBAL_TYPENAMES")==0){ length=GLOBAL_NTYPES; if (GLOBAL_TYPENAMES!=NULL){ if (verbose){ printf("warning! overwriting old %s\n",vname);}} GLOBAL_TYPENAMES = (char **) tcalloc(length,sizeof(char *));for (nv=0;nv<length;nv++){ GLOBAL_TYPENAMES[nv] = (char *) tcalloc(128,sizeof(char)); scanf("%[^,;]",GLOBAL_TYPENAMES[nv]); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_TYPENAMES[nv]);}}}
  if (strcmp(vname,"GLOBAL_NVARS")==0){ scanf("%d",&GLOBAL_NVARS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_NVARS);}}
  if (strcmp(vname,"GLOBAL_NCLUSTERS")==0){ scanf("%d",&GLOBAL_NCLUSTERS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_NCLUSTERS);}}
  if (strcmp(vname,"GLOBAL_NEURON_MODEL")==0){ scanf("%d",&GLOBAL_NEURON_MODEL); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_NEURON_MODEL);}}
  if (strcmp(vname,"GLOBAL_SNX_VERSION")==0){ scanf("%d",&GLOBAL_SNX_VERSION); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_SNX_VERSION);}}
  if (strcmp(vname,"GLOBAL_SNX_NSTATES_")==0){ length=GLOBAL_NVARS; if (GLOBAL_SNX_NSTATES_!=NULL){ tfree(GLOBAL_SNX_NSTATES_); GLOBAL_SNX_NSTATES_=NULL;} GLOBAL_SNX_NSTATES_ = (int *) tcalloc(length,sizeof(int)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_SNX_NSTATES_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_SNX_NSTATES_[nv]);}}}
  if (strcmp(vname,"GLOBAL_LINK_SPARSE_OR_DENSE")==0){ scanf("%d",&GLOBAL_LINK_SPARSE_OR_DENSE); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_LINK_SPARSE_OR_DENSE);}}
  if (strcmp(vname,"GLOBAL_CLUSTER_MAKE_OR_READ")==0){ scanf("%d",&GLOBAL_CLUSTER_MAKE_OR_READ); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_CLUSTER_MAKE_OR_READ);}}
  if (strcmp(vname,"GLOBAL_CLUSTER_READ_FILE")==0){ scanf("%[^,;]",GLOBAL_CLUSTER_READ_FILE); if (verbose){ printf("%s read to be %s\n",vname,GLOBAL_CLUSTER_READ_FILE);}}
  if (strcmp(vname,"GLOBAL_CLUSTER_PRA_")==0){ length=GLOBAL_NTYPES; if (GLOBAL_CLUSTER_PRA_!=NULL){ tfree(GLOBAL_CLUSTER_PRA_); GLOBAL_CLUSTER_PRA_=NULL;} GLOBAL_CLUSTER_PRA_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(GLOBAL_CLUSTER_PRA_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,GLOBAL_CLUSTER_PRA_[nv]);}}}
  if (strcmp(vname,"GLOBAL_REWEIGHT_STRENGTH")==0){ scanf("%lf",&GLOBAL_REWEIGHT_STRENGTH); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_REWEIGHT_STRENGTH);}}
  if (strcmp(vname,"GLOBAL_REWEIGHT_MULTIPLY_OR_ADD")==0){ scanf("%d",&GLOBAL_REWEIGHT_MULTIPLY_OR_ADD); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_REWEIGHT_MULTIPLY_OR_ADD);}}
  if (strcmp(vname,"GLOBAL_REWEIGHT_READ_FILE")==0){ scanf("%[^,;]",GLOBAL_REWEIGHT_READ_FILE); if (verbose){ printf("%s read to be %s\n",vname,GLOBAL_REWEIGHT_READ_FILE);}}
  if (strcmp(vname,"GLOBAL_VARNAMES")==0){ length=GLOBAL_NVARS; if (GLOBAL_VARNAMES!=NULL){ if (verbose){ printf("warning! overwriting old %s\n",vname);}} GLOBAL_VARNAMES = (char **) tcalloc(length,sizeof(char *));for (nv=0;nv<length;nv++){ GLOBAL_VARNAMES[nv] = (char *) tcalloc(128,sizeof(char)); scanf("%[^,;]",GLOBAL_VARNAMES[nv]); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_VARNAMES[nv]);}}}
  if (strcmp(vname,"GLOBAL_INDEXING_sra_LENGTH")==0){ scanf("%d",&GLOBAL_INDEXING_sra_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_INDEXING_sra_LENGTH);}}
  if (strcmp(vname,"GLOBAL_INDEXING_REFILE_sra")==0){ length=GLOBAL_NVARS; if (GLOBAL_INDEXING_REFILE_sra!=NULL){ tfree(GLOBAL_INDEXING_REFILE_sra); GLOBAL_INDEXING_REFILE_sra=NULL;} GLOBAL_INDEXING_REFILE_sra = (int *) tcalloc(length,sizeof(int)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_INDEXING_REFILE_sra[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_INDEXING_REFILE_sra[nv]);}}}
  if (strcmp(vname,"GLOBAL_INDEXING_CHECKOUT_sra")==0){ length=GLOBAL_INDEXING_sra_LENGTH; if (GLOBAL_INDEXING_CHECKOUT_sra!=NULL){ tfree(GLOBAL_INDEXING_CHECKOUT_sra); GLOBAL_INDEXING_CHECKOUT_sra=NULL;} GLOBAL_INDEXING_CHECKOUT_sra = (int *) tcalloc(length,sizeof(int)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_INDEXING_CHECKOUT_sra[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_INDEXING_CHECKOUT_sra[nv]);}}}
  if (strcmp(vname,"GLOBAL_LENGTHRA")==0){ length=GLOBAL_NTYPES; if (GLOBAL_LENGTHRA!=NULL){ tfree(GLOBAL_LENGTHRA); GLOBAL_LENGTHRA=NULL;} GLOBAL_LENGTHRA = (int *) tcalloc(length,sizeof(int)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_LENGTHRA[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_LENGTHRA[nv]);}}}
  if (strcmp(vname,"VOLTAGE_")==0){ length=GLOBAL_NVARS; if (VOLTAGE_!=NULL){ tfree(VOLTAGE_); VOLTAGE_=NULL;} VOLTAGE_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(VOLTAGE_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,VOLTAGE_[nv]);}}}
  if (strcmp(vname,"GLOBAL_VESICLE_DEPLETION")==0){ scanf("%lf",&GLOBAL_VESICLE_DEPLETION); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_VESICLE_DEPLETION);}}
  if (strcmp(vname,"VOLTAGE_THRESHOLD_S")==0){ scanf("%lf",&VOLTAGE_THRESHOLD_S); if (verbose){ printf("%s read to be %lf\n",vname,VOLTAGE_THRESHOLD_S);}}
  if (strcmp(vname,"VOLTAGE_THRESHOLD_D")==0){ scanf("%lf",&VOLTAGE_THRESHOLD_D); if (verbose){ printf("%s read to be %lf\n",vname,VOLTAGE_THRESHOLD_D);}}
  if (strcmp(vname,"TAU_REF")==0){ scanf("%lf",&TAU_REF); if (verbose){ printf("%s read to be %lf\n",vname,TAU_REF);}}
  if (strcmp(vname,"TAU_")==0){ length=GLOBAL_NVARS; if (TAU_!=NULL){ tfree(TAU_); TAU_=NULL;} TAU_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(TAU_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,TAU_[nv]);}}}
  if (strcmp(vname,"CONDUCTANCE_")==0){ length=GLOBAL_NVARS; if (CONDUCTANCE_!=NULL){ tfree(CONDUCTANCE_); CONDUCTANCE_=NULL;} CONDUCTANCE_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(CONDUCTANCE_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,CONDUCTANCE_[nv]);}}}
  if (strcmp(vname,"CONDUCTANCE_DS")==0){ scanf("%lf",&CONDUCTANCE_DS); if (verbose){ printf("%s read to be %lf\n",vname,CONDUCTANCE_DS);}}
  if (strcmp(vname,"CONDUCTANCE_SD")==0){ scanf("%lf",&CONDUCTANCE_SD); if (verbose){ printf("%s read to be %lf\n",vname,CONDUCTANCE_SD);}}
  if (strcmp(vname,"CURRENT_INJECTION_S")==0){ scanf("%lf",&CURRENT_INJECTION_S); if (verbose){ printf("%s read to be %lf\n",vname,CURRENT_INJECTION_S);}}
  if (strcmp(vname,"CURRENT_INJECTION_D")==0){ scanf("%lf",&CURRENT_INJECTION_D); if (verbose){ printf("%s read to be %lf\n",vname,CURRENT_INJECTION_D);}}
  if (strcmp(vname,"AUTAPSES_OFF")==0){ length=GLOBAL_NTYPES; if (AUTAPSES_OFF!=NULL){ tfree(AUTAPSES_OFF); AUTAPSES_OFF=NULL;} AUTAPSES_OFF = (int *) tcalloc(length,sizeof(int)); for (nv=0;nv<length;nv++){ scanf("%d",&(AUTAPSES_OFF[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,AUTAPSES_OFF[nv]);}}}
  if (strcmp(vname,"GLOBAL_CS_ORN_SCALE")==0){ scanf("%lf",&GLOBAL_CS_ORN_SCALE); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_CS_ORN_SCALE);}}
  if (strcmp(vname,"CS_ORN_")==0){ length=GLOBAL_NTYPES; if (CS_ORN_!=NULL){ tfree(CS_ORN_); CS_ORN_=NULL;} CS_ORN_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(CS_ORN_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,CS_ORN_[nv]);}}}
  if (strcmp(vname,"CS_ORN_mainak_stim_")==0){ length=GLOBAL_NTYPES; if (CS_ORN_mainak_stim_!=NULL){ tfree(CS_ORN_mainak_stim_); CS_ORN_mainak_stim_=NULL;} CS_ORN_mainak_stim_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(CS_ORN_mainak_stim_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,CS_ORN_mainak_stim_[nv]);}}}
  if (strcmp(vname,"CS_ORN_wilson_stim_")==0){ length=GLOBAL_NTYPES; if (CS_ORN_wilson_stim_!=NULL){ tfree(CS_ORN_wilson_stim_); CS_ORN_wilson_stim_=NULL;} CS_ORN_wilson_stim_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(CS_ORN_wilson_stim_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,CS_ORN_wilson_stim_[nv]);}}}
  if (strcmp(vname,"GLOBAL_CS_SCALE")==0){ scanf("%lf",&GLOBAL_CS_SCALE); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_CS_SCALE);}}
  if (strcmp(vname,"GLOBAL_CS_PRETYPE_SCALE_")==0){ length=GLOBAL_NTYPES; if (GLOBAL_CS_PRETYPE_SCALE_!=NULL){ tfree(GLOBAL_CS_PRETYPE_SCALE_); GLOBAL_CS_PRETYPE_SCALE_=NULL;} GLOBAL_CS_PRETYPE_SCALE_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(GLOBAL_CS_PRETYPE_SCALE_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,GLOBAL_CS_PRETYPE_SCALE_[nv]);}}}
  if (strcmp(vname,"GLOBAL_CS_POSTYPE_SCALE_")==0){ length=GLOBAL_NTYPES; if (GLOBAL_CS_POSTYPE_SCALE_!=NULL){ tfree(GLOBAL_CS_POSTYPE_SCALE_); GLOBAL_CS_POSTYPE_SCALE_=NULL;} GLOBAL_CS_POSTYPE_SCALE_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(GLOBAL_CS_POSTYPE_SCALE_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,GLOBAL_CS_POSTYPE_SCALE_[nv]);}}}
  if (strcmp(vname,"GLOBAL_CS_SRA_SCALE_")==0){ length=GLOBAL_INDEXING_sra_LENGTH; if (GLOBAL_CS_SRA_SCALE_!=NULL){ tfree(GLOBAL_CS_SRA_SCALE_); GLOBAL_CS_SRA_SCALE_=NULL;} GLOBAL_CS_SRA_SCALE_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(GLOBAL_CS_SRA_SCALE_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,GLOBAL_CS_SRA_SCALE_[nv]);}}}
  if (strcmp(vname,"CS__")==0){ length=GLOBAL_NTYPES*GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH; if (CS__!=NULL){ tfree(CS__); CS__=NULL;} CS__ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(CS__[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,CS__[nv]);}}}
  if (strcmp(vname,"SPARSE__")==0){ length=GLOBAL_NTYPES*GLOBAL_NTYPES; if (SPARSE__!=NULL){ tfree(SPARSE__); SPARSE__=NULL;} SPARSE__ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(SPARSE__[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,SPARSE__[nv]);}}}
  if (strcmp(vname,"SPARSEHIT__")==0){ length=GLOBAL_NTYPES*GLOBAL_NTYPES; if (SPARSEHIT__!=NULL){ tfree(SPARSEHIT__); SPARSEHIT__=NULL;} SPARSEHIT__ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(SPARSEHIT__[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,SPARSEHIT__[nv]);}}}
  if (strcmp(vname,"SPARSE_OGLI__")==0){ length=GLOBAL_NTYPES*GLOBAL_NTYPES; if (SPARSE_OGLI__!=NULL){ tfree(SPARSE_OGLI__); SPARSE_OGLI__=NULL;} SPARSE_OGLI__ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(SPARSE_OGLI__[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,SPARSE_OGLI__[nv]);}}}
  if (strcmp(vname,"SPARSEHIT_OGLI__")==0){ length=GLOBAL_NTYPES*GLOBAL_NTYPES; if (SPARSEHIT_OGLI__!=NULL){ tfree(SPARSEHIT_OGLI__); SPARSEHIT_OGLI__=NULL;} SPARSEHIT_OGLI__ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(SPARSEHIT_OGLI__[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,SPARSEHIT_OGLI__[nv]);}}}
  if (strcmp(vname,"GLOBAL_P_FAIL_SCALE")==0){ scanf("%lf",&GLOBAL_P_FAIL_SCALE); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_P_FAIL_SCALE);}}
  if (strcmp(vname,"P_FAIL__")==0){ length=GLOBAL_NTYPES*GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH; if (P_FAIL__!=NULL){ tfree(P_FAIL__); P_FAIL__=NULL;} P_FAIL__ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(P_FAIL__[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,P_FAIL__[nv]);}}}
  if (strcmp(vname,"GLOBAL_ODOR_BASE")==0){ scanf("%d",&GLOBAL_ODOR_BASE); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_ODOR_BASE);}}
  if (strcmp(vname,"GLOBAL_ODORra_START")==0){ length=GLOBAL_NTYPES; if (GLOBAL_ODORra_START!=NULL){ odortfree(GLOBAL_ODORra_START);} GLOBAL_ODORra_START=NULL; strra=(char **) tcalloc(length,sizeof(char *)); for (nv=0;nv<length;nv++){ strra[nv]=(char *)tcalloc(GLOBAL_LENGTHRA[nv]+2,sizeof(char)); scanf("%[^,;]",strra[nv]); if (nv<length-1){ scanf("%c",comma_vs_semicolon);}} GLOBAL_ODORra_START=odormake(length,GLOBAL_LENGTHRA,GLOBAL_ODOR_BASE,strra,NULL); for (nv=0;nv<length;nv++){ tfree(strra[nv]);strra[nv]=NULL;} tfree(strra);strra=NULL; if (verbose){ printf("%s read to be ",vname); odorfprintf(stdout,GLOBAL_ODORra_START);}}
  if (strcmp(vname,"INPUT_CONTRAST_START")==0){ scanf("%lf",&INPUT_CONTRAST_START); if (verbose){ printf("%s read to be %lf\n",vname,INPUT_CONTRAST_START);}}
  if (strcmp(vname,"INPUT_PULSE")==0){ scanf("%lf",&INPUT_PULSE); if (verbose){ printf("%s read to be %lf\n",vname,INPUT_PULSE);}}
  if (strcmp(vname,"ORN_BOTHER")==0){ scanf("%d",&ORN_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,ORN_BOTHER);}}
  if (strcmp(vname,"GLOBAL_ORN_TIMES_")==0){ length=4; if (GLOBAL_ORN_TIMES_!=NULL){ tfree(GLOBAL_ORN_TIMES_); GLOBAL_ORN_TIMES_=NULL;} GLOBAL_ORN_TIMES_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(GLOBAL_ORN_TIMES_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,GLOBAL_ORN_TIMES_[nv]);}}}
  if (strcmp(vname,"ORN_BACKRATE")==0){ scanf("%lf",&ORN_BACKRATE); if (verbose){ printf("%s read to be %lf\n",vname,ORN_BACKRATE);}}
  if (strcmp(vname,"SPIKETOL")==0){ scanf("%d",&SPIKETOL); if (verbose){ printf("%s read to be %d\n",vname,SPIKETOL);}}
  if (strcmp(vname,"GLOBAL_SPACE_SMOOTHER")==0){ scanf("%d",&GLOBAL_SPACE_SMOOTHER); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_SPACE_SMOOTHER);}}
  if (strcmp(vname,"GLOBAL_TI")==0){ scanf("%lf",&GLOBAL_TI); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_TI);}}
  if (strcmp(vname,"GLOBAL_TF")==0){ scanf("%lf",&GLOBAL_TF); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_TF);}}
  if (strcmp(vname,"GLOBAL_DTmax")==0){ scanf("%lf",&GLOBAL_DTmax); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_DTmax);}}
  if (strcmp(vname,"GLOBAL_dtadapt")==0){ scanf("%d",&GLOBAL_dtadapt); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_dtadapt);}}
  if (strcmp(vname,"FIDDLE_ROW")==0){ scanf("%d",&FIDDLE_ROW); if (verbose){ printf("%s read to be %d\n",vname,FIDDLE_ROW);}}
  if (strcmp(vname,"FIDDLE_COL")==0){ scanf("%d",&FIDDLE_COL); if (verbose){ printf("%s read to be %d\n",vname,FIDDLE_COL);}}
  if (strcmp(vname,"FIDDLE_PARAMETER")==0){ scanf("%d",&FIDDLE_PARAMETER); if (verbose){ printf("%s read to be %d\n",vname,FIDDLE_PARAMETER);}}
  if (strcmp(vname,"DRAW_FLAG")==0){ scanf("%d",&DRAW_FLAG); if (verbose){ printf("%s read to be %d\n",vname,DRAW_FLAG);}}
  if (strcmp(vname,"DRAW_FLAG2")==0){ scanf("%d",&DRAW_FLAG2); if (verbose){ printf("%s read to be %d\n",vname,DRAW_FLAG2);}}
  if (strcmp(vname,"STEPS_PER_DRAW")==0){ scanf("%d",&STEPS_PER_DRAW); if (verbose){ printf("%s read to be %d\n",vname,STEPS_PER_DRAW);}}
  if (strcmp(vname,"STD_VIEW")==0){ scanf("%lf",&STD_VIEW); if (verbose){ printf("%s read to be %lf\n",vname,STD_VIEW);}}
  if (strcmp(vname,"GRAYSCALE")==0){ scanf("%d",&GRAYSCALE); if (verbose){ printf("%s read to be %d\n",vname,GRAYSCALE);}}
  if (strcmp(vname,"SUPERGLOBAL_DRAW_FLAG")==0){ scanf("%d",&SUPERGLOBAL_DRAW_FLAG); if (verbose){ printf("%s read to be %d\n",vname,SUPERGLOBAL_DRAW_FLAG);}}
  if (strcmp(vname,"OUTPUT_DUMP_EVERY")==0){ scanf("%lf",&OUTPUT_DUMP_EVERY); if (verbose){ printf("%s read to be %lf\n",vname,OUTPUT_DUMP_EVERY);}}
  if (strcmp(vname,"POWER_BOTHER")==0){ scanf("%d",&POWER_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,POWER_BOTHER);}}
  if (strcmp(vname,"GLOBAL_POWER_LENGTH")==0){ scanf("%d",&GLOBAL_POWER_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_POWER_LENGTH);}}
  if (strcmp(vname,"GLOBAL_POWER_INDEXING_NTYPE_LENGTH")==0){ scanf("%d",&GLOBAL_POWER_INDEXING_NTYPE_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_POWER_INDEXING_NTYPE_LENGTH);}}
  if (strcmp(vname,"GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT")==0){ length=GLOBAL_POWER_INDEXING_NTYPE_LENGTH; if (GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT!=NULL){ tfree(GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT); GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT=NULL;} GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT = (int *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT[nv]);}}}
  if (strcmp(vname,"GLOBAL_POWER_INDEXING_NTYPE_REFILE")==0){ length=GLOBAL_NTYPES; if (GLOBAL_POWER_INDEXING_NTYPE_REFILE!=NULL){ tfree(GLOBAL_POWER_INDEXING_NTYPE_REFILE); GLOBAL_POWER_INDEXING_NTYPE_REFILE=NULL;} GLOBAL_POWER_INDEXING_NTYPE_REFILE = (int *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_POWER_INDEXING_NTYPE_REFILE[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_POWER_INDEXING_NTYPE_REFILE[nv]);}}}
  if (strcmp(vname,"GLOBAL_POWER_INDEXING_NVAR_LENGTH")==0){ scanf("%d",&GLOBAL_POWER_INDEXING_NVAR_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_POWER_INDEXING_NVAR_LENGTH);}}
  if (strcmp(vname,"GLOBAL_POWER_INDEXING_NVAR_CHECKOUT")==0){ length=GLOBAL_POWER_INDEXING_NVAR_LENGTH; if (GLOBAL_POWER_INDEXING_NVAR_CHECKOUT!=NULL){ tfree(GLOBAL_POWER_INDEXING_NVAR_CHECKOUT); GLOBAL_POWER_INDEXING_NVAR_CHECKOUT=NULL;} GLOBAL_POWER_INDEXING_NVAR_CHECKOUT = (int *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_POWER_INDEXING_NVAR_CHECKOUT[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_POWER_INDEXING_NVAR_CHECKOUT[nv]);}}}
  if (strcmp(vname,"GLOBAL_POWER_INDEXING_NVAR_REFILE")==0){ length=GLOBAL_NVARS; if (GLOBAL_POWER_INDEXING_NVAR_REFILE!=NULL){ tfree(GLOBAL_POWER_INDEXING_NVAR_REFILE); GLOBAL_POWER_INDEXING_NVAR_REFILE=NULL;} GLOBAL_POWER_INDEXING_NVAR_REFILE = (int *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_POWER_INDEXING_NVAR_REFILE[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_POWER_INDEXING_NVAR_REFILE[nv]);}}}
  if (strcmp(vname,"GLOBAL_POWER_TRAJECTORY_LENGTH")==0){ scanf("%d",&GLOBAL_POWER_TRAJECTORY_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_POWER_TRAJECTORY_LENGTH);}}
  if (strcmp(vname,"GLOBAL_POWER_WINDOW_LENGTH")==0){ scanf("%d",&GLOBAL_POWER_WINDOW_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_POWER_WINDOW_LENGTH);}}
  if (strcmp(vname,"GLOBAL_POWER_WINDOW_UPDATE_EVERY")==0){ scanf("%d",&GLOBAL_POWER_WINDOW_UPDATE_EVERY); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_POWER_WINDOW_UPDATE_EVERY);}}
  if (strcmp(vname,"GLOBAL_POWER_maxra_")==0){ length=GLOBAL_NVARS; if (GLOBAL_POWER_maxra_!=NULL){ tfree(GLOBAL_POWER_maxra_); GLOBAL_POWER_maxra_=NULL;} GLOBAL_POWER_maxra_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(GLOBAL_POWER_maxra_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,GLOBAL_POWER_maxra_[nv]);}}}
  if (strcmp(vname,"GLOBAL_POWER_minra_")==0){ length=GLOBAL_NVARS; if (GLOBAL_POWER_minra_!=NULL){ tfree(GLOBAL_POWER_minra_); GLOBAL_POWER_minra_=NULL;} GLOBAL_POWER_minra_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(GLOBAL_POWER_minra_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,GLOBAL_POWER_minra_[nv]);}}}
  if (strcmp(vname,"GLOBAL_POWER_CYCLE_BOTHER")==0){ scanf("%d",&GLOBAL_POWER_CYCLE_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_POWER_CYCLE_BOTHER);}}
  if (strcmp(vname,"GLOBAL_POWER_CORRELATION_BOTHER")==0){ scanf("%d",&GLOBAL_POWER_CORRELATION_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_POWER_CORRELATION_BOTHER);}}
  if (strcmp(vname,"RHO_BOTHER")==0){ scanf("%d",&RHO_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,RHO_BOTHER);}}
  if (strcmp(vname,"GLOBAL_RHO_LENGTH")==0){ scanf("%d",&GLOBAL_RHO_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_RHO_LENGTH);}}
  if (strcmp(vname,"GLOBAL_RHO_INDEXING_NVAR_LENGTH")==0){ scanf("%d",&GLOBAL_RHO_INDEXING_NVAR_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_RHO_INDEXING_NVAR_LENGTH);}}
  if (strcmp(vname,"GLOBAL_RHO_INDEXING_NVAR_CHECKOUT")==0){ length=GLOBAL_RHO_INDEXING_NVAR_LENGTH; if (GLOBAL_RHO_INDEXING_NVAR_CHECKOUT!=NULL){ tfree(GLOBAL_RHO_INDEXING_NVAR_CHECKOUT); GLOBAL_RHO_INDEXING_NVAR_CHECKOUT=NULL;} GLOBAL_RHO_INDEXING_NVAR_CHECKOUT = (int *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_RHO_INDEXING_NVAR_CHECKOUT[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_RHO_INDEXING_NVAR_CHECKOUT[nv]);}}}
  if (strcmp(vname,"GLOBAL_RHO_INDEXING_NVAR_REFILE")==0){ length=GLOBAL_NVARS; if (GLOBAL_RHO_INDEXING_NVAR_REFILE!=NULL){ tfree(GLOBAL_RHO_INDEXING_NVAR_REFILE); GLOBAL_RHO_INDEXING_NVAR_REFILE=NULL;} GLOBAL_RHO_INDEXING_NVAR_REFILE = (int *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_RHO_INDEXING_NVAR_REFILE[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_RHO_INDEXING_NVAR_REFILE[nv]);}}}
  if (strcmp(vname,"GLOBAL_RHO_NBINRA_")==0){ length=GLOBAL_NVARS; if (GLOBAL_RHO_NBINRA_!=NULL){ tfree(GLOBAL_RHO_NBINRA_); GLOBAL_RHO_NBINRA_=NULL;} GLOBAL_RHO_NBINRA_ = (int *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_RHO_NBINRA_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_RHO_NBINRA_[nv]);}}}
  if (strcmp(vname,"PTREE_BOTHER")==0){ scanf("%d",&PTREE_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,PTREE_BOTHER);}}
  if (strcmp(vname,"GLOBAL_PTREE_BITBYBIT")==0){ scanf("%d",&GLOBAL_PTREE_BITBYBIT); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_BITBYBIT);}}
  if (strcmp(vname,"GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER")==0){ scanf("%d",&GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER);}}
  if (strcmp(vname,"GLOBAL_PTREE_ZZZ")==0){ scanf("%d",&GLOBAL_PTREE_ZZZ); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_ZZZ);}}
  if (strcmp(vname,"GLOBAL_PTREE_NREGIONS")==0){ scanf("%d",&GLOBAL_PTREE_NREGIONS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_NREGIONS);}}
  if (strcmp(vname,"GLOBAL_PTREE_NLEGS")==0){ scanf("%d",&GLOBAL_PTREE_NLEGS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_NLEGS);}}
  if (strcmp(vname,"GLOBAL_PTREE_LEGTIME")==0){ scanf("%d",&GLOBAL_PTREE_LEGTIME); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_LEGTIME);}}
  if (strcmp(vname,"GLOBAL_PTREE_EVENT_THRESHOLD")==0){ scanf("%d",&GLOBAL_PTREE_EVENT_THRESHOLD); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_EVENT_THRESHOLD);}}
  if (strcmp(vname,"GLOBAL_PTREE_EVENT_WITHIN")==0){ scanf("%lf",&GLOBAL_PTREE_EVENT_WITHIN); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_PTREE_EVENT_WITHIN);}}
  if (strcmp(vname,"GLOBAL_PTREE_REGION_TYPE")==0){ scanf("%d",&GLOBAL_PTREE_REGION_TYPE); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_REGION_TYPE);}}
  if (strcmp(vname,"GLOBAL_PTREE_DUMP_TYPE")==0){ scanf("%d",&GLOBAL_PTREE_DUMP_TYPE); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_PTREE_DUMP_TYPE);}}
  if (strcmp(vname,"CLUSTERDATA_BOTHER")==0){ scanf("%d",&CLUSTERDATA_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,CLUSTERDATA_BOTHER);}}
  if (strcmp(vname,"SUITE_BOTHER")==0){ scanf("%d",&SUITE_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,SUITE_BOTHER);}}
  if (strcmp(vname,"SUITE_NODORS")==0){ scanf("%d",&SUITE_NODORS); if (verbose){ printf("%s read to be %d\n",vname,SUITE_NODORS);}}
  if (strcmp(vname,"SUITE_ODORra_BACON")==0){ length=GLOBAL_NTYPES; if (SUITE_ODORra_BACON!=NULL){ odortfree(SUITE_ODORra_BACON);} SUITE_ODORra_BACON=NULL; strra=(char **) tcalloc(length,sizeof(char *)); for (nv=0;nv<length;nv++){ strra[nv]=(char *)tcalloc(GLOBAL_LENGTHRA[nv]+2,sizeof(char)); scanf("%[^,;]",strra[nv]); if (nv<length-1){ scanf("%c",comma_vs_semicolon);}} SUITE_ODORra_BACON=odormake(length,GLOBAL_LENGTHRA,GLOBAL_ODOR_BASE,strra,NULL); for (nv=0;nv<length;nv++){ tfree(strra[nv]);strra[nv]=NULL;} tfree(strra);strra=NULL; if (verbose){ printf("%s read to be ",vname); odorfprintf(stdout,SUITE_ODORra_BACON);}}
  if (strcmp(vname,"SUITE_CS_ORN_SCALE")==0){ scanf("%lf",&SUITE_CS_ORN_SCALE); if (verbose){ printf("%s read to be %lf\n",vname,SUITE_CS_ORN_SCALE);}}
  if (strcmp(vname,"SUITE_INPUT_CONTRAST")==0){ scanf("%lf",&SUITE_INPUT_CONTRAST); if (verbose){ printf("%s read to be %lf\n",vname,SUITE_INPUT_CONTRAST);}}
  if (strcmp(vname,"SUITE_NCONCENTRATIONS")==0){ scanf("%d",&SUITE_NCONCENTRATIONS); if (verbose){ printf("%s read to be %d\n",vname,SUITE_NCONCENTRATIONS);}}
  if (strcmp(vname,"SUITE_NBICUCULLINES")==0){ scanf("%d",&SUITE_NBICUCULLINES); if (verbose){ printf("%s read to be %d\n",vname,SUITE_NBICUCULLINES);}}
  if (strcmp(vname,"SUITE_NINSTANCES")==0){ scanf("%d",&SUITE_NINSTANCES); if (verbose){ printf("%s read to be %d\n",vname,SUITE_NINSTANCES);}}
  if (strcmp(vname,"SUITE_POWER_BOTHER")==0){ scanf("%d",&SUITE_POWER_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,SUITE_POWER_BOTHER);}}
  if (strcmp(vname,"SUITE_PTREE_BOTHER")==0){ scanf("%d",&SUITE_PTREE_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,SUITE_PTREE_BOTHER);}}
  if (strcmp(vname,"SUITE_CLUSTERDATA_BOTHER")==0){ scanf("%d",&SUITE_CLUSTERDATA_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,SUITE_CLUSTERDATA_BOTHER);}}
  if (strcmp(vname,"SUITE_4_REWEIGHT_STRENGTH")==0){ scanf("%lf",&SUITE_4_REWEIGHT_STRENGTH); if (verbose){ printf("%s read to be %lf\n",vname,SUITE_4_REWEIGHT_STRENGTH);}}
  if (strcmp(vname,"SUITE_7_CLEANUP")==0){ scanf("%d",&SUITE_7_CLEANUP); if (verbose){ printf("%s read to be %d\n",vname,SUITE_7_CLEANUP);}}
  if (strcmp(vname,"SUITE_8_CLEANUP")==0){ scanf("%d",&SUITE_8_CLEANUP); if (verbose){ printf("%s read to be %d\n",vname,SUITE_8_CLEANUP);}}
  if (strcmp(vname,"SUITE_8_wilson_axon_stim")==0){ scanf("%lf",&SUITE_8_wilson_axon_stim); if (verbose){ printf("%s read to be %lf\n",vname,SUITE_8_wilson_axon_stim);}}
  if (strcmp(vname,"SUITE_NSECONDS")==0){ scanf("%d",&SUITE_NSECONDS); if (verbose){ printf("%s read to be %d\n",vname,SUITE_NSECONDS);}}
  if (strcmp(vname,"SUITE_DUMPEVERY")==0){ scanf("%d",&SUITE_DUMPEVERY); if (verbose){ printf("%s read to be %d\n",vname,SUITE_DUMPEVERY);}}
  if (strcmp(vname,"SUITE_BITBYBIT_RECORD")==0){ scanf("%d",&SUITE_BITBYBIT_RECORD); if (verbose){ printf("%s read to be %d\n",vname,SUITE_BITBYBIT_RECORD);}}
  if (strcmp(vname,"SUITE_BITBYBIT_REMOVE")==0){ scanf("%d",&SUITE_BITBYBIT_REMOVE); if (verbose){ printf("%s read to be %d\n",vname,SUITE_BITBYBIT_REMOVE);}}
  if (strcmp(vname,"SUITE_SINDEXMAX")==0){ scanf("%d",&SUITE_SINDEXMAX); if (verbose){ printf("%s read to be %d\n",vname,SUITE_SINDEXMAX);}}
  if (strcmp(vname,"SUITE_DINDEXMAX")==0){ scanf("%d",&SUITE_DINDEXMAX); if (verbose){ printf("%s read to be %d\n",vname,SUITE_DINDEXMAX);}}
  if (strcmp(vname,"SUITE_TINDEXMAX")==0){ scanf("%d",&SUITE_TINDEXMAX); if (verbose){ printf("%s read to be %d\n",vname,SUITE_TINDEXMAX);}}
  if (strcmp(vname,"LYAPUNOV_BOTHER")==0){ scanf("%d",&LYAPUNOV_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,LYAPUNOV_BOTHER);}}
  if (strcmp(vname,"CAICOR_BOTHER")==0){ scanf("%d",&CAICOR_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,CAICOR_BOTHER);}}
  if (strcmp(vname,"GLOBAL_CAICOR_NBINS")==0){ scanf("%d",&GLOBAL_CAICOR_NBINS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_CAICOR_NBINS);}}
  if (strcmp(vname,"HHLIB_BOTHER")==0){ scanf("%d",&HHLIB_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,HHLIB_BOTHER);}}
  if (strcmp(vname,"GLOBAL_HHLIB_MAKE_OR_READ")==0){ scanf("%d",&GLOBAL_HHLIB_MAKE_OR_READ); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_HHLIB_MAKE_OR_READ);}}
  if (strcmp(vname,"GLOBAL_HHLIB_READ_FILE")==0){ scanf("%[^,;]",GLOBAL_HHLIB_READ_FILE); if (verbose){ printf("%s read to be %s\n",vname,GLOBAL_HHLIB_READ_FILE);}}
  if (strcmp(vname,"GLOBAL_HHLIB_NBINS")==0){ scanf("%d",&GLOBAL_HHLIB_NBINS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_HHLIB_NBINS);}}
  if (strcmp(vname,"GLOBAL_HHLIB_USE_AFTER")==0){ scanf("%d",&GLOBAL_HHLIB_USE_AFTER); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_HHLIB_USE_AFTER);}}
  if (strcmp(vname,"GLOBAL_HHLIB_TAU_SKIP")==0){ scanf("%lf",&GLOBAL_HHLIB_TAU_SKIP); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_HHLIB_TAU_SKIP);}}
  if (strcmp(vname,"GLOBAL_HHLIB_INDEXING_TRIGGER")==0){ scanf("%d",&GLOBAL_HHLIB_INDEXING_TRIGGER); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_HHLIB_INDEXING_TRIGGER);}}
  if (strcmp(vname,"GLOBAL_HHLIB_LOGFLAGRA")==0){ length=GLOBAL_NVARS; if (GLOBAL_HHLIB_LOGFLAGRA!=NULL){ tfree(GLOBAL_HHLIB_LOGFLAGRA); GLOBAL_HHLIB_LOGFLAGRA=NULL;} GLOBAL_HHLIB_LOGFLAGRA = (int *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_HHLIB_LOGFLAGRA[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_HHLIB_LOGFLAGRA[nv]);}}}
  if (strcmp(vname,"GLOBAL_HHLIB_LIBUSEFLAGRA")==0){ length=GLOBAL_NVARS; if (GLOBAL_HHLIB_LIBUSEFLAGRA!=NULL){ tfree(GLOBAL_HHLIB_LIBUSEFLAGRA); GLOBAL_HHLIB_LIBUSEFLAGRA=NULL;} GLOBAL_HHLIB_LIBUSEFLAGRA = (int *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_HHLIB_LIBUSEFLAGRA[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_HHLIB_LIBUSEFLAGRA[nv]);}}}
  if (strcmp(vname,"GLOBAL_HHLIB_INDEXING_NVAR_LENGTH")==0){ scanf("%d",&GLOBAL_HHLIB_INDEXING_NVAR_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_HHLIB_INDEXING_NVAR_LENGTH);}}
  if (strcmp(vname,"GLOBAL_HHLIB_INDEXING_NVAR_CHECKOUT")==0){ length=GLOBAL_HHLIB_INDEXING_NVAR_LENGTH; if (GLOBAL_HHLIB_INDEXING_NVAR_CHECKOUT!=NULL){ tfree(GLOBAL_HHLIB_INDEXING_NVAR_CHECKOUT); GLOBAL_HHLIB_INDEXING_NVAR_CHECKOUT=NULL;} GLOBAL_HHLIB_INDEXING_NVAR_CHECKOUT = (int *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_HHLIB_INDEXING_NVAR_CHECKOUT[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_HHLIB_INDEXING_NVAR_CHECKOUT[nv]);}}}
  if (strcmp(vname,"GLOBAL_HHLIB_INDEXING_NVAR_REFILE")==0){ length=GLOBAL_NVARS; if (GLOBAL_HHLIB_INDEXING_NVAR_REFILE!=NULL){ tfree(GLOBAL_HHLIB_INDEXING_NVAR_REFILE); GLOBAL_HHLIB_INDEXING_NVAR_REFILE=NULL;} GLOBAL_HHLIB_INDEXING_NVAR_REFILE = (int *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_HHLIB_INDEXING_NVAR_REFILE[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_HHLIB_INDEXING_NVAR_REFILE[nv]);}}}
  if (strcmp(vname,"GLOBAL_HHLIB_maxra_")==0){ length=GLOBAL_NVARS; if (GLOBAL_HHLIB_maxra_!=NULL){ tfree(GLOBAL_HHLIB_maxra_); GLOBAL_HHLIB_maxra_=NULL;} GLOBAL_HHLIB_maxra_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(GLOBAL_HHLIB_maxra_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,GLOBAL_HHLIB_maxra_[nv]);}}}
  if (strcmp(vname,"GLOBAL_HHLIB_minra_")==0){ length=GLOBAL_NVARS; if (GLOBAL_HHLIB_minra_!=NULL){ tfree(GLOBAL_HHLIB_minra_); GLOBAL_HHLIB_minra_=NULL;} GLOBAL_HHLIB_minra_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(GLOBAL_HHLIB_minra_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,GLOBAL_HHLIB_minra_[nv]);}}}
  if (strcmp(vname,"SNXDATA_BOTHER")==0){ scanf("%d",&SNXDATA_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,SNXDATA_BOTHER);}}
  if (strcmp(vname,"GLOBAL_SNXDATA_LOOKBACK")==0){ scanf("%d",&GLOBAL_SNXDATA_LOOKBACK); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_SNXDATA_LOOKBACK);}}
  if (strcmp(vname,"ISI_BOTHER")==0){ scanf("%d",&ISI_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,ISI_BOTHER);}}
  if (strcmp(vname,"GLOBAL_ISI_MAX")==0){ scanf("%lf",&GLOBAL_ISI_MAX); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_ISI_MAX);}}
  if (strcmp(vname,"GLOBAL_ISI_MIN")==0){ scanf("%lf",&GLOBAL_ISI_MIN); if (verbose){ printf("%s read to be %lf\n",vname,GLOBAL_ISI_MIN);}}
  if (strcmp(vname,"GLOBAL_ISI_NBINS")==0){ scanf("%d",&GLOBAL_ISI_NBINS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_ISI_NBINS);}}
  if (strcmp(vname,"LATTICE3D_BOTHER")==0){ scanf("%d",&LATTICE3D_BOTHER); if (verbose){ printf("%s read to be %d\n",vname,LATTICE3D_BOTHER);}}
  if (strcmp(vname,"GLOBAL_LATTICE3D_NBINS")==0){ scanf("%d",&GLOBAL_LATTICE3D_NBINS); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_LATTICE3D_NBINS);}}
  if (strcmp(vname,"GLOBAL_LATTICE3D_LENGTH")==0){ scanf("%d",&GLOBAL_LATTICE3D_LENGTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_LATTICE3D_LENGTH);}}
  if (strcmp(vname,"GLOBAL_LATTICE3D_DEPTH")==0){ scanf("%d",&GLOBAL_LATTICE3D_DEPTH); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_LATTICE3D_DEPTH);}}
  if (strcmp(vname,"GLOBAL_LATTICE3D_TET_VS_CUBE")==0){ scanf("%d",&GLOBAL_LATTICE3D_TET_VS_CUBE); if (verbose){ printf("%s read to be %d\n",vname,GLOBAL_LATTICE3D_TET_VS_CUBE);}}
  if (strcmp(vname,"GLOBAL_LATTICE3D_CS__")==0){ length=GLOBAL_LATTICE3D_LENGTH*GLOBAL_LATTICE3D_LENGTH; if (GLOBAL_LATTICE3D_CS__!=NULL){ tfree(GLOBAL_LATTICE3D_CS__); GLOBAL_LATTICE3D_CS__=NULL;} GLOBAL_LATTICE3D_CS__ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(GLOBAL_LATTICE3D_CS__[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,GLOBAL_LATTICE3D_CS__[nv]);}}}
  if (strcmp(vname,"GLOBAL_LATTICE3D_CS_")==0){ length=GLOBAL_LATTICE3D_LENGTH; if (GLOBAL_LATTICE3D_CS_!=NULL){ tfree(GLOBAL_LATTICE3D_CS_); GLOBAL_LATTICE3D_CS_=NULL;} GLOBAL_LATTICE3D_CS_ = (double *) tcalloc(length,sizeof(double)); for (nv=0;nv<length;nv++){ scanf("%lf",&(GLOBAL_LATTICE3D_CS_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose){ printf("%s[%d] read to be %lf\n",vname,nv,GLOBAL_LATTICE3D_CS_[nv]);}}}
  if (strcmp(vname,"END")==0){ /* do nothing */ if (verbose){ printf("end of input reached\n");}}
/*   if (strcmp(vname,"yy")==0){ scanf("%zz",&yy); if (verbose){ printf("%s read to be %zz\n",vname,yy);}} */
}

void dumpoutput(char *filename)
{
  /* prints an output file */
  int verbose=0;
  char text[256];
  char vname[128];
  FILE *fp=NULL;
  int length=0,nv=0;
  if (verbose){ printf(" %% [entering dumpoutput]\n");}
  sprintf(text,"./%s%s",filename,GLOBAL_STRING);
  if (filename==NULL || (fp = fopen(text, "w")) == NULL){ printf("dumpoutput to stdout\n"); fp = stdout;}
  sprintf(vname,"GLOBAL_verbose"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_verbose);
  sprintf(vname,"ON_MY_COMPUTER"); fprintf(fp,"%s= %d;\n",vname,ON_MY_COMPUTER);
  sprintf(vname,"GLOBAL_STRING"); fprintf(fp,"%s= %s;\n",vname,GLOBAL_STRING);
  sprintf(vname,"GLOBAL_RECORD_NUMBER"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_RECORD_NUMBER);
  sprintf(vname,"GLOBAL_SPIKEINPUT_RSEED"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_SPIKEINPUT_RSEED);
  sprintf(vname,"RUN_DONE"); fprintf(fp,"%s= %d;\n",vname,RUN_DONE);
  sprintf(vname,"GLOBAL_NTYPES"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_NTYPES);
  sprintf(vname,"GLOBAL_TYPENAMES"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%s%s",GLOBAL_TYPENAMES[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_NVARS"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_NVARS);
  sprintf(vname,"GLOBAL_NCLUSTERS"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_NCLUSTERS);
  sprintf(vname,"GLOBAL_NEURON_MODEL"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_NEURON_MODEL);
  if (GLOBAL_NEURON_MODEL==7){
    sprintf(vname,"GLOBAL_SNX_VERSION"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_SNX_VERSION);
    sprintf(vname,"GLOBAL_SNX_NSTATES_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_SNX_NSTATES_[nv],nv<length-1?",":";");} fprintf(fp,"\n");}
  sprintf(vname,"GLOBAL_LINK_SPARSE_OR_DENSE"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_LINK_SPARSE_OR_DENSE);
  sprintf(vname,"GLOBAL_CLUSTER_MAKE_OR_READ"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_CLUSTER_MAKE_OR_READ);
  sprintf(vname,"GLOBAL_CLUSTER_READ_FILE"); fprintf(fp,"%s= %s;\n",vname,GLOBAL_CLUSTER_READ_FILE);
  sprintf(vname,"GLOBAL_CLUSTER_PRA_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",GLOBAL_CLUSTER_PRA_[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_REWEIGHT_STRENGTH"); fprintf(fp,"%s= %lg;\n",vname,GLOBAL_REWEIGHT_STRENGTH);
  sprintf(vname,"GLOBAL_REWEIGHT_MULTIPLY_OR_ADD"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_REWEIGHT_MULTIPLY_OR_ADD);
  sprintf(vname,"GLOBAL_REWEIGHT_READ_FILE"); fprintf(fp,"%s= %s;\n",vname,GLOBAL_REWEIGHT_READ_FILE);
  sprintf(vname,"GLOBAL_VARNAMES"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%s%s",GLOBAL_VARNAMES[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_INDEXING_sra_LENGTH"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_INDEXING_sra_LENGTH);
  sprintf(vname,"GLOBAL_INDEXING_REFILE_sra"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_INDEXING_REFILE_sra[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_INDEXING_CHECKOUT_sra"); fprintf(fp,"%s= ",vname); length=GLOBAL_INDEXING_sra_LENGTH; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_INDEXING_CHECKOUT_sra[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_LENGTHRA"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_LENGTHRA[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"VOLTAGE_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",VOLTAGE_[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_VESICLE_DEPLETION"); fprintf(fp,"%s= %lg;\n",vname,GLOBAL_VESICLE_DEPLETION);
  sprintf(vname,"VOLTAGE_THRESHOLD_S"); fprintf(fp,"%s= %lg;\n",vname,VOLTAGE_THRESHOLD_S);
  sprintf(vname,"VOLTAGE_THRESHOLD_D"); fprintf(fp,"%s= %lg;\n",vname,VOLTAGE_THRESHOLD_D);
  sprintf(vname,"TAU_REF"); fprintf(fp,"%s= %lg;\n",vname,TAU_REF);
  sprintf(vname,"CONDUCTANCE_DS"); fprintf(fp,"%s= %lg;\n",vname,CONDUCTANCE_DS);
  sprintf(vname,"CONDUCTANCE_SD"); fprintf(fp,"%s= %lg;\n",vname,CONDUCTANCE_SD);
  sprintf(vname,"CURRENT_INJECTION_S"); fprintf(fp,"%s= %lg;\n",vname,CURRENT_INJECTION_S);
  sprintf(vname,"CURRENT_INJECTION_D"); fprintf(fp,"%s= %lg;\n",vname,CURRENT_INJECTION_D);
  sprintf(vname,"TAU_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",TAU_[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"CONDUCTANCE_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",CONDUCTANCE_[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"AUTAPSES_OFF"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",AUTAPSES_OFF[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_CS_ORN_SCALE"); fprintf(fp,"%s= %lg;\n",vname,GLOBAL_CS_ORN_SCALE);
  sprintf(vname,"CS_ORN_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",CS_ORN_[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  switch (GLOBAL_NEURON_MODEL){
  case 5: sprintf(vname,"CS_ORN_mainak_stim_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",CS_ORN_mainak_stim_[nv],nv<length-1?",":";");} fprintf(fp,"\n"); break;
  case 6: sprintf(vname,"CS_ORN_wilson_stim_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",CS_ORN_wilson_stim_[nv],nv<length-1?",":";");} fprintf(fp,"\n"); break;
  default: break;}
  sprintf(vname,"GLOBAL_CS_SCALE"); fprintf(fp,"%s= %lg;\n",vname,GLOBAL_CS_SCALE);
  sprintf(vname,"GLOBAL_CS_PRETYPE_SCALE_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",GLOBAL_CS_PRETYPE_SCALE_[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_CS_POSTYPE_SCALE_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",GLOBAL_CS_POSTYPE_SCALE_[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_CS_SRA_SCALE_"); fprintf(fp,"%s= ",vname); length=GLOBAL_INDEXING_sra_LENGTH; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",GLOBAL_CS_SRA_SCALE_[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"CS__"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES*GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",CS__[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"SPARSE__"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES*GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",SPARSE__[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"SPARSEHIT__"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES*GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",SPARSEHIT__[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"SPARSE_OGLI__"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES*GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",SPARSE_OGLI__[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"SPARSEHIT_OGLI__"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES*GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",SPARSEHIT_OGLI__[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_P_FAIL_SCALE"); fprintf(fp,"%s= %lg;\n",vname,GLOBAL_P_FAIL_SCALE);
  sprintf(vname,"P_FAIL__"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES*GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",P_FAIL__[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_ODOR_BASE"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_ODOR_BASE);
  sprintf(vname,"GLOBAL_ODORra_START"); fprintf(fp,"%s= ",vname); odorfprintf(fp,GLOBAL_ODORra_START);
  sprintf(vname,"INPUT_CONTRAST_START"); fprintf(fp,"%s= %lg;\n",vname,INPUT_CONTRAST_START);
  sprintf(vname,"INPUT_PULSE"); fprintf(fp,"%s= %lg;\n",vname,INPUT_PULSE);
  sprintf(vname,"ORN_BOTHER"); fprintf(fp,"%s= %d;\n",vname,ORN_BOTHER);
  if (ORN_BOTHER){
    sprintf(vname,"GLOBAL_ORN_TIMES_"); fprintf(fp,"%s= ",vname); length=4; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",GLOBAL_ORN_TIMES_[nv],nv<length-1?",":";");} fprintf(fp,"\n");}
  sprintf(vname,"ORN_BACKRATE"); fprintf(fp,"%s= %lg;\n",vname,ORN_BACKRATE);
  sprintf(vname,"SPIKETOL"); fprintf(fp,"%s= %d;\n",vname,SPIKETOL);
  sprintf(vname,"GLOBAL_SPACE_SMOOTHER"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_SPACE_SMOOTHER);
  sprintf(vname,"GLOBAL_TI"); fprintf(fp,"%s= %lg;\n",vname,GLOBAL_TI);
  sprintf(vname,"GLOBAL_TF"); fprintf(fp,"%s= %lg;\n",vname,GLOBAL_TF);
  sprintf(vname,"GLOBAL_DTmax"); fprintf(fp,"%s= %lg;\n",vname,GLOBAL_DTmax);
  sprintf(vname,"GLOBAL_dtadapt"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_dtadapt);
  sprintf(vname,"FIDDLE_ROW"); fprintf(fp,"%s= %d;\n",vname,FIDDLE_ROW);
  sprintf(vname,"FIDDLE_COL"); fprintf(fp,"%s= %d;\n",vname,FIDDLE_COL);
  sprintf(vname,"FIDDLE_PARAMETER"); fprintf(fp,"%s= %d;\n",vname,FIDDLE_PARAMETER);
  sprintf(vname,"DRAW_FLAG"); fprintf(fp,"%s= %d;\n",vname,DRAW_FLAG);
  sprintf(vname,"DRAW_FLAG2"); fprintf(fp,"%s= %d;\n",vname,DRAW_FLAG2);
  sprintf(vname,"STEPS_PER_DRAW"); fprintf(fp,"%s= %d;\n",vname,STEPS_PER_DRAW);
  sprintf(vname,"STD_VIEW"); fprintf(fp,"%s= %lg;\n",vname,STD_VIEW);
  sprintf(vname,"GRAYSCALE"); fprintf(fp,"%s= %d;\n",vname,GRAYSCALE);
  sprintf(vname,"SUPERGLOBAL_DRAW_FLAG"); fprintf(fp,"%s= %d;\n",vname,SUPERGLOBAL_DRAW_FLAG);
  sprintf(vname,"OUTPUT_DUMP_EVERY"); fprintf(fp,"%s= %lg;\n",vname,OUTPUT_DUMP_EVERY);
  sprintf(vname,"POWER_BOTHER"); fprintf(fp,"%s= %d;\n",vname,POWER_BOTHER);
  sprintf(vname,"GLOBAL_POWER_LENGTH"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_POWER_LENGTH);
  sprintf(vname,"GLOBAL_POWER_INDEXING_NTYPE_LENGTH"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_POWER_INDEXING_NTYPE_LENGTH);
  sprintf(vname,"GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT"); fprintf(fp,"%s= ",vname); length=GLOBAL_POWER_INDEXING_NTYPE_LENGTH; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_POWER_INDEXING_NTYPE_REFILE"); fprintf(fp,"%s= ",vname); length=GLOBAL_NTYPES; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_POWER_INDEXING_NTYPE_REFILE[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_POWER_INDEXING_NVAR_LENGTH"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_POWER_INDEXING_NVAR_LENGTH);
  sprintf(vname,"GLOBAL_POWER_INDEXING_NVAR_CHECKOUT"); fprintf(fp,"%s= ",vname); length=GLOBAL_POWER_INDEXING_NVAR_LENGTH; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_POWER_INDEXING_NVAR_CHECKOUT[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_POWER_INDEXING_NVAR_REFILE"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_POWER_INDEXING_NVAR_REFILE[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_POWER_TRAJECTORY_LENGTH"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_POWER_TRAJECTORY_LENGTH);
  sprintf(vname,"GLOBAL_POWER_WINDOW_LENGTH"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_POWER_WINDOW_LENGTH);
  sprintf(vname,"GLOBAL_POWER_WINDOW_UPDATE_EVERY"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_POWER_WINDOW_UPDATE_EVERY);
  sprintf(vname,"GLOBAL_POWER_maxra_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",GLOBAL_POWER_maxra_[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_POWER_minra_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",GLOBAL_POWER_minra_[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_POWER_CYCLE_BOTHER"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_POWER_CYCLE_BOTHER);
  sprintf(vname,"GLOBAL_POWER_CORRELATION_BOTHER"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_POWER_CORRELATION_BOTHER);
  sprintf(vname,"RHO_BOTHER"); fprintf(fp,"%s= %d;\n",vname,RHO_BOTHER);
  sprintf(vname,"GLOBAL_RHO_LENGTH"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_RHO_LENGTH);
  sprintf(vname,"GLOBAL_RHO_INDEXING_NVAR_LENGTH"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_RHO_INDEXING_NVAR_LENGTH);
  sprintf(vname,"GLOBAL_RHO_INDEXING_NVAR_CHECKOUT"); fprintf(fp,"%s= ",vname); length=GLOBAL_RHO_INDEXING_NVAR_LENGTH; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_RHO_INDEXING_NVAR_CHECKOUT[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_RHO_INDEXING_NVAR_REFILE"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_RHO_INDEXING_NVAR_REFILE[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"GLOBAL_RHO_NBINRA_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_RHO_NBINRA_[nv],nv<length-1?",":";");} fprintf(fp,"\n");
  sprintf(vname,"PTREE_BOTHER"); fprintf(fp,"%s= %d;\n",vname,PTREE_BOTHER);
  sprintf(vname,"GLOBAL_PTREE_BITBYBIT"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_PTREE_BITBYBIT);
  sprintf(vname,"GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER);
  sprintf(vname,"GLOBAL_PTREE_ZZZ"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_PTREE_ZZZ);
  sprintf(vname,"GLOBAL_PTREE_NREGIONS"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_PTREE_NREGIONS);
  sprintf(vname,"GLOBAL_PTREE_NLEGS"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_PTREE_NLEGS);
  sprintf(vname,"GLOBAL_PTREE_LEGTIME"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_PTREE_LEGTIME);
  sprintf(vname,"GLOBAL_PTREE_EVENT_THRESHOLD"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_PTREE_EVENT_THRESHOLD);
  sprintf(vname,"GLOBAL_PTREE_EVENT_WITHIN"); fprintf(fp,"%s= %lg;\n",vname,GLOBAL_PTREE_EVENT_WITHIN);
  sprintf(vname,"GLOBAL_PTREE_REGION_TYPE"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_PTREE_REGION_TYPE);
  sprintf(vname,"GLOBAL_PTREE_DUMP_TYPE"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_PTREE_DUMP_TYPE);
  sprintf(vname,"CLUSTERDATA_BOTHER"); fprintf(fp,"%s= %d;\n",vname,CLUSTERDATA_BOTHER);
  sprintf(vname,"SUITE_BOTHER"); fprintf(fp,"%s= %d;\n",vname,SUITE_BOTHER);
  sprintf(vname,"SUITE_NODORS"); fprintf(fp,"%s= %d;\n",vname,SUITE_NODORS);
  if (SUITE_ODORra_BACON!=NULL){
    sprintf(vname,"SUITE_ODORra_BACON"); fprintf(fp,"%s= ",vname); odorfprintf(fp,SUITE_ODORra_BACON);}
  sprintf(vname,"SUITE_CS_ORN_SCALE"); fprintf(fp,"%s= %lg;\n",vname,SUITE_CS_ORN_SCALE);
  sprintf(vname,"SUITE_INPUT_CONTRAST"); fprintf(fp,"%s= %lg;\n",vname,SUITE_INPUT_CONTRAST);
  sprintf(vname,"SUITE_NCONCENTRATIONS"); fprintf(fp,"%s= %d;\n",vname,SUITE_NCONCENTRATIONS);
  sprintf(vname,"SUITE_NBICUCULLINES"); fprintf(fp,"%s= %d;\n",vname,SUITE_NBICUCULLINES);
  sprintf(vname,"SUITE_NINSTANCES"); fprintf(fp,"%s= %d;\n",vname,SUITE_NINSTANCES);
  sprintf(vname,"SUITE_POWER_BOTHER"); fprintf(fp,"%s= %d;\n",vname,SUITE_POWER_BOTHER);
  sprintf(vname,"SUITE_PTREE_BOTHER"); fprintf(fp,"%s= %d;\n",vname,SUITE_PTREE_BOTHER);
  sprintf(vname,"SUITE_CLUSTERDATA_BOTHER"); fprintf(fp,"%s= %d;\n",vname,SUITE_CLUSTERDATA_BOTHER);
  sprintf(vname,"SUITE_4_REWEIGHT_STRENGTH"); fprintf(fp,"%s= %lg;\n",vname,SUITE_4_REWEIGHT_STRENGTH);
  sprintf(vname,"SUITE_7_CLEANUP"); fprintf(fp,"%s= %d;\n",vname,SUITE_7_CLEANUP);
  sprintf(vname,"SUITE_8_CLEANUP"); fprintf(fp,"%s= %d;\n",vname,SUITE_8_CLEANUP);
  sprintf(vname,"SUITE_8_wilson_axon_stim"); fprintf(fp,"%s= %lg;\n",vname,SUITE_8_wilson_axon_stim);
  sprintf(vname,"SUITE_NSECONDS"); fprintf(fp,"%s= %d;\n",vname,SUITE_NSECONDS);
  sprintf(vname,"SUITE_DUMPEVERY"); fprintf(fp,"%s= %d;\n",vname,SUITE_DUMPEVERY);
  sprintf(vname,"SUITE_BITBYBIT_RECORD"); fprintf(fp,"%s= %d;\n",vname,SUITE_BITBYBIT_RECORD);
  sprintf(vname,"SUITE_BITBYBIT_REMOVE"); fprintf(fp,"%s= %d;\n",vname,SUITE_BITBYBIT_REMOVE);
  sprintf(vname,"SUITE_SINDEXMAX"); fprintf(fp,"%s= %d;\n",vname,SUITE_SINDEXMAX);
  sprintf(vname,"SUITE_DINDEXMAX"); fprintf(fp,"%s= %d;\n",vname,SUITE_DINDEXMAX);
  sprintf(vname,"SUITE_TINDEXMAX"); fprintf(fp,"%s= %d;\n",vname,SUITE_TINDEXMAX);
  sprintf(vname,"LYAPUNOV_BOTHER"); fprintf(fp,"%s= %d;\n",vname,LYAPUNOV_BOTHER);
  sprintf(vname,"CAICOR_BOTHER"); fprintf(fp,"%s= %d;\n",vname,CAICOR_BOTHER);
  sprintf(vname,"GLOBAL_CAICOR_NBINS"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_CAICOR_NBINS);
  sprintf(vname,"HHLIB_BOTHER"); fprintf(fp,"%s= %d;\n",vname,HHLIB_BOTHER);
  if (HHLIB_BOTHER){
    sprintf(vname,"GLOBAL_HHLIB_MAKE_OR_READ"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_HHLIB_MAKE_OR_READ);
    sprintf(vname,"GLOBAL_HHLIB_READ_FILE"); fprintf(fp,"%s= %s;\n",vname,GLOBAL_HHLIB_READ_FILE);
    sprintf(vname,"GLOBAL_HHLIB_NBINS"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_HHLIB_NBINS);
    sprintf(vname,"GLOBAL_HHLIB_USE_AFTER"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_HHLIB_USE_AFTER);
    sprintf(vname,"GLOBAL_HHLIB_TAU_SKIP"); fprintf(fp,"%s= %lg;\n",vname,GLOBAL_HHLIB_TAU_SKIP);
    sprintf(vname,"GLOBAL_HHLIB_INDEXING_TRIGGER"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_HHLIB_INDEXING_TRIGGER);
    sprintf(vname,"GLOBAL_HHLIB_LOGFLAGRA"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_HHLIB_LOGFLAGRA[nv],nv<length-1?",":";");} fprintf(fp,"\n");
    sprintf(vname,"GLOBAL_HHLIB_LIBUSEFLAGRA"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_HHLIB_LIBUSEFLAGRA[nv],nv<length-1?",":";");} fprintf(fp,"\n");
    sprintf(vname,"GLOBAL_HHLIB_INDEXING_NVAR_LENGTH"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_HHLIB_INDEXING_NVAR_LENGTH);
    sprintf(vname,"GLOBAL_HHLIB_INDEXING_NVAR_CHECKOUT"); fprintf(fp,"%s= ",vname); length=GLOBAL_HHLIB_INDEXING_NVAR_LENGTH; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_HHLIB_INDEXING_NVAR_CHECKOUT[nv],nv<length-1?",":";");} fprintf(fp,"\n");
    sprintf(vname,"GLOBAL_HHLIB_INDEXING_NVAR_REFILE"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%d%s",GLOBAL_HHLIB_INDEXING_NVAR_REFILE[nv],nv<length-1?",":";");} fprintf(fp,"\n");
    sprintf(vname,"GLOBAL_HHLIB_maxra_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",GLOBAL_HHLIB_maxra_[nv],nv<length-1?",":";");} fprintf(fp,"\n");
    sprintf(vname,"GLOBAL_HHLIB_minra_"); fprintf(fp,"%s= ",vname); length=GLOBAL_NVARS; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",GLOBAL_HHLIB_minra_[nv],nv<length-1?",":";");} fprintf(fp,"\n");}
  sprintf(vname,"SNXDATA_BOTHER"); fprintf(fp,"%s= %d;\n",vname,SNXDATA_BOTHER);
  sprintf(vname,"GLOBAL_SNXDATA_LOOKBACK"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_SNXDATA_LOOKBACK);
  sprintf(vname,"ISI_BOTHER"); fprintf(fp,"%s= %d;\n",vname,ISI_BOTHER);
  if (ISI_BOTHER){
    sprintf(vname,"GLOBAL_ISI_MAX"); fprintf(fp,"%s= %lg;\n",vname,GLOBAL_ISI_MAX);
    sprintf(vname,"GLOBAL_ISI_MIN"); fprintf(fp,"%s= %lg;\n",vname,GLOBAL_ISI_MIN);
    sprintf(vname,"GLOBAL_ISI_NBINS"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_ISI_NBINS);}
  sprintf(vname,"LATTICE3D_BOTHER"); fprintf(fp,"%s= %d;\n",vname,LATTICE3D_BOTHER);
  if (LATTICE3D_BOTHER){
    sprintf(vname,"GLOBAL_LATTICE3D_NBINS"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_LATTICE3D_NBINS);
    sprintf(vname,"GLOBAL_LATTICE3D_LENGTH"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_LATTICE3D_LENGTH);
    sprintf(vname,"GLOBAL_LATTICE3D_DEPTH"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_LATTICE3D_DEPTH);
    sprintf(vname,"GLOBAL_LATTICE3D_TET_VS_CUBE"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_LATTICE3D_TET_VS_CUBE);
    sprintf(vname,"GLOBAL_LATTICE3D_CS__"); fprintf(fp,"%s= ",vname); length=GLOBAL_LATTICE3D_LENGTH*GLOBAL_LATTICE3D_LENGTH; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",GLOBAL_LATTICE3D_CS__[nv],nv<length-1?",":";");} fprintf(fp,"\n");
    sprintf(vname,"GLOBAL_LATTICE3D_CS_"); fprintf(fp,"%s= ",vname); length=GLOBAL_LATTICE3D_LENGTH; for (nv=0;nv<length;nv++){ fprintf(fp,"%lg%s",GLOBAL_LATTICE3D_CS_[nv],nv<length-1?",":";");} fprintf(fp,"\n");}
  sprintf(vname,"END"); fprintf(fp,"END= 0");
  if (fp!=stdout){ fclose(fp);}
/*   sprintf(vname,"yy"); fprintf(fp,"%s= ",vname); length=ww; for (nv=0;nv<length;nv++){ fprintf(fp,"%zz%s",yy[nv],nv<length-1?",":";");} fprintf(fp,"\n"); */
/*   sprintf(vname,"yy"); fprintf(fp,"%s= %zz;\n",vname,yy); */
  if (verbose){ printf(" %% [finishing dumpoutput]\n");}
}

void sparse_link_dump_ascii(struct neuronarray *Nra,char *fgvn)
{
  /* assumes fgvn starts with "./" 
     format:
     Nra->ntypes;
     Nra->lengthra[0],Nra->lengthra[1],...,Nra->lengthra[Nra->ntypes-1];
     LINK= ntA,nrA,NUMBER_OF_CONNECTED_NEURONS,nt1,nr1,nt2,nr2,nt3,nr3,...,ntN,nrN;
     LINK= ntB,nrB,NUMBER_OF_CONNECTED_NEURONS,nt1,nr1,nt2,nr2,nt3,nr3,...,ntN,nrN;
     ... 
     END;
  */
  char filename[512];
  int nt1=0,nr1=0;
  FILE *fp=NULL;
  struct neuron *n1=NULL,*n2=NULL;
  struct llist *L=NULL;
  struct litem *l=NULL;
  if (fgvn==NULL){ sprintf(filename,"./sparse_link_dump_ascii_%s",GLOBAL_STRING_2);}
  else /* if (fgvn!=NULL) */{ sprintf(filename,"%s",fgvn);}
  if ((fp=fopen(filename,"w"))==NULL){ printf(" warning, cannot open %s in sparse_link_dump_ascii\n",filename); fp=stdout;}
  fprintf(fp,"%d;\n",Nra->ntypes);
  for (nr1=0;nr1<Nra->ntypes-1;nr1++){ fprintf(fp,"%d,",Nra->lengthra[nr1]);}  
  fprintf(fp,"%d;\n",Nra->lengthra[Nra->ntypes-1]);
  for (nt1=0;nt1<Nra->ntypes;nt1++){ for (nr1=0;nr1<Nra->lengthra[nt1];nr1++){
    n1 = nget(Nra,nt1,nr1);
    L=llistmake();
    llitem2llist(n1->sparse_link,L);
    if (L->length>=1){
      fprintf(fp,"LINK= %d,%d,%d,",n1->type,n1->index,L->length);
      l=L->first;
      while (l!=NULL && l->child!=NULL){
	n2 = (struct neuron *)l->item;
	fprintf(fp,"%d,%d,",n2->type,n2->index);
	l=l->child;}
      l=L->last;
      n2 = (struct neuron *)l->item;
      fprintf(fp,"%d,%d;\n",n2->type,n2->index);}
    llisttfree(L);L=NULL;}}
  fprintf(fp,"END;");
  if (fp!=stdout){ fclose(fp);}
}

void sparse_link_read_ascii(struct neuronarray *Nra,char *fgvn)
{
  /* reads an ascii file dumped by sparse_link_dump_ascii */
  int verbose=GLOBAL_verbose;
  char filename[512];
  FILE *fp=NULL;
  char vname[128],equals[2],space[1],semicolon[1],comma_vs_semicolon[1];
  int nr=0,nt=0,nl=0,nl1=0,nr2=0,nt2=0,temp=0;
  struct neuron *n1=NULL,*n2=NULL;
  if (verbose){ printf(" %% [entering sparse_link_read_ascii]\n");}
  if (fgvn==NULL){ sprintf(filename,"./sparse_link_dump_ascii_%s",GLOBAL_STRING_2);}
  else /* if (fgvn!=NULL) */{ sprintf(filename,"%s",fgvn);}
  if (verbose){ printf("filename: %s\n",filename);}
  if ((fp=fopen(filename,"r"))==NULL){ printf(" warning, cannot open %s in sparse_link_read_ascii\n",filename); fp=stdin; exit(0);}
  fscanf(fp,"%d",&temp); if (verbose){ printf("ntypes read as %d, ",temp);}
  if (Nra->ntypes!=temp){ if (verbose){ printf("wrong\n");} exit(0);}
  else if (Nra->ntypes==temp){ if (verbose){ printf("right\n");}}
  fscanf(fp,"%c",semicolon);fscanf(fp,"%c",space);
  for (nt=0;nt<Nra->ntypes;nt++){
    fscanf(fp,"%d",&temp); if (verbose){ printf("lengthra[%d] read as %d, ",nt,temp);}
    if (Nra->lengthra[nt]!=temp){ if (verbose){ printf("wrong\n");} exit(0);}
    else if (Nra->lengthra[nt]==temp){ if (verbose){ printf("right\n");}}
    fscanf(fp,"%c",semicolon);}
  fscanf(fp,"%c",space);
  for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
    n1=nget(Nra,nt,nr); 
    if (n1->sparse_link!=NULL){ llitemtfree(n1->sparse_link,NULL); n1->sparse_link=NULL;}
    n1->sparse_link = llitemmake(); n1->sparse_link->item=n1;}}
  do{
    fscanf(fp,"%[^=]",vname);
    if (verbose){ printf(" vname %s\n",vname);}
    if (strcmp(vname,"LINK")==0){
      fscanf(fp,"%s",equals);fscanf(fp,"%c",space);
      fscanf(fp,"%d",&nt);fscanf(fp,"%c",comma_vs_semicolon);
      fscanf(fp,"%d",&nr);fscanf(fp,"%c",comma_vs_semicolon);
      fscanf(fp,"%d",&nl);fscanf(fp,"%c",comma_vs_semicolon);
      if (verbose){ printf("reading nt %d nr %d nl %d\n",nt,nr,nl);}
      if (nl>0){
	n1 = nget(Nra,nt,nr);
	nl1=0;
	while (nl1<nl){
	  fscanf(fp,"%d",&nt2);fscanf(fp,"%c",comma_vs_semicolon);
	  fscanf(fp,"%d",&nr2);fscanf(fp,"%c",comma_vs_semicolon);
	  if (verbose){ printf("connecting to nt2 %d nr2 %d \n",nt2,nr2);}
	  n2 = nget(Nra,nt2,nr2);
	  if (llitemaddorfind(0,n1->sparse_link,n2,&void_compare)==NULL){ llitemaddorfind(1,n1->sparse_link,n2,&void_compare);}
	  nl1+=1;}
	llitembalance(n1->sparse_link); n1->sparse_link=llitemclimb(n1->sparse_link);
	//assert(strcmp(comma_vs_semicolon,";")==0);
	fscanf(fp,"%c",space);}}}
  while (strcmp(vname,"LINK")==0);
  if (fp!=stdin){ fclose(fp);}
  if (verbose){ printf(" %% [finished sparse_link_read_ascii]\n");}
}    
