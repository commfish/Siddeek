// WAG GKC assessment model using 1981-1984 retained catch and full data sets for 1985/86 to 2015/16 
// This program uses standardized observer CPUE indices for 1991/92-2015/16 (estimated by GLM) in the likelihood function
// It considers two selectivity and catchability estimation periods, pre-rationalization 1985/86-2004/05, and post-rationalization  2005/06-2015/16
//
DATA_SECTION
// 
// read the data file
//
  init_int commstyr       //commercial fishery start year (1985) for input catch at length
  init_int commendyr      //commercial fishery & observer data end year (2015) for input catch, bycatch, and CPUE
  init_int nbins          // number of length bins (17) for all types of data
//
// Catch and relative abundance information
//
  init_matrix  retdcatch(commstyr,commendyr,1,nbins)                 // retained catch by length and year
  init_matrix  totalcatch(commstyr+5,commendyr,1,nbins)              // total catch by length and year, no data available for 1993
//
// mid lengths and weight
//
  init_ivector length_bins(1,nbins)                                 // mid length of each bin in mm CL
  init_vector retweight(1,nbins)                                    // weight at mid length in gram
//
// Elapsed time  
//
  init_vector eltime1(commstyr-4,commendyr)                         // elapsed time from July 1st to mid fishing time
//
// Effective sample size

  init_vector retcatcheffsample(commstyr,commendyr)                 // input retained catch effective sample size
  init_vector totalcatcheffsample(commstyr+5,commendyr)             // input total catch effective sample size
  init_vector gdiscdcatcheffsample(commstyr+4,commendyr)            // input groundfish discarded catch effective sample size
// 
// CPUE indices
//  
      init_vector obslegalcpueindex(commstyr+6,commendyr)           // Observer legal size crab CPUE index
      init_vector obslegalcpueindex_var(commstyr+6,commendyr)       // Observer legal size crab standard error of CPUE index
//
// groundfish discarded catch by length
//
      init_matrix gdiscdcatch(commstyr+4,commendyr,1,nbins)         // groundfish discarded catch 
//
// Fitted CPUE
//
      init_vector NBpredictedobslegalcpue(commstyr+6,commendyr)    // Negative binomial fitted model predicted CPUE, not used in the minimization
//
// Annual retained catch in number of crabs for 1981/82 to 1984/85
//
     init_vector retdcatch1(commstyr-4,commstyr-1)    
//
// Observer pot sample size for 1990 to 2015. 
//
     init_vector potsampsize(commstyr+5,commendyr)   
// 
//##################
// Read tagging data for Andre's Tag function
   init_int MaxTimeAtLib
   init_int NtagData                         // total number of tag records
   init_matrix TagData(1,NtagData,1,5)        // Tag data
   init_int DoTagDiag                         // a signal for diagnostic calculation
//#######################
   init_int eof;
!! if(eof !=999){cout<<"data reading error"<< endl;exit(1);}
!! cout<<"End of reading data file success and the end of file code is:"<<eof<<endl;
// End of reading normal data file
// Read control file
!! ad_comm::change_datafile_name("Adak8515Sc1Base.ctl");      // open the control file
  init_number alphar                                           // male mean length of recruitment
  init_int recbin                                              // number of length bins for recruit distribution
  init_int sel_phase                                           // selectivity phase
   init_number lowF                                            // lower limit of the pot F boundary
   init_number Mphase                                           // M estimation phase, not used here
  init_vector like_wght(1,13)                                  // weights for various likelihoods
  init_number m_disc                                           // pot fishery handling mortality
  init_number m_gdisc                                          // groundfish discard mortality
  init_number ryear1                                           // start year for mean R calculation for F35 calculation
  init_number ryear2                                           // end year for mean R calculation for F35 calculation
  init_int endctl
!! if(endctl !=9999){cout<<"control file reading error"<< endl;exit(1);}
!! cout<<"End of reading control file"<<endctl<<endl;
!! cout<<alphar<<endl;
!! cout<<recbin<<endl;
!! cout<<sel_phase<<endl;
!! cout<<lowF<<endl;
!! cout<<Mphase<<endl;
!! cout<<like_wght<<endl;
!! cout<<m_disc<<endl;
!! cout<<m_gdisc<<endl;
!! cout<<ryear1<<endl;
!! cout<<ryear2<<endl;
//
  vector obsretdcatchN(commstyr,commendyr)                  // observed annual retained catch number
  vector obsretdcatchB(commstyr,commendyr)                  // observed annual retained catch biomass
   vector obstotalcatchN(commstyr+5,commendyr)              // observed annual discarded catch number
  vector obstotalcatchB(commstyr+5,commendyr)               // observed annual discarded catch biomass
  vector obsgdiscdcatchN(commstyr+4,commendyr)              // observed annual groundfish bycatch number
  vector obsgdiscdcatchB(commstyr+4,commendyr)              // observed annual groundfish bycatch biomass
  matrix obsretlencomp(commstyr,commendyr,1,nbins)          // observed retained length composition by year
  matrix obstotallencomp(commstyr+5,commendyr,1,nbins)      // observed total length composition by year
  matrix obsgdiscdlencomp(commstyr+4,commendyr,1,nbins)     // observed groundfish discard length composition by year
  int NTagDiag2
!! NTagDiag2 = 0;
//
INITIALIZATION_SECTION
//   M  0.18
      log_a  2.54                                          // linear molt incremant intercep parameter 
   G_b -9.29                                               // linear molt increment slope parameter
    log_aa -2.50                                           // molt probabilty logistic parameter 
   log_b  4.95                                            // molt probability logistic parameter 
   log_T04L50 4.83                                          // total selectivity log L50 for 1985-2004
   log_T12L50 4.92                                           // total selectivity log L50 for 2005-2015
   log_R04L50 4.91                                          // retention curve log L50 for 1985-2015
    log_betar -0.73                                        // recruit distribution log beta parameter
    logq2 -0.63                                            // pot fishery log catchability for 1985-2004
   logq3 -0.96                                             // pot fishery log catchability for 2005-2015
       log_mean_rec 0.90                                // log mean number of recruits
     log_mean_Fpot -1.07                                   // log mean pot fishery F
    log_mean_Fground -9.51                                   // log mean groundfish bycatch F
   prelegal_var 0.02                                        // additional CPUE index varaince
    stdx 3.68                                               // normal distribution standard deviation parameter for the growth matrix 
//
PARAMETER_SECTION
 init_bounded_dev_vector rec_dev(commstyr-24,commendyr+1,-5.,5.,2)             // recruit deviation vector
 init_bounded_dev_vector Fpot_dev(commstyr-4,commendyr,-5.,5.,2)                //pot fishery F deviation vector
 init_bounded_dev_vector Fground_dev(commstyr+4,commendyr,-10.,15.,2)          // Groundfish bycatch F deviation vector
//
   init_bounded_number log_a(1.,4.5,2)                                         //  boudary for log_a of G=a+bl 
   init_bounded_number G_b(-12.0,-5.0,2)                                     //  boudary for b of G=a+bl 
 init_bounded_number log_aa(-4.61,-1.39,2)                                     // boudary for log_aa of logistic molt probability
 init_bounded_number log_b(3.869,5.05,2)                                          // boudary for log_b of logistic molt probability
    init_bounded_number stdx(0.1,12.0,3)                                        //  boudary for standard deviation of normal probability for growth
//
   init_bounded_number log_T04delta(0.,4.4,sel_phase)                          // boudary for total selectivity log differencee of L95 and L50 for 1985-2004
   init_bounded_number log_T12delta(0.,4.4,sel_phase)                          // boudary for total selectivity log differencee of L95 and L50 for 2005-2015
   init_bounded_number log_R04delta(0.,4.4,sel_phase)                          // boudary for retention curve log differencee of L95 and L50 for 1985-2015
//
 init_bounded_number log_T04L50(4.0,5.0,sel_phase)                            // boudary for total selectivity log L50 for 1985-2004
 init_bounded_number log_T12L50(4.0,5.0,sel_phase)                            // boudary for total selectivity log L50 for 2005 -2015
 init_bounded_number log_R04L50(4.0,5.0,sel_phase+1)                         // boudary for retention curve log L50 for 1985-2015
 init_bounded_number log_betar(-12,12.,3)                                      // boudary for recruitment distribution parameter log betar
 init_bounded_number logq2(-9.,2.25,sel_phase+2)                               // boudary for pot fishery log catchability for 1985-2004
 init_bounded_number logq3(-9.,2.25,sel_phase+3)                               // boudary for pot fishery log catchability for 2005 - 2015
 init_bounded_number log_mean_rec(0.01,5.,1)                                   // boudary for log mean recruit abundance
 init_bounded_number log_mean_Fpot(-15.0,lowF,1)                                // boudary for log mean F in pot fishery
 init_bounded_number log_mean_Fground(-15.0,-1.6,1)                            // boudary for log mean F in groundfish fishery bycatch
//
//  init_bounded_number M(0.01,0.5,Mphase)   
  init_number M(-1)
!! M=0.225;                                                                    // fixed M
//
  init_bounded_number prelegal_var(0.0,0.15,6)                                // boudary for legal cpue index additional variance
//  
  init_bounded_number Ftemp(0.,0.75,7)                                        // boudary for Fofl iteration
//
// Standard deviation report
//
      sdreport_vector mat_biomass(commstyr-25,commendyr+1)
      sdreport_vector legal_biomass(commstyr-25,commendyr+1)
      sdreport_number retain_harvest
      sdreport_number predofl 
     sdreport_number refpredofl                                                    // F35 reference point
      sdreport_number refretain_harvest                                             // F35 reference point
//
// Likelihood profiles
//
//
   objective_function_value f
//
  number x1
  number surv1
// 
// survival2 and survival3 functions are needed for the intermediate equilibrium abundance calculation
// 
  vector surv2(commstyr-4,commendyr)
//
  vector surv3(commstyr-4,commendyr)
  number surv4
  vector surv5(commstyr-4,commendyr)
  number qtemp
    number temp1
   number temp2
   number temp3
   number temp4
   number temp5
   number temp6
   number temp7
  matrix newsh(commstyr,commendyr+1,1,nbins)
  matrix retsigmasquare(commstyr,commendyr,1,nbins)
  matrix totalsigmasquare(commstyr+5,commendyr,1,nbins)
  matrix gdiscdsigmasquare(commstyr+4,commendyr,1,nbins)
  vector retcatchsum(commstyr,commendyr)
  vector totalcatchsum(commstyr,commendyr)
  vector gdiscdcatchsum(commstyr,commendyr)
  matrix retcatchlen(commstyr,commendyr,1,nbins)
  matrix totalcatchlen(commstyr,commendyr,1,nbins)
  matrix gdiscdcatchlen(commstyr,commendyr,1,nbins)
  matrix predretcatchlen(commstyr,commendyr,1,nbins)
  matrix predtotalcatchlen(commstyr,commendyr,1,nbins)
  matrix predgdiscdcatchlen(commstyr,commendyr,1,nbins)
  matrix predtotalcpue(commstyr,commendyr,1,nbins)
  vector predttotalcpue(commstyr,commendyr)
  matrix predretcpue(commstyr,commendyr,1,nbins)
  vector predtretcpue(commstyr,commendyr)
  vector predtretcpue91_15(commstyr+6,commendyr)
  matrix predtotalcatchNtemp(commstyr-4,commendyr+1,1,nbins)
  matrix predtotalcatchN(commstyr-4,commendyr+1,1,nbins)
  matrix preddiscdcatchN(commstyr-4,commendyr+1,1,nbins)
  vector predtotalcatchB(commstyr-4,commendyr)
  matrix predretdcatchN(commstyr-4,commendyr+1,1,nbins)
  vector predtretdcatchN(commstyr-4,commstyr-1)
  vector predretdcatchB(commstyr-4,commendyr)
  matrix predgdiscdcatchN(commstyr-4,commendyr+1,1,nbins)
  vector predgdiscdcatchB(commstyr-4,commendyr)
  number preddiscdcatchBB
  number predgdiscdcatchBB
//
   number like_retdcatchB
  number like_totalcatchB
  number like_gdiscdcatchB
  number like_rec_dev
  number like_retlencomp
  number like_totallencomp
  number like_gdiscdlencomp
  number like_retcpue
  number like_F
  number like_Mpenalty
  number like_gF
   number like_finalF
   number like_LLyr
   number like_fpen
   number like_meanFpot
   number like_initialLF
   number like_catch1
   vector totalCBWght(commstyr+5,commendyr)
//
   vector predretcatcheffsample(commstyr,commendyr)
   vector predgdiscdcatcheffsample(commstyr+4,commendyr)
   vector predtotalcatcheffsample(commstyr+5,commendyr)
//#############
  number fpen1
  number fpen2
//#############
//  vector expalphaN(1,nbins)
//  vector first_propN(1,nbins)
//
// Equilibrium initial population number by size
//
  matrix NN(commstyr-25,commstyr,1,nbins)
  vector Neq(1,nbins)
  vector x(1,nbins)
//
   number avg_F
   vector F_pot(commstyr-4,commendyr)
   number final_F
   vector gfmort_tem(commstyr+4,commendyr)
   vector gfmort_final(commstyr-4,commendyr)
   number mean_gfmort
    matrix z(commstyr-4,commendyr+1,1,nbins)
   vector zz(1,nbins)
//
  vector totalselect(1,nbins)
  vector totalselect1(1,nbins)
  vector totalselect2(1,nbins)
  vector retselect(1,nbins)
  vector retselect1(1,nbins)
  vector retselect2(1,nbins)
  vector  gselect(1,nbins)
  vector mean_gselect(1,nbins)
  vector tselectivity(1,nbins)
//
  matrix ret_stdresid(commstyr,commendyr,1,nbins)
  matrix total_stdresid(commstyr+5,commendyr,1,nbins)
  matrix gdiscd_stdresid(commstyr+4,commendyr,1,nbins)
//
    vector legalabund(commstyr-25,commendyr+1)
//  vector legal_biomass(commstyr-25,commendyr+1)
  vector matabund(commstyr-25,commendyr+1)
//  vector mat_biomass(commstyr-25,commendyr+1)
  vector harvestrate(commstyr-4,commendyr)
//  number predofl
   vector explegalB(commstyr,commendyr+1)
   vector legalcatch(commstyr,commendyr)   
  matrix  matabundtemp(commstyr,commendyr+1,1,nbins)
//   number retain_harvest
  number total_harvest
  number meanmatbiomass
   number meanmatbiomasstemp
//
  vector RR(commstyr-24,commendyr+1)
  vector rec_len(1,nbins)
     number alpha_rec
 vector molt(1,nbins)
  number t5
 vector t11(1,60)
  vector tt1(1,nbins)
  matrix len_len(1,nbins,1,nbins)
//
  vector meanx(1,nbins)
 matrix tt(1,nbins,1,nbins)
//
//################## 
// Andre's tag module
   number TagLike
   matrix S(1,nbins,1,nbins);  
   3darray TransX(1,2,1,nbins,1,nbins)
   4darray TransX2(1,2,1,MaxTimeAtLib,1,nbins,1,nbins)
   4darray  PredTagTable(1,2,1,MaxTimeAtLib,1,2,1,nbins)
   matrix  TagDiag2(1,NtagData,1,9)
//
//#########################
//
// F35 calculation variables
//
  matrix ref_na(1,100,1,nbins)
  matrix ref_na_m(1,100,1,nbins)
  matrix reftotalcatchNtemp(1,100,1,nbins) 
  matrix refretdcatchN(1,100,1,nbins) 
  matrix refdiscdcatchN(1,100,1,nbins) 
  matrix reftotalcatchN(1,100,1,nbins) 
  matrix refgdiscdcatchN(1,100,1,nbins)
  vector ref_catch_tot(1,100)
  vector ref_catch_ret(1,100)
  vector ref_discdcatch(1,100)
  vector ref_gdiscdcatch(1,100)
   number ref_F
   number ref_Ftrawl
   number ref_mr
  vector ref_Ftot(1,nbins)
  vector ref_mbio(1,101)
  vector ref_mbio215(1,100)
  vector ref_totc(1,101)
  vector ref_retc(1,101)
    number eb35
   number b25
    number f35
    number f40
    number i35
    number i40
    number ofl_f
    number last_mmb
    number last_mmbtemp
//    number refretain_harvest   //F35
//    number refpredofl          //F35
    number reftotal_harvest
    number refdiscdcatchBB
    number refgdiscdcatchBB
//
// Fournier's suggestion to centralize for estimation of confounding parameters
    number jcenter
//#####
//
PRELIMINARY_CALCS_SECTION
//
  int t; int l; int ii;
//
// For centerizing lengths for growth increment distribution normal model
//
   jcenter=double(sum(length_bins))/double(nbins);
//
//
//
// change the total catch biomass weight by observer pot sample size with the maximum set at 250. Septmeber 2015 CPT suggestion
//
    for(t = commstyr+5; t <= commendyr; t++)
    {totalCBWght(t) = potsampsize(t)*like_wght(6)/max(potsampsize);}
//
//
// convert the catch and bycatch into millions of crabs. All output are in millions and metric tonnes
//
  retdcatch /= 1000000.;    
  totalcatch /= 1000000.;                                                       
  gdiscdcatch /= 1000000.;                                                       

//
// 1981 to 1984 yearly catches in millions of crabs
//
  retdcatch1 /= 1000000.;
//
// Observed retained,pot and groundfish discard, and total pot catch by year in millions of crabs/metric tons
//
    obsretdcatchN.initialize();
    obsretdcatchB.initialize();
      obsretdcatchN=rowsum(retdcatch);
      obsretdcatchB=retdcatch*retweight;
//
     obsgdiscdcatchN.initialize();
    obsgdiscdcatchB.initialize();
     for(t = commstyr+4; t <= commendyr; t++)
  {
    if(t==1993)continue;                                    // no data
     for(l = 1; l <= nbins; l++)
    {
     obsgdiscdcatchN(t) += gdiscdcatch(t,l);
     obsgdiscdcatchB(t) += gdiscdcatch(t,l)*retweight(l);  
    }
  }
//
         obstotalcatchN.initialize();
         obstotalcatchB.initialize();
       for(t = commstyr+5; t <= commendyr; t++)
    {
          for(l = 1; l <= nbins; l++)
         {
         obstotalcatchN(t) += totalcatch(t,l);
         obstotalcatchB(t) += totalcatch(t,l)*retweight(l);  
         }
     }
//
    obsretlencomp.initialize();
   for(t = commstyr; t <= commendyr; t++)
  {
   for(l = 1; l <= nbins; l++)
    {
    obsretlencomp(t,l)=retdcatch(t,l)/obsretdcatchN(t);
    }
   }
//
    obsgdiscdlencomp.initialize();
   for(t = commstyr+4; t <= commendyr; t++)
  {   
    if(t==1993)continue;                                   // no data
   for(l = 1; l <= nbins; l++)
    {
    obsgdiscdlencomp(t,l)=gdiscdcatch(t,l)/obsgdiscdcatchN(t);
    } 
  }
//
         obstotallencomp.initialize();
       for(t = commstyr+5; t <= commendyr; t++)
    {
           for(l = 1; l <= nbins; l++)
         {
          obstotallencomp(t,l)=totalcatch(t,l)/obstotalcatchN(t);
         }
    }
//
   cout<<" preliminary calculation completed" <<endl;
//
PROCEDURE_SECTION
//
// |----------------------------------------------------------------|
// MAIN PROGRAM FUNCTION CALL
// |----------------------------------------------------------------|
//  
   totalselectivity();                          // total selectivity
   retselectivity();                            // retention
   get_moltingp();                              // molt
   growth_matrix();                             // growth 
  TaggingLikelihood();                          // tagging amodule
   get_mortality();                             // mortality
   gselectivity();                              // groundfish discard selectivity
   gmortality();                                // groundfish bycatch fishing mortality
  recruit();                                   // recruit ditribution
   get_initial_equilbrium_LFQ();                // initial equilibrium size distribution 
   get_population_number();                     // population dynamics
    get_cpue();                                  // cpue
    get_lengthcomp();                            // length composition
  if(last_phase())calc_final_statistics();     // last phase of optimization
   if(sd_phase())                               // after completing the optimization
   {
   get_reference_points();                      // F35
   } 
   objective_function();    
//
// |--------------------------------------------------------------|
// |  END OF MAIN PROGRAM FUNCTION CALLS
// |--------------------------------------------------------------|
//
//***************************************************
//
FUNCTION calc_final_statistics
//
// Effective sample size calculation from McAllister and Ianelli variance formula
//
  int t; int i; int l;
  temp1=0.; temp2=0.; temp3=0.;temp4=0.; temp5=0.; temp6=0.;
  predretcatcheffsample.initialize();
  predgdiscdcatcheffsample.initialize();
  predtotalcatcheffsample.initialize();
  for(t = commstyr;t <= commendyr; t++)
    {
     for(l = 1; l <= nbins; l++)
      {
       temp1 +=predretcatchlen(t,l)*(1.0-predretcatchlen(t,l));
       temp2 += square(obsretlencomp(t,l)-predretcatchlen(t,l));
      }
      if(temp2>0.)predretcatcheffsample(t)= temp1/temp2;
    }
//
       for(t = commstyr+4;t <= commendyr; t++)
    {
       if(t==1991)continue;
     for(l = 1; l <= nbins; l++)
         {
          temp3 +=predgdiscdcatchlen(t,l)*(1.0-predgdiscdcatchlen(t,l));
          temp4 += square(obsgdiscdlencomp(t,l)-predgdiscdcatchlen(t,l));
         }
      if(temp4>0.)predgdiscdcatcheffsample(t)= temp3/temp4;
    }
//
   for(t = commstyr+5; t <= commendyr; t++)
    {
      if(t==1993)continue;
        for(l = 1; l <= nbins; l++)
       {
        temp5 +=predtotalcatchlen(t,l)*(1.0-predtotalcatchlen(t,l));
        temp6 += square(obstotallencomp(t,l)-predtotalcatchlen(t,l));
       }
     if(temp6>0.)predtotalcatcheffsample(t)=temp5/temp6;
   }
//
// Predicted legal male abundance in number (millions) and weight (t) by year at survey time
//
     explegalB.initialize();
     for(t = commstyr;t <= commendyr+1; t++)
    {
     for(l = 8; l <= nbins; l++)        // bin number 8 starts from 136 mm CL, middle length 138 mm CL
      {
      legalabund(t) += newsh(t,l);
      legal_biomass(t) += newsh(t,l)*retweight(l);
      if(t<2005)                                                     // total selectivity are in two blocks, pre- and post rationalization time periods
      {explegalB(t) += newsh(t,l)*totalselect1(l)*retweight(l);}  
       else
      {explegalB(t) += newsh(t,l)*totalselect2(l)*retweight(l);}
      }
    }
//
//Predicted pot fishery total harvest rate based on predicted legal size catch divided by legal abundance
//
   for(t = commstyr; t <= commendyr; t++)
    {
       if(legalabund(t)>0.)
         {harvestrate(t) = legalcatch(t)/legalabund(t);}  
       else
         {harvestrate(t) =0.;}
    }
//
// Predicted mature male (>=121 mmCL) abundance in number (millions) and weight (t) by year on next 15 February upto
// final data avialability year, ie mature abundance estimates are on 86Feb15 .... 2016Feb15
//
    for(t = commstyr;t <= commendyr; t++)
    {
        for(l = 5; l <= nbins; l++)        // bin number 5 starts from 121 mm CL (50% maturity length)
      {
       matabund(t) += newsh(t,l)*surv4-(predtotalcatchN(t,l)+predgdiscdcatchN(t,l))*surv5(t);
       mat_biomass(t) += (newsh(t,l)*surv4-(predtotalcatchN(t,l)+predgdiscdcatchN(t,l))*surv5(t))*retweight(l);
      }
    }
         meanmatbiomasstemp=0.0;
    for(t = commstyr;t <= commendyr; t++)
    meanmatbiomasstemp +=mat_biomass(t);
    meanmatbiomass=meanmatbiomasstemp/(commendyr-commstyr+1);          // mean mature biomass for 1986Feb15 to 2016Feb15
//
// Predicted mature male (>=121 mmCL) abundance in number (millions) and weight (t) on next year 15 February for one step forward projection
//
   for(l = 1; l <= nbins; l++)
    {
    zz(l)=0.8*Ftemp*totalselect2(l)*retselect2(l)+0.2*Ftemp*totalselect2(l)+0.65*gfmort_final(commendyr)*gselect(l);   // updated total mortality formula
    predtotalcatchNtemp(commendyr+1,l) = newsh(commendyr+1,l)*surv3(commendyr)*(1.-exp(-zz(l)))*Ftemp*totalselect2(l)/zz(l);   
     predretdcatchN(commendyr+1,l) = predtotalcatchNtemp(commendyr+1,l)*retselect2(l);
     preddiscdcatchN(commendyr+1,l) = m_disc*predtotalcatchNtemp(commendyr+1,l)*(1.0-retselect2(l));
    predtotalcatchN(commendyr+1,l) = predretdcatchN(commendyr+1,l)+preddiscdcatchN(commendyr+1,l);
     predgdiscdcatchN(commendyr+1,l) = m_gdisc*newsh(commendyr+1,l)*surv3(commendyr)*(1.-exp(-zz(l)))*gfmort_final(commendyr)*gselect(l)/zz(l);   
    }
//
        for(l = 5; l <= nbins; l++)         
      {
        matabundtemp(commendyr+1,l)= newsh(commendyr+1,l)*surv4-(predtotalcatchN(commendyr+1,l)+predgdiscdcatchN(commendyr+1,l))*surv5(commendyr);
       matabund(commendyr+1) +=      matabundtemp(commendyr+1,l);
       mat_biomass(commendyr+1) +=  matabundtemp(commendyr+1,l)*retweight(l);
      }
//
// predicted F_OFL by Tier 4 formula, M multiplier is set to 1
//
   if(mat_biomass(commendyr+1)>=meanmatbiomass)final_F = M;
   if(mat_biomass(commendyr+1)<meanmatbiomass && mat_biomass(commendyr+1)>0.25*meanmatbiomass)
      {final_F = M*(mat_biomass(commendyr+1)/meanmatbiomass-0.1)/0.9;}
     if(mat_biomass(commendyr+1)<= 0.25*meanmatbiomass)final_F=0.;
//
//  Use the F to get the predicted harvest for 2015/16:
//
      retain_harvest=0.;
    total_harvest=0.;
   preddiscdcatchBB=0.;
    predgdiscdcatchBB=0.;
    predofl=0.;
     for(l = 1; l <= nbins; l++)
      {
        zz(l)=0.8*final_F*totalselect2(l)*retselect2(l)+0.2*final_F*totalselect2(l)+0.65*gfmort_final(commendyr)*gselect(l);  
    total_harvest +=  retweight(l)*newsh(commendyr+1,l)*surv3(commendyr)*(1.-exp(-zz(l)))*final_F*totalselect2(l)/zz(l);
    retain_harvest += retweight(l)*(newsh(commendyr+1,l))*surv3(commendyr)*(1.-exp(-zz(l)))*final_F*totalselect2(l)*retselect2(l)/zz(l);
    predgdiscdcatchBB += retweight(l)*m_gdisc*(newsh(commendyr+1,l))*surv3(commendyr)*(1.-exp(-zz(l)))*gfmort_final(commendyr)*gselect(l)/zz(l);
      }
    preddiscdcatchBB = m_disc*(total_harvest-retain_harvest);
//
   predofl=retain_harvest+predgdiscdcatchBB+preddiscdcatchBB;
//
//
// Tier 3  F35 Calculation module
//
FUNCTION get_reference_points
    cout<<"start reference points"<<endl;
   dvariable mr;
   int i, j, l, k, ll;
// 
// reference recruitment
//
   ref_mr=0.;
   ref_Ftrawl=0.;
   mr = 0.;
  for(i=ryear1; i<=ryear2; i++)
  {
     mr += mfexp(log_mean_rec+rec_dev(i));
  }
  ref_mr = mr/(ryear2-ryear1+1.0);           //mean recruitment in millions for B35 estimation
//
//  Ten-year mean reference trawl fishing mortality
//
  for (i = commendyr-9; i<=commendyr; i++)
  {
    ref_Ftrawl +=gfmort_final(i);
  }
// 
   ref_Ftrawl=ref_Ftrawl/10.0;
//
// per recruit calculation
//
  for (j = 1; j<= 101; j++)   
 {
     ref_F = 0.01*j-0.01;
//
     ref_catch_tot.initialize();
     ref_catch_ret.initialize();
     ref_gdiscdcatch.initialize();
     ref_na.initialize();
     ref_na_m.initialize();
// 
//initial year
//
   ref_na_m(1)  = newsh(commendyr);    
   ref_na(1)=ref_na_m(1);
//
// Now add Recruits
//    
   for (i=2;i<=100;i++)     
     {
      ref_na_m(i) += 1000000.0 *rec_len;  
     }
//
//numbers at length for reference point estimation
//
   for (i=1;i<100;i++)   
   {
// Predicted catch for substraction
        for(l = 1; l <= nbins; l++)
         {
       ref_Ftot(l)=0.8*ref_F*totalselect2(l)*retselect2(l)+0.2*ref_F*totalselect2(l)+0.65*ref_Ftrawl*gselect(l);                 
          reftotalcatchNtemp(i,l) = ref_na(i,l)*surv3(commendyr)*(1.-mfexp(-ref_Ftot(l)))*ref_F*totalselect2(l)/ ref_Ftot(l);  
          refretdcatchN(i,l) = reftotalcatchNtemp(i,l)*retselect2(l);
          refdiscdcatchN(i,l) = m_disc*reftotalcatchNtemp(i,l)*(1.0-retselect2(l));
          reftotalcatchN(i,l) = refretdcatchN(i,l)+refdiscdcatchN(i,l);
          refgdiscdcatchN(i,l) = m_gdisc*ref_na(i,l)*surv3(commendyr)*(1.-mfexp(-ref_Ftot(l)))*ref_Ftrawl*gselect(l)/ ref_Ftot(l);  
         }
//      
// 
       for(l = 1; l <= nbins; l++)
        {
         for(ll = 1; ll <= l; ll++)
         {
       dvariable ref_tmp = (ref_na(i,ll)*surv1 - (reftotalcatchN(i,ll)+refgdiscdcatchN(i,ll))*surv2(commendyr))*len_len(ll,l);
      ref_na_m(i+1,l) += ref_tmp;      
         }
        }
         ref_na(i+1) = ref_na_m(i+1);
////
         ref_mbio215(i)=0.0;
     for(l = 5; l <= nbins; l++)        // bin number 5 starts from 121 mm CL, middle length 123 mm CL
      {
       ref_mbio215(i) += (ref_na(i,l)*surv4-(reftotalcatchN(i,l)+refgdiscdcatchN(i,l))*surv5(commendyr))*retweight(l);    
      }
////
//#####
    for(l = 1; l <= nbins; l++)
     {
      ref_Ftot(l)=0.8*ref_F*totalselect2(l)*retselect2(l)+0.2*ref_F*totalselect2(l)+0.65*ref_Ftrawl*gselect(l);             
      reftotalcatchNtemp(i,l) = ref_na(i,l)*surv3(commendyr)*(1.-mfexp(-ref_Ftot(l)))*ref_F*totalselect2(l)/ref_Ftot(l);  
      ref_catch_tot(i)  += retweight(l)*reftotalcatchNtemp(i,l);
      ref_catch_ret(i)  += retweight(l)*reftotalcatchNtemp(i,l)*retselect2(l);
      ref_gdiscdcatch(i)  += retweight(l)*m_gdisc*(ref_na(i,l))*surv3(commendyr)*(1.-mfexp(-ref_Ftot(l)))*ref_Ftrawl*gselect(l)/ref_Ftot(l);
     }
      ref_discdcatch(i) = m_disc*(ref_catch_tot(i)-ref_catch_ret(i));
  }
   ref_mbio(j) = ref_mbio215(99)/1000.0;               //MMB/R, kg/R
   ref_totc(j) = ref_catch_tot(99)/1000.0;             // total catch /R
   ref_retc(j) = ref_catch_ret(99)/1000.0;             // retained catch /R
 }
//
    i35 = 0;
    i40 = 0;
  for (j = 1; j<= 101; j++)
  {
    if (i35 < 1.0)
    {
      if (ref_mbio(j) <= 0.35*ref_mbio(1))
      {
         f35 = 0.01*j-0.01;
         eb35 = ref_mbio(j)*ref_mr/1000.0;      // in t
         b25 = 0.25*eb35;
         i35 = 2.0;
      }
    }
    if (i40 < 1.0)
    {
      if (ref_mbio(j) <= 0.40*ref_mbio(1))
      {
         f40 = 0.01*j-0.01;
         i40 = 2.0;
      }
    }
  }
//
// Final year OFL and MMB estimation
//
       last_mmbtemp=0.0;
         ofl_f = f35;
     for(l = 1; l <= nbins; l++)
       {
    ref_Ftot(l)=0.8*ofl_f*totalselect2(l)*retselect2(l)+0.2*ofl_f*totalselect2(l)+0.65*ref_Ftrawl*gselect(l);                 
        predtotalcatchNtemp(commendyr,l) = newsh(commendyr,l)*surv3(commendyr)*(1.-mfexp(-ref_Ftot(l)))*ofl_f*totalselect2(l)/ref_Ftot(l);  
        predretdcatchN(commendyr,l) = predtotalcatchNtemp(commendyr,l)*retselect2(l);
        preddiscdcatchN(commendyr,l) = m_disc*predtotalcatchNtemp(commendyr,l)*(1.0-retselect2(l));
        predtotalcatchN(commendyr,l) = predretdcatchN(commendyr,l)+preddiscdcatchN(commendyr,l);
        predgdiscdcatchN(commendyr,l) = m_gdisc*newsh(commendyr,l)*surv3(commendyr)*(1.-mfexp(-ref_Ftot(l)))*ref_Ftrawl*gselect(l)/ref_Ftot(l);  
       }
//
         for(l = 5; l <= nbins; l++)        // bin number 5 starts from 121 mm CL, middle length 123 mm CL
         {
          matabundtemp(commendyr,l)= newsh(commendyr,l)*surv4-(predtotalcatchN(commendyr,l)+predgdiscdcatchN(commendyr,l))*surv5(commendyr);
          last_mmbtemp += matabundtemp(commendyr,l)*retweight(l); 
         }
       last_mmb =  last_mmbtemp;                            
//
      if (last_mmb < eb35)
  {
     for (k = 1; k<10; k++)   //  this loop is the iteration to find the ofl_f when B<B35%, 10 iterations are sufficient
     {
        last_mmbtemp=0.;
        for(l = 1; l <= nbins; l++)
        {
      ref_Ftot(l)=0.8*ofl_f*totalselect2(l)*retselect2(l)+0.2*ofl_f*totalselect2(l)+0.65*ref_Ftrawl*gselect(l);                 
        predtotalcatchNtemp(commendyr,l) = newsh(commendyr,l)*surv3(commendyr)*(1.-mfexp(-ref_Ftot(l)))*ofl_f*totalselect2(l)/ref_Ftot(l);  
        predretdcatchN(commendyr,l) = predtotalcatchNtemp(commendyr,l)*retselect2(l);
        preddiscdcatchN(commendyr,l) = m_disc*predtotalcatchNtemp(commendyr,l)*(1.0-retselect2(l));
        predtotalcatchN(commendyr,l) = predretdcatchN(commendyr,l)+preddiscdcatchN(commendyr,l);
        predgdiscdcatchN(commendyr,l) = m_gdisc*newsh(commendyr,l)*surv3(commendyr)*(1.-mfexp(-ref_Ftot(l)))*ref_Ftrawl*gselect(l)/ref_Ftot(l);  
        }
//
         for(l = 5; l <= nbins; l++)        // bin number 5 starts from 121 mm CL, middle length 123 mm CL
         {
          matabundtemp(commendyr,l)= newsh(commendyr,l)*surv4-(predtotalcatchN(commendyr,l)+predgdiscdcatchN(commendyr,l))*surv5(commendyr);
          last_mmbtemp += matabundtemp(commendyr,l)*retweight(l); 
         }   
            last_mmb =  last_mmbtemp; 
//
      if (last_mmb < b25)
        {
           ofl_f = 0.0;
        }
        else
        {
            ofl_f = f35*(last_mmb/eb35-0.1)/0.9;
        }
     }
  }
//  
// Use the F to get the predicted harvest for 2015/16:
//
    refretain_harvest=0.;
    reftotal_harvest=0.;
    refdiscdcatchBB=0.;
    refgdiscdcatchBB=0.;
   for(l = 1; l <= nbins; l++)
      {
    ref_Ftot(l)=0.8*ofl_f*totalselect2(l)*retselect2(l)+0.2*ofl_f*totalselect2(l)+0.65*ref_Ftrawl*gselect(l);                 
    reftotal_harvest +=  retweight(l)*newsh(commendyr+1,l)*surv3(commendyr)*(1.-mfexp(-ref_Ftot(l)))*ofl_f*totalselect2(l)/ref_Ftot(l);
    refretain_harvest += retweight(l)*(newsh(commendyr+1,l))*surv3(commendyr)*(1.-mfexp(-ref_Ftot(l)))*ofl_f*totalselect2(l)*retselect2(l)/ref_Ftot(l);
    refgdiscdcatchBB += retweight(l)*m_gdisc*(newsh(commendyr+1,l))*surv3(commendyr)*(1.-mfexp(-ref_Ftot(l)))*ref_Ftrawl*gselect(l)/ref_Ftot(l);
      }
    refdiscdcatchBB = m_disc*(reftotal_harvest-refretain_harvest);
////
   refpredofl=refretain_harvest+refgdiscdcatchBB+refdiscdcatchBB;
// 
//####
  ref_mbio/=1000.0;
  ref_totc/=1000.0;
  ref_retc/=1000.0;
//####
  ofstream report2("Adak8515Sc1BaseF35.out");
  report2 <<"Reference MMB per recruit in t on 2/15 as F = 0.00, 0.01, ... 1.0"<<endl;
  report2 << ref_mbio<<endl;
  report2 <<"Reference total catch per recruit in t as F = 0.00, 0.01, ... 1.0"<<endl;
  report2 << ref_totc<<endl;
  report2 <<"Reference retained catch per recruit in t as F = 0.00, 0.01, ... 1.0"<<endl;
  report2 << ref_retc<<endl;
  report2 <<"F35: "<<f35<<"  B35 in t: "<<eb35<<" Mean R in millions: "<<ref_mr<<endl;
  report2 <<"F40: "<<f40<<endl;
  report2 <<"ref_Ftrawl = (mean of 10 years trawl F)  "<<endl;
  report2 <<ref_Ftrawl<<endl;
   report2 <<"Total catch OFL in t=  "<<refpredofl<<endl;
   report2 <<"Retained catch OFL in t =  "<<refretain_harvest<<endl;
   report2 <<"Male pot bycatch OFL in t=  "<<refdiscdcatchBB<<endl;
   report2 <<"Trawl bycatch OFL  in t=  "<<refgdiscdcatchBB<<endl;
   report2 <<"OFL F =  "<<ofl_f<<endl;
   report2 <<"MMB at terminal year in t =  "<<last_mmb<<endl;
//*****************************************
  cout<<"  end reference points calculation "<<endl; 
//*****************************************
//*****************************************
// Two sets of parameters for fishery selectivity for pre- and post rationalization periods
//
FUNCTION totalselectivity
   int t;
   int j;
   for(t = commstyr;t <= commendyr; t++)
   {
     for(j = 1; j <= nbins; j++)
     {
     if(t<= 2004)totalselect1(j)=1./(1.+mfexp(-log(19.)*(double(length_bins(j))-mfexp(log_T04L50))/(mfexp(log_T04delta))));
     if(t>2004)totalselect2(j)=1./(1.+mfexp(-log(19.)*(double(length_bins(j))-mfexp(log_T12L50))/(mfexp(log_T12delta))));
     }
   }
//
FUNCTION retselectivity
// a single retention curve for both periods
//
 int t;
 int j;
 for(t = commstyr;t <= commendyr; t++)
 {
  for(j = 1; j <= nbins; j++)
   {
    if(t<= 2004)retselect1(j)=1./(1.+mfexp(-log(19.)*(double(length_bins(j))-mfexp(log_R04L50))/(mfexp(log_R04delta))));
    if(t> 2004)retselect2(j)=1./(1.+mfexp(-log(19.)*(double(length_bins(j))-mfexp(log_R04L50))/(mfexp(log_R04delta))));
    }
  }
//
FUNCTION gselectivity
   int j;
  for(j = 1; j <= nbins; j++)
   {
//  gselect(j)=1./(1.+mfexp(-log(19.)*(double(length_bins(j))-mfexp(log_gL50))/(mfexp(log_gLdelta))));
     gselect(j)=1.;                                                                              // gf selectivity is set to one
   }
//
FUNCTION gmortality
    int t;
    gfmort_tem.initialize();
   mean_gfmort.initialize();
 for(t = commstyr+4;t <= commendyr; t++)
    {
      gfmort_tem(t)=mfexp(log_mean_Fground+Fground_dev(t));
     }
    temp7=0.;
  for(t = commstyr+4;t <= commstyr+9; t++)
     if(t==1993)continue;
       {temp7 += gfmort_tem(t);}
      mean_gfmort = temp7/((commstyr+9)-(commstyr+4));   // 5 year average
//
    gfmort_final.initialize();
 for(t = commstyr-4;t <= commendyr; t++)
  {
    if(t<(commstyr+4))
     {gfmort_final(t)=mean_gfmort;}
    else
    {gfmort_final(t)=mfexp(log_mean_Fground+Fground_dev(t));}
     if(t==1993)gfmort_final(t)=0.0;
  }
//********************////////////////////////////////////////////////////////////////////
//
FUNCTION recruit
////
  int ilen, il2, j, ii;
  dvariable t2, t3, t4;
  rec_len.initialize();
  t11.initialize();
//
// The gamma function is truncated and numerically simplified.
//
       alpha_rec = mfexp(log_betar);
      t4=alphar/alpha_rec;                               // alphar is the mean recruit length
      ii = (recbin)*5;
      ilen = 5;
      t2 = 0.0;
      t3 = (t4-1.0)*log(alphar)-t4-t4*log(alpha_rec);
      for (j = 1; j<=ii;j++)
      {
         t5 = ((j-1)+100.0)*1.0;
        t11(j) = (t4-1.0)*log(t5)-t5/alpha_rec-t4*log(alpha_rec)-t3;
        if (t11(j) < -1.0e30) t11(j)=-1.0e30;
        t11(j) = mfexp(t11(j));
        t2 = t2 + t11(j);
      }
      for (j=1;j<=recbin-0;j++)
      {
        for (il2=1;il2<=ilen;il2++) rec_len(j) += t11((j-1)*ilen+il2)/t2;
      }
//*******************
FUNCTION get_moltingp
//
//assuming a declining logistic function
//
    int j;
    for(j = 1; j <= nbins; j++)
     {
        molt(j)=1./(1.+mfexp(mfexp(log_aa)*(double(length_bins(j))-mfexp(log_b))));
     }
//
FUNCTION growth_matrix
//
//estimate growth matrix with the normal probability model
//
   int i, j;
   len_len.initialize();
   tt.initialize();
   for (j=1;j<=nbins;j++)
     {
//#####
      meanx(j) = mfexp(log_a)+G_b*((j*5.0+98.-jcenter)/jcenter);    // centralized length to reduce the confounding between log_a and G_b
//#####
     }
// 
  for (j=1;j< nbins;j++)
   {
      for (i=j;i<=nbins;i++)
       {
        tt(j,i)=((i-j+1)*5.0-2.5-meanx(j))/stdx;
          if(i==j)
           {len_len(j,i) = cumd_norm(tt(j,i));}
          else if(i<nbins)
           {len_len(j,i) = cumd_norm(tt(j,i))-cumd_norm(tt(j,i-1));}
          else
           {len_len(j,i) = 1.0-cumd_norm(tt(j,i-1));}
       }
    }
    len_len(nbins,nbins) = 1.;
//
// Regulate the growth matrix by molt probability
//
   for (j=1;j<=nbins;j++)
   {
      for (i=j;i<=nbins;i++)
      {
        len_len(j,i) = len_len(j,i)*molt(j);
        if (i==j) len_len(j,i) = len_len(j,i)+(1.0-molt(j));
      }
   }
//
FUNCTION get_mortality
  int t;
  F_pot.initialize();
  surv2.initialize();
  surv3.initialize();
  surv5.initialize();
  surv1 = exp(-M);
  surv4=exp(-0.627397*M);  // survival on Feb 15 of crabs appeared on July 1
//#####
//
// surv 2, 3, and 5 calculation for use during the 1981/82 to 1983/84 period
//
   for(t = commstyr-4; t <= commendyr; t++)    
   {
    F_pot(t)=mfexp(log_mean_Fpot+Fpot_dev(t));
    surv3(t) =exp(-eltime1(t)*M);
     surv2(t) = exp((eltime1(t)-1.)*M);
    surv5(t) = exp((eltime1(t)-0.627397)*M);
   }
//
// Tagging likelihood
// Andre's module
// Tagged crabs released in 1991 1997 2000 2003 and 2006 and recoveries obtained up to 2012
//****
// 
//Andre's tag function
//
FUNCTION TaggingLikelihood
 int Itag, Isex, RelSize, RecSize, TimeAtLib,Isize, Jsize,Ksize;
 int Iyear, Ipnt;
 dvariable MeanLenRec, Total, Total1, Slope, Intc, Pred, VarLenRec, Cnt, wt; 
 dvariable Bot, Top, MeanPred, MeanObs, VarPred, EffN;
 dvar_vector Preds(1,nbins);
 dvar_vector Obs(1,nbins),selexS(1,nbins);
 dvar_matrix TempX(1,nbins,1,nbins), TempY(1,nbins,1,nbins), PredTag(1,NtagData,0,3), PredProb(1,NtagData,1,nbins);
 dvar_matrix  tgrowth(1,nbins,1,nbins);
// dvar_3darray TransX(1,2,1,nbins,1,nbins);
// dvar_4darray TransX2(1,2,1,MaxTimeAtLib,1,nbins,1,nbins), PredTagTable(1,2,1,MaxTimeAtLib,1,2,1,nbins);
//
   tgrowth=trans(len_len);                 // transpose the growth matrix
   for (Isize=1;Isize<=nbins;Isize++)
        for (Jsize=1;Jsize<=nbins;Jsize++)
  TransX(1,Isize,Jsize)=tgrowth(Isize,Jsize);
   selexS=totalselect1;
//
 TransX2.initialize();
// Set up (one-year-ahead)
        for (Isize=1;Isize<=nbins;Isize++)
        for (Jsize=1;Jsize<=nbins;Jsize++)
        {
         TempX(Isize,Jsize) = TransX(1,Isize,Jsize);
         TransX2(1,1,Isize,Jsize) +=TempX(Isize,Jsize);   
        }
//
// Project ahead
     for (Iyear=2;Iyear<=6;Iyear++)
      {
       TempY.initialize();
         for (Isize=1;Isize<=nbins;Isize++)
         for (Jsize=1;Jsize<=nbins;Jsize++)
         for (Ksize=1;Ksize<=nbins;Ksize++)
            TempY(Isize,Jsize) +=  TransX(1,Isize,Ksize)*TempX(Ksize,Jsize);
            TempX = TempY;  
         for (Isize=1;Isize<=nbins;Isize++)
         for (Jsize=1;Jsize<=nbins;Jsize++)
         TransX2(1,Iyear,Isize,Jsize) += TempX(Isize,Jsize);   
      }
//
// Likelihood
//
 TagLike = 0;
 for (Itag=1;Itag<=NtagData;Itag++)
  {
   Isex = TagData(Itag,1);
   RelSize = TagData(Itag,2);
   RecSize = TagData(Itag,3);
   TimeAtLib = TagData(Itag,4);
//#####
// Changing weights for the tagging component
//
  if(TimeAtLib == 1)
    {wt=1.0;}  
   else if(TimeAtLib == 2)
     {wt=1.0;} 
   else if(TimeAtLib == 3)
      {wt=1.0;}  
  else if(TimeAtLib == 4)
      {wt=1.0;}  
 else if(TimeAtLib == 5)
      {wt=1.0;}   
   else
      {wt=1.0;}
//
//#####
   Total = 0.;
   for (Isize=RelSize;Isize<=nbins;Isize++)
  
    Total += selexS(Isize)*TransX2(Isex,TimeAtLib,Isize,RelSize);
   Pred = selexS(RecSize)*TransX2(Isex,TimeAtLib,RecSize,RelSize)/Total+0.00001;
   TagLike = TagLike + -1.0*wt*float(TagData(Itag,5))*log(Pred);   // wt is included
  }
 if (DoTagDiag == 1)
  {
   cout << "Doing" << endl;
   PredTagTable.initialize();
   PredProb.initialize();
   for (Itag=1;Itag<=NtagData;Itag++)
    {
     Isex = TagData(Itag,1);
     RelSize = TagData(Itag,2);
     RecSize = TagData(Itag,3);
     TimeAtLib = TagData(Itag,4);
     Preds.initialize();
     Total = 0;
     for (Isize=RelSize;Isize<=nbins;Isize++)
      Total += selexS(Isize)*TransX2(Isex,TimeAtLib,Isize,RelSize);
     for (Isize=RelSize;Isize<=nbins;Isize++)
      {
       Preds(Isize) = selexS(Isize)*TransX2(Isex,TimeAtLib,Isize,RelSize)/Total;
       PredProb(Itag,Isize) = Preds(Isize);
      } 
    PredTag(Itag,0) = length_bins(RelSize);
     PredTag(Itag,1) = Preds(RecSize);
     PredTag(Itag,2) = length_bins(RecSize);
     MeanLenRec = 0;
     for (Isize=1;Isize<=nbins;Isize++)
      MeanLenRec += length_bins(Isize)*Preds(Isize);
     PredTag(Itag,3) = MeanLenRec;
   PredTagTable(Isex,TimeAtLib,1,RecSize) += TagData(Itag,5);
   for (Isize=1;Isize<=nbins;Isize++)
      PredTagTable(Isex,TimeAtLib,2,Isize) += TagData(Itag,5)*Preds(Isize);
    }  
   Ipnt = 0;
   for (TimeAtLib=1;TimeAtLib<=6;TimeAtLib++)  
    {
     for (RelSize=1;RelSize<nbins;RelSize++)   
      {
       Preds.initialize();
       Obs.initialize();
       Cnt = 0.0;
       for (Itag=1;Itag<=NtagData;Itag++)
        if (TagData(Itag,4) == TimeAtLib & TagData(Itag,2) == RelSize)
         {
          cout << TagData(Itag) << endl;
          Cnt += float(TagData(Itag,5));
          Preds = PredProb(Itag);
          Obs(TagData(Itag,3)) += float(TagData(Itag,5));
         }
       if (Cnt > 0.)                   
        {
         Ipnt += 1;
        cout << Cnt << " " << Obs << endl;  
         Bot = 0.; Top = 0.; MeanPred = 0.; MeanObs = 0.; VarPred = 0.;
        for (RecSize=RelSize;RecSize<=nbins;RecSize++)
          {
           Obs(RecSize) /= Cnt;
           Top += Preds(RecSize)*(1.0-Preds(RecSize));
           Bot += square(Obs(RecSize) - Preds(RecSize));
           MeanPred += Preds(RecSize)*length_bins(RecSize);
           MeanObs += Obs(RecSize)*length_bins(RecSize);
           VarPred += Preds(RecSize)*length_bins(RecSize)*length_bins(RecSize);
          }
         VarPred = (VarPred-MeanPred*MeanPred)/Cnt; 
         EffN = Top/Bot; 
         TagDiag2(Ipnt,1) = TimeAtLib;
         TagDiag2(Ipnt,2) = RelSize;
         TagDiag2(Ipnt,3) = length_bins(RelSize);
         TagDiag2(Ipnt,4) = Cnt;
         TagDiag2(Ipnt,5) = MeanPred;
         TagDiag2(Ipnt,6) = MeanObs;
         TagDiag2(Ipnt,7) = sqrt(VarPred);
         TagDiag2(Ipnt,8) = EffN/Cnt;
         TagDiag2(Ipnt,9) = (MeanObs-MeanPred)/sqrt(VarPred);
        } 
      }    
     }
    NTagDiag2 = Ipnt; 
  }
//****
//
// This is the equilibrium calculation module
//
FUNCTION get_initial_equilbrium_LFQ;
  int t, l, ll;
   dvar_vector rt(1,nbins);
   dmatrix S(1,nbins,1,nbins);
   dmatrix Id=identity_matrix(1,nbins);
   dvar_vector x(1,nbins);
   dvar_matrix At(1,nbins,1,nbins);
   dvar_matrix A(1,nbins,1,nbins);
   dvar_vector rt_temp(1,nbins);
   dvar_vector rt_mean(1,nbins);  
   Neq.initialize();
   NN.initialize();
   S.initialize();
   predtotalcatchNtemp.initialize();
   predtotalcatchN.initialize();
   predretdcatchN.initialize();
   predretdcatchB.initialize();
   predtretdcatchN.initialize();
   predtotalcatchB.initialize();
   preddiscdcatchN.initialize();
   predgdiscdcatchN.initialize();
   predgdiscdcatchB.initialize();
//
// mean recruits from 1996 to 2015 estimates
//    
    rt_temp.initialize();
   rt_mean.initialize();
  for(l=1;l<=nbins;l++)
   {
    for(t = commstyr+11;t <= commendyr; t++)
    {
     rt_temp(l) +=mfexp(log_mean_rec+rec_dev(t))*rec_len(l);
    }
     rt_mean(l)=rt_temp(l)/((commendyr)-(commstyr+11)+1); 
   }
//
   rt = rt_mean;
// 
    for(int i=1;i<=nbins; i++)
     {S(i,i)=value(mfexp(-M));}
     At=trans(len_len);
     A= Id-At*S;
    Neq = solve(A,rt);
//
     totalselect=totalselect1;
     retselect=retselect1;
//
       RR.initialize();
//
//// If equilibrium pop starts in 1960
////
  for(t = commstyr-25; t <= commstyr-5; t++)
   {
        RR(t+1)=mfexp(log_mean_rec+rec_dev(t+1));
          if(t == commstyr-25)
           {
             for(l = 1; l <= nbins; l++)
             {
              NN(t,l) = Neq(l);
             }
           }
      for(l = 1; l <= nbins; l++)
       {
          for(ll = 1; ll <= l; ll++)
           {
            NN(t+1,l) += NN(t,ll)*surv1*len_len(ll,l);
           }
        NN(t+1,l) += RR(t+1)*rec_len(l);
       }
   }
////###
  for(t = commstyr-4; t <= commstyr-1; t++)
   {
        RR(t+1)=mfexp(log_mean_rec+rec_dev(t+1));
//#####
       for(l = 1; l <= nbins; l++)
        {
        z(t,l)=0.8*F_pot(t)*totalselect(l)*retselect(l)+0.2*F_pot(t)*totalselect(l)+0.65*gfmort_final(t)*gselect(l);   // updated total mortality
        predtotalcatchNtemp(t,l) = NN(t,l)*surv3(t)*(1.-exp(-z(t,l)))*F_pot(t)*totalselect(l)/z(t,l);                  // total extraction by the pot gear
        predretdcatchN(t,l) = predtotalcatchNtemp(t,l)*retselect(l);                                                   //  retention of crab   
        predtretdcatchN(t) += predretdcatchN(t,l);                                                                    // total predicted retained number of crabs for the likelihood input
        preddiscdcatchN(t,l) = m_disc*predtotalcatchNtemp(t,l)*(1.0-retselect(l));                                    // directed pot fishery discard with handling mortality applied
        predtotalcatchN(t,l) = predretdcatchN(t,l)+preddiscdcatchN(t,l);                                              // pot fishery total catch (dead)
        predtotalcatchB(t) += predtotalcatchNtemp(t,l)*retweight(l);                                                  // total extraction in weight by the pot gear for likelihood input
        predretdcatchB(t) += predretdcatchN(t,l)*retweight(l);                                                            
        predgdiscdcatchN(t,l) = m_gdisc*NN(t,l)*surv3(t)*(1.-exp(-z(t,l)))*gfmort_final(t)*gselect(l)/z(t,l);         // groundfish bycatch with average groundfish mortality applied
        predgdiscdcatchB(t) += predgdiscdcatchN(t,l)*retweight(l);                                                       
        }
        for(l = 1; l <= nbins; l++)
        {
         for(ll = 1; ll <= l; ll++)
          {
           NN(t+1,l) += (NN(t,ll)*surv1 - (predtotalcatchN(t,ll)+predgdiscdcatchN(t,ll))*surv2(t))*len_len(ll,l);
          }
         NN(t+1,l) += RR(t+1)*rec_len(l);
        }
   }
//
//#######
//
// Predicted legal male abundance in number (millions) and weight (t) by year at survey time
//
     legalabund.initialize();
     legal_biomass.initialize();
   for(t = commstyr-25; t <= commstyr-1; t++)
    {
     for(l = 8; l <= nbins; l++)        // bin number 8 starts from 136 mm CL, middle length 138 mm CL
      {
      legalabund(t) += NN(t,l);
      legal_biomass(t) += NN(t,l)*retweight(l);
      }
    }
//
// Predicted mature male (>=121 mmCL) abundance in number (millions) and weight (t) by year on next 15 February upto
// final data avialability year, ie mature abundance estimates are on 86Feb15 .... 2016Feb15
//
       matabund.initialize();
     mat_biomass.initialize();
    for(t = commstyr-25;t <= commstyr-5; t++)
    {
        for(l = 5; l <= nbins; l++)        // bin number 5 starts from 121 mm CL (50% maturity length)
      {
       matabund(t) += NN(t,l)*surv4;
       mat_biomass(t) += NN(t,l)*surv4*retweight(l);
      }
    }
//
    for(t = commstyr-4;t <= commstyr-1; t++)
    {
        for(l = 5; l <= nbins; l++)        // bin number 5 starts from 121 mm CL (50% maturity length)
      {
       matabund(t) += NN(t,l)*surv4-(predtotalcatchN(t,l)+predgdiscdcatchN(t,l))*surv5(t);
       mat_biomass(t) += (NN(t,l)*surv4-(predtotalcatchN(t,l)+predgdiscdcatchN(t,l))*surv5(t))*retweight(l);
      }
     }
//#######   
//*******************************************   
//
FUNCTION get_population_number
  int t, i, l, j, ll;
//
////#############################################################Revised
   newsh.initialize();                  
   fpen2=0.0;
  for(t = commstyr; t <= commendyr; t++)
  {
  predtotalcatchB(t) =0.;
  predretdcatchB(t) = 0.; 
  predgdiscdcatchB(t) = 0.;
  legalcatch(t)=0.;
  for(l = 1; l <= nbins; l++)
  { 
    predtotalcatchNtemp(t,l) = 0.;                  
    predretdcatchN(t,l) = 0.;                                                      
    preddiscdcatchN(t,l) = 0.;                                        
    predtotalcatchN(t,l) = 0.;                                                  
    predgdiscdcatchN(t,l) = 0.;          
   }                                                   
  }
// 
//
// Population dynamics
//
  for(t = commstyr; t <= commendyr; t++)
 {
    if(t<2005)
    {totalselect=totalselect1;
     retselect=retselect1;}
    else
    {totalselect=totalselect2;
     retselect=retselect2;}
   RR(t+1)=mfexp(log_mean_rec+rec_dev(t+1));
  if(t==commendyr)RR(commendyr+1)=(RR(commendyr-2)+ RR(commendyr-1)+RR(commendyr))/3.;
//
   if(t == commstyr)
   {
   for(l = 1; l <= nbins; l++)
    {
        newsh(t,l) = NN(t,l);      // equilirium estimate of initial stock size is input here (start year 1985)
    }
   }
//
// Predicted retained and discard catch
//
   for(l = 1; l <= nbins; l++)
    {
    z(t,l)=0.8*F_pot(t)*totalselect(l)*retselect(l)+0.2*F_pot(t)*totalselect(l)+0.65*gfmort_final(t)*gselect(l);      // updated total mortality
    predtotalcatchNtemp(t,l) = newsh(t,l)*surv3(t)*(1.-exp(-z(t,l)))*F_pot(t)*totalselect(l)/z(t,l);                  // total extraction by the pot gear
    predretdcatchN(t,l) = predtotalcatchNtemp(t,l)*retselect(l);                                                      //  retention of crab   
    preddiscdcatchN(t,l) = m_disc*predtotalcatchNtemp(t,l)*(1.0-retselect(l));                                        // directed pot fishery discard with handling mortality applied
    predtotalcatchN(t,l) = predretdcatchN(t,l)+preddiscdcatchN(t,l);                                                  // pot fishery total catch (dead)
    predtotalcatchB(t) += predtotalcatchNtemp(t,l)*retweight(l);                                                      // total extraction in weight by the pot gear for likelihood input
    predretdcatchB(t) += predretdcatchN(t,l)*retweight(l);                                                            
    predgdiscdcatchN(t,l) = m_gdisc*newsh(t,l)*surv3(t)*(1.-exp(-z(t,l)))*gfmort_final(t)*gselect(l)/z(t,l);          // groundfish bycatch with average groundfish mortality applied
    predgdiscdcatchB(t) += predgdiscdcatchN(t,l)*retweight(l);                                                       
   if(l>7)legalcatch(t) += predtotalcatchN(t,l);
   }

//
     for(l = 1; l <= nbins; l++)
     {
    for(ll = 1; ll <= l; ll++)
       {
       dvariable yyy = (newsh(t,ll)*surv1 - (predtotalcatchN(t,ll)+predgdiscdcatchN(t,ll))*surv2(t))*len_len(ll,l);
      newsh(t+1,l) += posfun(yyy,0.000001,fpen2);                                                                     // posfunction flag
       }
    newsh(t+1,l) += RR(t+1)*rec_len(l);
     }
 }
FUNCTION get_cpue
//
//  Predicted total, retained and discard CPUE in the fishery by length and year
//
  int t; int l;
  dvariable xxx;
   predtotalcpue.initialize();
   predretcpue.initialize();
   fpen1=0.0;
   for(t = commstyr; t <= commendyr; t++)
  {
     if(t <= 2004)
    {qtemp = logq2;
     totalselect=totalselect1;
     retselect=retselect1;}
    else
   {qtemp = logq3;
    totalselect=totalselect2;
     retselect=retselect2;}
   for(l = 1; l <= nbins; l++)
     {
   xxx = mfexp(qtemp)*totalselect(l)*(newsh(t,l) - 0.5*(predtotalcatchN(t,l)+predgdiscdcatchN(t,l)))*surv3(t); // number of crabs
   predtotalcpue(t,l)= posfun(xxx,0.000001,fpen1);                                                             // posfunction flag
   predretcpue(t,l)   =retselect(l)*predtotalcpue(t,l);
     }
  }
     predttotalcpue.initialize();
     predtretcpue.initialize();
   for(t = commstyr; t <= commendyr; t++)
   {
    for(l = 1; l <= nbins; l++)
     {
      predttotalcpue(t) +=predtotalcpue(t,l);
      predtretcpue(t) +=predretcpue(t,l);
     }
   }
//
//  Predicted annual retained cpue for 1991 to 2015
//
       for(t = commstyr+6; t <= commendyr; t++)
        {predtretcpue91_15(t)=predtretcpue(t);}
//
FUNCTION get_lengthcomp
 int t; int l;
//
//   Predicted retained and total catch length composition by year in the fishery
//
      retcatchsum.initialize();
      totalcatchsum.initialize();
      gdiscdcatchsum.initialize();
      retcatchlen.initialize();
      totalcatchlen.initialize();
      gdiscdcatchlen.initialize();
//
    for(t = commstyr; t <= commendyr; t++)
    {
       for(l = 1; l <= nbins; l++)       
      {
       retcatchlen(t,l) = predretdcatchN(t,l);
       retcatchsum(t)  += predretdcatchN(t,l);
       totalcatchlen(t,l) = predtotalcatchNtemp(t,l); 
       totalcatchsum(t) +=  predtotalcatchNtemp(t,l); 
       gdiscdcatchlen(t,l) = predgdiscdcatchN(t,l);
       gdiscdcatchsum(t) +=  predgdiscdcatchN(t,l);
      }
    }
      predretcatchlen.initialize();
      predtotalcatchlen.initialize();
      predgdiscdcatchlen.initialize();
      retsigmasquare.initialize();
   for(t = commstyr; t <= commendyr; t++)
    {
        for(l = 1; l <= nbins; l++)   
        {
        if(retcatchsum(t)>0.)predretcatchlen(t,l) = retcatchlen(t,l)/retcatchsum(t);
        if(totalcatchsum(t)>0.)predtotalcatchlen(t,l) = totalcatchlen(t,l)/totalcatchsum(t);
        if(gdiscdcatchsum(t)>0.)predgdiscdcatchlen(t,l) = gdiscdcatchlen(t,l)/gdiscdcatchsum(t);
        }
    }
//
//   Length composition variance estimation using multinomial formula for retained catch:
//
     retsigmasquare.initialize();
     for(t = commstyr; t <= commendyr; t++)
    {
        for(l = 1; l <= nbins; l++)   
        {
        retsigmasquare(t,l)= ((1.0-obsretlencomp(t,l))*obsretlencomp(t,l)+0.1/17.0)/(retcatcheffsample(t));
        }
    }
//
//   Length composition variance estimation using multinomial formula for groundfish discard:
//
         gdiscdsigmasquare.initialize();
    for(t = commstyr+4; t <= commendyr; t++)
     {
       if(t==1993)continue;
    for(l = 1; l <= nbins; l++)
      {
       gdiscdsigmasquare(t,l)= ((1.0-obsgdiscdlencomp(t,l))*obsgdiscdlencomp(t,l)+0.1/17.0)/(gdiscdcatcheffsample(t));
      }
     }
//
//   Length composition variance estimation using multinomial formula for total catch:
//
     totalsigmasquare.initialize();
     for(t = commstyr+5; t <= commendyr; t++)
    {
          for(l = 1; l <= nbins; l++)
       {
        totalsigmasquare(t,l)= ((1.0-obstotallencomp(t,l))*obstotallencomp(t,l)+0.1/17.0)/(totalcatcheffsample(t));
       }
    }
//
// standardized residuals of length compositions of retained catch for bubble plot
//
     ret_stdresid.initialize();
     for(t = commstyr; t <= commendyr; t++)
    {
        for(l = 1; l <= nbins; l++)
        {
          ret_stdresid(t,l)=(obsretlencomp(t,l)-predretcatchlen(t,l))/sqrt(2.*retsigmasquare(t,l));

     }
    }
//
// standardized residuals of length compositions of groundfish discard catch for bubble plot
//
    gdiscd_stdresid.initialize();
  for(t = commstyr+4; t <= commendyr; t++)
     {
       if(t==1993)continue;
    for(l = 1; l <= nbins; l++)
      {
        gdiscd_stdresid(t,l)=(obsgdiscdlencomp(t,l)-predgdiscdcatchlen(t,l))/sqrt(2.*gdiscdsigmasquare(t,l));
     }
    }
//
// standardized residuals of length compositions of total catch for bubble plot
//
    total_stdresid.initialize();
     for(t = commstyr+5; t <= commendyr; t++)
    {
         for(l = 1; l <= nbins; l++)
       {
          total_stdresid(t,l)=(obstotallencomp(t,l)-predtotalcatchlen(t,l))/sqrt(2.*totalsigmasquare(t,l));

       }
    }
//
//########################
FUNCTION objective_function
  int t;  int j; int i;
//
//Robust liklihood approach for length composition fit
//
 like_retlencomp=0.;
  for(t = commstyr; t <= commendyr; t++)
  {
 for(j = 6; j <= nbins; j++)  // j=6 start from 128 since mostly zero entries before this size bin for retained catch
    {
  like_retlencomp  += like_wght(1)*(0.5*log(6.286*retsigmasquare(t,j))-(log(mfexp(-1.*square(obsretlencomp(t,j)-predretcatchlen(t,j))/(2.*retsigmasquare(t,j)))+0.01)));
    }
  }
//
  like_totallencomp=0.;
   for(t = commstyr+5; t <= commendyr; t++)
  {
 for(j = 1; j <= nbins; j++)
   {
  like_totallencomp += like_wght(2)*(0.5*log(6.286*totalsigmasquare(t,j))-(log(mfexp(-1.*square(obstotallencomp(t,j)-predtotalcatchlen(t,j))/(2.*totalsigmasquare(t,j)))+0.01)));
   }
  }
//
// Ignored groundfish length composition in the fit
//
  like_gdiscdlencomp=0.;
//    for(t = commstyr+4; t <= commendyr; t++)
//  {
//    if(t== 1993)continue;          // no effective sample size in 1993
//    for(j = 1; j <= nbins; j++)      
//   {
//  like_gdiscdlencomp += like_wght(3)*(0.5*log(6.286*gdiscdsigmasquare(t,j))-(log(mfexp(-1.*square(obsgdiscdlencomp(t,j)-predgdiscdcatchlen(t,j))/(2.*gdiscdsigmasquare(t,j)))+0.01)));
//   }
//  }
//
// CPUE likelihood
//
   like_retcpue=0.;
   for(t = commstyr+6; t <= commendyr; t++)    // 1991 to 2015
    {
     like_retcpue += like_wght(4)*(0.5*(log(6.286*(obslegalcpueindex_var(t)+prelegal_var))+square(log(obslegalcpueindex(t)+0.001)-log(predtretcpue(t)+0.001))/(obslegalcpueindex_var(t)+prelegal_var)));
    }
//
//  Catch and bycatch biomass likelihood
//
   like_retdcatchB=0.;
   for(t = commstyr; t <= commendyr; t++)
  {
   like_retdcatchB += like_wght(5)*square(log(obsretdcatchB(t)+0.001)-log(predretdcatchB(t)+0.001));
  }
//
   like_totalcatchB=0.;
   for(t = commstyr+5; t <= commendyr; t++)               
  {
     like_totalcatchB += totalCBWght(t)*square(log(obstotalcatchB(t)+0.001)-log(predtotalcatchB(t)+0.001));  // total catch biomass is weighted by scaled number of pots
  }
//
   like_gdiscdcatchB=0.;
   for(t = commstyr+4; t <= commendyr; t++)
  {
    if(t==1993)continue;                  // no data
   like_gdiscdcatchB += like_wght(7)*square(log(obsgdiscdcatchB(t)+0.001)-log(predgdiscdcatchB(t)+0.001));
   }
//
// Recruit deviation penalty
//
    like_rec_dev=0.;
  if(active(rec_dev))
 {
  for(t = commstyr-24; t <= commendyr+1; t++)
     {
      like_rec_dev += like_wght(8)*square(rec_dev(t));
     }
  }
//
//
// Pot F dev penalty
//
  like_F=0.;
  if(active(Fpot_dev))
  {
   if(current_phase()>=sel_phase)
    {
     for(t = commstyr-4; t <= commendyr; t++)
       {
        like_F += 0.001*square(Fpot_dev(t));         
       }
    }
   else
    {
     for(t = commstyr-4; t <= commendyr; t++)
       {
        like_F += like_wght(9)*square(Fpot_dev(t));
       }
    }
  }
//
// Groundfish bycatch F dev penalty
//
   like_gF=0.;
  if(active(Fground_dev))
  {
    if(current_phase()>=sel_phase)
     {
          for(t = commstyr+4; t <= commendyr; t++)
       { 
         like_gF += 0.001*square(Fground_dev(t));     
       }
     }
  else
    {
     for(t = commstyr+4; t <= commendyr; t++)
     {
      
     like_gF += like_wght(10)*square(Fground_dev(t));
     }
    }
  }
//
     like_Mpenalty= 0.;              // likelihood not implemented
//    like_Mpenalty= like_wght(11)*square(log(M)-log(0.18));    
//
// Tagging likelihood
//
    like_LLyr=like_wght(12)*TagLike;
//
//
   like_finalF =0.;
  if(last_phase())
   {like_finalF = square(Ftemp-final_F);}
//
  like_fpen= 0.;
  like_fpen = 1000.*(fpen1+fpen2);
//
// 1981 to 1984 observed vs predicted catch numbers
// 
  like_catch1=0.;
  for(t = commstyr-4; t <= commstyr-1; t++)
   {like_catch1 += like_wght(13)*square(predtretdcatchN(t)-retdcatch1(t));}
//
  f = like_retlencomp+like_totallencomp+like_gdiscdlencomp+like_retcpue+like_retdcatchB+like_totalcatchB
    +like_gdiscdcatchB+like_rec_dev+like_F+like_gF+like_Mpenalty+like_LLyr+like_finalF+like_catch1;
//
    if(like_fpen==0. && last_phase())exit(1);
     f +=like_fpen;
//
//
//
REPORT_SECTION
 cout<<endl<<endl<<"Completed phase: "<<current_phase()<<endl<<endl<<endl;
//
  report<<"predicted pot fishery annual total harvest rate"<<endl;
  report<<harvestrate<<endl;
//
  report<<"predicted annual total cpue"<<endl;
  report<<predttotalcpue<<endl;
//
  report<<"observed annual retained cpue index"<<endl;
  report<<obslegalcpueindex<<endl;
  report<<"predicted annual retained cpue index"<<endl;
  report<<predtretcpue91_15<<endl;
//
  report<<"observed retained catch proportion by length and year"<<endl;
  report<<obsretlencomp<<endl;
  report<<"Predicted retained catch proprtion by length and year"<<endl;
  report<<predretcatchlen<<endl;
//
  report<<"observed total catch proportion by length and year"<<endl;
  report<<obstotallencomp<<endl;
  report<<"Predicted total catch proportion by length and year"<<endl;
  report<<predtotalcatchlen<<endl;
//
  report<<"observed groundfish discard catch proportion by length and year"<<endl;
  report<<obsgdiscdlencomp<<endl;
  report<<"Predicted groundfish discard catch proportion by length and year"<<endl;
  report<<predgdiscdcatchlen<<endl;
//
  report<<"predicted annual recruitment"<<endl;
  report<<mfexp(log_mean_rec+rec_dev)<<endl;
//
  report<<"predicted annual pot fishery total F"<<endl;
  report<<mfexp(log_mean_Fpot+Fpot_dev)<<endl;
//
  report<<"predicted annual groundfish fishery F"<<endl;
  report<<gfmort_final<<endl;
//
  report<<"predicted annual legal male abundance in millions 1960 to 2015."<<endl;
  report<<legalabund<<endl;
//
  report<<"predicted annual legal male biomass in t from 1960 to 2015"<<endl;
  report<<legal_biomass<<endl;
//
  report<<"predicted annual mature male abundance in millions 1961Feb to 2016Feb."<<endl;
  report<<matabund<<endl;
//
  report<<"predicted annual mature male biomass in t 1961Feb to 2016Feb."<<endl;
  report<<mat_biomass<<endl;
//
  report<<"predicted mean mature biomass based on 1986 Feb to 2016 Feb mature biomass"<<endl;
  report<<meanmatbiomass<<endl;
//
  report<<"predicted Full Selection F based on 1986 Feb to 2016 Feb mean mature biomass in the CR"<<endl;
  report<<final_F<<endl;
//
  report<<"Predicted 2016/17 retained harvest in t based on 1986 Feb to 2016 Feb mean mature biomass in the CR"<<endl;
  report<<retain_harvest<<endl;
//
  report<<"Predicted 2016/17 pot fishery discard removal in for hm=0.2"<<endl;
  report<<preddiscdcatchBB<<endl;
//
  report<<"Predicted 2016/17 groundfish fishery discard removal in t for hm=0.65"<<endl;
  report<<predgdiscdcatchBB<<endl;
//
  report<<"predicted total selectivity during pre rationalization"<<endl;  
  report<<totalselect1<<endl;
//
  report<<"predicted retention curve during pre rationalization"<<endl;  
  report<<retselect1<<endl;
//
  report<<"predicted total selectivity during post rationalization"<<endl;
  report<<totalselect2<<endl;
//
  report<<"predicted retention curve during post rationalization"<<endl;  
  report<<retselect2<<endl;   
//
  report<<"predicted groundfish fishery selectivity throught all years"<<endl;
  report<<gselect<<endl;
//
  report<<"predicted size transition matrix"<<endl;
  report<<len_len<<endl;
//
  report<<"predicted recruitment distribution"<<endl;
  report<<rec_len<<endl;
//
// Negtive Log Likelihood values:
//
  report<<"like_retlencomp"<<"   "<<like_retlencomp<<endl;
//
  report<<"like_totallencomp"<<"   "<<like_totallencomp<<endl;
//
  report<<"like_gdiscdlencomp"<<"   "<<like_gdiscdlencomp<<endl;
//
  report<<"like_retcpue"<<"   "<<like_retcpue<<endl;
//
  report<<"like_retdcatchB"<<"   "<<like_retdcatchB<<endl;
//
  report<<"like_totalcatchB"<<"   "<<like_totalcatchB<<endl;
//
  report<<"like_gdiscdcatchB"<<"   "<<like_gdiscdcatchB<<endl;
//
  report<<"like_rec_dev"<<"   "<<like_rec_dev<<endl;
//
  report<<"like_F"<<"   "<<like_F<<endl;
//
  report<<"like_gF"<<"   "<<like_gF<<endl;
//
  report<<"like_LLyr "<<"   "<<like_LLyr<<endl;
//
  report<<"like_finalF"<<"   "<<like_finalF<<endl;
//
  report<<"like_fpen"<<"   "<<like_fpen<<endl;
//
  report<<"sum of Ret&DiscardLFQ Negative log likelihood"<<"   "<<like_retlencomp+like_totallencomp+like_gdiscdlencomp<<endl;
//
  report<<"Total Negative log_likelihood"<<"   "<<f<<endl;
//
  report<<"newshell abundance by size by year"<<endl;
  report<<newsh<<endl;
//
  report<<"predicted equilibrium initial size composition"<<endl;
  report<<Neq<<endl;
//
  report<<"predicted retained catch t by year"<<endl;
  report<<predretdcatchB<<endl;
//
  report<<"observed retained catch t by year"<<endl;
  report<<obsretdcatchB<<endl;
//
  report<<"predicted pot total catch t by year"<<endl;
  report<<predtotalcatchB<<endl;
//
  report<<"observed pot total catch t by year"<<endl;
  report<<obstotalcatchB<<endl;
//
  report<<"predicted groundfish discard catch t by year"<<endl;
  report<<predgdiscdcatchB<<endl;
//
  report<<"observed groundfish discard catch t by year"<<endl;
  report<<obsgdiscdcatchB<<endl;
//
  report<<"Retained standardized residuals"<<endl;
  report<<ret_stdresid<<endl;
//
  report<<"Pot total standardized residuals"<<endl;
  report<<total_stdresid<<endl;
//
  report<<"Groundfish Discard standardized residuals"<<endl;
  report<<gdiscd_stdresid<<endl;
//
  report<< "predicted molt probability  "<< endl;
  report<<molt<<endl;
//
  report<<"observed retained catch effective sample"<<endl;
  report<<retcatcheffsample<<endl;
//
  report<<"predicted retained catch effective sample"<<endl;
  report<<predretcatcheffsample<<endl;
//
  report<<"observed total catch effective sample"<<endl;
  report<<totalcatcheffsample<<endl;
//
  report<<"predicted total catch effective sample"<<endl;
  report<<predtotalcatcheffsample<<endl;
//
  report<<"observed groudfish discard catch effective sample"<<endl;
  report<<gdiscdcatcheffsample<<endl;
//
  report<<"predicted groudfish discard effective sample"<<endl;
  report<<predgdiscdcatcheffsample<<endl;
// 
   report<<"NB Predicted Observer Legal CPUE"<<endl;
   report<<NBpredictedobslegalcpue<<endl; 
//
   report<<"1981 to 1984 catch likelihood"<<endl;
   report<<like_catch1 <<endl; 
//
   report<<"observed 1981 to 1984 catch numbers "<<endl;
   report<<retdcatch1 <<endl; 
//
//
    cout<<endl<<endl;
    int II;
    report<<" Yr, RelL, RelS, Cnt, MeanPrd, MeanObs, sq(VarPrd), Effn/Cnt,  (Obs-Prd)/sq(Var)"<<endl;
    for (II=1;II<=NTagDiag2;II++)
      report<<TagDiag2(II) <<endl;
//
  report<<"Table output  "<<endl;
  report<< PredTagTable<<endl; 
//
// Another report file for R
//
  ofstream report1("Adak8515Sc1Base.out");
  report1<<obsretlencomp<<endl;     //observed retained catch proportion by length and year
  report1<<predretcatchlen<<endl;   //Predicted retained catch proprtion by length and year
  report1<<obstotallencomp<<endl;   // observed total catch proportion by length and year
  report1<<predtotalcatchlen<<endl;  //Predicted total catch proportion by length and year
  report1<<obsgdiscdlencomp<<endl;   //observed groundfish discard catch proportion by length and year
  report1<<predgdiscdcatchlen<<endl;  //Predicted groundfish discard catch proportion by length and year
//
//
   report<<"likelihood M "<<endl;
   report<<like_Mpenalty <<endl; 
////
  report<<"predicted retained catch in millions from 1981...."<<endl;
   report<<predtretdcatchN<<endl;
//
//
GLOBALS_SECTION
//
RUNTIME_SECTION
//one number for each phase, if more phases then uses the last number
  maximum_function_evaluations 50000,50000,50000,60000,1000000,100000,100000
  convergence_criteria 1e-5,1e-5,1e-5,1e-6,1e-6,1e-6,1e-6
//
TOP_OF_MAIN_SECTION
  arrmblsize = 10000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(10000000); // this may be incorrect in
// the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000000);
//
