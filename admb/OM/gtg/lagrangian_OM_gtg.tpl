//  ******************************************************************
//  Lagragian movement model
//  
//  Created by Catarina Wor on 2015-04-13.
//  Copyright (c) 2015. All rights reserved.
//  Comments:
//  ******************************************************************


DATA_SECTION

	// Read new seed in
	int seed;

	LOC_CALCS
		ifstream ifs( "seed.txt" ); // if this file is available
		ifs>>seed; //read in the seed
		seed += 1; // add 1 to the seed
		ofstream ofs( "seed.txt" ); //put out to seed.txt
		ofs<<seed<<endl; //the new value of the seed
		cout<<"seed"<<seed<<endl;
	END_CALCS

	//======================
	// model dimensions
	//======================

	init_int syr;
	init_int nyr;
	init_int rep_yr;
	init_int proj_yr;
	init_int sage;
	init_int nage;
	init_int smon;
	init_int nmon;
	init_int sarea;
	init_int narea;

	init_int ngroup;

	init_int nlen;
	init_int inilen;
	init_int lenstp;


	init_int nations;
	init_vector border(1,nations-1);

	init_int fisharea;
	init_vector fishbound(1,fisharea-1);


	//======================
	//model parameters -- fixed
	//======================
	init_number Ro;
	init_number h;
	init_number m;
	init_number fe;
	init_number q;
	init_number fbeta;
	init_number sigR;		
	
	init_number tau_c; 
	init_number tau_l; 
	init_number tau_survey;		
	init_number tau_surveyl;		
	init_number err;

	//======================
	//length parameters -- fixed
	//======================
	init_vector Linf(1,ngroup);
	init_number vbk;
	init_number to;
	
	init_number cvl;
	init_number awt;
	init_number bwt;



	//======================
	//Survey input 
	//======================
	init_int surv_mon;
	init_int surv_nobs;
	init_vector surv_yrs(1,surv_nobs);


	//=================================
	//parameters in estimation model 
	//
	//=================================
	init_number mo; 		
	init_number maxPos50;	
	init_number maxPossd;
	init_number cvPos;				
	//init_vector maxPossd(1,ngroup);			
	//init_vector cvPos(1,ngroup); 

	//init_vector wa(sage,nage);
	init_vector fa(sage,nage);
	init_vector va(sage,nage);
	init_vector minPos(sage,nage);


	init_number Fmult;
	init_matrix pTotEffyear(1,fisharea,syr,nyr);
	init_matrix TotEffmonth(1,fisharea,smon,nmon);

	init_vector effPwr(sarea,narea);

	
	

	init_number rhoB;
	

	init_int eof;
	
	

	LOC_CALCS
		
		if( eof != 999 )
		{
			cout<<"pTotEffyear "<<pTotEffyear<<endl;
			cout<<"effPwr "<<effPwr<<endl;
			cout<<"Error reading data.\n Fix it."<<endl;
			cout<< "eof is: "<<eof<<endl;
			ad_exit(1);
		}

	END_CALCS


	!! ad_comm::change_datafile_name("wt.dat");
	

	init_vector wt(rep_yr+1,proj_yr);


	init_int eoffw;

	LOC_CALCS
		
		if( eoffw != 999 )
		{
			cout<<"wt "<<wt<<endl;
			cout<<"Error reading recruitment time series.\n Fix it."<<endl;
			cout<< "eoffw is: "<<eoffw<<endl;
			ad_exit(1);
		}

	END_CALCS
	




	!! ad_comm::change_datafile_name("HCR.dat");
	init_int satype;

	init_vector nationTACprop(1,nations);

	init_int hcr;

	init_number slope_hcr;

	init_number intercept_hcr;

	init_int eoff;

	LOC_CALCS
		
		if( eoff != 999 )
		{
			cout<<"slope_hcr "<<slope_hcr<<endl;
			cout<<"Error reading HCR controls.\n Fix it."<<endl;
			cout<< "eoff is: "<<eoff<<endl;
			ad_exit(1);
		}

	END_CALCS
	

	!! ad_comm::change_datafile_name("catlim_hist.dat");
	
	init_matrix catlim_hist(rep_yr+1,nyr,1,nations);


	init_int eofct;

	LOC_CALCS
		
		if( eofct != 999 )
		{
			cout<<"catlim_hist "<<catlim_hist<<endl;
			cout<<"Error reading HCR controls.\n Fix it."<<endl;
			cout<< "eofct is: "<<eofct<<endl;
			ad_exit(1);
		}

	END_CALCS
	



	//===================================
	// accessory quantities and counters
	//===================================
	
	int nrep_yr;
	int nproj_yr;
	int ntstp;
	int ststp;
	int ptstp;
	int tmon;

   	vector age(sage,nage);
   	vector len(1,nlen);
    vector areas(sarea,narea);
   	ivector fishingr(1,fisharea);
   	ivector nationareas(1,nations);
   	
   	vector epsilon(1,surv_nobs);
	//vector wt(syr,proj_yr);
	vector gt(nyr,proj_yr);
	vector ot(nyr,proj_yr);
	vector vt(syr,proj_yr);
   	vector wx(syr,proj_yr);


   	int n_rg;
   	!!n_rg  = (narea-sarea+1) * ngroup;
   	
   	ivector   n_area(1,n_rg);
	ivector  n_group(1,n_rg);

   	imatrix  pntr_rg(sarea,narea,1,ngroup);
   		
   	LOC_CALCS	


		int ih,r,g;
		ih = 0;

		for(r=sarea; r<=narea; r++)
		{
			for(g=1; g<=ngroup; g++)
			{
				ih ++;
				pntr_rg(r,g) = ih;
				n_area(ih)  = r;
				n_group(ih) = g;
			}
		}
		
		nrep_yr = rep_yr-syr+1;
       	nproj_yr = proj_yr-syr+1;

		tmon = nmon-smon+1;
		ststp =	tmon * nrep_yr;
		ntstp = tmon * (nyr-syr+1);
		ptstp = tmon * (proj_yr-syr+1);
		
		age.fill_seqadd(sage,1);
		len.fill_seqadd(inilen,lenstp);
		areas.fill_seqadd(sarea,1);
		nationareas.initialize();
		fishingr.initialize();

		dvector natmp1(1,fisharea);
		dvector natmp2(1,nations);
			
		
		//fishingr
		natmp1(1)=sarea;

		for(int n=1; n<=fisharea-1; n++)
		{
			natmp1(n+1)=fishbound(n);

			for(int a=sarea;a<=narea;a++)
			{
				if(areas(a)>=natmp1(n)&areas(a)<fishbound(n))
				{
					fishingr(n)++;
				}
			}
		}

		fishingr(fisharea)= narea-sarea+1 - sum(fishingr(1,fisharea-1));

		//nationareas
		natmp2(1)=sarea;

		for(int n=1; n<=nations-1; n++)
		{
			natmp2(n+1)=border(n);

			for(int a=sarea;a<=narea;a++)
			{
				if(areas(a)>=natmp2(n)&areas(a)<border(n))
				{
					nationareas(n)++;
				}
			}
		}

		nationareas(nations)=narea-sarea+1 - sum(nationareas(1,nations-1));
			
		random_number_generator rng(seed);
		//wt.fill_randn(rng);
		//wt*=sigR;

		vt.fill_randn(rng);
		vt*=0.1;

		wx.fill_randn(rng);
		wx*=0.08;

		gt.fill_randn(rng);
		gt*=0.3;

		epsilon.fill_randn(rng);
		epsilon*=0.1;

		ot(nyr)=gt(nyr);
		for(int i=nyr+1; i<=proj_yr; i++){
			ot(i)=rhoB*ot(i-1)+gt(i);

		}
		
	END_CALCS

		
	 	ivector indyr(1,ptstp);
 		ivector indmonth(1,ptstp);
		ivector indnatarea(sarea,narea);
		ivector indfisharea(sarea,narea);
		ivector pcat(1,fisharea);
		vector ntmp(1,fisharea+1);
		vector ntmp1(1,nations+1);

		int tot_pcat;

		imatrix pmat(1,fisharea,1,ptstp);
		matrix TotEffyear(1,fisharea,syr,nyr);
		//matrix TotEffyear_rep(1,fisharea,rep_yr+1,proj_yr);

		number delta;
		vector Xini(1,ngroup);
		
	

       LOC_CALCS

      		int aa =0;
       			
       		for(int y=syr;y<=proj_yr;y++)		
       		{
       			for(int ii=smon;ii<=nmon;ii++)
       			{
       				aa++;
       				indyr(aa) = y;
       				indmonth(aa) = ii;
       			}
       		} 

       		for(int n=1;n<=fisharea;n++)
       		{
       			//TotEffyear(n)(syr,nyr) = Fmult* pTotEffyear(n)(syr,nyr);
       			TotEffyear(n)(syr,nyr) =Fmult* pTotEffyear(n)(syr,nyr); //* exp(vt(syr,nyr))
       		}
       			
   			

       		//calc for indfisharea
       		ntmp(1) = sarea;
       			
       		for(int n=1;n<=fisharea;n++)
       		{
       			ntmp(n+1)= ntmp(n)+fishingr(n);	
       			indfisharea(ntmp(n),ntmp(n+1)-1)=n;
       		}

       		ntmp(fisharea+1) = narea;
       		indfisharea(narea) = fisharea;

       			//cout<<"Ok after indfisharea calcs"<< endl;

       		//calc for indnatarea
       		ntmp1(1) = sarea;
       			
       		for(int n=1;n<=nations;n++)
       		{
       			ntmp1(n+1)= ntmp1(n)+nationareas(n);	
       			indnatarea(ntmp1(n),ntmp1(n+1)-1)=n;
       		}

       		ntmp1(nations+1) = narea;
       		indnatarea(narea) = nations;
       			

       		pcat.initialize();
       		pmat.initialize();
       		for(int n=1;n<=fisharea;n++)
       		{
       			//for(int ii=rep_yr+1;ii<=nyr;ii++)
       			//{
       			//	TotEffyear_rep(n)(ii)= TotEffyear(n)(ii);
       			//}

       			//for(int ia=nyr+1;ia<=proj_yr;ia++)
       			//{
       			//	TotEffyear_rep(n)(ia)= TotEffyear(n)(nyr);
       			//}

       			for(int i=nrep_yr*tmon+1;i<=ptstp;i++)
       			{
       				if(TotEffmonth(n)(indmonth(i))>0.0)
       				{
       					pcat(n)++;
       					pmat(n)(i)=1;
       				} 							
       			}
       		}

   				
       		tot_pcat=sum(pcat);  

       		//calc_propg
       		delta = 2.0*3.0/ngroup;	

       		for(int g=1; g<=ngroup ;g++)
			{
				Xini(g) = delta*((g)-ngroup/2.0);
			}
			//cout<<"Xini is "<<Xini<<endl;	

			ofstream nfs("catlim.txt");
			nfs<<"#catlim" << endl << "0 0" <<endl;
			
       			
	END_CALCS

	


PARAMETER_SECTION

	objective_function_value no_f;

	//init_number nop;

	//derived quantities
	number kappa;
	number phiE;
	number So;
	number SBo;
	number beta;
	number Bo;


	number m_tsp;
	
	vector lxo(sage,nage);
	vector wa(sage,nage);
	vector za(sage,nage);
	vector SB(1,ptstp);
	vector tB(1,ptstp);
	vector survB(1,surv_nobs);
	
	vector maxPos(sage,nage);
	vector varPos(sage,nage);
	vector varPosg(sage,nage);
	vector prop_ng(1,ngroup);
	vector yYieldtotal(syr,proj_yr);

	vector ytB(syr,proj_yr);
	vector ySB(syr,proj_yr);
	vector Yieldtotal(1,ptstp);

	vector totcatchLen(rep_yr+1,proj_yr);
	vector Ut(rep_yr+1,proj_yr);
	//vector UtLen(rep_yr+1,proj_yr);
	vector spr(syr,proj_yr);
	vector phie(syr,proj_yr);
	vector msy(syr,proj_yr);
	vector fmsy(syr,proj_yr);
	vector ftarget(nyr,proj_yr);
	vector seen_bt(nyr,proj_yr);
	vector seen_sb(nyr,proj_yr);

	vector catlim(1,nations);


	//matrix Catage(1,ntstp,sage,nage);
	matrix Catlen(1,ptstp,1,nlen);
	//matrix Fatage(1,ntstp,sage,nage);
 	//matrix meanPosX(1,ntstp,sage-1,nage);
 	matrix tVBarea(1,ptstp,sarea,narea);
	matrix totVBnation(1,ptstp,1,nations);
	matrix VulBage(1,ptstp,sage,nage);
	matrix Effarea(1,ptstp,sarea,narea);
 	matrix ywa(syr,proj_yr,sage,nage);

 	matrix YieldNat(1,ptstp,1,nations);
	matrix yYieldtotalage(syr,proj_yr,sage,nage);
 	matrix Watage_comm(rep_yr+1,proj_yr,sage,nage);
 	//matrix tot_comm_obsCatage(rep_yr+1,proj_yr,sage,nage);
 	matrix comm_obsCatage(rep_yr+1,proj_yr,sage,nage);
 	matrix comm_obsCatLen(rep_yr+1,proj_yr,1,nlen);
 	matrix surv_obsCatage(1,surv_nobs,sage,nage);
 	//matrix surv_obsCatLen(1,surv_nobs,1,nlen);
 	matrix yNage(syr,proj_yr,sage,nage);
 	matrix yFatage(syr,proj_yr,sage,nage);
 	matrix seltotal(syr,proj_yr,sage,nage);
 	matrix yCatchtotalage(syr,proj_yr,sage,nage);
 	matrix yCatchtotalLen(syr,proj_yr,1,nlen);
 	matrix yYieldNat(syr,proj_yr,1,nations);
 	matrix histTAC(nyr,proj_yr,1,nations);

 	3darray selfisharea(syr,proj_yr,1,fisharea,sage-2,nage);
 	3darray selnation(syr,proj_yr,1,nations,sage-2,nage); 	
 		
	3darray Nage(1,ngroup,1,ptstp,sage-2,nage); 
 	3darray VulB(1,ngroup,1,ptstp,sage-2,nage);
 	3darray Lage(1,ngroup,smon,nmon,sage-2,nage); 
 	3darray Wage(1,ngroup,smon,nmon,sage-2,nage); 
 	//3darray CatageG(1,ngroup,1,ntstp,sage,nage);
 	//3darray Lstd(1,ngroup,1,ntstp,sage,nage); 
 	
	
 	

 	//3darray totB(1,ngroup,1,ntstp,sage,nage);
 	
 	3darray PosX(1,ngroup,smon,nmon,sage-1,nage);
 	//3darray NAreaLen(1,ntstp,sarea,narea,1,nlen);
 	//3darray CatchAreaLen(1,ptstp,sarea,narea,1,nlen);
 	3darray CatchNatAge(1,ptstp,1,nations,sage-2,nage);
 	3darray CatchFGAge(1,ptstp,1,fisharea,sage-2,nage);
 	3darray CatchNatLen(1,ptstp,1,nations,1-2,nlen);
 	3darray CatchFGLen(1,ptstp,1,fisharea,1-2,nlen);
 	//3darray EffNatAge(1,fisharea,1,ntstp,sage-2,nage);
 	// yYieldNatAge(syr,proj_yr,1,fisharea,sage-2,nage);
 	//3darray yCatchStateAge(syr,proj_yr,1,nations,sage-2,nage);

 	3darray yCatchNatLen(syr,proj_yr,1,nations,1-2,nlen);
 	3darray yCatchFGLen(syr,proj_yr,1,fisharea,1-2,nlen);
 	3darray yCatchNatAge(syr,proj_yr,1,nations,sage-2,nage);
 	3darray yCatchFGAge(syr,proj_yr,1,fisharea,sage-2,nage);
 	
 	3darray NAreaAgeG(1,n_rg,1,ptstp,sage-3,nage);
 	3darray NAreaLenG(1,n_rg,1,ptstp,1-3,nlen);
 	3darray CatchAreaAgeG(1,n_rg,1,ptstp,sage-3,nage);
 	3darray CatchAreaLenG(1,n_rg,1,ptstp,1-3,nlen);
 	3darray propVBarea(1,n_rg,1,ptstp,sage-3,nage);
 	4darray P_al(1,ngroup,smon,nmon,sage,nage,1-3,nlen); 

 	matrix obsCatchFGAge(1,tot_pcat,sage-3,nage);

PRELIMINARY_CALCS_SECTION


	//cout<<"wt"<<wt(rep_yr+1,nyr)<<endl;
	//exit(1);
	read_catlim();
	initialization();

	calc_InitPos_gtg();
	incidence_functions();

	//cout<<"aqui??"<<endl;

	calc_first_year();

	move_grow_die();

	exit(1);
	
	
	//run_projections();
	
	//output_true();
	
	//STUFF FROM LAGRANGIAN SIMEVAL
	//output_pin();
	//output_dat();
	//output_gtgdat();


	//save_OMrep();
	
	

PROCEDURE_SECTION



FUNCTION double cnorm(const dvariable x, const dvariable& mu, const dvariable& sd)

	double rst;
	double stx;

	stx = value((x-mu)/sd);

	rst = cumd_norm(stx);
	

	return(rst);


FUNCTION dvar_vector cnorm2(const double& x, const dvar_vector& mu, const dvar_vector& sd)

	dvar_vector rst(sage,nage);
	dvar_vector stx(sage,nage);

	for(int a= sage; a<= nage; a++)
	{
		stx(a) = (x-mu( a ))/sd( a );
		rst(a) = cumd_norm(stx( a ));
	}

	return(rst);


FUNCTION calc_InitPos_gtg
	
	dvar_vector inimu(sage,nage);
	dvector Xtemp(1,ngroup);

	dvariable dif;

	maxPos.initialize();
	calcmaxpos(0.0);
	
	varPos=maxPos*cvPos;
	varPosg=sqrt((varPos*varPos)/(ngroup*ngroup*4));

	
	
	//Xtemp=value(Xini);

	//cout<<"Xtemp"<<Xtemp<<endl;
	
	dif=(Xini(2)-Xini(1))/2.;

	inimu.initialize();

	prop_ng(1) = cumd_norm(value(Xini(1)+dif))-cumd_norm(value(Xini(1)-dif));
	
	//cout<<"dif"<<dif<<endl;


	for(int g=2; g<=ngroup ;g++)
	{
		prop_ng(g) =cumd_norm(value(Xini(g)+dif))-cumd_norm(value(Xini(g)-dif));
	}

	//prop_ng = prop_ng/sum(prop_ng);
	//cout<<"sum(prop_ng)"<<sum(prop_ng)<<endl;
	prop_ng = prop_ng/value(sum(prop_ng));

	//cout<<"prop_ng"<<prop_ng<<endl;

	dvar_vector meanPosX(sage,nage);
	 meanPosX.initialize();

	 for(int m=smon;m<=nmon;m++){
			meanPosX(sage,nage) = minPos + (maxPos - minPos) * (0.5+0.5*sin(m*PI/6. - mo*PI/6. -PI/2.)); 
		for(int g=1; g<=ngroup ;g++)
		{
			PosX(g)(m)(sage-1) = g;
			PosX(g)(m)(sage,nage) =  Xini(g)*varPos+meanPosX(sage,nage);
	 	
	 	}
	}


	cout<<"ok after calc_InitPos_gtg"<<endl;



FUNCTION incidence_functions

	
	lxo(sage)=1.;
	for(int a = sage+1; a<= nage; a++){
		lxo(a) = lxo(a-1)*mfexp(-m);
	}	
	lxo(nage) += lxo(nage)* mfexp(-m); 


	for(int im=smon;im<=nmon;im++){
		calc_length_comps(im);
	}

	kappa 	= 4*h/(1-h);


	for(int g=1;g<=ngroup;g++){
		phiE += elem_prod(lxo,fa)*Wage(g)(1)(sage,nage)*prop_ng(g);
		Bo 	+= lxo*Ro* Wage(g)(1)(sage,nage)*prop_ng(g);
	}

	So 		= kappa/(phiE);
	SBo 	= phiE*Ro;
	beta 	= (kappa-1)/SBo;


	//cout<<"kappa "<<kappa<<endl;
	//cout<<"SBo "<<SBo <<endl;
	//cout<<"lxo "<<lxo <<endl;
	//cout<<"wa "<<wa <<endl;
	//cout<<"Wage(g)(1)(sage,nage) "<< Wage(1)(1)(sage,nage) <<endl;
	

	m_tsp = m/tmon;

	//cout<<"m_tsp"<<m_tsp<<endl;
	za 	= m_tsp+va*fe;

	cout<<"Ok after incidence_functions"<<endl;  

FUNCTION void calc_numbers_at_age(const int& ii, const dvariable& expwt)

	dvar_vector propBarea(sarea,narea);
	//dvar_vector tnage(sage,nage);

	//tnage.initialize();
	
	
			
			
	switch (indmonth(ii)){

		
        case 1:  
        	
        	//cout<< "prop_ng "<<sum(prop_ng)<<endl;
        	//cout<< "expwt "<<expwt<<endl;
        	
        	//exit(1);

            for(int g=1;g<=ngroup;g++)
			{
				Nage(g)(ii)(sage-2)=ii;
				Nage(g)(ii)(sage-1)=g;	
				
            	Nage(g)(ii)(sage) = (So*SB(ii-nmon)/(1.+beta*SB(ii-nmon)))*mfexp(expwt- (sigR*sigR/2.))* prop_ng(g);//*1.0/dg;

            	//(4.*h*Ro*SB(ii-nmon))/(SBo*(1-h)+SB(ii-nmon)*(5*h-1))*mfexp(expwt- (sigR*sigR/2))* prop_ng(g);
            	//
            	//cout<<"chegou aqui? "<<endl;

            	for(int a = sage+1;a<=nage;a++)
            	{
					propBarea.initialize();
			
					for(int rr =sarea; rr<=narea; rr++)
					{		
						propBarea(rr) = cnorm(areas(rr)+0.5,PosX(g)(nmon)(a-1),varPosg(a-1))-cnorm(areas(rr)-0.5,PosX(g)(nmon)(a-1),varPosg(a-1));	
					}

					//dvariable psarea;
					//psarea=cnorm(areas(narea)+0.5,PosX(g)(nmon)(a-1),varPosg(a-1))-cnorm(areas(sarea)-0.5,PosX(g)(nmon)(a-1),varPosg(a-1));
						
					Nage(g)(ii)(a) = Nage(g)(ii-1)(a-1)*propBarea*mfexp(-(m_tsp+q*Effarea(ii-1)*va(a-1))) +
								     Nage(g)(ii-1)(a-1)*(1.0-sum(propBarea))*mfexp(-(m_tsp));
				}
				
				Nage(g)(ii)(nage) += Nage(g)(ii-1)(nage)*propBarea*mfexp(-(m_tsp+q*Effarea(ii-1)*va(nage))) +
								     Nage(g)(ii-1)(nage)*(1.0-sum(propBarea))*mfexp(-(m_tsp));

				//sum(elem_div(elem_prod(Nage(g)(ii-1)(nage-1)*propBarea,mfexp(-(m_tsp+q*Effarea(ii-1)*va(nage-1)))),
            	//				    (1.0-mfexp(-(m_tsp+q*Effarea(ii-1)*va(nage))))))+
            	//				    (Nage(g)(ii-1)(nage-1)*(1.0-sum(propBarea))*mfexp(-m_tsp))/(1.-mfexp(-m_tsp));
            	
            	//cout<<"propBarea)"<<sum(propBarea)<<endl;
            	yNage(indyr(ii))(sage,nage) += Nage(g)(ii)(sage,nage);
				
            	VulB(g)(ii)(sage-2) = ii;
        		VulB(g)(ii)(sage-1) = g;
	
				VulB(g)(ii)(sage,nage) = elem_prod(elem_prod(Nage(g)(ii)(sage,nage),va),Wage(g)(indmonth(ii))(sage,nage));
				VulBage(ii)(sage,nage) += VulB(g)(ii)(sage,nage);

				//tnage(sage,nage) += Nage(g)(ii)(sage,nage);
				SB(ii) += elem_prod(Nage(g)(ii)(sage,nage),fa)*Wage(g)(indmonth(ii))(sage,nage);
				tB(ii) += Nage(g)(ii)(sage,nage)*Wage(g)(indmonth(ii))(sage,nage);
					//totB(g)(ii)(sage,nage) = elem_prod(Nage(g)(ii)(sage,nage),wa(sage,nage));
			

			}

			ySB(indyr(ii))= SB(ii); // multiplied by two so it is both males and females
            ytB(indyr(ii))= tB(ii);
            calc_empirical_wa(indyr(ii));

            //cout<<"Nage(g)(ii)(sage) "<<Nage(g)(ii)(sage)<<endl;
            //exit(1);

        break;
            	
        default: 

            for(int g=1;g<=ngroup;g++)
			{
				Nage(g)(ii)(sage-2)=ii;
				Nage(g)(ii)(sage-1)=g;	

            	for(int a = sage;a<=nage;a++)
            	{
            		propBarea.initialize();
			
					for(int rr =sarea; rr<=narea; rr++)
					{		
						propBarea(rr) = cnorm(areas(rr)+0.5,PosX(g)(indmonth(ii)-1)(a),varPosg(a))-cnorm(areas(rr)-0.5,PosX(g)(indmonth(ii)-1)(a),varPosg(a));	
					}

					//cout<<"aqui?? "<<endl
						
            		Nage(g)(ii)(a) = Nage(g)(ii-1)(a)*propBarea*mfexp(-(m_tsp+q*Effarea(ii-1)*va(a)))+
            						 Nage(g)(ii-1)(a)*(1.0-sum(propBarea))*mfexp(-(m_tsp));
            	}

            	VulB(g)(ii)(sage-2) = ii;
        		VulB(g)(ii)(sage-1) = g;
	
				VulB(g)(ii)(sage,nage) = elem_prod(elem_prod(Nage(g)(ii)(sage,nage),va),Wage(g)(indmonth(ii))(sage,nage));
				VulBage(ii)(sage,nage) += VulB(g)(ii)(sage,nage);

				//tnage(sage,nage) += Nage(g)(ii)(sage,nage);
				SB(ii) += elem_prod(Nage(g)(ii)(sage,nage),fa)*Wage(g)(indmonth(ii))(sage,nage);
				tB(ii) += Nage(g)(ii)(sage,nage)*Wage(g)(indmonth(ii))(sage,nage);
				//totB(g)(ii)(sage,nage) = elem_prod(Nage(g)(ii)(sage,nage),wa(sage,nage));
			            
        	}

        break;    	
        		
	}


	
FUNCTION void calc_effarea(const int& ii,const int& ia, const dvar_vector& ctlim)

	dvar_vector tmp1(sarea,narea);
	dvar_vector tmp2(sarea,narea);
	tmp1.initialize();
	tmp2.initialize();

	for(int n=1; n<=nations;n++){
       totVBnation(ii,n) = sum(pow(tVBarea(ii)(ntmp1(n),ntmp1(n+1)-1.0)+1e-30,fbeta));		 
	}

	//cout<<"ctlim "<<ctlim <<endl;

	for(int r= sarea; r<=narea; r++)
	{
		if(yYieldNat(indyr(ii))(indnatarea(r))<ctlim(indnatarea(r))){

			double lmc;

			if(indmonth(ii)>1){

				lmc = value(yYieldNat(indyr(ii-1))(indnatarea(r))/ctlim(indnatarea(r)));			
			}else{
				lmc = 0.;
			}

			if(lmc>.85){

				tmp1(r)= (pow(tVBarea(ii)(r)+1e-30,fbeta)/(totVBnation(ii)(indnatarea(r)))) * effPwr(r);
				tmp2(r) = tmp1(r)*TotEffyear(indfisharea(r))(indyr(ia));
				Effarea(ii)(r) = tmp2(r)*TotEffmonth(indfisharea(r))(indmonth(ii))* (1.-lmc);
			}else{

		//	if(lmc>1.){
		//		tmp1(r)= (pow(tVBarea(ii)(r)+1e-30,fbeta)/(totVBnation(ii)(indnatarea(r)))) * effPwr(r);
		//		tmp2(r) = tmp1(r)*TotEffyear(indfisharea(r))(indyr(ia));
		//		Effarea(ii)(r) = tmp2(r)*TotEffmonth(indfisharea(r))(indmonth(ii))* (lmc-1.);
				//Effarea(ii)(r) = 0.0;
		//		cout<<"catlim(indnatarea(r)) "<<catlim(indnatarea(r))<<endl;
		//		cout<<"lmc "<<lmc<<endl;
		//		cout<<"indmonth "<<indmonth(ii)<<endl;
		//		cout<<"yYieldNat(indyr(ii-1) "<<yYieldNat(indyr(ii-1))<<endl;

		//	}else{

				tmp1(r)= (pow(tVBarea(ii)(r)+1e-30,fbeta)/(totVBnation(ii)(indnatarea(r)))) * effPwr(r);
				tmp2(r) = tmp1(r)*TotEffyear(indfisharea(r))(indyr(ia));
				Effarea(ii)(r) = tmp2(r)*TotEffmonth(indfisharea(r))(indmonth(ii));
			//	cout<<"veio pra ca? "<<endl;
			}


			//cout<<"yYieldNat(indyr(ii))(indnatarea(r)) "<<yYieldNat(indyr(ii))(indnatarea(r)) <<endl;
			//cout<<"ctlim(indnatarea(r))"<<ctlim(indnatarea(r)) <<endl;
		}else{

					Effarea(ii)(r) = 0.0;
					
				//cout<<"closed fisheries year is"<<indyr(ii)<<"mes"<<indmonth(ii)<<"pais"<<indnatarea(r) <<endl;
		}
	}
		
	//cout<<"Ok after calc_effarea"<<endl;


FUNCTION void calc_position(const int& ii)

	dvar_vector meanPosX(sage,nage);
	meanPosX.initialize();


	meanPosX(sage,nage) = minPos + (maxPos - minPos) * (0.5+0.5*sin(indmonth(ii)*PI/6. - mo*PI/6. -PI/2.)); 


	int g, r, ig, b,a;

	for(ig=1;ig<=n_rg ;ig++)
	{
		
		g = n_group(ig);
		r = n_area(ig);

		NAreaAgeG(ig)(ii)(sage-3)= ii;
		NAreaAgeG(ig)(ii)(sage-2)= g;
		NAreaAgeG(ig)(ii)(sage-1)= r;

		NAreaLenG(ig)(ii)(1-3) = ii;
		NAreaLenG(ig)(ii)(1-2) = g;
		NAreaLenG(ig)(ii)(1-1) = r;



		PosX(g)(indmonth(ii))(sage-1) = g;
		PosX(g)(indmonth(ii))(sage,nage) =  Xini(g)*varPos+meanPosX(sage,nage);
		

		NAreaAgeG(ig)(ii)(sage,nage) = elem_prod(Nage(g)(ii)(sage,nage),(cnorm2(areas(r)+0.5,PosX(g)(indmonth(ii))(sage,nage),varPosg)-cnorm2(areas(r)-0.5,PosX(g)(indmonth(ii))(sage,nage),varPosg)));
	
		propVBarea(ig)(ii)(sage-3) = ii;
		propVBarea(ig)(ii)(sage-2) = g;
		propVBarea(ig)(ii)(sage-1) = r;
		propVBarea(ig)(ii)(sage,nage) =  elem_prod(VulB(g)(ii)(sage,nage), (cnorm2(areas(r)+0.5,PosX(g)(indmonth(ii))(sage,nage),varPosg)-cnorm2(areas(r)-0.5,PosX(g)(indmonth(ii))(sage,nage),varPosg)));
			
		tVBarea(ii)(r) += sum(propVBarea(ig)(ii)(sage,nage));		

		//for(b=1;b<=nlen;b++)
		//{
		//		NAreaLen(ii)(r)(b) += elem_prod(Nage(g)(ii)(sage,nage),(cnorm(areas(r)+0.5,PosX(g)(ii),varPosg)-cnorm(areas(r)-0.5,PosX(g)(ii),varPosg)))*P_la(g)(indmonth(ii))(b)(sage,nage);
		//}

		for(a = sage; a<=nage;a++)
		{	
				//cout<<"sage "<<sage<<endl;
				
				//cout<<"NAreaAgeG(ig)(ii)(a) "<<NAreaAgeG(ig)(ii)(a)<<endl;
				//cout<<"P_al*NAreaAgeG(ig)(ii)(a)"<<NAreaAgeG(ig)(ii)(a)*P_al(g)(indmonth(ii))( a )( 1,nlen )<<endl;
				
				NAreaLenG(ig)(ii)(1,nlen) += NAreaAgeG(ig)(ii)(a)*P_al(g)(indmonth(ii))( a )( 1,nlen );
	
		}


	}




	//cout<<"Ok after calc_position"<<endl;


FUNCTION void calc_catage(const int& ii)

		int ig,igg,r,rr,gg,g,n;
			
		dvar_vector tmpCAAG(sage,nage);
		dvar_vector yYN(1,nations);

		tmpCAAG.initialize();
		yYN.initialize();

		yYN =	yYieldNat(indyr(ii));
		//cout<<" indyr(ii)"<<indyr(ii)<<endl;
		//cout<<" yYieldNat(indyr(ii))"<<yYieldNat(indyr(ii))<<endl;

	for(int rr=sarea;rr<=narea;rr++)
	{
			dvar_vector tmpc1(sage,nage);
			dvar_vector tmpc2(sage,nage);
		
			tmpc1.initialize();
			tmpc2.initialize();
		

			tmpc1 =elem_div(q*Effarea(ii)(rr)*va,q*Effarea(ii)(rr)*va+m_tsp);
			tmpc2 =elem_prod(tmpc1,(1-mfexp(-(q*Effarea(ii)(rr)*va+m_tsp))));


		for(int gg=1;gg<=ngroup;gg++)
		{
		//gg = n_group(igg);
		//rr = n_area(igg);



			
		//cout<<" tmpc1 "<<tmpc1<<endl;
		//cout<<" tmpc2 "<<tmpc2<<endl;
		//cout<<" NAreaAgeG(ig)(comm_obsCatageii)(sage,nage) "<<NAreaAgeG(ig)(ii)(sage,nage)<<endl;

			tmpCAAG(sage,nage) = elem_prod(tmpc2,NAreaAgeG(pntr_rg(rr,gg))(ii)(sage,nage));		
		
			yYN(indnatarea(rr)) += tmpCAAG(sage,nage)*Wage(gg)(indmonth(ii))(sage,nage);					
			
			
		}

		if(yYN(indnatarea(rr))>.9*catlim(indnatarea(rr))){
				//cout<<"yYN(indnatarea(rr))"<<(indnatarea(rr))<<" "<< yYN(indnatarea(rr))<<endl;	
		
			//if(yYN(indnatarea(rr))<1.3*catlim(indnatarea(rr))){
				Effarea(ii)(rr) *= 0.01;
				//Effarea(ii)(rr) *= 0.1;
			//}else{	
				//cout<<" aqui"<<endl;
				//cout<<"closed"<< indyr(ii) <<"nation is"<< indnatarea(rr)<<endl;
			//	Effarea(ii)(rr) = 0;
			//}
		}

		
	}

    













	for(int ig=1;ig<=n_rg;ig++)
	{
		g = n_group(ig);
		r = n_area(ig);


		dvar_vector tmpc1(sage,nage);
		dvar_vector tmpc2(sage,nage);

		tmpc1.initialize();
		tmpc2.initialize();

		CatchAreaAgeG(ig)(ii)(sage-3)= ii;
		CatchAreaAgeG(ig)(ii)(sage-2)= g;
		CatchAreaAgeG(ig)(ii)(sage-1)= r;

		CatchNatAge(ii)(indnatarea(r))(sage-2) = ii;
		CatchNatAge(ii)(indnatarea(r))(sage-1) = indnatarea(r);

		CatchFGAge(ii)(indfisharea(r))(sage-2) = ii;
		CatchFGAge(ii)(indfisharea(r))(sage-1) = indfisharea(r);

		yCatchFGAge(indyr(ii))(indfisharea(r))(sage-2) = indyr(ii);
		yCatchFGAge(indyr(ii))(indfisharea(r))(sage-1) = indfisharea(r);

		yCatchNatAge(indyr(ii))(indnatarea(r))(sage-2) = indyr(ii);
		yCatchNatAge(indyr(ii))(indnatarea(r))(sage-1) = indfisharea(r);



		tmpc1 =elem_div(q*Effarea(ii)(r)*va,q*Effarea(ii)(r)*va+m_tsp);
		tmpc2 =elem_prod(tmpc1,(1-mfexp(-(q*Effarea(ii)(r)*va+m_tsp))));

		//cout<<" tmpc1 "<<tmpc1<<endl;
		//cout<<" tmpc2 "<<tmpc2<<endl;
		//cout<<" NAreaAgeG(ig)(comm_obsCatageii)(sage,nage) "<<NAreaAgeG(ig)(ii)(sage,nage)<<endl;

		CatchAreaAgeG(ig)(ii)(sage,nage) = elem_prod(tmpc2,NAreaAgeG(ig)(ii)(sage,nage));		
		CatchNatAge(ii)(indnatarea(r))(sage,nage) += CatchAreaAgeG(ig)(ii)(sage,nage);//+0.0000001;
		CatchFGAge(ii)(indfisharea(r))(sage,nage) += CatchAreaAgeG(ig)(ii)(sage,nage);//+0.0000001;
		yCatchFGAge(indyr(ii))(indfisharea(r))(sage,nage) += CatchAreaAgeG(ig)(ii)(sage,nage);//+0.0000001;
		yCatchNatAge(indyr(ii))(indnatarea(r))(sage,nage) += CatchAreaAgeG(ig)(ii)(sage,nage);//+0.0000001;
		

		yYieldtotal(indyr(ii)) += CatchAreaAgeG(ig)(ii)(sage,nage)*Wage(g)(indmonth(ii))(sage,nage);
		Yieldtotal(ii) += CatchAreaAgeG(ig)(ii)(sage,nage)*Wage(g)(indmonth(ii))(sage,nage);
		yYieldtotalage(indyr(ii))(sage,nage) += elem_prod(CatchAreaAgeG(ig)(ii)(sage,nage),Wage(g)(indmonth(ii))(sage,nage));
		YieldNat(ii)(indnatarea(r)) += CatchAreaAgeG(ig)(ii)(sage,nage)*Wage(g)(indmonth(ii))(sage,nage);					
		yYieldNat(indyr(ii))(indnatarea(r)) += CatchAreaAgeG(ig)(ii)(sage,nage)*Wage(g)(indmonth(ii))(sage,nage);			

		
	}



	for(int n=1;n<=nations;n++)
	{
		yCatchtotalage(indyr(ii))(sage,nage) += CatchNatAge(ii)(n)(sage,nage);
		
	}
	
			
			
			
			//totcat(sage,nage) += CatageG(g)(ii+(im-1))(sage,nage);

		

		//cout<<"CatchNatAge(ii) is "<<CatchNatAge(ii)<<endl;
		
		
		//cout<<"Ok after calc_catage"<<endl;

FUNCTION void calc_catlen(const int& ii)
	
			int ig,r,g,b,a;
			
	for(int ig=1;ig<=n_rg;ig++)
	{
				
		g = n_group(ig);
		r = n_area(ig);

		
		CatchAreaLenG(ig)(ii)(1-3)= ii;
		CatchAreaLenG(ig)(ii)(1-2)= g;
		CatchAreaLenG(ig)(ii)(1-1)= r;
				
		CatchNatLen(ii)(indnatarea(r))(-1) = ii;
		CatchNatLen(ii)(indnatarea(r))(0) = indnatarea(r);

		CatchFGLen(ii)(indfisharea(r))(-1) = ii;
		CatchFGLen(ii)(indfisharea(r))(0) = indfisharea(r);
				
		for(a = sage; a<=nage;a++)
		{
								
			CatchAreaLenG(ig)(ii)(1,nlen) += CatchAreaAgeG(ig)(ii)(a)*P_al(g)(indmonth(ii))( a )( 1,nlen );
		}
		
		CatchNatLen(ii)(indnatarea(r))(1,nlen) += CatchAreaLenG(ig)(ii)(1,nlen);
		CatchFGLen(ii)(indfisharea(r))(1,nlen) += CatchAreaLenG(ig)(ii)(1,nlen);
				
		Catlen(ii)(1,nlen) += CatchAreaLenG(ig)(ii)(1,nlen);
				
		yCatchFGLen(indyr(ii))(indfisharea(r))(-1) = indyr(ii);
		yCatchFGLen(indyr(ii))(indfisharea(r))(0) = indfisharea(r);
		yCatchFGLen(indyr(ii))(indfisharea(r))(1,nlen) += CatchAreaLenG(ig)(ii)(1,nlen);			
				
		yCatchNatLen(indyr(ii))(indnatarea(r))(1-2) = indyr(ii);
		yCatchNatLen(indyr(ii))(indnatarea(r))(1-1) = indnatarea(r);
		yCatchNatLen(indyr(ii))(indnatarea(r))(1,nlen) += CatchAreaLenG(ig)(ii)(1,nlen);
				
		yCatchtotalLen(indyr(ii))(1,nlen) += CatchAreaLenG(ig)(ii)(1,nlen);
				
	}

			
		//cout<<"Ok after calc_catLen"<<endl;


	
FUNCTION initialization
			
	Nage.initialize();
 	VulB.initialize();
 	VulBage.initialize();

 	PosX.initialize();
 	yYieldtotal.initialize();
 	Yieldtotal.initialize();


 	CatchNatAge.initialize();
 	CatchFGAge.initialize();
 	NAreaAgeG.initialize();
 	NAreaLenG.initialize();
 	CatchAreaAgeG.initialize();
 	CatchAreaLenG.initialize();
 	propVBarea.initialize();

 	
 	tB.initialize();
 	tVBarea.initialize();
 	
 	totVBnation.initialize();
	
	Effarea.initialize();
 	


 	//NAreaLen.initialize();

 	//CatageG.initialize();
 	//Catage.initialize();
 	yNage.initialize();
 	yFatage.initialize();
 	ywa.initialize();
 	wa.initialize();
 	//CatchAreaLen.initialize();
 	yCatchtotalage.initialize();
 	yCatchtotalLen.initialize();
 	yYieldNat.initialize();
 	YieldNat.initialize();
 	yCatchNatLen.initialize();
 	yCatchFGLen.initialize();


FUNCTION calc_first_year  

	

 	for(int g=1;g<=ngroup;g++)
	{
		Nage(g)(1)(sage-1)=g;
		Nage(g)(1)(sage-2)=1;
	
			
		Nage(g)(1)(sage) = Ro*prop_ng(g);
		//(So*SBo/(1.+beta*SBo))* prop_ng(g);

        for(int a=sage+1 ; a <= nage ; a++)
		{
			Nage(g)(1)(a) = Nage(g)(1)(a-1) * mfexp(-m);
		}
		Nage(g)(1)(nage) += (Nage(g)(1)(nage)*mfexp(-m));
		//Nage(g)(1)(nage) /= (1.-mfexp(-m));
      		          	
        VulB(g)(1)(sage-2) = 1;
        VulB(g)(1)(sage-1) = g;
		VulB(g)(1)(sage,nage) = elem_prod(elem_prod(Nage(g)(1)(sage,nage),va),Wage(g)(1)(sage,nage));
	
		SB(1) += elem_prod(Nage(g)(1)(sage,nage),fa)*Wage(g)(1)(sage,nage);
		tB(1) += Nage(g)(1)(sage,nage)*Wage(g)(1)(sage,nage);
		yNage(indyr(1))(sage,nage) += Nage(g)(1)(sage,nage);
		
		
	}

	//cout<<"yNage(indyr(1))(sage,nage)"<<yNage(indyr(1))(sage,nage)<<endl;


	ytB(indyr(1))= tB(1);
	ySB(indyr(1))= SB(1);
	calc_empirical_wa(indyr(1));
	
	
	calc_position(1);


	dvar_vector tmp1(sarea,narea);
	dvar_vector tmp2(sarea,narea);
	tmp1.initialize();
	tmp2.initialize();

	for(int n=1; n<=nations;n++){
       totVBnation(1,n) = sum(pow(tVBarea(1)(ntmp1(n),ntmp1(n+1)-1.0)+1e-30,fbeta));		 
	}

	for(int r= sarea; r<=narea; r++)
	{
		
		tmp1(r)= (pow(tVBarea(1)(r)+1e-30,fbeta)/(totVBnation(1)(indnatarea(r)))) * effPwr(r);
		tmp2(r) = tmp1(r)*TotEffyear(indfisharea(r))(1);
		Effarea(1)(r) = tmp2(r)*TotEffmonth(indfisharea(r))(1);
	}

	//calc_effarea(1,1,catlim);
	
	

	cout<<"Ok after calc_first_year "<<endl;
 

FUNCTION move_grow_die

		//if varps and mu are time varying need to move these to inside the loop
		//calc_InitPos_gtg();
		
		//maxPos.initialize();		
		//calcmaxpos();
		
		//varPos=maxPos*cvPos;
		//varPosg=sqrt((varPos*varPos)/(ngroup*ngroup*4));

	dvar_vector catl(1,nations);
	catl.fill("{100000,100000}");



	for(int ie=2;ie<=ststp;ie++)
	{
		//cout<<"ie is "<<indyr(ie)<<" " << ie<<endl;
		calcmaxpos(0.0);
		//calcmaxpos(wx(syr));
		calc_numbers_at_age(ie,0.0 );	//wt(indyr(ie))	
		//calc_cumulative_numbers(ie);
		calc_position(ie);
		
		calc_effarea(ie,ie, catlim);
		calc_catage(ie);
		calc_catlen(ie);
		

	}
	//exit(1);

	int p;
 	p=1;

 	int svyr;
 	svyr = 1;

	for(int i=ststp+1;i<=ntstp;i++){


		maxPos.initialize();
		calcmaxpos(wx(indyr(i)));
		calc_numbers_at_age(i,wt(indyr(i)));


		if(indmonth(i)==smon){

			calc_hist_catlim(indyr(i));
			read_catlim();
			//cout<< "year is"<< indyr(i)<<endl;
			//cout<< "catlim is"<< catlim <<endl;
			//cout<< "hiscatch is "<< catlim_hist(indyr(i))<<endl;
			

		}


		//calc_cumulative_numbers(i);
	
		//cout<<"indyr(i) "<<indyr(i)<<endl;
		//cout<<"(i) "<<(i)<<endl;
		//cout<<"indmonth(i) "<<indmonth(i)<<endl;
		
		calc_position(i);
		
	
		calc_effarea(i,i,catlim);	
		calc_catage(i);
		calc_catlen(i);	

		




		for(int n=1;n<=fisharea;n++)
		{					
			if(TotEffmonth(n)(indmonth(i))>0)
       		{
					clean_catage(i,p,n);
					p++;
       		}	
		}


		
		if(indmonth(i)==nmon){
			//survey calculations
			
			if(surv_yrs(svyr)==indyr(i)){
				
				survey_data(svyr);
				svyr +=1;
			}

			catage_comm(indyr(i));

			catlen_comm(indyr(i));
			calc_wt_comm(indyr(i));
			calc_selectivity(indyr(i));
			calc_msy(indyr(i));
			//cout<< "yYieldNat is"<< yYieldNat(indyr(i)) <<endl;
		
		}
	}

	//calc_known_catlim(nyr);
	//read_catlim();
	//histTAC(nyr)=catlim;


	run_projections(p,svyr);
	 //exit(1);


	cout<<"Ok after move_grow_die"<<endl;


FUNCTION void clean_catage(const int& ii,const int& pp,const int& nn)
	
	
				dvector pa(nage,sage);
       			pa.initialize();		

       			obsCatchFGAge(pp)(sage-3) = ii;
       			obsCatchFGAge(pp)(sage-2) = indmonth(ii);
				obsCatchFGAge(pp)(sage-1) = nn;
				
				pa = value(CatchFGAge(ii)(nn)(sage,nage)/(sum(CatchFGAge(ii)(nn)(sage,nage))))+0.0000001;
				
			 	obsCatchFGAge(pp)(sage,nage) = rmvlogistic(pa,tau_c,seed+ii);
				
				//cout<<" obsCatchNatAge(pp)(sage,nage) "<<obsCatchNatAge(pp)(sage,nage)<<endl;

	//cout<<"Ok after clean_catage"<<endl;
	


FUNCTION dvar_vector calcmaxpos(const dvariable& expwx)
	
		 					
		//maxPos(sage,nage) = (1./(1.+mfexp(-(age-maxPos50)/maxPossd)))*mfexp(expwx-0.1*0.1/2);
		maxPos(sage,nage) = (1./(1.+mfexp(-(age-maxPos50)/maxPossd)))*mfexp(expwx);//*mfexp(expwx)
		//maxPos(sage,nage) = (1./(1.+mfexp(-(age-maxPos50)/maxPossd)));
		maxPos(sage,nage) *= (narea-minPos(sage));
		maxPos(sage,nage) += minPos(sage);		
	
		return(maxPos);
		
	//cout<<"Ok after calcmaxpos"<<endl;



FUNCTION void catage_comm(const int& ii)

	dvector tmpca(sage,nage);
    tmpca.initialize();
	
	for(int i=(ii-syr+1)*tmon-tmon+1;i<=(ii-syr+1)*tmon;i++)
	{
		
		for(int n=1 ; n<=fisharea;n++)
		{
							
			if(TotEffmonth(n)(indmonth(i))>0)
       		{

       			tmpca(sage,nage) += value(CatchFGAge(i)(n)(sage,nage));
       		}
       		
		}
	}
			dvector pa(sage,nage);
       		pa.initialize();
      
    	//cout<<"tmpca(sage,nage) "<<sum(tmpca(sage,nage))<<endl;
		pa = ((tmpca)/(sum(tmpca)))+0.0000001;			
		//pa = value((tot_comm_obsCatage(ii)(sage,nage))/(sum(tot_comm_obsCatage(ii)(sage,nage))));			
		
		comm_obsCatage(ii)(sage,nage) = rmvlogistic(pa,tau_c,seed+ii);	
		Ut(ii) = Yieldtotal(ii)/tB((ii-syr+1)*tmon-(nmon-smon));
		//ytB(ii) = tB((ii-syr+1)*tmon-(nmon-smon));
	
		


FUNCTION void catlen_comm(const int& ii)

	///olha esse
	//for(int i=rep_yr*nmon+1;i<=ntstp;i++)
	//{

		dvector tempcl(1,nlen);
		 tempcl.initialize();

	for(int i=(ii-syr+1)*tmon-tmon+1;i<=(ii-syr+1)*tmon;i++)
	{
		
		for(int n=1;n<=fisharea;n++)
		{
							
			if(TotEffmonth(n)(indmonth(i))>0.0)
       		{

       			tempcl(1,nlen) += value(CatchFGLen(i)(n)(1,nlen));
       			//tempcl(1,nlen) += value(Catlen(i)(1,nlen));
				
				
       		
				//cout<<"CatchFGLen(i)("<<n<<")(1,nlen) "<<CatchFGLen(i)(n)(1,nlen)<<endl;
       		}
       		
		}
	}

	//cout<<"ii is  "<<ii<<endl;
	//cout<<"tempcl(1,nlen) "<< sum(tempcl(1,nlen))<<endl;


			dvector pl(1,nlen);
       		pl.initialize();
      
		totcatchLen(ii) = 	sum(tempcl(1,nlen));
    
		pl = ((tempcl(1,nlen))/(sum(tempcl(1,nlen))))+0.0000001;			
		//pa = value((tot_comm_obsCatage(ii)(sage,nage))/(sum(tot_comm_obsCatage(ii)(sage,nage))));			
		
		//cout<<"pl "<<pl<<endl;
		//cout<<"totcatchLen(ii) "<<totcatchLen(ii)<<endl;

		comm_obsCatLen(ii)(1,nlen) = rmvlogistic(pl,tau_l,seed+ii);	
		//comm_obsCatLen(ii)(1,nlen) = pl;		
		comm_obsCatLen(ii)(1,nlen) *= totcatchLen(ii);	
		//UtLen(ii) = totcatchLen(ii)/tB((ii-syr+1)*tmon-(nmon-smon));
		//cout<<"comm_obsCatLen"<<endl<<comm_obsCatLen<<endl;
		//exit(1);


FUNCTION void calc_wt_comm(const int& ii)
   

   	Watage_comm(ii)(sage,nage) = yYieldtotalage(ii)/yYieldtotal(ii);
   	
   	
   
FUNCTION  void survey_data(const int& ii)
			
			int ind_sv;
			ind_sv = (surv_yrs(ii)-syr+1)*(tmon)-(nmon-surv_mon);
			
			if(indmonth(ind_sv)==surv_mon)
	       	{	
	       		survB(ii)=sum(VulBage(ind_sv)(sage,nage))* mfexp(epsilon(ii));
	       		
	
	       		dvector pp(sage,nage);
	       		pp.initialize();
	       		for(int g=1;g<=ngroup;g++)
				{    
					//pp += value(CatchNatAge(ind_sv)(n)(sage,nage)); 			
					pp +=  value(Nage(g)(ind_sv)(sage,nage));
					//pp += value((CatchNatAge(ind_sv)(n)(sage,nage))/(sum(CatchNatAge(ind_sv)(n)(sage,nage))+0.01))+0.000001;
				}	
				dvector ppp(nage,sage);
	       		ppp.initialize();
	       		ppp = (pp)/(sum(pp));
				surv_obsCatage(ii)(sage,nage) = rmvlogistic(ppp,tau_survey,seed+ii);
       		}	
	
			
    //cout<<"Ok after survey_data"<< endl;

	//FUNCTION  void survey_dataLen(const int& ii)
	//		
	//		int ind_sv;
	//		ind_sv = (surv_yrs(ii)-syr+1)*(tmon)-(nmon-surv_mon);
	//		
	//		
	//       		survB(ii)=sum(VulBage(ind_sv)(sage,nage))* mfexp(epsilon(ii));
	//       		
	//
	//       		dvector ppl(1,nlen);
	//       		ppl.initialize();
	//       		for(int n=1;n<=fisharea;n++)
	//			{    
	//				ppl += value(CatchNatLen(ind_sv)(n)(1,nlen)); 			
	//				//pp += value((CatchNatAge(ind_sv)(n)(sage,nage))/(sum(CatchNatAge(ind_sv)(n)(sage,nage))+0.01))+0.000001;
	//			}	
	//			dvector pppl(1,nlen);
	//       		pppl.initialize();
	//       		pppl = ((ppl)/(sum(ppl)+0.01))+0.000001;
	//			surv_obsCatLen(ii)(sage,nage) = rmvlogistic(pppl,tau_surveyl,seed+ii);
	//      	
	
			
    //cout<<"Ok after survey_dataLen"<< endl;


FUNCTION void calc_length_comps(const int& mi)
	
	int g,a, b ;

	//cout<<"ii is "<<ii<<endl;
	wa.initialize();


	for(int g=1;g<=ngroup;g++)
	{
		dvar_vector ag(sage,nage);
		
		ag.fill_seqadd(sage,1);
		ag += ((mi-1)/12.);
			
		
		//cout<< "mi is "<< ((mi-1)/12.)<< endl;
		//cout<< "ag is "<< ag<< endl;
		//ag = age;
		Lage(g)(mi)(sage-2) = mi;
		Lage(g)(mi)(sage-1) = g;

		Wage(g)(mi)(sage-2) = mi;
		Wage(g)(mi)(sage-1) = g;

		Lage(g)(mi)(sage,nage) = Linf(g)*(1.- mfexp(-vbk*(ag-to)));
		Wage(g)(mi)(sage,nage) = awt*pow(Lage(g)(mi)(sage,nage),bwt);

		
		//Lstd(g)(mi)(sage,nage) = Lage(g)(mi)(sage,nage)*cvl;  		  //std for length at age


		dvar_vector z1(1,nlen); 				// intermediate steps for calculating proportion of age at length
		dvar_vector z2(1,nlen); 				// intermediate steps for calculating proportion of age at length
		dvar_vector std(sage,nage); 					// std for length at age curve

		std = Lage(g)(mi)(sage,nage)*cvl;

	

		 //Calculate proportion of length at age class

		//cout<<"std "<< endl<< std<<endl;
		//cout<<"Lage(g)(ii)(sage,nage) "<< endl<< Lage(g)(mi)(sage,nage)<<endl;

 		for( a = sage; a <= nage; a++ )
		{
			z1 = (( len - lenstp * 0.5 )-Lage(g)(mi)(a))/std( a );
			z2 = (( len + lenstp * 0.5 )-Lage(g)(mi)(a))/std( a );
			
		//cout<<"lenstp "<< endl<< lenstp<<endl;
		//cout<<"len "<< endl<< len<<endl;
			

		//cout<<"z1 "<< endl<< z1<<endl;
		//cout<<"z2 "<< endl<< z2<<endl;

			for( b=1; b<= nlen; b++ )
			{
				P_al(g)(mi)( a )( 0 )=b;
				P_al(g)(mi)( a )( 1 -2)=mi;
				P_al(g)(mi)( a )( 1 -3)=g;
				P_al(g)(mi)( a )( b )=cumd_norm( z2( b ))-cumd_norm( z1( b )); // calculates the proportion of a given age given your length
			
				
			}
		
			// plus length group
			P_al(g)(mi)( a)( nlen )= 1.00 -cumd_norm( z1( nlen ));
		
		}

		//cout<<"sum(P_la(g)(mi)( nlen )( sage,nage ))"<<sum(P_la(g)(indmonth(ii))( b )( sage,nage ))<<endl;
		
	
		wa(sage,nage) += Wage(g)(1)(sage,nage)*prop_ng(g);
		//	wa(sage,nage) = Wage(1)(1)(sage,nage);

		
		

		
	}
	//exit(1);

	//cout<<"Wage(g)(10)(sage,nage)"<<Wage(10)(1)(sage,nage)<<endl;
	//cout<<"prop_ng(g) "<<sum(prop_ng)<<endl;
	
	


FUNCTION void calc_selectivity(const int& ii)

	int n, nn;

	double sumcat;
	sumcat=value(sum(yCatchtotalage(ii)(sage,nage)));

	if(sumcat>0){

		for(n=1;n<=fisharea;n++){


			
			selfisharea(ii)(n)(sage-1)=ii;
			selfisharea(ii)(n)(sage-2)=n;
			selfisharea(ii)(n)(sage,nage) = elem_div(yCatchFGAge(ii)(n)(sage,nage),yNage(ii)(sage,nage))/max(elem_div(yCatchFGAge(ii)(n)(sage,nage),yNage(ii)(sage,nage)));

		}
	
		for(nn=1;nn<=nations;nn++){
			
			selnation(ii)(nn)(sage-1)=ii;
			selnation(ii)(nn)(sage-2)=nn;
			selnation(ii)(nn)(sage,nage) = elem_div(yCatchNatAge(ii)(nn)(sage,nage),yNage(ii)(sage,nage))/max(elem_div(yCatchNatAge(ii)(nn)(sage,nage),yNage(ii)(sage,nage)));	

		}

		seltotal(ii)(sage,nage) = elem_div(yCatchtotalage(ii)(sage,nage),yNage(ii)(sage,nage))/max(elem_div(yCatchtotalage(ii)(sage,nage),yNage(ii)(sage,nage)));
	}else{
		seltotal(ii)(sage,nage)=seltotal(ii-1)(sage,nage);
	}
		//cout<<"seltotal(ii)(sage,nage) "<<seltotal(ii)<<endl;
		//cout<<"yCatchtotalage(ii)(sage,nage) "<<yCatchtotalage(ii)(sage,nage)<<endl;


FUNCTION void calc_empirical_wa(const int& ii)


 		ywa(ii)(sage,nage).initialize();
		dvar_vector ypropg(1,ngroup);
		for(int g=1; g<=ngroup; g++){
			ypropg(g) =  sum(Nage(g)((ii-1)*tmon+1)(sage,nage))/sum(yNage(ii)(sage,nage));
			//cout<<"Nage(g)((ii-1)*tmon+1)(sage,nage) "<<Nage(g)((ii-1)*tmon+1)(sage,nage)<<endl;
			//cout<<"yNage(ii)(sage,nage) "<<yNage(ii)(sage,nage)<<endl;
			//cout<<"ypropg(g) "<<ypropg(g)<<endl;
			ywa(ii)(sage,nage) += Wage(g)(1)(sage,nage)*ypropg(g);
		}


FUNCTION calc_spr	

	int i, ii, a, g;

	
	dvector lz(sage,nage);
	
	
	for(ii=syr; ii<=proj_yr;ii++){
	
		yFatage(ii)(sage,nage) = elem_div(yCatchtotalage(ii)(sage,nage),yNage(ii)(sage,nage));

		lz.initialize();
		lz(sage) = 1.;
		for(a=sage+1; a<=nage;a++){
			lz(a)= value(lz(a-1)*mfexp(-m-yFatage(ii)(a-1)));
		}
		lz(nage) += value(lz(nage)*mfexp(-m-yFatage(ii)(nage)));

		//lz(nage) /= value(1.-mfexp(-m-yFatage(ii)(nage)));


		phie(ii)=elem_prod(lz,fa)*ywa(ii);
		spr(ii)=phie(ii)/phiE;

	
	}

	cout<<"Ok after calculating SPR"<<endl;






FUNCTION void run_stock_assessment(const int& ii,const int& svy)
			
		cout<<"running stock assessment"<<endl;
	
		switch (satype) {
	            case 1: 

	            //Length-SRA assessment 
	            
	            write_LSRA_data_file(ii,svy);
	
	            write_LSRA_ctl_file(ii,svy);

	
	            #if defined __APPLE__ || defined __linux
	             //cout<<m_est_fmsy<<endl;
	            	system("cd /Users/catarinawor/Documents/length_SRA/examples/closedloop/&& make clean && make");
	           
	            #endif
	
				//exit(1);
	
	            //#if defined _WIN32 || defined _WIN64
	            //system("lagrangian_SS.exe");
	            //#endif
	
	            //calc_catlim();
				//read_catlim();
	
				//system("cd ../../../R/read_mse && make");
	
	            break;
	
	            case 2:
	
	            write_iscam_data_file(ii,svy);
	            write_iscam_ctl_file(ii,svy);
	
	
	            #if defined __APPLE__ || defined __linux
	            // cout<<m_est_fmsy<<endl;
	            system("cd /Users/catarinawor/Documents/iSCAM/examples/hakelag/DATA && make clean && make");// && "exec make");
	           
	            #endif
	
	            //calc_catlim();
				//read_catlim();
				
	
				//system("cd ../../../R/read_mse && make readRdattwo");
	
	
	            //need to figure out how to write the winows version of this
	            //#if defined _WIN32 || defined _WIN64
	
	            //system("lagrangian_SS.exe");
	
	            //#endif
	
	
	
	             break;
	
	         }



FUNCTION calc_catlim

	
	//variables

	double fspr_SA;
	dvar_vector seltotal_SA(sage,nage);
	dvar_vector yNage_SA(sage,nage);
	double Bo_SA;
	double ytB_SA;
	

	cifstream ifs_cip("../TAC_input.dat");
	
    ifs_cip >> fspr_SA;
    ifs_cip >> seltotal_SA;
    ifs_cip >> yNage_SA;
    ifs_cip >> Bo_SA;
    ifs_cip >> ytB_SA;

    
    //calculate new catlim
	double BBo;
	dvariable TAC;
	//dvector TACnation(1,nations);


	BBo = ytB_SA/Bo_SA;


	if(BBo>0.4){

		TAC = elem_prod(elem_div(fspr_SA*seltotal_SA,fspr_SA*seltotal_SA+m),elem_prod(yNage_SA,(1-mfexp(-fspr_SA*seltotal_SA-m))))*wa; 
		
		

	}else{
		
		if(BBo>0.1){

			double y;

			y = (BBo-0.1)*100/(0.4-0.1);

			TAC = y*elem_prod(elem_div(fspr_SA*seltotal_SA,fspr_SA*seltotal_SA+m),elem_prod(yNage_SA,(1-mfexp(-fspr_SA*seltotal_SA-m))))*wa; 

			

		}else{

			TAC = 0;

		}
	}

	//cout<< "TAC" << endl << TAC<< endl;
	catlim = TAC * nationTACprop;
	cout<< "catlim" << endl << catlim<< endl;
	//cout<< "seltotal_SA" << endl << seltotal_SA<< endl;

	ofstream afs("catlim.txt");
	afs<<"#catlim" << endl << catlim <<endl;




FUNCTION void calc_msy(const int& ii )

  //This function calculates MSY in the lazy and slow way. 
 

	dvector ftest(1,700);
	ftest.fill_seqadd(0,0.001);
 	int k, kk ;
	int NF=size_count(ftest);
	
	//selage= seltotal(ii);
	dvector ye(1,NF);
	ye.initialize();

	//calc_selectivity(ii);

		
	for(k=1; k<=NF; k++)
	{
		dvector lz(sage,nage);
		dvector qka(sage,nage);
		lz.initialize();
		qka.initialize();
			
		dvariable phiq;
		dvariable phiz;
		dvariable req;

		phiq.initialize();
		phiz.initialize();
		req.initialize();
			
		lz.initialize();
		lz(sage) = 1.;

		for(int a=sage+1; a<=nage;a++){
			lz(a)= value(lz(a-1)*mfexp(-m-seltotal(ii-1)(a)*ftest(k)));
		}
		lz(nage) += value(lz(nage)*mfexp(-m-seltotal(ii-1)(nage)*ftest(k)));

		//lz(nage) /= value(1.-mfexp(-m-seltotal(ii-1)(nage)*ftest(k)));
		//cout<<"lz"<<lz<<endl;
		//cout<<"lxo"<<lxo<<endl;
 		
			phiz= elem_prod(lz,fa)*ywa(ii)/2.;
			
			qka = value(elem_div(elem_prod(elem_prod(ywa(ii-1),seltotal(ii-1)),(1.-mfexp(-m-seltotal(ii-1)*ftest(k)))),mfexp(-m-seltotal(ii-1)*ftest(k))));
			
			phiq = qka*lz;
			
			req = Ro*(kappa-phiE/phiz)/(kappa-1.);
			
			ye(k)= value(ftest(k)*req*phiq);

			//cout<<"ywa "<<ywa<<endl;
			//cout<<"kappa "<<kappa<<endl;
			//cout<<"phiE "<<phiE<<endl;
			//cout<<"phiE/phiz "<<phiE/phiz<<endl;
			//cout<<"req "<<req<<endl;
			//cout<<"ftest(k) "<<ftest(k)<<endl;
			//cout<<"ye(k) "<<ye(k)<<endl;

		}
		
		//cout<<"yield"<<ye<<endl;
		msy(ii)= max(ye);
		double mtest;	
		for(kk=1; kk<=NF; kk++)
		{
			mtest=ye(kk);
				
			if(mtest==msy(ii)){
				fmsy(ii)=ftest(kk);
			} 
		}

		//cout<<"seltotal(ii) "<<seltotal(ii)<<endl;
		//exit(1);
	

FUNCTION void calc_Fspr_target(const int& ii, double target )

  //This function calculates MSY in the lazy and slow way. 
 

	dvector ftest(1,800);
	ftest.fill_seqadd(0,0.001);
 	int k, kk, a;
	int NF=size_count(ftest);


	//selage= seltotal(ii);
	dvector tmp_phie(1,NF);
	dvector tmp_target(1,NF);
	tmp_phie.initialize();

		
	for(k=1; k<=NF; k++)
	{

		dvector lzt(sage,nage);
		lzt.initialize();

		lzt(sage)=1.0;
		
		for( a=sage+1; a<=nage; a++){
			lzt(a) = lzt(a-1) *value(mfexp(-m-seltotal(ii-1)(a)*ftest(k)));
		}
		lzt(nage) += value(lzt(nage)* mfexp(-m-seltotal(ii-1)(nage)*ftest(k)));
		

		//lzt(nage) /= 1.0 - value(mfexp(-m-seltotal(ii-1)(nage)*ftest(k)));
		
		//cout<<"ftest(k) "<<ftest(k)<<endl;
		//cout<<"seltotal(ii-1) "<<seltotal(ii-1)<<endl;
		//cout<<"lzt "<<lzt<<endl;

		tmp_phie(k) = value(elem_prod(lzt,fa)*ywa(ii));
		//cout<<"tmp_phie(k) "<<tmp_phie(k)<<endl;

		tmp_target(k) = fabs(value(tmp_phie(k)/phiE - target));
		//cout<<"tmp_target(kk) "<<tmp_target(k)<<endl;
		//cout<<"tmp_phie(k)/phiE "<<tmp_phie(k)/phiE <<endl;

	}
		
		double ttest;

		ttest =  min(tmp_target);

		for(kk=1; kk<=NF; kk++)
		{
			if(tmp_target(kk)==ttest){
				ftarget(ii)=ftest(kk);
			} 
		}
				
		cout<<"ftarget(ii) "<<ftarget(ii)<<endl;

		//exit(1);





FUNCTION void calc_known_catlim(const int& ii)

	
	//variables

	double fspr_true;

	
	double Bo_true;
	double SBo_true;
	//double vBo_true;
	double seen_bt;
	double seen_sbt;
	double seen_vbt;
	

	//vBo_true= value(elem_prod(lxo,seltotal(ii-1))*Ro* ywa(1));
	
	//cout<<"vBo_true "<< vBo_true<<endl;
	//Bo_true=value(Bo);
	SBo_true=value(SBo);
	
	seen_sbt = value(ySB(indyr(ii))*mfexp(ot(ii)));
	seen_bt = value(ytB(indyr(ii))*mfexp(ot(ii)));
	//seen_vbt = value((seltotal(ii-1)*elem_prod(yNage(ii),ywa(ii)))*mfexp(ot(ii)));

    //calculate new catlim
	double BBo;
	double SBSBo;
	dvariable TAC;
	//dvector TACnation(1,nations);


	SBSBo = seen_sbt/SBo_true;
	BBo = seen_bt/Bo_true;

	switch(hcr){

	case 1:
	//40:10 harvest control rule	

		if(SBSBo>0.4){
			//cout<<"ftarget(ii)"<<ftarget(ii)<<endl;

		TAC = (1.-mfexp(-ftarget(ii)*seltotal(ii-1)))*(elem_prod(yNage(ii),ywa(ii))*mfexp(ot(ii)));


			//TAC = value((seltotal(ii-1)*elem_prod(yNage(ii),ywa(ii)))*mfexp(ot(ii)))
			//TAC = (1.-mfexp(-ftarget(ii)))*seltotal(ii)*elem_prod(yNage(ii),ywa(ii));
			
			//TAC = elem_prod(elem_div(ftarget(ii)*seltotal(ii-1),ftarget(ii)*seltotal(ii-1)+m),elem_prod(yNage(ii),(1-mfexp(-ftarget(ii)*seltotal(ii-1)-m))))*ywa(ii)*mfexp(ot(ii)); 
			
			//cout<<"all good?? "<<endl;
		 	//cout<<"BBo "<<BBo<<endl;


		}else{
		
			if(SBSBo>0.1){

				double y;

				//y = (BBo-0.1)/(0.4-0.1);

				y = (seen_sbt-0.1*SBo_true)*((0.4/seen_sbt)/(0.4-0.1));
				
				//TAC = elem_prod(elem_div(ftarget(ii)*seltotal(ii-1),ftarget(ii)*seltotal(ii-1)+m),elem_prod(yNage(ii),(1-mfexp(-ftarget(ii)*seltotal(ii-1)-m))))*ywa(ii)*mfexp(ot(ii))*y; 
				//TAC = (1.-mfexp(-ftarget(ii)))*seen_vbt*y;
				TAC = (1.-mfexp(-ftarget(ii)*seltotal(ii-1)))*(elem_prod(yNage(ii),ywa(ii))*mfexp(ot(ii)))*y;
				//cout<<"ftarget(ii) "<<ftarget(ii)<<endl;
				//cout<<"seltotal(ii)"<<seltotal(ii)<<endl;

				//cout<<"BBo "<<BBo<<endl;
				//cout<<"recovery TAC is "<<TAC<<endl;

				
			}else{

				
				TAC = 0;

			}
		}


		if(TAC>600){
			TAC = 600.;
		}

		//cout<<"ii "<<ii<<endl<<"BBo "<<BBo<<endl;
		//cout<<"TAC "<<TAC<<endl;
			

		break;


	case 2:
	//Slope-intercept linear HCR 	

		double yint;

		//yint = - slope_hcr*intercept_hcr*SBo_true;


		if(SBSBo>intercept_hcr){

			//TAC= slope_hcr*seen_vbt;
			//TAC= slope_hcr*(seen_bt- Bo_true*intercept_hcr);
			TAC= slope_hcr* seen_bt*((seen_sbt- SBo_true*intercept_hcr)/SBo_true);


		}else{

			TAC= 0.0;

		}
		

		if(TAC>600){
			TAC = 600.;
		}

		if(TAC<0){
			TAC =0.;
		}


		break;


	}

	//cout<< "SBSBo " << SBSBo<< endl;
	//cout<< "TAC "<< ii << endl << TAC<< endl;
	catlim = TAC * nationTACprop;
	//cout<< "catlim" << endl << catlim<< endl;
	//cout<< "seltotal_SA" << endl << seltotal_SA<< endl;

	ofstream afs("catlim.txt");
	afs<<"#catlim" << endl << catlim <<endl;


FUNCTION void calc_hist_catlim(const int& ii)

	
	ofstream afs("catlim.txt");
	afs<<"#catlim" << endl << catlim_hist(ii) <<endl;



FUNCTION void run_projections(const int& pi,const int& svi)

	int svyr;
	svyr=svi;

	int p;
	p=pi;


	
	for(int ib=syr;ib<=nyr;ib++){
		if(ib==surv_yrs(svyr)) svyr++;
	}

	//catage_comm(nyr);
	
	//run_stock_assessment(nyr,svyr-1);
	

	cout<<"começa a projecao"<<endl;
	//exit(1);


	for(int ii=ntstp+1;ii<=ptstp;ii++)
	{ 

		maxPos.initialize();		
		calcmaxpos(wx(indyr(ii)));
		calc_numbers_at_age(ii,wt(indyr(ii)));



		if(indmonth(ii)==smon){
			
				
			calc_msy(indyr(ii));
			calc_Fspr_target(indyr(ii), 0.40);
			calc_known_catlim(indyr(ii));
			read_catlim();
			histTAC(indyr(ii))=catlim;
		}
		calc_position(ii);
		
		calc_effarea(ii,ntstp,catlim);
		calc_catage(ii);
		calc_catlen(ii);





		//cout<<"ii is"<<ii<<endl;
		
		for(int n=1;n<=fisharea;n++)
		{					
			if(TotEffmonth(n)(indmonth(ii))>0)
       		{
					clean_catage(ii,p,n);
					p++;
       		}	
		}
		
		

		if(indmonth(ii)==nmon){
			//survey calculations

			if(surv_yrs(svyr)==indyr(ii)){
				survey_data(svyr);
				svyr +=1;
			}
			
			catage_comm(indyr(ii));
			catlen_comm(indyr(ii));
			calc_wt_comm(indyr(ii));
			
			calc_selectivity(indyr(ii));
			//cout<<"yYieldNat(indyr(ii))"<<yYieldNat(indyr(ii))<<endl;
	
			

			
			//run_stock_assessment(indyr(ii),svyr-1);

			//exit(1);
		}
		

		
		
	}

	output_CL();
	output_true();
	//calc_spr();



	//FUNCTION save_OMrep

	//system("cd ../../../R/read_mse && make readROM");
	// Terminal year of projection.
   

	
FUNCTION read_catlim

	/* This Function reads in the catch limits in the text file
	This is usually based on the TAC. 
 	*/

	cifstream ifs_clm("catlim.txt");

	// Terminal year of projection.
    ifs_clm >> catlim;
    




FUNCTION output_true
	
	ofstream ofs("lagrangian_OM_gtg.rep");

	ofs<<"OM type" << endl << "gtg" <<endl;
	ofs<<"seed" << endl << seed <<endl;
	ofs<<"mo" << endl << mo <<endl;
	ofs<<"tau_c" << endl << tau_c<<endl;
	ofs<<"maxPos50" << endl << maxPos50 <<endl;
	ofs<<"maxPossd" << endl << maxPossd <<endl;
	ofs<<"cvPos" << endl << cvPos <<endl;
	ofs<<"fbeta " << endl << fbeta <<endl;
	ofs<<"kappa " << endl << kappa <<endl;
	ofs<<"Fmult" << endl <<  Fmult <<endl; //* exp(vt(syr,nyr)-(0.1*0.1/2))
	ofs<<"syr" << endl << syr <<endl;
	ofs<<"nyr" << endl << nyr <<endl;
	ofs<<"rep_yr" << endl << rep_yr <<endl;
	ofs<<"proj_yr" << endl << proj_yr <<endl;
	ofs<<"sage" << endl << sage <<endl;
	ofs<<"nage" << endl << nage <<endl;
	ofs<<"ngroup" << endl << ngroup <<endl;
	ofs<<"prop_ng" << endl << prop_ng <<endl;
	ofs<<"smon" << endl << smon <<endl;
	ofs<<"nmon" << endl << nmon <<endl;
	ofs<<"sarea" << endl << sarea <<endl;
	ofs<<"narea" << endl << narea <<endl;
	ofs<<"Ro " << endl << Ro <<endl;
	ofs<<"SBo " << endl << SBo <<endl;
	ofs<<"So " << endl << So <<endl;
	ofs<<"h " << endl << h <<endl;
	ofs<<"m " << endl << m <<endl;
	ofs<<"fisharea" << endl << fisharea <<endl;
	ofs<<"nations" << endl << nations <<endl;
	ofs<<"maxPos" << endl << maxPos <<endl;
	ofs<<"minPos" << endl << minPos <<endl;
	ofs<<"varPos" << endl << varPos <<endl;
	//ofs<<"PosX" << endl << PosX <<endl;	
	ofs<<"wa" << endl << wa <<endl;	
	ofs<<"ywa" << endl << ywa <<endl;	
	ofs<<"Wage" << endl << Wage <<endl;	
	ofs<<"Lage" << endl << Lage <<endl;	
	ofs<<"SB" << endl << SB <<endl;
	ofs<<"tB" << endl << tB <<endl;
	ofs<<"Ut" << endl << Ut <<endl;
	ofs<<"ytB" << endl << ytB <<endl;
	ofs<<"survB" << endl << survB <<endl;
	ofs<<"VulB" << endl << VulB <<endl;
	ofs<<"Nage" << endl << Nage <<endl;
	ofs<<"Effarea"<< endl << Effarea <<endl;
	ofs<<"comm_obsCatage"<< endl << comm_obsCatage <<endl;
	ofs<<"comm_obsCatLen"<< endl << comm_obsCatLen <<endl;
	ofs<<"surv_obsCatage"<< endl << surv_obsCatage <<endl;
	ofs<<"totVBnation" << endl << totVBnation <<endl;
	ofs<<"CatchNatAge"<< endl << CatchNatAge<<endl;
	ofs<<"indyr"<< endl << indyr<<endl;
	ofs<<"indmonth"<< endl << indmonth<<endl;
	ofs<<"indnatarea"<< endl << indnatarea<<endl;
	ofs<<"propVBarea"<< endl << propVBarea <<endl;
	ofs<<"selfisharea"<< endl << selfisharea <<endl;
	ofs<<"selnation"<< endl << selnation <<endl;
	ofs<<"seltotal"<< endl << seltotal <<endl;
	ofs<<"yCatchtotalage"<< endl << yCatchtotalage <<endl;
	ofs<<"yYieldNat"<< endl << yYieldNat <<endl;
	//ofs<<"yCatchStateAge"<< endl << yCatchStateAge <<endl;
	ofs<<"yNage"<< endl << yNage <<endl;
	//ofs<<"Fatage"<< endl << Fatage <<endl;
	ofs<<"yFatage"<< endl << yFatage <<endl;
	ofs<<"Yieldtotal"<<endl<<Yieldtotal<<endl;
	ofs<<"spr"<< endl <<spr <<endl;
	ofs<<"phie"<< endl <<phie <<endl;
	ofs<<"phiE"<< endl <<phiE <<endl;
	ofs<<"P_al"<< endl <<P_al <<endl;
	ofs<<"fmsy"<< endl <<fmsy <<endl;
	ofs<<"msy"<< endl <<msy <<endl;
	ofs<<"histTAC"<< endl <<histTAC <<endl;
	ofs<<"intercept_hcr"<< endl <<intercept_hcr <<endl;
	ofs<<"slope_hcr"<< endl <<slope_hcr <<endl;


	cout<<"OK after output_true"<<endl;


FUNCTION output_CL
	
	ofstream ofs("CL_gtg.rep");

	ofs<<"OM type" << endl << "gtg" <<endl;
	ofs<<"seed" << endl << seed <<endl;
	ofs<<"mo" << endl << mo <<endl;
	//ofs<<"tau_c" << endl << tau_c<<endl;
	ofs<<"maxPos50" << endl << maxPos50 <<endl;
	ofs<<"maxPossd" << endl << maxPossd <<endl;
	ofs<<"cvPos" << endl << cvPos <<endl;
	ofs<<"fbeta " << endl << fbeta <<endl;
	ofs<<"kappa " << endl << kappa <<endl;
	ofs<<"Fmult" << endl <<  Fmult <<endl; //* exp(vt(syr,nyr)-(0.1*0.1/2))
	ofs<<"syr" << endl << syr <<endl;
	ofs<<"nyr" << endl << nyr <<endl;
	ofs<<"rep_yr" << endl << rep_yr <<endl;
	ofs<<"proj_yr" << endl << proj_yr <<endl;
	ofs<<"sage" << endl << sage <<endl;
	ofs<<"nage" << endl << nage <<endl;
	ofs<<"ngroup" << endl << ngroup <<endl;
	ofs<<"prop_ng" << endl << prop_ng <<endl;
	ofs<<"smon" << endl << smon <<endl;
	ofs<<"nmon" << endl << nmon <<endl;
	ofs<<"sarea" << endl << sarea <<endl;
	ofs<<"narea" << endl << narea <<endl;
	ofs<<"Ro " << endl << Ro <<endl;
	ofs<<"SBo " << endl << SBo <<endl;
	ofs<<"Bo " << endl << Bo <<endl;
	ofs<<"So " << endl << So <<endl;
	ofs<<"h " << endl << h <<endl;
	ofs<<"m " << endl << m <<endl;
	ofs<<"fisharea" << endl << fisharea <<endl;
	ofs<<"nations" << endl << nations <<endl;
	ofs<<"maxPos" << endl << maxPos <<endl;
	ofs<<"minPos" << endl << minPos <<endl;
	ofs<<"varPos" << endl << varPos <<endl;
	//ofs<<"PosX" << endl << PosX <<endl;	
	ofs<<"wa" << endl << wa <<endl;	
	ofs<<"ywa" << endl << ywa <<endl;	
	//ofs<<"Wage" << endl << Wage <<endl;	
	//ofs<<"Lage" << endl << Lage <<endl;	
	ofs<<"SB" << endl << SB <<endl;
	ofs<<"tB" << endl << tB <<endl;
	ofs<<"Ut" << endl << Ut <<endl;
	ofs<<"ytB" << endl << ytB <<endl;
	ofs<<"ySB" << endl << ySB <<endl;
	ofs<<"survB" << endl << survB <<endl;
	//ofs<<"VulB" << endl << VulB <<endl;
	//ofs<<"Nage" << endl << Nage <<endl;
	//ofs<<"Effarea"<< endl << Effarea <<endl;
	ofs<<"comm_obsCatage"<< endl << comm_obsCatage <<endl;
	ofs<<"comm_obsCatLen"<< endl << comm_obsCatLen <<endl;
	ofs<<"surv_obsCatage"<< endl << surv_obsCatage <<endl;
	//ofs<<"totVBnation" << endl << totVBnation <<endl;
	//ofs<<"CatchNatAge"<< endl << CatchNatAge<<endl;
	ofs<<"indyr"<< endl << indyr<<endl;
	ofs<<"indmonth"<< endl << indmonth<<endl;
	ofs<<"indnatarea"<< endl << indnatarea<<endl;
	ofs<<"tVBarea"<< endl << tVBarea <<endl;
	//ofs<<"selfisharea"<< endl << selfisharea <<endl;
	ofs<<"selnation"<< endl << selnation <<endl;
	ofs<<"seltotal"<< endl << seltotal <<endl;
	ofs<<"yCatchtotalage"<< endl << yCatchtotalage <<endl;
	ofs<<"yYieldNat"<< endl << yYieldNat <<endl;
	//ofs<<"yCatchStateAge"<< endl << yCatchStateAge <<endl;
	ofs<<"yNage"<< endl << yNage <<endl;
	//ofs<<"Fatage"<< endl << Fatage <<endl;
	ofs<<"yFatage"<< endl << yFatage <<endl;
	ofs<<"Yieldtotal"<<endl<<Yieldtotal<<endl;
	ofs<<"spr"<< endl <<spr <<endl;
	ofs<<"phie"<< endl <<phie <<endl;
	ofs<<"phiE"<< endl <<phiE <<endl;
	ofs<<"P_al"<< endl <<P_al <<endl;
	ofs<<"fmsy"<< endl <<fmsy <<endl;
	ofs<<"msy"<< endl <<msy <<endl;
	ofs<<"histTAC"<< endl <<histTAC <<endl;
	ofs<<"intercept_hcr"<< endl <<intercept_hcr <<endl;
	ofs<<"slope_hcr"<< endl <<slope_hcr <<endl;


	cout<<"OK after output_CL"<<endl;


	
FUNCTION output_pin

	//Generate initial values as random picks from a set of likely guesses  

	
	random_number_generator rngmo(seed);
	random_number_generator rngcvPos(seed);
	random_number_generator rngmaxPos50(seed);
	random_number_generator rngmaxPossd(seed);
	random_number_generator rngFmult(seed);
	
	double tmp_mo;

	dvector guess_mo(1,6);
	dvector guess_cvPos(1,9);
	dvector guess_maxPos50(1,10);
	dvector guess_maxPossd(1,6);
	dvector guess_Fmult(1,10);


	guess_mo.fill_seqadd(0.5,0.5);
	guess_cvPos.fill_seqadd(0.06,0.01);
	guess_maxPos50.fill_seqadd(1.5,0.5);
	guess_maxPossd.fill_seqadd(1.0,0.25);
	guess_Fmult.fill_seqadd(2.6,0.1);


	
	tmp_mo 		= ceil(randu(rngmo)*(mo)+mo);
		
	ofstream ifs("../../mov_est/gtg/lagrangian_est_gtg.pin");

	ifs<<"#log_mo \n "  << log(guess_mo(ceil(randu(rngmo)*5)))<<endl;
	//ifs<<"#log_mo \n "  << log(tmp_mo) <<endl;
	//ifs<<"#mo \n "  << log(mo) <<endl;
	ifs<<"#cvPos \n" << log(guess_cvPos(ceil(randu(rngcvPos)*8))) <<endl;	
	//ifs<<"#cvPos \n" << log(cvPos) <<endl;	
	ifs<<"# maxPos50 \n" << log(guess_maxPos50(ceil(randu(rngmaxPos50)*9))) <<endl;
	//ifs<<"# maxPos50 \n" << log(maxPos50) <<endl;
	ifs<<"# maxPossd \n"<< log(guess_maxPossd(ceil(randu(rngmaxPossd)*5))) <<endl;
	//ifs<<"# maxPossd \n"<< log(maxPossd) <<endl;
	//ifs<<"# Fmult \n"<< log(Fmult) <<endl;
	ifs<<"# Fmult \n" << log(guess_Fmult(ceil(randu(rngFmult)*9))) <<endl;
	ifs<<"#wt \n" << wt(rep_yr+1,nyr)*err <<endl;
	//ifs<<"#wx \n" << wx(rep_yr+1,nyr)*0.0 <<endl;


FUNCTION output_dat

	//ofstream afs("lagrangian_est_gtg.dat");
	ofstream afs("../../mov_est/simple/lagrangian_est.dat");
	afs<<"# syr " << endl << rep_yr+1 <<endl;
	afs<<"# nyr " << endl << nyr <<endl;
	afs<<"# sage " << endl << sage <<endl;
	afs<<"# nage " << endl << nage <<endl;
	afs<<"# smon " << endl << smon <<endl;
	afs<<"# nmon " << endl << nmon <<endl;
	afs<<"# sarea " << endl << sarea <<endl;
	afs<<"# narea " << endl << narea <<endl;
	//afs<<"# ngroup " << endl << ngroup <<endl; //gtg
	afs<<"# fisharea " << endl << fisharea <<endl;
	afs<<"# fishbound " << endl << fishbound <<endl;	
	afs<<"# nations " << endl << nations <<endl;
	afs<<"# border " << endl << border <<endl;	
	afs<<"# Ro " << endl << Ro <<endl;
	afs<<"# h " << endl << h <<endl;
	afs<<"# m " << endl << m <<endl;
	afs<<"# fe " << endl << fe <<endl;
	afs<<"# q " << endl << q <<endl;
	afs<<"# fbeta " << endl << fbeta <<endl;
	afs<<"# sigR " << endl << sigR <<endl;	
	afs<<"# weight at age " << endl << wa <<endl;
	afs<<"# fecundity at age " << endl << fa <<endl;
	afs<<"# vulnerability at age " << endl << va << endl;
	afs<<"# minPos "<< endl << minPos <<endl;
	//afs<<"# Fmult "<< endl << Fmult <<endl;
	afs<<"# Total effort by country and year " << endl << trans(trans(pTotEffyear).sub(rep_yr+1,nyr))  <<endl;
	afs<<"# Total effort by country and month " << endl << TotEffmonth <<endl;
	afs<<"# effPwr"<< endl << effPwr <<endl;	
	afs<<"# dMinP " << endl << 0.1e-13<<endl;
	afs<<"# tstp month area catage " << endl << obsCatchFGAge <<endl;	
	afs<<"# yYieldtotal"<< endl << yYieldtotal(rep_yr+1,nyr) <<endl;
	afs<<"# eof " << endl << 999 <<endl;


FUNCTION output_gtgdat

	
	ofstream afs("../../mov_est/gtg/lagrangian_est_gtg.dat");
	afs<<"# syr " << endl << rep_yr+1 <<endl;
	afs<<"# nyr " << endl << nyr <<endl;
	afs<<"# sage " << endl << sage <<endl;
	afs<<"# nage " << endl << nage <<endl;
	afs<<"# smon " << endl << smon <<endl;
	afs<<"# nmon " << endl << nmon <<endl;
	afs<<"# sarea " << endl << sarea <<endl;
	afs<<"# narea " << endl << narea <<endl;
	afs<<"# ngroup " << endl << ngroup <<endl; //gtg
	afs<<"# nations " << endl << nations <<endl;
	afs<<"# border " << endl << border <<endl;	
	afs<<"# fisharea " << endl << fisharea <<endl;
	afs<<"# fishbound " << endl << fishbound <<endl;	
	afs<<"# Ro " << endl << Ro <<endl;
	afs<<"# h " << endl << h <<endl;
	afs<<"# m " << endl << m <<endl;
	afs<<"# fe " << endl << fe <<endl;
	afs<<"# q " << endl << q <<endl;
	afs<<"# fbeta " << endl << fbeta <<endl;
	afs<<"# sigR " << endl << sigR <<endl;	
	afs<<"# weight at age " << endl << wa <<endl;
	afs<<"# fecundity at age " << endl << fa <<endl;
	afs<<"# vulnerability at age " << endl << va <<endl;
	afs<<"# minPos "<< endl << minPos <<endl;
	//afs<<"# Fmult "<< endl << Fmult <<endl;
	afs<<"# Total effort by country and year " << endl << trans(trans(pTotEffyear).sub(rep_yr+1,nyr)) <<endl;
	afs<<"# Total effort by country and month " << endl << TotEffmonth <<endl;
	afs<<"# effPwr"<< endl << effPwr <<endl;	
	afs<<"# dMinP " << endl << 0.1e-9<<endl;
	afs<<"# tstp month area catage " << endl << obsCatchFGAge <<endl;	
	afs<<"# yYieldtotal"<< endl << yYieldtotal(rep_yr+1,nyr) <<endl;
	afs<<"# eof " << endl << 999 <<endl;


FUNCTION void write_LSRA_ctl_file(const int& ii,const int& svy)

	dvector wt_init(sage,nage-1);
	dvector wt_yrs(rep_yr+1,ii);

	wt_init.initialize();
	wt_yrs.initialize();
	ofstream lfs("/Users/catarinawor/Documents/length_SRA/examples/closedloop/hakelag.ctl");
	lfs<<"## ------------------------------------------------------------------------------------ ##"<< endl;
	lfs<<"## CONTROL FILE TEMPLATE                                                                ##"<< endl;
	lfs<<"## ------------------------------------------------------------------------------------ ##"<< endl;
	lfs<<"## ------------------------------------------------------------------------------------ ##"<< endl;
	lfs<<"## CONTROLS FOR LEADING PARAMETERS                                         ##"<< endl;
	lfs<<"##  Prior descriptions:                                                    ##"<< endl;
	lfs<<"##                      -0 uniform      (0,0)                              ##"<< endl;
	lfs<<"##                      -1 normal       (p1=mu,p2=sig)                     ##"<< endl;
	lfs<<"##                      -2 lognormal    (p1=log(mu),p2=sig)                ##"<< endl;
	lfs<<"##                      -3 beta         (p1=alpha,p2=beta)                 ##"<< endl;
	lfs<<"##                      -4 gamma        (p1=alpha,p2=beta)                 ##"<< endl;
	lfs<<"## ——————————————————————————————————————————————————————————————————————— ##"<< endl;
	lfs<<"## npar 															"<< endl;		
	lfs<<"4 																"<< endl;
	lfs<<"## ival    	lb      ub        phz     prior   p1      p2      #parameter   ##"<< endl;
	lfs<<"## ——————————————————————————————————————————————————————————————————————— ##"<< endl;
	//lfs<<  log(Ro)  <<"\t"<< log(Ro*0.01) <<"\t"<< log(Ro*100) <<"\t"<<" 1" <<"\t"<<"0"	<<"\t"<< log(Ro*0.01) <<"\t"<< log(Ro*100)  <<"\t"<< "## Ro"<< endl;
	lfs<<  log(Ro)  <<"\t"<< log(Ro*0.001) <<"\t"<< log(Ro*1000) <<"\t"<<" 1" <<"\t"<<"0"	<<"\t"<< log(Ro*0.001) <<"\t"<< log(Ro*1000)  <<"\t"<< "## Ro"<< endl;
	lfs<<  log(yNage(rep_yr+1)(sage)*0.8)  <<"\t"<< log(Ro*0.001) <<"\t"<< log(Ro*1000) <<"\t"<<" 1" <<"\t"<<"0"	<<"\t"<< log(Ro*0.001) <<"\t"<< log(Ro*1000)  <<"\t"<< "## R_init "<< endl;
	//lfs<<  log(Ro*0.8)  <<"\t"<< log(Ro*0.01) <<"\t"<< log(Ro*100) <<"\t"<<" 1" <<"\t"<<"0"	<<"\t"<< log(Ro*0.01) <<"\t"<< log(Ro*100)  <<"\t"<< "## R_init "<< endl;	
	lfs<<  log(kappa*0.9)  <<"\t"<< log(kappa*0.1) <<"\t"<< log(kappa*10) <<"\t"<<" 1" <<"\t"<<"1"	<<"\t"<< log(kappa) <<"\t"<< 0.9  <<"\t"<< "## log_reck"<< endl;
	lfs<<  log(1.1)  <<"\t"<< log(0.1) <<"\t"<< log(3.0) <<"\t"<<" -1" <<"\t"<<"0"	<<"\t"<< log(0.1) <<"\t"<< log(3.0)  <<"\t"<< "## "<< endl;
	lfs<<"## ——————————————————————————————————————————————————————————————————————— ##"<< endl;
	lfs<<"##initial values for recruitment deviations ##"<< endl;
	lfs<<"# wt "<< endl;
	lfs<< wt_init << endl;
	lfs<< wt_yrs << endl;
	lfs<<"#log(q) prior - same codes as above"<< endl;
	lfs<<"#prior   p1      p2 "<< endl;
	lfs<<"1		   0		 0.1"<< endl;
	lfs<<"#closedloop"<<endl<< 1 << endl;
	lfs<<"##initial values for recruitment deviations in first year ##"<< endl;
	lfs<<"999"<< endl;


FUNCTION void write_LSRA_data_file(const int& ii,const int& svy)

	
	ofstream rfs("/Users/catarinawor/Documents/length_SRA/examples/closedloop/hakelag.dat");
	rfs<<"# DATA FILE FOR Length-SRA  " << endl;
	rfs<<"## ------------------------------------------------------------------------- ##"<< endl;
	rfs<<"# syr "<< endl;
	rfs<<rep_yr+1<< endl;
	rfs<<"# eyr "<< endl;
	rfs<< ii << endl;
	rfs<<"# sage "<< endl;
	rfs<<sage<< endl;
	rfs<<"# nage "<< endl;
	rfs<< nage << endl;
	rfs<<"# slen "<< endl;
	rfs<< inilen << endl;
	rfs<<"# nlen "<< endl;
	rfs<< nlen << endl;
	rfs<<"# lstp "<< endl;
	rfs<< lenstp << endl;
	rfs<<"# SR function "<< endl;
	rfs<< 1 << endl;
	rfs<<"# m "<< endl;
	rfs<< m << endl;
	rfs<<"# alw "<< endl;
	rfs<< awt << endl;
	rfs<< "# blw" << endl;
	rfs<< bwt << endl;
	rfs<< "# w@age " << endl;
	rfs<< ywa(ii)(sage,nage) << endl;
	rfs<< "# wasw " << endl;
	rfs<< "0 " << endl;
	rfs<< "# vul " << endl;
	rfs<< va << endl;
	rfs<< "# fec " << endl;
	rfs<< elem_prod(fa,ywa(ii)(sage,nage)) << endl;
	rfs<< "# nyt " << endl;
	rfs<< svy << endl;
	rfs<< "# iyr  " << endl;
	rfs<< surv_yrs(1,svy) << endl;
	rfs<< "# it  " << endl;
	rfs<< survB(1,svy) << endl;
	rfs<< "# Clt   " << endl;
	rfs<< (comm_obsCatLen).sub(rep_yr+1,ii) << endl;
	rfs<< "# linf"  <<endl;    
	rfs<< mean(Linf) <<endl;     
	rfs<< "# vbk"  <<endl;    
	rfs<< vbk <<endl;       
	rfs<< "# to"  <<endl;    
	rfs<< to <<endl;   
	rfs<< "# cvl"  <<endl;    
	rfs<< cvl <<endl;   
	rfs<< "# cv_it"  <<endl;    
	rfs<< "0.1" <<endl;  
	rfs<< "# sigVul "  <<endl;    
	rfs<< "100" <<endl;   
	rfs<< "#uinit"  <<endl;    
	rfs<< "0.1" <<endl;   
	rfs<< "# eof " <<endl;                              
	rfs<< "999 " <<endl;   



FUNCTION void write_iscam_ctl_file(const int& ii,const int& svy)

	ofstream mfs("/Users/catarinawor/Documents/iSCAM/examples/hakelag/DATA/hakelag.ctl");
	mfs<<"## ------------------------------------------------------------------------------------ ##"<< endl;
	mfs<<"## CONTROL FILE TEMPLATE                                                                ##"<< endl;
	mfs<<"## ------------------------------------------------------------------------------------ ##"<< endl;
	mfs<<"## ------------------------------------------------------------------------------------ ##"<< endl;
	mfs<<"## CONTROLS FOR LEADING PARAMETERS                                                      ##"<< endl;
	mfs<<"##  Prior descriptions:                                                                 ##"<< endl;
	mfs<<"##                      -0 uniform      (0,0)                                           ##"<< endl;
	mfs<<"##                      -1 normal       (p1=mu,p2=sig)                                  ##"<< endl;
	mfs<<"##                      -2 lognormal    (p1=log(mu),p2=sig)                             ##"<< endl;
	mfs<<"##                      -3 beta         (p1=alpha,p2=beta)                              ##"<< endl;
	mfs<<"##                      -4 gamma        (p1=alpha,p2=beta)                              ##"<< endl;	
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## npar"<<endl<< "7"<< endl;
	mfs<<"## ival         lb      ub      phz     prior   p1      p2        #parameter            ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<< 1.2  		 <<"\t"<< -1.0 <<"\t"<< 4  	 <<"\t"<< 1  <<"\t"<< 1  <<"\t"<<  1.2 	 <<"\t"<<  0.4   <<"\t"<<"#log_ro   ##"<<endl;
   	mfs<< 0.862 	 <<"\t"<<  0.2 <<"\t"<< 1.0  <<"\t"<< 2  <<"\t"<< 3  <<"\t"<<  9.909 <<"\t"<< 2.959  <<"\t"<<"#steepness    ##"<<endl;
   	mfs<< -1.537117  <<"\t"<< -4.0 <<"\t"<< 0.0  <<"\t"<< -3 <<"\t"<< 1  <<"\t"<< -1.609 <<"\t"<< 0.1    <<"\t"<<"#log_m g&b    ##"<<endl;
   	mfs<< 2.0  		 <<"\t"<< -4.0 <<"\t"<< 10.0 <<"\t"<< 1  <<"\t"<< 0  <<"\t"<< -3.0 	 <<"\t"<< 10.0   <<"\t"<<"#log_avgrec   ##"<<endl;
   	mfs<< 1.0  		 <<"\t"<< -4.0 <<"\t"<< 10.0 <<"\t"<< 1  <<"\t"<< 0  <<"\t"<< -3.0   <<"\t"<< 10.0   <<"\t"<<"#log_recinit  ##"<<endl;
   	mfs<< 0.03764649 <<"\t"<< 0.01 <<"\t"<< 0.99 <<"\t"<< 4  <<"\t"<< 3  <<"\t"<<  5.98  <<"\t"<< 155.98 <<"\t"<<"#rho          ##"<<endl;
   	mfs<< 10.4909967 <<"\t"<< 0.01 <<"\t"<< 10.0 <<"\t"<< 3  <<"\t"<< 4  <<"\t"<<  4.01  <<"\t"<< 2.01   <<"\t"<<"#sigma_r      ##"<<endl;
	//mfs<< 0.9805792  <<"\t"<<  -3 <<"\t"<< 5.0  <<"\t"<< 1  <<"\t"<< 1  <<"\t"<<  1.2 	 <<"\t"<<  0.4   <<"\t"<<"#log_ro   ##"<<endl;
   	//mfs<< 0.862 	 <<"\t"<<  0.2 <<"\t"<< 1.0  <<"\t"<< 2  <<"\t"<< 3  <<"\t"<<  9.909 <<"\t"<< 2.959  <<"\t"<<"#steepness    ##"<<endl;
   	//mfs<< -1.537117  <<"\t"<< -4.0 <<"\t"<< 0.0  <<"\t"<< -3 <<"\t"<< 1  <<"\t"<< -1.609 <<"\t"<< 0.1    <<"\t"<<"#log_m g&b    ##"<<endl;
   	//mfs<< 2.0  		 <<"\t"<< -4.0 <<"\t"<< 10.0 <<"\t"<< 4  <<"\t"<< 0  <<"\t"<< -3.0 	 <<"\t"<< 10.0   <<"\t"<<"#log_avgrec   ##"<<endl;
   	//mfs<< 1.0  		 <<"\t"<< -4.0 <<"\t"<< 10.0 <<"\t"<< 4  <<"\t"<< 0  <<"\t"<< -3.0   <<"\t"<< 10.0   <<"\t"<<"#log_recinit  ##"<<endl;
   	//mfs<< 0.03764649 <<"\t"<< 0.01 <<"\t"<< 0.99 <<"\t"<< 4  <<"\t"<< 3  <<"\t"<<  5.98  <<"\t"<< 155.98 <<"\t"<<"#rho          ##"<<endl;
   	//mfs<< 10.4909967 <<"\t"<< 0.01 <<"\t"<< 10.0 <<"\t"<< 1  <<"\t"<< 4  <<"\t"<<  4.01  <<"\t"<< 2.01   <<"\t"<<"#sigma_r      ##"<<endl;
   	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## CONTROL PARAMETERS FOR AGE/SIZE COMPOSITION DATA FOR na_gears                        ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## Likelihood type for each gear:														  ##"<< endl;
	mfs<<"##     -1 : multivariate logistic (dmvlogistic) 										  ##"<< endl;
	mfs<<"##     -2 : multinomial, sample size based on input data 								  ##"<< endl;
	mfs<<"##     -3 : logistic_normal, no autocorrelation, AR1, AR2. 							  ##"<< endl;
	mfs<<"##     -4 : logistic_normal, AR1 														  ##"<< endl;
	mfs<<"##     -5 : logistic_normal, AR2 														  ##"<< endl;
	mfs<<"##     -6 : multinomial with estimated effective samples size (log_age_tau2 phz >0) 	  ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## Number of columns == na_gears." << endl;
	mfs<<" 1 "	 <<"\t"<< " 2 " 	<<"\t"<< "## • Gear Index"<<endl;
   	mfs<<" 1 "	 <<"\t"<< " 1 "  	<<"\t"<< "## • Likelihood type"<<endl;
   	mfs<<"0.000" <<"\t"<< "0.000"	<<"\t"<< "## • Minimum proportion for aggregation & tail compression"<<endl;
   	mfs<<"0.0000"<<"\t"<< "0.0000" 	<<"\t"<< "## • Small constant to add to comps & renormalize"<<endl;
   	mfs<<" 1 "	 <<"\t"<< " 1 " 	<<"\t"<< "## • phase for log_age_tau2 estimation. "<<endl;
   	mfs<<" 2 "	 <<"\t"<< " 2 " 	<<"\t"<< "## • phase for phi1 estimation: bounded (-1,1) AR1 "<<endl;
   	mfs<<"-2 "	 <<"\t"<< "-2 " 	<<"\t"<< "## • phase for phi2 estimation: bounded (0,1)  AR2 "<<endl;
   	mfs<<"-2 "	 <<"\t"<< "-2 " 	<<"\t"<< "## • phase for degrees of freedom for student T. "<<endl;
   	mfs<<"-12345" <<"\t"<<  "## • int check (-12345)"<<endl;
   	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## SELECTIVITY CONTROLS                                                                 ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## - Each gear must have at least one selectivity and retention curve.					"<< endl;
	mfs<<"## - Use a -ve phase with the gear index to mirror another gear.  Note 					"<< endl;
	mfs<<"##   that if you mirror another gear, it must have the same sel type and                  "<< endl;
	mfs<<"##   age and year nodes so that the arrays are the same shape & block years.				"<< endl;
	mfs<<"##																						"<< endl;
	mfs<<"## • Index       = gear index for selectivity curve.										"<< endl;
	mfs<<"## • sel_type    = type of selectivity function (see Legend). 							"<< endl;
	mfs<<"## • sel_mu      = mean age/length 50% selectivity. 										"<< endl;
	mfs<<"## • sel_sd      = std in 50% SELECTIVITY 												"<< endl;
	mfs<<"## • sex_dep     = 0 -> no;  1 -> offset for sex 2. 										"<< endl;
	mfs<<"## • size_nodes  = # of nodes for age/size cubic spline. 									"<< endl;
	mfs<<"## • year_nodes  = # of nodes for time varying bicubic spline. 							"<< endl;
	mfs<<"## • phz_mirror  = phase of estimation (-ve phase to mirror selextivity index) 			"<< endl;
	mfs<<"## • lam1        = penalty weight for 2nd differences (w = 1/(2•sig^2)). 					"<< endl;
	mfs<<"## • lam2        = penalty weight for dome-shaped selectivity. 							"<< endl;
	mfs<<"## • lam3        = penalty weight for time-varying selectivity. 							"<< endl;
	mfs<<"## • start_block = year index for first year of selectivity curve.						"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## sel_nBlocks    ret_nBlocks      ## Gear Description.						"<< endl;
	mfs<<" 2 " <<"\t"<< " 1 " <<"\t"<< "## Commercial retained 					"<< endl;
    mfs<<" 1 " <<"\t"<< " 0 " <<"\t"<< "## survey retained    					"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;					
	mfs<<"## Selectivity P(capture of all size/age)													"<< endl;
	mfs<<"## slx_dControls																			"<< endl;				
	mfs<<"## • index for sex (0=both, 1=female, 2=male)												"<< endl;
	mfs<<"##        sel   sel  sel       age    year  phz or                    start  end        ##"<< endl;
	mfs<<"## Index  type  mu   sd   sex  nodes  nodes mirror lam1  lam2  lam3 | block  block      ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
    mfs<<"1" <<"\t"<<"4" <<"\t"<<"3.5" <<"\t"<<"0.45" <<"\t"<<"0" <<"\t"<<"7" <<"\t"<<"1" <<"\t"<<"2" <<"\t"<<"150.0" <<"\t"<<"1.0" <<"\t"<<"1.0" <<"\t"<<rep_yr <<"\t"<<"1999"<<endl;
    mfs<<"1" <<"\t"<<"4" <<"\t"<<"3.5" <<"\t"<<"0.45" <<"\t"<<"0" <<"\t"<<"7" <<"\t"<<"1" <<"\t"<<"2" <<"\t"<<"150.0" <<"\t"<<"1.0" <<"\t"<<"1.0" <<"\t"<<"1999" <<"\t"<<ii<<endl;
    mfs<<"2" <<"\t"<<"2" <<"\t"<<"1.5" <<"\t"<<"0.45" <<"\t"<<"0" <<"\t"<<"4" <<"\t"<<"5" <<"\t"<<"2" <<"\t"<<"200.0" <<"\t"<<"50.0"<<"\t"<<"50.0"<<"\t"<<rep_yr <<"\t"<<ii<<endl;
    mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## Retention P(retaining size/age)													 	"<< endl;														
	mfs<<"## ret_dControls																			"<< endl;
	mfs<<"## • index for sex (0=both, 1=female, 2=male)												"<< endl;
	mfs<<"##        sel   sel  sel       age    year  phz or                    start  end        ##"<< endl;
	mfs<<"## Index  type  mu   sd   sex  nodes  nodes mirror lam1  lam2  lam3 | block  block      ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
  	mfs<<"-1"<<"\t"<<"2"<<"\t"<<"10"<<"\t"<<"2.0"<<"\t"<<"1"<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"-2"<<"\t"<<"0.0"<<"\t"<<"0.0"<<"\t"<<"0.0"<<"\t"<<syr<<"\t"<<ii <<endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## LEGEND FOR SELECTIVITY TYPES (sel_type)                                              ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## sel  | No.        | 																	"<< endl;
	mfs<<"## type | parameters | Description 														"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"##   1  |   2        | • Logistic curve with mean and standard dev at age = p(50%) 		"<< endl;
	mfs<<"##   2  | (nages-1)  | • Age-specific selectivity coefficients for (sage):(nage-1) 		"<< endl;
	mfs<<"##   3  | age_nodes  | • Age-specific coefficients based on cubic-spline interpolation 	"<< endl;
	mfs<<"##   4  | n*age_nodes| • Annual age-specific coeffs using cubic-spline interpolation 		"<< endl;
	mfs<<"##   5  | nyr*nage   | • Bicubic spline interpolation over time & age. 					"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## TO BE DEPRECATED															 "<< endl;	
	mfs<<"## ------------------------------------------------------------------------- ##"<< endl;
	mfs<<"## SELECTIVITY PARAMETERS Columns for gear                                   ##"<< endl;
	mfs<<"## OPTIONS FOR SELECTIVITY (isel_type):                                      ##"<< endl;
	mfs<<"##      1) logistic selectivity parameters                                   ##"<< endl;
	mfs<<"##      2) selectivity coefficients                                          ##"<< endl;
	mfs<<"##      3) a constant cubic spline with age-nodes                            ##"<< endl;
	mfs<<"##      4) a time varying cubic spline with age-nodes                        ##"<< endl;
	mfs<<"##      5) a time varying bicubic spline with age & year nodes.              ##"<< endl;
	mfs<<"##      6) fixed logistic (set isel_type=6, and estimation phase to -1)      ##"<< endl;
	mfs<<"##      7) logistic function of body weight.                                 ##"<< endl;
	mfs<<"##      8) logistic with weight deviations (3 parameters)                    ##"<< endl;
	mfs<<"##      11) logistic selectivity with 2 parameters based on mean length      ##"<< endl;
	mfs<<"##      12) length-based selectivity coefficients with spline interpolation  ##"<< endl;
	mfs<<"##      sig=0.05 0.10 0.15 0.20 0.30 0.40 0.50                               ##"<< endl;
	mfs<<"##      wt =200. 50.0 22.2 12.5 5.56 3.12 2.00                               ##"<< endl;
	mfs<<"## ------------------------------------------------------------------------- ##"<< endl;
  	mfs<<"5   "<<"\t"<<" 2   "<<"\t"<<"        # 1  -selectivity type ivector(isel_type) for gear"<< endl;
  	mfs<<"2.5 "<<"\t"<<" 2.5 "<<"\t"<<"        # 2  -Age/length at 50% selectivity (logistic)"<< endl;
  	mfs<<"0.45"<<"\t"<<" 0.45"<<"\t"<<"        # 3  -STD at 50% selectivity (logistic)"<< endl;
  	mfs<<"4   "<<"\t"<<" 5   "<<"\t"<<"        # 4  -No. of age nodes for each gear (0=ignore)"<< endl;
  	mfs<<"5   "<<"\t"<<" 5   "<<"\t"<<"        # 5  -No. of year nodes for 2d spline(0=ignore)"<< endl;
  	mfs<<"1   "<<"\t"<<" 1   "<<"\t"<<"        # 6  -Phase of estimation (-1 for fixed)"<< endl;
  	mfs<<"150 "<<"\t"<<" 150 "<<"\t"<<"     # 7  -Penalty wt for 2nd differences w=1/(2*sig^2)"<< endl;
  	mfs<<"2.00"<<"\t"<<" 2.00"<<"\t"<<"  # 8  -Penalty wt for dome-shaped w=1/(2*sig^2)"<< endl;
  	mfs<<"12.5"<<"\t"<<" 12.5"<<"\t"<<"    # 9  -Penalty wt for time-varying selectivity"<< endl;
  	mfs<<"1   "<<"\t"<<" 1   "<<"\t"<<"       # 10 -n_sel_blocks (number of selex blocks)"<< endl;
	mfs<<"## ------------------------------------------------------------------------- ##"<< endl;
	mfs<<"## Start year of each time block: 1 row for each gear 						 "<< endl;
	mfs<<rep_yr<< endl;
	mfs<<rep_yr<< endl;
	mfs<<"##"<< endl;
	mfs<<"##"<< endl;
	mfs<<"##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## TIME VARYING NATURAL MORTALIIY RATES                                                 ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## TYPE: "<< endl;
	mfs<<"##      0 = constant natural mortality "<< endl;
	mfs<<"##      1 = Random walk (deviates constrained by variance in M) "<< endl;
	mfs<<"##      2 = Cubic Spline (deviates constrined by nodes & node-placement)"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ## "<< endl;
    mfs<<"0"<< endl;
	mfs<<"## Phase of estimation"<< endl;
	mfs<<"  -3"<< endl;
	mfs<<"## STDEV in m_dev for Random walk"<< endl;
	mfs<<"  0.01"<< endl;
	mfs<<"## Number of nodes for cubic spline"<< endl;
	mfs<<"  0"<< endl;
	mfs<<"## Year position of the knots (vector must be equal to the number of nodes)"<< endl;
  	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## ABUNDANCE OBSERVATION MODELS 															"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## QTYPE:"<< endl;
	mfs<<"##    0 = FIXED SURVEY Q (specify log(mean) for prior log(mean))"<< endl;
	mfs<<"##    1 = CONSTANT Q     (use MLE for q and optional informative prior)"<< endl;
	mfs<<"##    2 = RANDOM WALK Q  (use prior mean & sd for penalized random walk)"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"1             # -number of surveys (n_it_nobs) 											"<< endl;
	mfs<<" 2            # -QTYPE (see legend above)													"<< endl;
	mfs<<" 0            # -prior log(mean)															"<< endl;
	mfs<<" 0.1            # -prior sd (set to 0 for uniformative prior)								"<< endl;
	mfs<<" 1            # -Estimation Phase 														"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## OTHER MISCELANEOUS CONTROLS                                                          ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"0           # 1  -verbose ADMB output (0=off, 1=on)"<< endl;
	mfs<<"1           # 2  -recruitment model (1=beverton-holt, 2=ricker)"<< endl;
	mfs<<"0.100       # 3  -std in observed catches in first phase."<< endl;
	mfs<<"0.0707      # 4  -std in observed catches in last phase."<< endl;
	mfs<<"0           # 5  -Assume unfished in first year (0=FALSE, 1=TRUE)"<< endl;
	mfs<<"0.00        # 6  -Minimum proportion to consider in age-proportions for dmvlogistic"<< endl;
	mfs<<"0.20        # 7  -Mean fishing mortality for regularizing the estimates of Ft"<< endl;
	mfs<<"0.10        # 8  -std in mean fishing mortality in first phase"<< endl;
	mfs<<"2.00        # 9  -std in mean fishing mortality in last phase"<< endl;
	mfs<<"-3          # 10 -DEPRECATED phase for estimating m_deviations (use -1 to turn off mdevs)"<< endl;
	mfs<<"0.1         # 11 -DEPRECATED std in deviations for natural mortality"<< endl;
	mfs<<"12          # 12 -DEPRECATED number of estimated nodes for deviations in natural mortality"<< endl;
	mfs<<"0.00        # 13 -fraction of total mortality that takes place prior to spawning"<< endl;
	mfs<<"0           # 14 -number of prospective years to add to syr. "<< endl;
	mfs<<"0           # 15 -switch for IFD distribution in selectivity simulations"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## MARKER FOR END OF CONTROL FILE (eofc)"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"999"<< endl;


FUNCTION void write_iscam_data_file(const int& ii,const int& svy)

	

	ofstream ufs("/Users/catarinawor/Documents/iSCAM/examples/hakelag/DATA/hakelag.dat");
	ufs<<"# DATA FILE FOR iSCAM  " << endl;
	ufs<<"## ------------------------------------------------------------------------- ##"<< endl;
	ufs<<"## MODEL DIMENSIONS  " << endl;
	ufs<<"## ------------------------------------------------------------------------- ##"<< endl;
	ufs<< 1  	<< "\t" << "# -number of areas            (narea) "  <<endl;
	ufs<< 1 	<< "\t" << "# -number of groups or stocks (ngroup)"  <<endl;
	ufs<< 1  	<< "\t" << "# -number of sexes            (nsex)"  <<endl;
	ufs<< rep_yr+1 	<< "\t" << "# -first year of data         (syr)"  <<endl;
	ufs<< ii 	<< "\t" << "# -last year of data          (nyr)"  <<endl;
	ufs<< sage 	<< "\t" << "# -age of youngest age class  (sage)"  <<endl;
	ufs<< nage 	<< "\t" << "# -age of plus group          (nage)"  <<endl;
	ufs<< 2 	<< "\t" << "# -number of gears            (ngear)"  <<endl;
	ufs<<"## ------------------------------------------------------------------------- ##"<< endl;
	ufs<<"## Allocation for each gear in (ngear), use 0 for survey gears. " << endl;
	ufs<<"## ------------------------------------------------------------------------- ##"<< endl;
	ufs<< "1 0"  <<endl;
	ufs<<"## ------------------------------------------------------------------------- ##"<< endl;
	ufs<<"## Age-schedule and population parameters                                    ##" << endl;
	ufs<<"## ------------------------------------------------------------------------- ##"<< endl;
	ufs<<"## Need one value for each area and sex.                                    " << endl;
	ufs<< "53.2  "  	<< "\t" << "# -asymptotic length (linf)  not updated         (nage)"  <<endl;
	ufs<< "0.3   "  	<< "\t" << "#-brody growth coefficient (k) not updated"  <<endl;
	ufs<< "-0.5  "   	<< "\t" << "#-theoretical age at zero length (to)"  <<endl;
	ufs<< "5e-6  "   	<< "\t" << "#-scaler in length-weight allometry"  <<endl;
	ufs<< "3.0  "		<< "\t" << "#-power parameter in length-weight allometry"  <<endl;
	ufs<< "3.45  "  	<< "\t" << "#-age at 50% maturity (approx with log(3.0)/k)"  <<endl;
	ufs<< "0.35  "  	<< "\t" << "#-std at 50% maturity (CV ~ 0.1)"  <<endl;
	ufs<< "1"  	<< "\t" << "#flag mat vec  (if nmat==1) then read this vector, else if nmat==0, ignore it."  <<endl;			
	ufs<< fa 	<< "\t" << "#mat vec"  <<endl;		
	ufs<<"## ------------------------------------------------------------------------- ##"<< endl;
	ufs<<"## Aging Error vectors (mean,sd,sage,nage)                                   ##" << endl;
	ufs<<"## ------------------------------------------------------------------------- ##"<< endl;
	ufs<< 1 	<< "\t" << "# Number of ageing error_definitions"  <<endl;		
	ufs<< "# 1 - No error"  <<endl;	
	ufs<< age  <<endl;		
	ufs<< "0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01"  <<endl;		
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<<"## TIME SERIES data                                						   ##" << endl;
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<<"## Observed catch from all gears, areas, and sex                             ##" << endl;
	ufs<<"## sex: 1=female, 2=male, 0=asexual                         				   ##" << endl;
	ufs<<"##               1 = catch in weight                                         ##" << endl;
	ufs<<"##               2 = catch in numbers                                        ##" << endl;
	ufs<<"##               3 = catch in spawn (roe)                                    ##" << endl;
	ufs<<"## n_ct_obs" << endl;
	ufs<<  ii -rep_yr <<endl;
	ufs<<"## Year gear area group sex type value" << endl;
	for(int i=rep_yr+1;i<=ii;i++){
		//   year 		gear 		area 	  group       sex       type   		value             se
		ufs<<  i <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 0 <<"\t"<< 1 <<"\t"<<yYieldtotal(i)<<"\t"<<0.1<<endl;
	}
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<<"## ABUNDANCE INDICES -A RAGGED ARRAY: (1,nit,1,nit_nobs,1,5)                 ##" << endl;
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<< 1 	<< "\t" << "# Number of abundance series         int(nit)"  <<endl;
	ufs<< svy 	<< "\t" << "# Number of observations in series   ivector(nit_nobs(1,nit))"  <<endl;
	ufs<< 2 	<< "\t" << "#Survey type (see key below)        ivector(survey_type(1,nit))"  <<endl;
	ufs<<"## 1 = survey is proportional to vulnerable numbers" << endl;
	ufs<<"## 2 = survey is proportional to vulnerable biomass" << endl;
	ufs<<"## 3 = survey is proportional to spawning biomass (e.g., a spawn survey)" << endl;
	ufs<<"## iyr     index(it) gear area group sex log(SE) log(PE)   timing" << endl;
	for(int l=1;l<=svy;l++){
		//  	iyr     		index(it) 		  gear       area 		group 		sex 	log(SE)  		log(PE)   	timing
		ufs<<surv_yrs(l) <<"\t"<< survB(l) <<"\t"<< 2 <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 0 <<"\t"<< 0.2<<"\t"<< 0.01<<"\t"<< 0.5 <<endl;
	}
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<<"## AGE COMPOSITION DATA (ROW YEAR, COL=AGE) Ragged object                    ##" << endl;
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<< 2 	<< "\t" << "# Number of gears with age-comps int(na_gears)"  <<endl;
	ufs<< ii-rep_yr << "\t" << svy	<< "\t" << "# Number of rows in the matrix   ivector(na_gears)"  <<endl;
	ufs<< sage << "\t" << sage	<< "\t" << "## ivector(na_gears) of youngest age-class"  <<endl;
	ufs<< nage << "\t" << nage	<< "\t" << "## ivector(na_gears) of oldest age-class + group"  <<endl;
	ufs<< 10 << "\t" << 10	<< "\t" << "## effective sample size for multinomial"  <<endl;
	ufs<< 1	<< "\t" << 1 << "\t" << "## Age composition flag (1==age comps, 0==length comps)"  <<endl;
	ufs<<"## year gear area group sex age_err | data columns (numbers or proportions)" << endl;
	for(int i=rep_yr+1;i<=ii;i++){
		//  year 	  gear 		area 		group 		sex 	age_err | data columns (numbers or proportions)
		ufs<<i <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 0 <<"\t"<< 1 <<"\t"<< comm_obsCatage(i)(sage,nage) <<endl;
	}
	for(int l=1;l<=svy;l++){	
		//  	year 	  		gear 	area 		group 		sex 	age_err | data columns (numbers or proportions)
		ufs<<surv_yrs(l)<<"\t"<< 2 <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 0 <<"\t"<< 1 <<"\t"<< surv_obsCatage(l)(sage,nage) <<endl;
	}	
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<<"## EMPIRICAL WEIGHT-AT-AGE DATA                                              ##" << endl;
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<<"## Number of weight-at-age tables (n_wt_tab)                                 ##" << endl;
	ufs<< 1 <<endl;
	ufs<<"## Number of rows in each weight-at-age table vector(n_wt_obs), use -99 if NA ##" << endl;
	ufs<< ii-rep_yr <<endl;
	ufs<<"## year gear area stock sex |age columns (sage, nage) of weight at age data   ##" << endl;
	for(int i=rep_yr+1;i<=ii;i++){
		ufs<<i <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 0 <<"\t"<< Watage_comm(i)(sage,nage) <<endl;
	}
	ufs<<"##1975 2     1     1    0  0.0550 0.1575 0.2987 0.3658 0.6143 0.6306 0.7873 0.8738 0.9678 0.9075 0.9700 1.6933 1.5000 1.9000 1.9555 2.7445 2.7445 2.7445 2.7445 2.7445 2.7445" << endl;
 	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<<"## MARKER FOR END OF DATA FILE (eof)                                         ##" << endl;
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<< 999 <<endl;                                
	

REPORT_SECTION

	//REPORT(Effarea);



TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 10000000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e8);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e8);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 

GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;

	#undef TINY
	#define TINY 1.e-08

	#include <admodel.h>
	#include <time.h>
	#include <statsLib.h>
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;

