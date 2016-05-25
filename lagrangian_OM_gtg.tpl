
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
		seed += 10; // add 1 to the seed
		ofstream ofs( "seed.txt" ); //put out to seed.txt
		ofs<<seed<<endl; //the new value of the seed
	END_CALCS

	//======================
	// model dimensions
	//======================

	init_int syr;
	init_int nyr;
	init_int rep_yr;
	init_int sage;
	init_int nage;
	init_int smon;
	init_int nmon;
	init_int sarea;
	init_int narea;

	init_int ngroup;


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
	init_number sigR;		
	
	init_number tau_c; 		
	init_number err;

	init_number fbeta;
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

	init_vector wa(sage,nage);
	init_vector fa(sage,nage);
	init_vector va(sage,nage);
	init_vector minPos(sage,nage);

						

	init_matrix TotEffyear(1,fisharea,syr,nyr);
	init_matrix TotEffmonth(1,fisharea,smon,nmon);

	init_vector effPwr(sarea,narea);

	init_int eof;
	
	

	LOC_CALCS
		
		if( eof != 999 )
		{
			cout<<"TotEffyear "<<TotEffyear<<endl;
			cout<<"effPwr "<<effPwr<<endl;
			cout<<"Error reading data.\n Fix it."<<endl;
			cout<< "eof is: "<<eof<<endl;
			ad_exit(1);
		}

	END_CALCS
	

	//===================================
	// accessory quantities and counters
	//===================================
	
	int ntstp;

   	vector age(sage,nage);
   	vector areas(sarea,narea);
   	ivector fishingr(1,fisharea);
   	ivector nationareas(1,nations);
   	vector wt(syr,nyr);

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
		
		ntstp = (nmon-smon+1) * (nyr-syr+1);
		age.fill_seqadd(sage,1);
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
			wt.fill_randn(rng);
			wt*=sigR;



	END_CALCS

	 	ivector indyr(1,ntstp);
 		ivector indmonth(1,ntstp);
		ivector indnatarea(sarea,narea);
		ivector indfisharea(sarea,narea);
		ivector pcat(1,fisharea);
		vector ntmp(1,fisharea+1);
		vector ntmp1(1,nations+1);

		int tot_pcat;

		matrix TotEffyear_rep(1,fisharea,rep_yr+1,nyr);

		number delta;
		vector Xini(1,ngroup);
		
	

       LOC_CALCS		
       			int aa =0;
       			
       			for(int y=syr;y<=nyr;y++)		
       			{
       				for(int ii=smon;ii<=nmon;ii++)
       				{
       					aa++;
       					indyr(aa) = y;
       					indmonth(aa) = ii;
       				}
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

       			cout<<"Ok after indfisharea calcs"<< endl;

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
       			for(int n=1;n<=fisharea;n++)
       			{
       				for(int ii=rep_yr+1;ii<=nyr;ii++)
       				{
       					TotEffyear_rep(n)(ii)= TotEffyear(n)(ii);
       				}

       				for(int i=rep_yr*nmon+1;i<=ntstp;i++)
       				{

       					if(TotEffmonth(n)(indmonth(i))>0.0)
       					{
       						pcat(n)++;
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

			
       			
	END_CALCS


PARAMETER_SECTION

	objective_function_value no_f;

	//derived quantities
	number kappa;
	number phie;
	number So;
	number Bo
	number beta;

	number m_tsp;
	
	vector lxo(sage,nage);
	vector za(sage,nage);
	vector SB(1,ntstp);
	vector tB(1,ntstp);
	
	vector maxPos(sage,nage);
	vector varPos(sage,nage);
	vector varPosg(sage,nage);

	vector prop_ng(1,ngroup);

	//matrix maxPos(1,ngroup,sage,nage);
	//matrix varPos(1,ngroup,sage,nage);


	matrix NationVulB(1,ntstp,1,nations);
	
 	matrix Effarea(1,ntstp,sarea,narea);
 	matrix tnage(1,ntstp,sage,nage);
 	matrix tVBarea(1,ntstp,sarea,narea);
 	matrix meanPosX(1,ntstp,sage-1,nage);
 	
 	 	
 	
	3darray Nage(1,ngroup,1,ntstp,sage-2,nage); 
 	3darray VulB(1,ngroup,1,ntstp,sage-2,nage);
	
 	3darray VBarea(1,ngroup,1,ntstp,sarea,narea);
 	3darray totB(1,ngroup,1,ntstp,sage,nage);
 	
 	3darray PosX(1,ngroup,1,ntstp,sage-1,nage);
 	3darray NAreaAge(1,ntstp,sarea,narea,sage,nage);
 	3darray CatchAreaAge(1,ntstp,sarea,narea,sage,nage);
 	3darray CatchNatAge(1,ntstp,1,fisharea,sage-2,nage);
 	3darray EffNatAge(1,fisharea,1,ntstp,sage-2,nage);
 	

 	4darray propVBarea(1,ntstp,1,ngroup,sarea,narea,sage-3,nage);

 	matrix obsCatchNatAge(1,tot_pcat,sage-3,nage);

PRELIMINARY_CALCS_SECTION

	incidence_functions();
	initialization();
	move_grow_die();
	clean_catage();
	
	output_true();
	output_dat();
	output_pin();
	
	exit(1);



PROCEDURE_SECTION

FUNCTION dvar_vector cnorm(const double& x, const dvar_vector& mu, const dvar_vector& sd)

	dvar_vector rst(sage,nage);
	dvar_vector stx(sage,nage);

	for(int a= sage; a<= nage; a++)
	{
		stx(a) = (x-mu( a ))/sd( a );
		rst(a) = cumd_norm(stx( a ));
	}

	return(rst);

FUNCTION dvar_vector calc_InitPos_gtg()
	
	dvector inimu(sage,nage);
	dvector Xtemp(1,ngroup);
	
	Xtemp(1)=value(Xini(1)*varPos(sage));
	inimu.initialize();

	prop_ng(1) = cnorm(Xtemp(1)+0.5,inimu,varPos)(sage)-cnorm(Xtemp(1)-0.5,inimu,varPos)(sage);
	for(int g=2; g<=ngroup ;g++)
	{
		Xtemp(g)=value(Xini(g)*varPos(sage));
		prop_ng(g) = cnorm(Xtemp(g)+0.5,inimu,varPos)(sage)-cnorm(Xtemp(g)-0.5,inimu,varPos)(sage);
	}

	prop_ng = prop_ng/sum(prop_ng);
	
	return(prop_ng);


FUNCTION incidence_functions

	
	lxo = mfexp(-m*age);
	lxo(nage) /= 1. - mfexp(-m); 

	kappa 	= 4*h/(1-h);
	phie	= lxo*fa;
	So 		= kappa/phie;
	Bo 		= kappa/So*Ro;
	beta 	= (kappa-1)/Bo;


	m_tsp = m/nmon;
	za 	= m_tsp+va*fe;

	maxPos.initialize();
	calcmaxpos();
	varPos=maxPos*cvPos;
	varPosg=sqrt((varPos*varPos)/(ngroup*ngroup*4));

	//cout<<"Ok after incidence_functions"<<endl;  

FUNCTION void calc_numbers_at_age(const int& ii)

	dvar_vector propBarea(sarea,narea);
	
	for(int g=1;g<=ngroup;g++)
	{
			Nage(g)(ii)(sage-2)=ii;
			Nage(g)(ii)(sage-1)=g;	
			
			
			switch (indmonth(ii)){
            	
            	case 1:  
            		
            		Nage(g)(ii)(sage) =((So*SB(ii-nmon)/(1.+beta*SB(ii-nmon)))*mfexp(wt(indyr(ii))*err))* prop_ng(g);//*1.0/dg;

            		for(int a = sage+1;a<=nage;a++)
            		{
            		
						propBarea.initialize();
			
						for(int rr =sarea; rr<=narea; rr++)
						{		
							propBarea(rr) = (cnorm(areas(rr)+0.5,PosX(g)(ii-1),varPosg)-cnorm(areas(rr)-0.5,PosX(g)(ii-1),varPosg))(a-sage);	
						}

						Nage(g)(ii)(a) = Nage(g)(ii-1)(a-1)*propBarea*mfexp(-(m_tsp+q*Effarea(ii-1)*va(a-1))) +
								  Nage(g)(ii-1)(a-1)*(1.0-sum(propBarea))*mfexp(-(m_tsp));
					}

					Nage(g)(ii)(nage) = sum(elem_div(elem_prod((Nage(g)(ii-1)(nage-1)*propBarea),mfexp(-(m_tsp+q*Effarea(ii-1)*va(nage-1)))),
            					(1.-mfexp(-(m_tsp+q*Effarea(ii-1)*va(nage))))))+
            					(Nage(g)(ii-1)(nage-1)*(1.0-sum(propBarea))*mfexp(-m_tsp))/(1.-mfexp(-m_tsp));
            	
            		//cout<<"Nage(g)(ii)(sage) "<<Nage(g)(ii)(sage)<<endl;
            		//exit(1);

            	break;
            	
            	default: 

            		for(int a = sage;a<=nage;a++)
            		{
            			propBarea.initialize();
			
						for(int rr =sarea; rr<=narea; rr++)
						{		
							propBarea(rr) = (cnorm(areas(rr)+0.5,PosX(g)(ii-1),varPosg)-cnorm(areas(rr)-0.5,PosX(g)(ii-1),varPosg))(a-sage+1);	
						}

            			Nage(g)(ii)(a) = Nage(g)(ii-1)(a)*propBarea*mfexp(-(m_tsp+q*Effarea(ii-1)*va(a)))+
            						  	 Nage(g)(ii-1)(a)*(1.0-sum(propBarea))*mfexp(-(m_tsp));
            		}
			
				
            
            	break;
        	}
            	

        	VulB(g)(ii)(sage-2) = ii;
        	VulB(g)(ii)(sage-1) = g;
	
			VulB(g)(ii)(sage,nage) = elem_prod(elem_prod(Nage(g)(ii)(sage,nage),va),wa);
	
			tnage(ii)(sage,nage) += Nage(g)(ii)(sage,nage);
			tB(ii) += Nage(g)(ii)(sage,nage)*wa;
			totB(g)(ii)(sage,nage) = elem_prod(Nage(g)(ii)(sage,nage),wa);
			


	}
	
	SB(ii) = elem_prod(tnage(ii),fa)*wa/2.0;


	//cout<<"Ok after calc_numbers_at_age"<<endl;
	
FUNCTION void calc_effarea(const int& ii)

	dvar_vector tmp1(sarea,narea);
	dvar_vector tmp2(sarea,narea);
	tmp1.initialize();
	tmp2.initialize();

	for(int n=1; n<=nations;n++){
       	NationVulB(ii,n) = sum(tVBarea(ii)(ntmp1(n),ntmp1(n+1)-1));		 
	}

	for(int r= sarea; r<=narea; r++)
	{
		tmp1(r)= pow((tVBarea(ii)(r)/((NationVulB(ii)(indnatarea(r)))+0.01))+0.0000001,fbeta) * effPwr(r);
		tmp2(r) = tmp1(r)*TotEffyear(indfisharea(r))(indyr(ii));
		Effarea(ii)(r) = tmp2(r)*TotEffmonth(indfisharea(r))(indmonth(ii));
	}
	//cout<<"Ok after calc_effarea"<<endl;


FUNCTION void calc_position(const int& ii)


	meanPosX(ii)(sage,nage) = minPos + (maxPos - minPos) * (0.5+0.5*sin(indmonth(ii)*PI/6 - mo*PI/6)); 



	int g, r, ig;

	for(ig=1;ig<=n_rg ;ig++)
	{
		
		g = n_group(ig);
		r = n_area(ig);

		
		PosX(g)(ii)(sage-1) = g;
		PosX(g)(ii)(sage,nage) =  Xini(g)*varPos+meanPosX(ii)(sage,nage);
		

		
		VBarea(g)(ii)(r) = VulB(g)(ii)(sage,nage) * (cnorm(areas(r)+0.5,PosX(g)(ii),varPosg)-cnorm(areas(r)-0.5,PosX(g)(ii),varPosg));
		NAreaAge(ii)(r) += elem_prod(Nage(g)(ii)(sage,nage),(cnorm(areas(r)+0.5,PosX(g)(ii),varPosg)-cnorm(areas(r)-0.5,PosX(g)(ii),varPosg)));
		tVBarea(ii)(r) += VBarea(g)(ii)(r);		

		propVBarea(ii)(g)(r)(sage-3) = ii;
		propVBarea(ii)(g)(r)(sage-2) = g;
		propVBarea(ii)(g)(r)(sage-1) = r;
		propVBarea(ii)(g)(r)(sage,nage) =  elem_prod(VulB(g)(ii)(sage,nage), (cnorm(areas(r)+0.5,PosX(g)(ii),varPosg)-cnorm(areas(r)-0.5,PosX(g)(ii),varPosg)));
			
		
	
	}

	//cout<<"Ok after calc_position"<<endl;


FUNCTION void calc_catage(const int& ii)

		int a,r;
		
		for(int r=sarea;r<=narea;r++)
		{

			CatchNatAge(ii)(indfisharea(r))(sage-2) = ii;
			CatchNatAge(ii)(indfisharea(r))(sage-1) = indfisharea(r);
			
			for(int a = sage; a<=nage;a++)
			{

				CatchAreaAge(ii)(r)(a) = q*Effarea(ii)(r)*va(a)/(q*Effarea(ii)(r)*va(a)+m_tsp)*(1-mfexp(-(q*Effarea(ii)(r)*va(a)+m_tsp)))*NAreaAge(ii)(r)(a);
				CatchNatAge(ii)(indfisharea(r))(a) += CatchAreaAge(ii)(r)(a);

			}

		}
		//cout<<"Ok after calc_catage"<<endl;
	
FUNCTION initialization
			
	NAreaAge.initialize();
 	CatchAreaAge.initialize();
 	CatchNatAge.initialize();
 	Nage.initialize();

 	calc_InitPos_gtg();

 	for(int g=1;g<=ngroup;g++)
	{
		Nage(g)(1)(sage-1)=g;
		Nage(g)(1)(sage-2)=1;
	
			
		Nage(g)(1)(sage) =(So*Bo/(1+beta*Bo))* prop_ng(g);

        for(int a=sage+1 ; a <= nage ; a++)
		{
			Nage(g)(1)(a) = Nage(g)(1)(a-1) * mfexp(-za(a-1));
		}
		Nage(g)(1)(nage) /= (1.-mfexp(-za(nage)));

		VulB(g)(1)(sage-2) = 1;         		          	
        VulB(g)(1)(sage-1) = g;
		VulB(g)(1)(sage,nage) = elem_prod(elem_prod(Nage(g)(1)(sage,nage),va),wa);
	
		tnage(1)(sage,nage) += Nage(g)(1)(sage,nage);
		tB(1) += Nage(g)(1)(sage,nage)*wa;
		totB(g)(1)(sage,nage) = elem_prod(Nage(g)(1)(sage,nage),wa);
	}

	
	SB(1) = elem_prod(tnage(1),fa)*wa/2.0;

	
	calc_position(1);
	
	calc_effarea(1);
	
	calc_catage(1);

	//cout<<"Ok after initialization"<<endl;
 

FUNCTION move_grow_die

	for(int i=2;i<=ntstp;i++)
	{
		calc_InitPos_gtg();
		calc_numbers_at_age(i);		
		
		maxPos.initialize();		
		calcmaxpos();
		varPos=maxPos*cvPos;
		varPosg=sqrt((varPos*varPos)/(ngroup*ngroup*4));

		
		calc_position(i);
		calc_effarea(i);
		calc_catage(i);
		//cout<<"i is "<<i<<endl;

	}

	//cout<<"Ok after move_grow_die"<<endl;

FUNCTION clean_catage
	
	int p;
	double tiny=1.e-10;
       
	p=1;
	for(int i=rep_yr*nmon+1;i<=ntstp;i++)
	{
		for(int n=1;n<=fisharea;n++)
		{
							
			if((TotEffmonth(n)(indmonth(i)))>0)
       		{
       			
       			dvector pa(nage,sage);
       			pa.initialize();

       			obsCatchNatAge(p)(sage-3) = i;
       			obsCatchNatAge(p)(sage-2) = indmonth(i);
				obsCatchNatAge(p)(sage-1) = n;
				//if(sum(CatchNatAge(i)(n)(sage,nage))>0.0000){
					pa = value((CatchNatAge(i)(n)(sage,nage))/sum(CatchNatAge(i)(n)(sage,nage)+0.01))+0.000001;
					obsCatchNatAge(p)(sage,nage) = rmvlogistic(pa,tau_c,seed+i);
					p++;	
				//}
       			
				//cout<<"obsCatchNatAge(p)(sage,nage)"<<obsCatchNatAge(p)(sage,nage)<<endl;
       		}	
		}
	}

	//cout<<"Ok after clean_catage"<<endl;
	
FUNCTION dvar_vector calcmaxpos()
	
		 					
		maxPos(sage,nage) = 1./(1.+mfexp(-(age-maxPos50)/maxPossd));
		maxPos(sage,nage) *= (narea-sarea);
		maxPos(sage,nage) += sarea;		
	
		return(maxPos);
		
	//cout<<"Ok after calcmaxpos"<<endl;
	


FUNCTION output_true
	
	ofstream ofs("lagrangian_OM_gtg.rep");

	ofs<<"mo" << endl << mo <<endl;
	ofs<<"tau_c" << endl << tau_c<<endl;
	ofs<<"maxPos50" << endl << maxPos50 <<endl;
	ofs<<"maxPossd" << endl << maxPossd <<endl;
	ofs<<"cvPos" << endl << cvPos <<endl;
	ofs<<"seed"<< endl << seed <<endl;
	ofs<<"syr" << endl << syr <<endl;
	ofs<<"nyr" << endl << nyr <<endl;
	ofs<<"rep_yr" << endl << rep_yr <<endl;
	ofs<<"sage" << endl << sage <<endl;
	ofs<<"nage" << endl << nage <<endl;
	ofs<<"ngroup" << endl << ngroup <<endl;
	ofs<<"prop_ng" << endl << prop_ng <<endl;
	ofs<<"smon" << endl << smon <<endl;
	ofs<<"nmon" << endl << nmon <<endl;
	ofs<<"sarea" << endl << sarea <<endl;
	ofs<<"narea" << endl << narea <<endl;
	ofs<<"nations" << endl << nations <<endl;
	ofs<<"maxPos" << endl << maxPos <<endl;
	ofs<<"minPos" << endl << minPos <<endl;
	ofs<<"varPos" << endl << varPos <<endl;
	ofs<<"PosX" << endl << PosX <<endl;	
	ofs<<"SB" << endl << SB <<endl;
	ofs<<"VulB" << endl << VulB <<endl;
	ofs<<"Nage" << endl << Nage <<endl;
	ofs<<"VBarea" << endl << VBarea <<endl;
	ofs<<"Effarea"<< endl << Effarea <<endl;
	ofs<<"EffNatAge"<< endl << EffNatAge<<endl;
	ofs<<"CatchNatAge"<< endl << CatchNatAge<<endl;
	ofs<<"CatchAreaAge"<< endl << CatchAreaAge<<endl;
	ofs<<"indyr"<< endl << indyr<<endl;
	ofs<<"indmonth"<< endl << indmonth<<endl;
	ofs<<"indnatarea"<< endl << indnatarea<<endl;
	//ofs<<"propVBarea"<< endl << propVBarea <<endl;
	ofs<<"propVBarea"<< endl << propVBarea <<endl;



	
FUNCTION output_pin

	//Generate initial values as random picks from a set of likely guesses  

	
	random_number_generator rngmo(seed);
	random_number_generator rngcvPos(seed);
	random_number_generator rngmaxPos50(seed);
	random_number_generator rngmaxPossd(seed);
	
	double tmp_mo;

	
	dvector guess_cvPos(1,6);
	dvector guess_maxPos50(1,10);
	dvector guess_maxPossd(1,8);

	
	guess_cvPos.fill_seqadd(0.05,0.05);
	guess_maxPos50.fill_seqadd(3,0.5);
	guess_maxPossd.fill_seqadd(0.5,0.5);

	
	tmp_mo 		= ceil(randu(rngmo)*(mo+3));
		
	ofstream ifs("lagrangian_est.pin");

	ifs<<"#log_mo \n "  << log(tmp_mo) <<endl;
	//ifs<<"#log_mo \n "  << log(mo) <<endl;
	ifs<<"#cvPos \n" << log(guess_cvPos(ceil(randu(rngcvPos)*5))) <<endl;	
	//ifs<<"#cvPos \n" << log(cvPos) <<endl;	
	ifs<<"# maxPos50 \n" << log(guess_maxPos50(ceil(randu(rngmaxPos50)*9))) <<endl;
	//ifs<<"# maxPos50 \n" << log(maxPos50) <<endl;
	ifs<<"# maxPossd \n"<< log(guess_maxPossd(ceil(randu(rngmaxPossd)*7))) <<endl;
	//ifs<<"# maxPossd \n"<< log(maxPossd) <<endl;
	ifs<<"#wt \n" << wt(rep_yr+1,nyr)*err <<endl;


FUNCTION output_dat

	//ofstream afs("lagrangian_est_gtg.dat");
	ofstream afs("lagrangian_est.dat");
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
	afs<<"# vulnerability at age " << endl << va <<endl;
	afs<<"# minPos "<< endl << minPos <<endl;
	afs<<"# Total effort by country and year " << endl << TotEffyear_rep <<endl;
	afs<<"# Total effort by country and month " << endl << TotEffmonth <<endl;
	afs<<"# effPwr"<< endl << effPwr <<endl;	
	afs<<"# dMinP " << endl << 0.1e-13<<endl;
	afs<<"# tstp month area catage " << endl << obsCatchNatAge <<endl;	
	afs<<"# eof " << endl << 999 <<endl;

	//afs<<"# syr " << endl << rep_yr+1 <<endl;
	//afs<<"# nyr " << endl << nyr <<endl;
	//afs<<"# sage " << endl << sage <<endl;
	//afs<<"# nage " << endl << nage <<endl;
	//afs<<"# smon " << endl << smon <<endl;
	//afs<<"# nmon " << endl << nmon <<endl;
	//afs<<"# sarea " << endl << sarea <<endl;
	//afs<<"# narea " << endl << narea <<endl;
	//afs<<"# nations " << endl << nations <<endl;
	//afs<<"# border " << endl << border <<endl;
	//afs<<"# ngroup "<< endl << ngroup <<endl;
	//afs<<"# prop_ng" << endl << prop_ng <<endl;
	//afs<<"# Ro " << endl << Ro <<endl;
	//afs<<"# h " << endl << h <<endl;
	//afs<<"# m " << endl << m <<endl;
	//afs<<"# fe " << endl << fe <<endl;
	//afs<<"# q " << endl << q <<endl;
	//afs<<"# sigR " << endl << sigR <<endl;
	//afs<<"# weight at age " << endl << wa <<endl;
	//afs<<"# fecundity at age " << endl << fa <<endl;
	//afs<<"# vulnerability at age " << endl << va <<endl;
	//afs<<"# minPos "<< endl << minPos <<endl;
	//afs<<"# Total effort by country and year " << endl << TotEffyear_rep <<endl;
	//afs<<"# Total effort by country and month " << endl << TotEffmonth <<endl;
	//afs<<"# effPwr"<< endl << effPwr <<endl;
	//afs<<"# dMinP " << endl << 0.1e-15 <<endl;
	//afs<<"# tstp month area catage " << endl << obsCatchNatAge <<endl;	
	//afs<<"# eof " << endl << 999 <<endl;
	

REPORT_SECTION

	REPORT(VBarea);
	REPORT(Effarea);



TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 900000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
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

