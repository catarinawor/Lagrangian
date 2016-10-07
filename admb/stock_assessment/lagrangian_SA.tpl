//  ******************************************************************
//  Lagragian movement model
//  
//  Created by Catarina Wor on 2015-04-13.
//  Copyright (c) 2015. All rights reserved.
//  Comments:
//  ******************************************************************


DATA_SECTION

	//======================
	// model dimensions
	//======================

	init_int syr;
	init_int nyr;
	init_int sage;
	init_int nage;
	init_int smon;
	init_int nmon;
	init_int sarea;
	init_int narea;

	init_int fisharea;
	init_vector fishbound(1,fisharea-1);
	
	init_int nations;
	init_vector border(1,nations-1);



	//======================
	//model parameters -- some are now estimated
	//======================

	//init_number Ro;//now estimated
	//init_number h;
	init_number m;
	init_number fe;
	init_number q;
	init_number fbeta;
	init_number sigR;
	

	init_vector wa(sage,nage);
	init_vector fa(sage,nage);
	init_vector va(sage,nage);
	init_vector minPos(sage,nage);
	
	init_number Fmult;
	init_matrix pTotEffyear(1,fisharea,syr,nyr);
	init_matrix TotEffmonth(1,fisharea,smon,nmon);

	init_vector effPwr(sarea,narea);
	
	//======================
	//control parameters
	//======================

	init_number dMinP;  //minimum value for multivariare logistic likelihood
	
	//===================================
	// accessory quantities and counters
	//===================================
	
	int ntstp;

   	vector age(sage,nage);
   	vector areas(sarea,narea);
   	ivector fishingr(1,fisharea);
   	ivector nationareas(1,nations);
   	
   	//int n_rg;
   	//!!n_rg  = (narea-sarea+1) * ngroup;
   	
   	//ivector   n_area(1,n_rg);
	//ivector  n_group(1,n_rg);

   	//imatrix  pntr_rg(sarea,narea,1,ngroup);
   		
   	LOC_CALCS

		//int ih,r,g;
		//ih = 0;

		//for(r=sarea; r<=narea; r++)
		//{
		//	for(g=1; g<=ngroup; g++)
		//	{
		//		ih ++;
		//		pntr_rg(r,g) = ih;
		//		n_area(ih)  = r;
		//		n_group(ih) = g;
		//	}
		//}

			ntstp = (nmon-smon+1) * (nyr-(syr-20)+1);

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

			
			

	END_CALCS

	 	ivector indyr(1,ntstp);
 		ivector indmonth(1,ntstp);
		ivector indnatarea(sarea,narea);
		ivector indfisharea(sarea,narea);
		ivector pcat(1,fisharea);
		vector ntmp(1,fisharea+1);
		vector ntmp1(1,nations+1);

		int tot_pcat;
		int itsp;
		

       LOC_CALCS

       			itsp = (nmon-smon +1)*20+1;		
       			int aa =0;
       			
       			for(int y=syr-20;y<=nyr;y++)		
       			{
       				for(int ii=smon;ii<=nmon;ii++)
       				{
       					aa++;
       					indyr(aa) = y;
       					indmonth(aa) = ii;
       				}
       			}

       			ntmp(1) = sarea;
       			
       			for(int n=1;n<=fisharea;n++)
       			{
       				ntmp(n+1)= ntmp(n)+fishingr(n);	
       				indfisharea(ntmp(n),ntmp(n+1)-1)=n;
       			}

       			ntmp(fisharea+1) = narea;
       			indfisharea(narea) = fisharea;

       			


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
       				for(int i=itsp;i<=ntstp;i++)
       				{
       					if(TotEffmonth(n)(indmonth(i))>0.0)
       					{
       						pcat(n)++;
       					}				
       				}
       			}

       			tot_pcat=sum(pcat);
       		
	END_CALCS

	init_matrix obsCatchNatAge(1,tot_pcat,sage-3,nage);

	init_number nItNobs;
	init_matrix obsIt(1,nItNobs,1,5);

	init_int eof;
	
	
	LOC_CALCS
		
		if( eof != 999 )
		{
			
			cout<< "effPwr is: "<<effPwr<<endl;
			cout<< "nItNobs is: "<<nItNobs<<endl;
			cout<< "Error reading data.\n Fix it."<<endl;
			cout<< "eof is: "<<eof<<endl;
			ad_exit(1);
		}

	END_CALCS

	//vector tau_c(1,nations);
	number tau_c;

	!! ad_comm::change_datafile_name("lagrangian_SA.ctl");
	
	// |---------------------------------------------------------------------------------|
	// | Control File - parameter intitial values and prior info
	// |---------------------------------------------------------------------------------|
	// | ilvec[1] 		-> fisheries catge data (1,fisharea)
	// | ilvec[2]       -> number of survey years       (1,nItNobs)
		

	init_int npar;
	init_matrix theta_control(1,npar,1,7);

	vector   theta_ival(1,npar);
	vector     theta_lb(1,npar);
	vector     theta_ub(1,npar);
	ivector   theta_phz(1,npar);
	ivector theta_prior(1,npar);
	
	LOC_CALCS
		
		theta_ival  = column(theta_control,1);
		theta_lb    = column(theta_control,2);
		theta_ub    = column(theta_control,3);
		theta_phz   = ivector(column(theta_control,4));
		theta_prior = ivector(column(theta_control,5));
		
	END_CALCS

	init_vector wt_ival(syr,nyr);
	


	// |---------------------------------------------------------------------------------|
	// | VECTOR DIMENSIONS FOR NEGATIVE LOG LIKELIHOODS
	// |---------------------------------------------------------------------------------|
	// | ilvec[1] 		-> fisheries catge data (1,fisharea)
	// | ilvec[2]       -> number of survey years       (1,nItNobs)
	
	
	ivector ilvec(1,3);

	!! ilvec(1) = fisharea;
	!! ilvec(2) = nItNobs;
	!! ilvec(3) = 1;

INITIALIZATION_SECTION
 	 theta theta_ival;
 	 wt wt_ival;

PARAMETER_SECTION

	//=================================
	//estimable parameters
	//=================================

	init_bounded_vector_vector theta(1,npar,1,1,theta_lb,theta_ub,theta_phz);

	//init_bounded_number log_mo(theta_control(1,2),theta_control(1,3),theta_control(1,4));
	
	//init_bounded_number log_cvPos(theta_control(2,2),theta_control(2,3),theta_control(2,4));
	//init_bounded_number log_maxPos50(theta_control(3,2),theta_control(3,3),theta_control(3,4));
	//init_bounded_number log_maxPossd(theta_control(4,2),theta_control(4,3),theta_control(4,4));




	//=================================
	//Population dynamics estimated parameters
	//=================================

	//init_bounded_number log_Ro(theta_control(5,2),theta_control(5,3),theta_control(5,4));
	//init_bounded_number h(theta_control(6,2),theta_control(6,3),theta_control(6,4),theta_control(6,5));

	//init_number log_mo;
	//init_number log_cvPos7
	//init_number log_maxPos50;
	//init_number log_maxPossd;

	init_bounded_number_vector wt(syr,nyr,-3.0,3.0,4);

	objective_function_value f;

	//=================================
	//derived quantities
	//=================================

	number kappa;
	
	number phiE;
	number So;
	number Bo;
	number Ro;
	number beta;
	number h;
	number log_avgrec;
	number sigma_r;

	//number Rbar;

	number m_tsp;
	
	number mo;
	number maxPos50;
	number maxPossd;
	number cvPos;

	number surv_q;

	number spr_opt;
	number fspr;	

	//vector wt(syr,nyr);

	vector lxo(sage,nage);
	vector za(sage,nage);
	vector SB(1,ntstp);
	vector tB(1,ntstp);
	vector ytB(syr-20,nyr);
	vector phie(syr,nyr);
	vector spr(syr,nyr);

	vector maxPos(sage,nage);
	vector varPos(sage,nage);

	vector SSB(syr,nyr);
	vector log_rt(syr,nyr);
	vector delta(syr,nyr);
	vector catlim(1,nations);
	
	vector it_hat(1,nItNobs);
	vector epsilon(1,nItNobs);
	//vector maxPos(sage,nage);
	//vector varPos(sage,nage);
	
	
	vector npvec(1,npar);

	matrix nlvec(1,3,1,ilvec);

	matrix NationVulB(1,ntstp,1,nations);
	matrix Nage(1,ntstp,sage,nage);
 	matrix VulB(1,ntstp,sage,nage);
 	matrix PosX(1,ntstp,sage,nage);
 	matrix Effage(1,ntstp,sage,nage);
 	matrix VBarea(1,ntstp,sarea,narea);
 	matrix Fatage(1,ntstp,sage,nage);
 	matrix Catage(1,ntstp,sage,nage);

 	matrix yFatage(syr-20,nyr,sage,nage);
 	matrix yNage(syr-20,nyr,sage,nage);
 	matrix seltotal(syr-20,nyr,sage,nage);
 	matrix yCatchtotalage(syr-20,nyr,sage,nage);


 	
 	//matrix propVBarea(1,ntstp,sarea,narea);
 	matrix Effarea(1,ntstp,sarea,narea);
 	matrix TotEffyear(1,fisharea,syr,nyr);
		
	
 	3darray NAreaAge(1,ntstp,sarea,narea,sage,nage);	
 	3darray CatchAreaAge(itsp,ntstp,sarea,narea,sage,nage);
 	3darray CatchNatAge(itsp,ntstp,1,fisharea,sage-2,nage);
 	3darray EffNatAge(1,fisharea,itsp,ntstp,sage-2,nage);
 	3darray selfisharea(syr,nyr,1,fisharea,sage-2,nage);
 	3darray selnation(syr,nyr,1,nations,sage-2,nage) 

 	3darray yCatchNatAge(syr,nyr,1,fisharea,sage,nage);
 	3darray yCatchStateAge(syr,nyr,1,nations,sage,nage);
 	
 	matrix predCatchNatAge(1,tot_pcat,sage-3,nage);

PRELIMINARY_CALCS_SECTION
		

PROCEDURE_SECTION

	incidence_functions();


	initialization();
	burn_in();
	
	move_grow_die();
	clean_catage();
	calc_selectivity();
	
	calcStockRecruitment();
	calc_obj_func();

FUNCTION dvar_vector cnorm(const double& x, const dvar_vector& mu, const dvar_vector& sd)

	dvar_vector rst(sage,nage);
	dvar_vector stx(sage,nage);


	for(int a= sage; a<= nage; a++)
	{
		stx(a) = (x-mu( a ))/sd( a );
		rst(a) = cumd_norm(stx( a ));
	}

	return(rst);

FUNCTION void calc_numbers_at_age(const int& ii, const dvariable& expwt )

		dvar_vector propBarea(sarea,narea);
		
		switch (indmonth(ii)) {
            case 1:           	
            	
            	//g_rt(indyr(ii)) = log_avgrec+ expwt;
            	Nage(ii)(sage) = mfexp(log_avgrec+ expwt);
            	//Nage(ii)(sage) = (So*SB(ii-nmon)/(1.+beta*SB(ii-nmon)))*(mfexp(expwt));

            	for(int a = sage+1;a<=nage;a++)
            	{
            		
					propBarea.initialize();
			
					for(int rr =sarea; rr<=narea; rr++)
					{		
						propBarea(rr) = (cnorm(areas(rr)+0.5,PosX(ii-1),varPos)-cnorm(areas(rr)-0.5,PosX(ii-1),varPos))(a-sage);	
					}
            	
					Nage(ii)(a) = (Nage(ii-1)(a-1)*propBarea)*mfexp(-(m_tsp+q*Effarea(ii-1)*va(a-1))) +
								  Nage(ii-1)(a-1)*(1-sum(propBarea))*mfexp(-(m_tsp));

            		//Nage(ii)(a) = Nage(ii-1)(a-1)*mfexp(-(m_tsp+q*Effage(ii-1)(a-1)*va(a-1)));
            	}
      
            	Nage(ii)(nage) = sum(elem_div(elem_prod((Nage(ii-1)(nage-1)*propBarea),mfexp(-(m_tsp+q*Effarea(ii-1)*va(nage-1)))),
            					(1.-mfexp(-(m_tsp+q*Effarea(ii-1)*va(nage))))))+
            					(Nage(ii-1)(nage-1)*(1.0-sum(propBarea))*mfexp(-m_tsp))/(1.-mfexp(-m_tsp));

            	yNage(indyr(ii))(sage,nage) = Nage(ii)(sage,nage);
            	ytB(indyr(ii)) = Nage(ii)*wa;
            	

            	break;
            	
            default: 
            
            	for(int a = sage;a<=nage;a++)
            	{		
					propBarea.initialize();
			
					for(int rr =sarea; rr<=narea; rr++)
					{		
						propBarea(rr) = (cnorm(areas(rr)+0.5,PosX(ii-1),varPos)-cnorm(areas(rr)-0.5,PosX(ii-1),varPos))(a-sage+1);	
					}

            		Nage(ii)(a) = Nage(ii-1)(a)*propBarea*mfexp(-(m_tsp+q*Effarea(ii-1)*va(a)))+
            					 Nage(ii-1)(a)*(1.0-sum(propBarea))*mfexp(-(m_tsp));
            	}


            	break;
        }
		
		VulB(ii) = elem_prod(elem_prod(Nage(ii),va),wa);
		SB(ii) = elem_prod(Nage(ii),fa)*wa/2.0;

		//cout<<"OK after calc_numbers_at_age"<<endl;		

FUNCTION void calc_effarea(const int& ii,const int& ie)

		dvar_vector tmp1(sarea,narea);
		dvar_vector tmp2(sarea,narea);
		tmp1.initialize();
		tmp2.initialize();

		

		for(int n=1; n<=nations;n++){
       		NationVulB(ii,n) = sum(pow(VBarea(ii)(ntmp1(n),ntmp1(n+1)-1.0)+0.00001,fbeta));
       	}



		for(int rr= sarea; rr<=narea; rr++)
		{
			
			//tmp1(rr)= pow((VBarea(ii)(rr)/(NationVulB(ii)(indnatarea(rr))+0.01))+0.0000001,fbeta) * effPwr(rr);	
			tmp1(rr)= (pow(VBarea(ii)(rr)+0.00001,fbeta)/(NationVulB(ii)(indnatarea(rr))+0.01)) * effPwr(rr);	
			tmp2(rr) = tmp1(rr)*TotEffyear(indfisharea(rr))(indyr(ie));
			Effarea(ii)(rr) = tmp2(rr)*TotEffmonth(indfisharea(rr))(indmonth(ii));

		}

		//cout<<"OK after calc_effarea"<<endl;		

FUNCTION void calc_position(const int& ii)

	varPos = maxPos*cvPos;

	PosX(ii) = minPos + (maxPos - minPos) * (0.5+0.5*sin(indmonth(ii)*PI/6 - mo*PI/6)); 

	

		for(int r = sarea;r <= narea;r++)
		{
			VBarea(ii)(r) = VulB(ii)*(cnorm(areas(r)+0.5,PosX(ii),varPos)-cnorm(areas(r)-0.5,PosX(ii),varPos));
			NAreaAge(ii)(r) = elem_prod(Nage(ii)(sage,nage),(cnorm(areas(r)+0.5,PosX(ii),varPos)-cnorm(areas(r)-0.5,PosX(ii),varPos)));
		}

		//cout<<"OK after calc_position"<<endl;


FUNCTION void calc_catage(const int& ii)

		int a,r;
		for(r = sarea;r <= narea;r++)
		{
			CatchNatAge(ii)(indfisharea(r))(sage-2) = ii;
			CatchNatAge(ii)(indfisharea(r))(sage-1) = indfisharea(r);

			for(a = sage; a<=nage;a++)
			{
				CatchAreaAge(ii)(r)(a) = q*Effarea(ii)(r)*va(a)/(q*Effarea(ii)(r)*va(a)+m_tsp)*(1-mfexp(-(q*Effarea(ii)(r)*va(a)+m_tsp)))*NAreaAge(ii)(r)(a);
				CatchNatAge(ii)(indfisharea(r))(a)+= CatchAreaAge(ii)(r)(a);
			}

			Catage(ii)(sage,nage) += CatchAreaAge(ii)(r)(sage,nage);
			yCatchNatAge(indyr(ii))(indfisharea(r))(sage,nage) += CatchAreaAge(ii)(r)(sage,nage);			
			yCatchStateAge(indyr(ii))(indnatarea(r))(sage,nage) += CatchAreaAge(ii)(r)(sage,nage);
			yCatchtotalage(indyr(ii))(sage,nage) += CatchAreaAge(ii)(r)(sage,nage);



		}

		//cout<<"OK after calc_catage"<<endl;



FUNCTION incidence_functions
	
	maxPos.initialize();

	lxo(sage)=1.;
	for(int a = sage+1; a<= nage; a++){
		lxo(a) = lxo(a-1)*mfexp(-m);
	}	
	lxo(nage) /= (1. - mfexp(-m)); 

	Ro = mfexp(theta(5)(1));
	h = theta(6)(1);
	log_avgrec = theta(7)(1);
	sigma_r = theta(8)(1);

	kappa 	= 4*h/(1-h);
	phiE	= elem_prod(lxo,fa)*wa;
	So 		= kappa/phiE;
	Bo 		= kappa/So*Ro;
	beta 	= (kappa-1)/Bo;

	m_tsp 	= m/nmon;
	za 		= m_tsp+va*fe;

	maxPos50 = mfexp(theta(3)(1));
	maxPossd = mfexp(theta(4)(1));
	cvPos 	 = mfexp(theta(2)(1));
	mo 	= mfexp(theta(1)(1));

	for(int n=1;n<=fisharea;n++)
    {
    	TotEffyear(n)(syr,nyr) = Fmult* pTotEffyear(n)(syr,nyr);
    }

    //cout<<"OK after incidence_functions"<<endl;

	
	

FUNCTION initialization
	

	NAreaAge.initialize();
 	CatchAreaAge.initialize();
 	CatchNatAge.initialize();
 	Nage.initialize();



	Nage(1,1) = So*Bo/(1+beta*Bo);

	for(int i=sage+1 ; i <= nage ; i++)
	{
		Nage(1,i) = Nage(1,i-1) * mfexp(-za(i-1));
	}
	Nage(1)(nage) /= (1.-mfexp(-za(nage)));


	yNage(syr-20)(sage,nage)= Nage(1)(sage,nage);
	
	

	VulB(1) = elem_prod(elem_prod(Nage(1),va),wa);


	SB(1) = elem_prod(Nage(1),fa)*wa/2;


	maxPos.initialize();
	tB(1) = Nage(1)*wa;
	ytB(syr) = Nage(1)*wa;
	
	calcmaxpos();
	
	calc_position(1);
	calc_effarea(1,itsp);

	//cout<<"OK after initialization"<<endl;	
	

FUNCTION burn_in



	for(int i=2;i<=itsp-1;i++)
	{
		
		calc_numbers_at_age(i,0.0);	
		
		maxPos.initialize();	
		tB(i) = Nage(i)*wa;		
		calcmaxpos();

		calc_position(i);	
		calc_effarea(i,itsp);	
		
	}
 	//cout<<"OK after burn_in"<<endl;	


FUNCTION move_grow_die


	for(int i=itsp;i<=ntstp;i++)
	{
		
		calc_numbers_at_age(i,wt(indyr(i)));
		if(indmonth(i)==1){
			SSB(indyr(i))  = elem_prod(Nage(i),fa)*wa/2.0;	
		}
		
		maxPos.initialize();	
		tB(i) = Nage(i)*wa;		
		calcmaxpos();	
		calc_position(i);
		calc_effarea(i,i);
		
		calc_catage(i);

		//cout<< "i is "<<i << endl;
		
	} 

FUNCTION calcStockRecruitment

	delta.initialize();
	dvar_vector tmp_rt(syr,nyr);
	
	//for(int i=syr;i<=nyr;i++)
	//{
		tmp_rt = elem_div((So*SSB),(1.+beta*SSB));
	//}
	

	delta = (log_avgrec+wt)-log(tmp_rt)+0.5*sigma_r*sigma_r;

FUNCTION clean_catage

	int p;
       
	p=1;
	for(int i=itsp;i<=ntstp;i++)
	{
		for(int n=1;n<=fisharea;n++)
		{					
			if(TotEffmonth(n)(indmonth(i))>0)
       		{
       			predCatchNatAge(p)(sage-3) = i;
       			predCatchNatAge(p)(sage-2) = indmonth(i);
				predCatchNatAge(p)(sage-1) = n;
       			predCatchNatAge(p)(sage,nage) = CatchNatAge(i)(n)(sage,nage);
				 p++;
       		}	
		}
	}

	//cout<<"OK after clean_catage"<<endl;


FUNCTION calc_survey_ll

	epsilon.initialize();
	dvector surv_yrs(1,nItNobs);
	dvar_vector Best(1,nItNobs);
	dvar_vector survBobs(1,nItNobs);
	dvar_vector wt(1,nItNobs);
	dvar_vector zt(1,nItNobs);
	dvar_vector (1,nItNobs);
	dvariable zbar;

	

	for(int i=1;i<=nItNobs;i++)
	{	
		surv_yrs(i) =  obsIt(i)(1);
		int ind_sv;
		ind_sv = surv_yrs(i)*(nmon-smon+1)-(nmon-obsIt(i)(5));

		if(indmonth(ind_sv)==obsIt(i)(5))
       	{	
       		Best(i) = sum(VulB(ind_sv)(sage,nage));
       		survBobs(i) = obsIt(i)(2);
       		wt(i) = obsIt(i)(2);
       	}
    }

    
    
    wt = 1.0/square(exp(wt));
    wt = wt/sum(wt);
    
    zt = log(survBobs);

	zt -= log(it_hat(1,nItNobs));
	zbar = zt*wt;
	surv_q = mfexp(zbar);
   
	it_hat(1,nItNobs) = surv_q * Best(1,nItNobs);
	zt  -= zbar;

	// Standardized observation error residuals.
	epsilon(1,nItNobs) = elem_div(zt,survBobs);




FUNCTION calc_obj_func

	//double tau_c;

	nlvec.initialize();
	npvec.initialize();
	
	
	
		for(int n = 1; n<=fisharea;n++)
		{
			dvar_matrix nu(1,pcat(n),sage,nage);
			dmatrix O(1,pcat(n),sage,nage);
			dvar_matrix P(1,pcat(n),sage,nage);
			
			O.initialize();
			P.initialize();
			nu.initialize();


			
			for(int i=1; i<=pcat(n);i++)
			{
				int ii;
				ii = sum(pcat(1,n))-pcat(n)+i;
				O(i) = (obsCatchNatAge(ii)(sage,nage))/sum(obsCatchNatAge(ii)(sage,nage));
				P(i) = ((predCatchNatAge(ii)(sage,nage))/(sum(predCatchNatAge(ii)(sage,nage))+0.01))+0.000001;
				

				
				
			}			

			nlvec(1)(n) =  dmvlogistic(O,P,nu,tau_c,dMinP);

		}
		
		for(int i=1;i<=nItNobs;i++){

			nlvec(2)(i) = dnorm(epsilon(i),0.0,1.0);
		}


	if( active(theta(5)) || active(theta(6)) )
	{
		
		nlvec(3) = dnorm(delta,sigma_r);
	}
		//cout<<"maxPos50 is "<<maxPos50<<endl;
		//cout<<"maxPossd is "<<maxPossd<<endl;
		//cout<<"cvPos is "<<cvPos<<endl;
		//cout<<"mo is "<<mo<<endl;
	//exit(1);

	//Priors 
	//prior for h
	for(int i=1;i<=npar;i++)
	{
		switch(theta_prior(i))
		{
				case 1:		//normal
					npvec(i) = dnorm(theta(i,1),theta_control(i,6),theta_control(i,7));
					break;
					
				case 2:		//lognormal CHANGED RF found an error in dlnorm prior. rev 116
					npvec(i) = dlnorm(theta(i,1),theta_control(i,6),theta_control(i,7));
					break;
					
				case 3:		//beta distribution (0-1 scale)
					double lb,ub;
					lb=theta_lb(i);
					ub=theta_ub(i);
					npvec(i) = dbeta((theta(i,1)-lb)/(ub-lb),theta_control(i,6),theta_control(i,7));
					break;
					
				case 4:		//gamma distribution
					npvec(i) = dgamma(theta(i,1),theta_control(i,6),theta_control(i,7));
					break;
					
				default:	//uniform density
					npvec(i) = log(1./(theta_control(i,3)-theta_control(i,2)));
					break;
			}
	}

	

	
	f=sum(nlvec)+sum(npvec);

	//cout<<"ok after calc_obj_func"<<endl;
	//f=sum(nlvec)/10000;
	


FUNCTION dvar_vector calcmaxpos()

	maxPos(sage,nage) = 1./(1.+mfexp(-(age-maxPos50)/maxPossd));
	maxPos(sage,nage) *= (narea-sarea);
	maxPos(sage,nage) += sarea;

			
	return(maxPos);


FUNCTION calc_spr
	/*** @brief Calculate SPR in the last year
	 * @details SPR=phi.e/phi.E 
	 * TODO read in SPR target
	 */

	 int i, ii, a;

	dvector lz(sage,nage);

	for(i=1; i<=ntstp;i++){
		Fatage(i)(sage,nage)=elem_div(Catage(i)(sage,nage),Nage(i)(sage,nage));
	}
	
	
	for(ii=syr; ii<=nyr;ii++){

		yFatage(ii)(sage,nage) = elem_div(yCatchtotalage(ii)(sage,nage),yNage(ii)(sage,nage));
			
		lz.initialize();
		lz(sage) = 1.;


		for(a=sage+1; a<=nage;a++){
			lz(a)= value(lz(a-1)*mfexp(-m-yFatage(ii)(a-1)));
		}
		lz(nage) /=  value(1.-mfexp(-m-yFatage(ii)(nage)));

		
		phie(ii)=elem_prod(lz,fa)*wa;
		spr(ii)=phie(ii)/phiE;


	}


FUNCTION calc_spr_optim	


	dvector fbars(1,4001);
	fbars.fill_seqadd(0.000,0.001);

	int  it, itt, a;

	int NF=size_count(fbars);

	dvector allspr(1,NF);
	dvector allphie(1,NF);
	dvector diffspr(1,NF);

	dvector     lz(sage,nage);
	dvector 	fage(sage,nage);
	dvector 	tmpsel(sage,nage);

	dvariable sol;


	for(it=1;it<=NF;it++)
	{
		tmpsel = value(seltotal(nyr)(sage,nage));
		fage = fbars(it)*tmpsel;

		lz.initialize();
				
		lz(sage) = 1.;

		for(a=sage+1; a<=nage;a++){			
			lz(a) = (lz(a-1)*mfexp(-m-fage(a-1)));
		}
		lz(nage) /=  (1.-mfexp(-m-fage(a-1)));

		allphie(it)=elem_prod(lz,fa)*wa;
		allspr(it)= value(allphie(it)/phiE);
		diffspr(it)= (allspr(it)-0.40)*(allspr(it)-0.40);


	}


	sol=min(diffspr);

	for(itt=1; itt<=NF; itt++)
	{
		if(sol==diffspr(itt)){
			spr_opt = allspr(itt);
			fspr = fbars(itt);
		} 
	}	
	

	ofstream ofs("spr.rep");


	ofs<<"diffspr" << endl << diffspr <<endl;
	ofs<<"allspr" << endl << allspr <<endl;
	ofs<<"allphie" << endl << allphie <<endl;
	ofs<<"phiE" << endl << phiE <<endl;
	
	




FUNCTION calc_selectivity

	int n, nn, i;

	for(i=syr; i<=nyr;i++){
		
		for(n=1;n<=nations;n++){
			
			selnation(i)(n)(sage-1)=i;
			selnation(i)(n)(sage-2)=n;
			selnation(i)(n)(sage,nage) = elem_div(yCatchStateAge(i)(n)(sage,nage),yNage(i)(sage,nage))/max(elem_div(yCatchStateAge(i)(n)(sage,nage),yNage(i)(sage,nage)));

		}
	
		for(nn=1;nn<=fisharea;nn++){
			
			selfisharea(i)(nn)(sage-1)=i;
			selfisharea(i)(nn)(sage-2)=nn;
			selfisharea(i)(nn)(sage,nage) = elem_div(yCatchNatAge(i)(nn)(sage,nage),yNage(i)(sage,nage))/max(elem_div(yCatchNatAge(i)(nn)(sage,nage),yNage(i)(sage,nage)));	

		}

		seltotal(i)(sage,nage) = elem_div(yCatchtotalage(i)(sage,nage),yNage(i)(sage,nage))/max(elem_div(yCatchtotalage(i)(sage,nage),yNage(i)(sage,nage)));
	}

FUNCTION TAC_input
	
	ofstream ofs("../OM/TAC_input.dat");
	ofs<<"#fspr" << endl << fspr <<endl;
	ofs<<"#seltotal" << endl << seltotal(nyr)(sage,nage) <<endl;	
	ofs<<"#yNage" << endl << yNage(nyr)(sage,nage) <<endl;	
	ofs<<"#Bo" << endl << Bo<<endl;
	ofs<<"#ytB" << endl << ytB(nyr) <<endl;
	



REPORT_SECTION

	calc_spr_optim();
	TAC_input();

	
	REPORT(mo);
	REPORT(maxPos50);
	REPORT(maxPossd);
	REPORT(cvPos);
	REPORT(h);
	REPORT(Ro);
	REPORT(syr);
	REPORT(nyr);
	REPORT(sage);
	REPORT(nage);
	REPORT(smon);
	REPORT(nmon);
	REPORT(sarea);
	REPORT(narea);
	REPORT(nations);
	REPORT(maxPos);
	REPORT(minPos);
	REPORT(varPos);
	REPORT(SB);
	REPORT(Nage);
	REPORT(VBarea);
	REPORT(Effage);
	REPORT(Effarea);
	REPORT(EffNatAge);
	REPORT(CatchNatAge);
	REPORT(CatchAreaAge);
	REPORT(indyr);
	REPORT(indmonth);
	REPORT(indnatarea);
	REPORT(fspr);
	REPORT(spr_opt);

	




TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
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

