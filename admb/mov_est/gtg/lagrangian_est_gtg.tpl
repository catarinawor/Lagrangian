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
	init_number fbeta;
	init_number sigR;

	init_vector wa(sage,nage);
	init_vector fa(sage,nage);
	init_vector va(sage,nage);
	init_vector minPos(sage,nage);
	
	//init_number Fmult;
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
	int tmon;
	int itsp;

   	vector age(sage,nage);
   	vector areas(sarea,narea);
   	ivector fishingr(1,fisharea);
   	ivector nationareas(1,nations);
   	
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

		tmon = nmon-smon+1;
		ntstp = tmon * (nyr-(syr-30)+1);
		itsp = tmon*30+1;

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

		imatrix pmat(1,fisharea,1,ntstp);
		//matrix TotEffyear(1,fisharea,syr,nyr);
		
		

		number delta;
		vector Xini(1,ngroup);

		

       LOC_CALCS

       				
       		int aa =0;
       			
       		for(int y=syr-30;y<=nyr;y++)		
       		{
       			for(int ii=smon;ii<=nmon;ii++)
       			{
       				aa++;
       				indyr(aa) = y;
       				indmonth(aa) = ii;
       			}
       		}

       		//for(int n=1;n<=fisharea;n++)
       		//{
       		//	TotEffyear(n)(syr,nyr) = Fmult* pTotEffyear(n)(syr,nyr);
       		//}      			


       		//calc for indfisharea
       		ntmp(1) = sarea;
       			
       		for(int n=1;n<=fisharea;n++)
       		{
       			ntmp(n+1)= ntmp(n)+fishingr(n);	
       			indfisharea(ntmp(n),ntmp(n+1)-1)=n;
       		}

       		ntmp(fisharea+1) = narea;
       		indfisharea(narea) = fisharea;

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
       			for(int i=itsp;i<=ntstp;i++)
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
			//cout<<"chegou aqui??"<<endl;
			//exit(1);
       		
	END_CALCS

	init_matrix obsCatchNatAge(1,tot_pcat,sage-3,nage);


	init_vector yYieldtotalobs(syr,nyr);

	init_int eof;
	
	
	LOC_CALCS
		
		if( eof != 999 )
		{
			cout<< "Error reading data.\n Fix it."<<endl;
			cout<< "minPos is: "<<minPos<<endl;
			cout<< "eof is: "<<eof<<endl;
			ad_exit(1);
		}

	END_CALCS

	//vector tau_c(1,nations);
	number tau_c;


PARAMETER_SECTION

	//=================================
	//estimable parameters
	//=================================

	//init_bounded_number mo(0.5,6,1);
	init_bounded_number log_mo(-0.6931472,1.25);
	init_bounded_number log_cvPos(-3,-0.7,1);
	init_bounded_number log_maxPos50(0.00,2.1,1);
	init_bounded_number log_maxPossd(-0.7,1.609438,1);
	init_bounded_number log_Fmult(-2.3,2.3,2);
	init_vector wt(syr,nyr,-1);

	objective_function_value f;

	//=================================
	//derived quantities
	//=================================

	
	number kappa;
	number phiE;
	number So;
	number Bo;
	number beta;

	number m_tsp;
	
	number mo;
	number maxPos50;
	number maxPossd;
	number cvPos;
	number Fmult;


	//vector wt(syr,nyr);

	vector lxo(sage,nage);
	vector za(sage,nage);
	vector SB(1,ntstp);
	//vector tB(1,ntstp);

	vector maxPos(sage,nage);
	vector varPos(sage,nage);
	vector varPosg(sage,nage);
	vector prop_ng(1,ngroup);
	vector yYieldtotal(syr,nyr);

	//matrix SBg(1,ntstp,1,ngroup);
	//matrix tBg(1,ntstp,1,ngroup);
	
	//matrix maxPos(1,ngroup,sage,nage);
	//matrix varPos(1,ngroup,sage,nage);

	//vector totcatch(syr,nyr);
	//vector totcatchLen(rep_yr+1,proj_yr);
	//vector Ut(syr,nyr);

	//matrix NationVulB(1,ntstp,1,nations);
	
	//matrix Catage(1,ntstp,sage,nage);
	//matrix Catlen(1,ntstp,1,nlen);
	
 	matrix tVBarea(1,ntstp,sarea,narea);
	matrix totVBnation(1,ntstp,1,nations);
	
	matrix Effarea(1,ntstp,sarea,narea);
	matrix TotEffyear(1,fisharea,syr,nyr);
 	

	
	
	
	//vector npvec(1,1);

	3darray Nage(1,ngroup,1,ntstp,sage-2,nage); 
 	3darray VulB(1,ngroup,1,ntstp,sage-2,nage);
 	//3darray CatageG(1,ngroup,itsp,ntstp,sage,nage);
	
 	3darray VBarea(1,ngroup,1,ntstp,sarea,narea);
 	//3darray totB(1,ngroup,1,ntstp,sage,nage);
	3darray PosX(1,ngroup,smon,nmon,sage-1,nage);
 	
 	
	
 	//3darray NAreaAge(1,ntstp,sarea,narea,sage,nage);	
 	//3darray CatchAreaAge(itsp,ntstp,sarea,narea,sage,nage);
 	3darray CatchNatAge(itsp,ntstp,1,fisharea,sage-2,nage);
 	//3darray EffNatAge(1,fisharea,itsp,ntstp,sage-2,nage);
 	
 	3darray NAreaAgeG(1,n_rg,1,ntstp,sage-3,nage);
 	3darray CatchAreaAgeG(1,n_rg,itsp,ntstp,sage-3,nage);
 	3darray propVBarea(1,n_rg,1,ntstp,sage-3,nage);

 	matrix predCatchNatAge(1,tot_pcat,sage-3,nage);

PRELIMINARY_CALCS_SECTION
		

PROCEDURE_SECTION

	

	initialization();
	incidence_functions();
	calc_InitPos_gtg();
	
	calc_first_year();
	burn_in();
	move_grow_die();


	
	calc_obj_func();

FUNCTION double cnorm(const double& x, const dvariable& mu, const dvariable& sd)

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

FUNCTION  calc_InitPos_gtg
	
	dvar_vector inimu(sage,nage);
	dvector Xtemp(1,ngroup);

	dvariable dif;



	maxPos.initialize();
	calcmaxpos();
	
	varPos=maxPos*cvPos;
	varPosg=sqrt((varPos*varPos)/(ngroup*ngroup*4.));

	
	
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
	

	//cout<<"Ok after calc_InitPos_gtg"<<endl;




	
	

FUNCTION initialization
			
	Nage.initialize();
 	VulB.initialize();
 	yYieldtotal.initialize();
 	//CatageG.initialize();
 	
 	VBarea.initialize();
 	PosX.initialize();


 	CatchNatAge.initialize();
 	NAreaAgeG.initialize();
 	CatchAreaAgeG.initialize();
 	propVBarea.initialize();
 	
 	predCatchNatAge.initialize();
 	//Catage.initialize();
 	tVBarea.initialize();
	
	totVBnation.initialize();
	
	Effarea.initialize();
 	

 	//cout<<"Ok after initilization"<<endl;
	
 	


FUNCTION incidence_functions
	
	lxo(sage)=1.;

	for(int a = sage+1; a<= nage; a++){
		lxo(a) = lxo(a-1)*mfexp(-m);
	}	
	lxo(nage) /= (1. - mfexp(-m)); 


	kappa 	= 4*h/(1-h);
	phiE	= elem_prod(lxo,fa)*wa;
	So 		= kappa/phiE;
	Bo 		= kappa/So*Ro;
	beta 	= (kappa-1)/Bo;


	m_tsp = m/nmon;
	za 	= m_tsp+va*fe;

	maxPos50 = mfexp(log_maxPos50);
	maxPossd = mfexp(log_maxPossd);
	cvPos 	 = mfexp(log_cvPos);
	mo 	= mfexp(log_mo);
	Fmult = mfexp(log_Fmult);
	
	for(int n=1;n<=fisharea;n++)
    {
    	TotEffyear(n)(syr,nyr) = Fmult* pTotEffyear(n)(syr,nyr);
    }  
	
	//cout<<"Ok after incidence_functions"<<endl;
	

	

FUNCTION void calc_numbers_at_age(const int& ii, const dvariable& expwt )
	
	dvar_vector propBarea(sarea,narea);
	dvar_vector tnage(sage,nage);

	tnage.initialize();
	
	
	for(int g=1;g<=ngroup;g++)
	{
			Nage(g)(ii)(sage-2)=ii;
			Nage(g)(ii)(sage-1)=g;	
			
			
			switch (indmonth(ii)){
            	
            	case 1:  
            		
            		Nage(g)(ii)(sage) =(So*SB(ii-nmon)/(1.+beta*SB(ii-nmon)))*mfexp(expwt)* prop_ng(g);//*1.0/dg;

            		for(int a = sage+1;a<nage;a++)
            		{
            		
						propBarea.initialize();
			
						for(int rr =sarea; rr<=narea; rr++)
						{		
							propBarea(rr) = cnorm(areas(rr)+0.5,PosX(g)(nmon)(a-1),varPosg(a-1))-cnorm(areas(rr)-0.5,PosX(g)(nmon)(a-1),varPosg(a-1));	
						}

						
						Nage(g)(ii)(a) = Nage(g)(ii-1)(a-1)*propBarea*mfexp(-(m_tsp+q*Effarea(ii-1)*va(a-1))) +
								  Nage(g)(ii-1)(a-1)*(1.0-sum(propBarea))*mfexp(-(m_tsp));
					}

				

					Nage(g)(ii)(nage) = sum(elem_div(elem_prod((Nage(g)(ii-1)(nage-1)*propBarea),mfexp(-(m_tsp+q*Effarea(ii-1)*va(nage-1)))),
            					(1.0-mfexp(-(m_tsp+q*Effarea(ii-1)*va(nage))))))+
            					(Nage(g)(ii-1)(nage-1)*(1.0-sum(propBarea))*mfexp(-m_tsp))/(1.-mfexp(-m_tsp));
            	
            		//exit(1);

            	break;
            	
            	default: 

            		for(int a = sage;a<=nage;a++)
            		{
            			propBarea.initialize();
			
						for(int rr =sarea; rr<=narea; rr++)
						{		
							propBarea(rr) = cnorm(areas(rr)+0.5,PosX(g)(indmonth(ii)-1)(a),varPosg(a))-cnorm(areas(rr)-0.5,PosX(g)(indmonth(ii)-1)(a),varPosg(a));	
						}

            			Nage(g)(ii)(a) = Nage(g)(ii-1)(a)*propBarea*mfexp(-(m_tsp+q*Effarea(ii-1)*va(a)))+
            						  	 Nage(g)(ii-1)(a)*(1.0-sum(propBarea))*mfexp(-(m_tsp));
            		}
			    
            	break;
        	}
            	
        	VulB(g)(ii)(sage-2) = ii;
        	VulB(g)(ii)(sage-1) = g;
			VulB(g)(ii)(sage,nage) = elem_prod(elem_prod(Nage(g)(ii)(sage,nage),va),wa);
		

			tnage(sage,nage) += Nage(g)(ii)(sage,nage);
			//tB(ii) += Nage(g)(ii)(sage,nage)*wa;
			//totB(g)(ii)(sage,nage) = elem_prod(Nage(g)(ii)(sage,nage),wa);
			


	}
	
	SB(ii) = elem_prod(tnage,fa)*wa/2.0;


	//cout<<"Ok after calc_numbers_at_age in year "<< indyr(ii)<<endl;
	//exit(1);


FUNCTION void calc_effarea(const int& ii,const int& ie)

	dvar_vector tmp1(sarea,narea);
	dvar_vector tmp2(sarea,narea);
	tmp1.initialize();
	tmp2.initialize();

	for(int n=1; n<=nations;n++){
       totVBnation(ii,n) = sum(pow(tVBarea(ii)(ntmp1(n),ntmp1(n+1)-1.0)+0.00001,fbeta));		 
	}

	for(int r= sarea; r<=narea; r++)
	{
			tmp1(r)= (pow(tVBarea(ii)(r)+0.00001,fbeta)/(totVBnation(ii)(indnatarea(r))+0.01)) * effPwr(r);
			tmp2(r) = tmp1(r)*TotEffyear(indfisharea(r))(indyr(ie));
			Effarea(ii)(r) = tmp2(r)*TotEffmonth(indfisharea(r))(indmonth(ii));
	}
	
	//cout<<"Ok after calc_effarea "<<endl;



FUNCTION void calc_position(const int& ii)


	//cout<<"aqui??"<<endl;
	dvar_vector meanPosX(sage,nage);
	meanPosX.initialize();

	meanPosX(sage,nage) = minPos + (maxPos - minPos) * (0.5+0.5*sin(indmonth(ii)*PI/6. - mo*PI/6. -PI/2.)); 


	int g, r, ig;

	for(ig=1;ig<=n_rg ;ig++)
	{
		
		g = n_group(ig);
		r = n_area(ig);

		NAreaAgeG(ig)(ii)(sage-3)= ii;
		NAreaAgeG(ig)(ii)(sage-2)= g;
		NAreaAgeG(ig)(ii)(sage-1)= r;

		
		
		
		VBarea(g)(ii)(r) = VulB(g)(ii)(sage,nage) * (cnorm2(areas(r)+0.5,PosX(g)(indmonth(ii)),varPosg)-cnorm2(areas(r)-0.5,PosX(g)(indmonth(ii)),varPosg));
		NAreaAgeG(ig)(ii)(sage,nage) = elem_prod(Nage(g)(ii)(sage,nage),(cnorm2(areas(r)+0.5,PosX(g)(indmonth(ii)),varPosg)-cnorm2(areas(r)-0.5,PosX(g)(indmonth(ii)),varPosg)));
	
		//NAreaAge(ii)(r) += elem_prod(Nage(g)(ii)(sage,nage),(cnorm(areas(r)+0.5,PosX(g)(ii),varPosg)-cnorm(areas(r)-0.5,PosX(g)(ii),varPosg)));
		tVBarea(ii)(r) += VBarea(g)(ii)(r);		


		propVBarea(ig)(ii)(sage-3) = ii;
		propVBarea(ig)(ii)(sage-2) = g;
		propVBarea(ig)(ii)(sage-1) = r;
		propVBarea(ig)(ii)(sage,nage) =  elem_prod(VulB(g)(ii)(sage,nage), (cnorm2(areas(r)+0.5,PosX(g)(indmonth(ii)),varPosg)-cnorm2(areas(r)-0.5,PosX(g)(indmonth(ii)),varPosg)));
			
		
	
	}



	//cout<<"Ok after calc_position"<<endl;



FUNCTION void calc_catage(const int& ii)

		//need to work on diminishing the dimensions here.

		int ig,r,g;
		
		

		

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

		CatchNatAge(ii)(indfisharea(r))(sage-2) = ii;
		CatchNatAge(ii)(indfisharea(r))(sage-1) = indfisharea(r);


		tmpc1 =elem_div(q*Effarea(ii)(r)*va,q*Effarea(ii)(r)*va+m_tsp);
		tmpc2 =elem_prod(tmpc1,(1-mfexp(-(q*Effarea(ii)(r)*va+m_tsp))));

		CatchAreaAgeG(ig)(ii)(sage,nage) = elem_prod(tmpc2,NAreaAgeG(ig)(ii)(sage,nage));		
		CatchNatAge(ii)(indfisharea(r))(sage,nage) += CatchAreaAgeG(ig)(ii)(sage,nage);

		yYieldtotal(indyr(ii)) += CatchNatAge(ii)(indfisharea(r))(sage,nage)*wa;
		
		//cout<<"CatchNatAge(ii)(indfisharea(r))(a)"<<CatchNatAge(ii)(indfisharea(r))(sage,nage)<<endl;
		
		//Catage(ii)(sage,nage) += CatchAreaAgeG(ig)(ii)(sage,nage);
		}

		
			
		//for(a = sage; a<=nage;a++)
		//{
		//	for(int r=sarea;r<=narea;r++)
		//	{		
		//		CatchAreaAge(ii)(r)(a) = q*Effarea(ii)(r)*va(a)/(q*Effarea(ii)(r)*va(a)+m_tsp)*(1-mfexp(-(q*Effarea(ii)(r)*va(a)+m_tsp)))*NAreaAge(ii)(r)(a);
		//		CatchNatAge(ii)(indfisharea(r))(a) += CatchAreaAge(ii)(r)(a);
		//
		//	}

		//}	
		//cout<<"Ok after calc_catage"<<endl;




	

FUNCTION calc_first_year
	
	
 	//calc_InitPos_gtg();

 	dvar_vector tnage(sage,nage);
 	tnage.initialize();

 	for(int g=1;g<=ngroup;g++)
	{
		Nage(g)(1)(sage-1)=g;
		Nage(g)(1)(sage-2)=1;
	
			
		Nage(g)(1)(sage,nage) =(So*Bo/(1.+beta*Bo))* prop_ng(g);


        for(int a=sage+1 ; a <= nage ; a++)
		{
			Nage(g)(1)(a) = Nage(g)(1)(a-1) * mfexp(-za(a-1));
		}
		Nage(g)(1)(nage) /= (1.-mfexp(-za(nage)));

            		          	
        VulB(g)(1)(sage-2) = 1;
        VulB(g)(1)(sage-1) = g;
		VulB(g)(1)(sage,nage) = elem_prod(elem_prod(Nage(g)(1)(sage,nage),va),wa);
	
		tnage(sage,nage) += Nage(g)(1)(sage,nage);
		//tB(1) += Nage(g)(1)(sage,nage)*wa;
		//totB(g)(1)(sage,nage) = elem_prod(Nage(g)(1)(sage,nage),wa);
	}

	


	
	SB(1) = elem_prod(tnage,fa)*wa/2.0;

	
	calc_position(1);
	
	calc_effarea(1,itsp);
	
	//cout<<"Ok after calc_first_year"<<endl;

FUNCTION burn_in

	

	for(int i=2;i<=itsp-1;i++)
	{
		//calc_InitPos_gtg();
		calc_numbers_at_age(i,0.0);	


		calc_position(i);	
		calc_effarea(i,itsp);	
	

		//cout<<"i is "<<i<<endl;
	}
	
 
 	


FUNCTION move_grow_die

	int p;
       
	p=1;

	

	for(int i=itsp;i<=ntstp;i++)
	{
		//calc_InitPos_gtg();
		calc_numbers_at_age(i,wt(indyr(i)));	
		//cout<< "i is "<<i << endl;
	

		calc_position(i);
		calc_effarea(i,i);
		
		calc_catage(i);
		//clean_catage(i);

		for(int n=1;n<=fisharea;n++)
		{					
			if(TotEffmonth(n)(indmonth(i))>0)
       		{
					clean_catage(i,p,n);
					p++;
       		}	
		}
		//cout<< "i is "<<i << endl;

		
		
	} 

	//cout<<"aqui??"<<endl;

	






FUNCTION void clean_catage(const int& ii,const int& pp,const int& nn)
	
	

       			predCatchNatAge(pp)(sage-3) = ii;
       			predCatchNatAge(pp)(sage-2) = indmonth(ii);
				predCatchNatAge(pp)(sage-1) = nn;
       			predCatchNatAge(pp)(sage,nage) = CatchNatAge(ii)(nn)(sage,nage);
				
	

	//cout<<"Ok after clean_catage"<<endl;


	//int p;
    //   
	//p=1;
	//for(int i=itsp;i<=ntstp;i++)
	//{
	//	for(int n=1;n<=fisharea;n++)
	//	{					
	//		if(TotEffmonth(n)(indmonth(i))>0)
    //   		{
    //   			predCatchNatAge(p)(sage-3) = i;
    //   			predCatchNatAge(p)(sage-2) = indmonth(i);
	//			predCatchNatAge(p)(sage-1) = n;
    //   			predCatchNatAge(p)(sage,nage) = CatchNatAge(i)(n)(sage,nage);
	//			 p++;
    //   		}	
	//	}
	//}



FUNCTION calc_obj_func


	dvar_vector nlvec(1,fisharea);

	nlvec.initialize();
	
	
		for(int n = 1; n<=fisharea;n++)
		{

			//cout<<"n"<<n<<endl;
			//cout<<"pcat(n)"<<pcat<<endl;

			dvar_matrix nu(1,pcat(n),sage,nage);
			dmatrix O(1,pcat(n),sage,nage);
			dvar_matrix P(1,pcat(n),sage,nage);
			
			O.initialize();
			P.initialize();
			nu.initialize();
			
			int ii;
			ii = sum(pcat(1,n))-pcat(n)+1;

		

			for(int i=1; i<=pcat(n);i++)
			{
			
				//cout<<"ii"<<ii<<endl;
					
				//O(i) = (obsCatchNatAge(ii)(sage,nage)+0.1e-30)/sum(obsCatchNatAge(ii)(sage,nage)+0.1e-5);
				//P(i) = (predCatchNatAge(ii)(sage,nage)+0.1e-30)/sum(predCatchNatAge(ii)(sage,nage)+0.1e-5);
				O(i) = (obsCatchNatAge(ii)(sage,nage))/sum(obsCatchNatAge(ii)(sage,nage));
				P(i) = (predCatchNatAge(ii)(sage,nage)/sum(predCatchNatAge(ii)(sage,nage)))+0.00000001;
				ii++;
				//cout<<"O("<<i<<")"<<O(i)<<endl;
				//cout<<"P("<<i<<")"<<P(i)<<endl;
			}			

			//cout<<"dmvlogistic(O,P,nu,tau_c(n),dMinP) is "<< dmvlogistic(O,P,nu,tau_c,dMinP)<<endl;
			//exit(1);									
			//nlvec(n) =  dmvlogistic(O,P,nu,tau_c(n),dMinP);
			//cout<<"nlvec"<<nlvec<<endl;
			//cout<<"nlvec"<<n<<dmvlogistic(O,P,nu,tau_c,dMinP)<<endl;



			nlvec(n) =  dmvlogistic(O,P,nu,tau_c,dMinP);
			//cout<<"chegou??"<<endl;

		}

		dvar_vector eta(syr,nyr);
		eta.initialize();
		dvar_vector nlcat(1,1);
		nlcat.initialize();
	
		for(int i = syr; i<= nyr; i++)
		{
			eta(i)=(log(yYieldtotalobs(i)) - log(yYieldtotal(i)));
		
		}
			
	if(last_phase()){
		nlcat(1) = norm2(eta);
	}else{
		nlcat(1) = 0.0;
	}
	//cout<<"Nage is"<<Nage<<endl;
	//output_true();
	//exit(1);

	cout<<"nlvec is"<<nlvec<<endl;
	//f=sum(nlvec)+sum(npvec);
	f=sum(nlvec)/100000+sum(nlcat);

	cout<<"nlcat is"<<nlcat<<endl;

	cout<<"f is"<<f<<endl;
	cout<<"maxPos50 is "<<maxPos50<<endl;
	cout<<"maxPossd is "<<maxPossd<<endl;
	cout<<"cvPos is "<<cvPos<<endl;
	cout<<"mo is "<<mo<<endl;
	cout<<"Fmult is "<<Fmult<<endl;
	
	//output_true();
	//exit(1);



FUNCTION dvar_vector calcmaxpos()
	
		maxPos(sage,nage) = 1./(1.+mfexp(-(age-maxPos50)/maxPossd));
		maxPos(sage,nage) *= (narea-minPos(sage));
		maxPos(sage,nage) += minPos(sage);		


		return(maxPos);

FUNCTION output_true
	
	ofstream ofs("firstrun.rep");

	ofs<<"mo" << endl << mo <<endl;
	ofs<<"tau_c" << endl << tau_c<<endl;
	ofs<<"maxPos50" << endl << maxPos50 <<endl;
	ofs<<"maxPossd" << endl << maxPossd <<endl;
	ofs<<"cvPos" << endl << cvPos <<endl;
	ofs<<"syr" << endl << syr <<endl;
	ofs<<"nyr" << endl << nyr <<endl;
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
	//ofs<<"PosX" << endl << PosX <<endl;	
	ofs<<"SB" << endl << SB <<endl;
	ofs<<"VulB" << endl << VulB <<endl;
	ofs<<"Nage" << endl << Nage <<endl;
	ofs<<"VBarea" << endl << VBarea <<endl;
	ofs<<"Effarea"<< endl << Effarea <<endl;
	ofs<<"CatchNatAge"<< endl << CatchNatAge<<endl;
	//ofs<<"CatchAreaAge"<< endl << CatchAreaAge<<endl;
	ofs<<"indyr"<< endl << indyr<<endl;
	ofs<<"indmonth"<< endl << indmonth<<endl;
	ofs<<"indnatarea"<< endl << indnatarea<<endl;
	//ofs<<"propVBarea"<< endl << propVBarea <<endl;
	ofs<<"propVBarea"<< endl << propVBarea <<endl;
	ofs<<"yYieldtotalobs"<< endl << yYieldtotalobs <<endl;
	ofs<<"yYieldtotal"<< endl << yYieldtotal <<endl;
	




REPORT_SECTION
	
	REPORT(mo);
	REPORT(maxPos50);
	REPORT(maxPossd);
	REPORT(cvPos);
	REPORT(Fmult);
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
	REPORT(Effarea);
	REPORT(CatchNatAge);
	//REPORT(CatchAreaAge);
	REPORT(indyr);
	REPORT(indmonth);
	REPORT(indnatarea);
	REPORT(yYieldtotalobs);
	REPORT(yYieldtotal);
	




TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 10000000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e10);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e10);
	gradient_structure::set_MAX_NVAR_OFFSET(50000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(50000);
	gradient_structure::set_MAX_DLINKS(40000);
 

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

