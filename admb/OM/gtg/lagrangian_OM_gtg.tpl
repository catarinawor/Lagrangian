
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
		ifstream ifs( "../../seed.txt" ); // if this file is available
		ifs>>seed; //read in the seed
		seed += 10; // add 1 to the seed
		ofstream ofs( "../../seed.txt" ); //put out to seed.txt
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
	init_number tau_survey;		
	init_number err;

	//======================
	//length parameters -- fixed
	//======================
	init_vector Linf(1,ngroup);
	init_number vbk;
	init_number to;


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

	init_vector wa(sage,nage);
	init_vector fa(sage,nage);
	init_vector va(sage,nage);
	init_vector minPos(sage,nage);

						
	init_number Fmult;
	init_matrix pTotEffyear(1,fisharea,syr,nyr);
	init_matrix TotEffmonth(1,fisharea,smon,nmon);

	init_vector effPwr(sarea,narea);

	init_vector wt(rep_yr+1,proj_yr);

	init_int satype;

	init_vector nationTACprop(1,nations);

	init_int eof;
	
	

	LOC_CALCS
		
		if( eof != 999 )
		{
			cout<<"Linf "<<Linf<<endl;
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
	int ststp;
	int tmon;

   	vector age(sage,nage);
   	vector areas(sarea,narea);
   	ivector fishingr(1,fisharea);
   	ivector nationareas(1,nations);
   	
   	vector epsilon(1,surv_nobs);

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
		ststp =	tmon * (nyr-syr+1);
		ntstp = tmon * (proj_yr-syr+1);
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
			//wt.fill_randn(rng);
			//wt*=sigR;

			epsilon.fill_randn(rng);
			epsilon*=0.2;

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
		matrix TotEffyear(1,fisharea,syr,nyr);
		matrix TotEffyear_rep(1,fisharea,rep_yr+1,proj_yr);

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
       				TotEffyear(n)(syr,nyr) = Fmult* pTotEffyear(n)(syr,nyr);
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
       			pmat.initialize();
       			for(int n=1;n<=fisharea;n++)
       			{
       				for(int ii=rep_yr+1;ii<=nyr;ii++)
       				{
       					TotEffyear_rep(n)(ii)= TotEffyear(n)(ii);
       				}

       				for(int ia=nyr+1;ia<=proj_yr;ia++)
       				{
       					TotEffyear_rep(n)(ia)= TotEffyear(n)(nyr);
       				}

       				for(int i=rep_yr*nmon+1;i<=ntstp;i++)
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

			ofstream nfs("catlim.txt");
			nfs<<"catlim" << endl << "100000 100000" <<endl;
			
       			
	END_CALCS

	


PARAMETER_SECTION

	objective_function_value no_f;

	//derived quantities
	number kappa;
	number So;
	number Bo
	number beta;
	number phiE;

	number m_tsp;
	
	vector lxo(sage,nage);
	vector za(sage,nage);
	vector SB(1,ntstp);
	vector tB(1,ntstp);
	vector survB(1,surv_nobs);
	
	vector maxPos(sage,nage);
	vector varPos(sage,nage);
	vector varPosg(sage,nage);
	vector prop_ng(1,ngroup);

	vector ytB(rep_yr+1,proj_yr);
	vector totcatch(rep_yr+1,proj_yr);
	vector Ut(rep_yr+1,proj_yr);
	vector spr(syr,proj_yr);
	vector phie(syr,proj_yr);

	vector catlim(1,nations);


	matrix Catage(1,ntstp,sage,nage);
	matrix Fatage(1,ntstp,sage,nage);
 	matrix meanPosX(1,ntstp,sage-1,nage);
 	matrix tVBarea(1,ntstp,sarea,narea);
	matrix totVBnation(1,ntstp,1,nations);
	matrix VulBage(1,ntstp,sage,nage);
	matrix Effarea(1,ntstp,sarea,narea);
 	matrix tnage(1,ntstp,sage,nage);

 	matrix tot_comm_obsCatage(rep_yr+1,proj_yr,sage,nage);
 	matrix comm_obsCatage(rep_yr+1,proj_yr,sage,nage);
 	matrix surv_obsCatage(1,surv_nobs,sage,nage);
 	matrix yNage(syr,proj_yr,sage,nage);
 	matrix yFatage(syr,proj_yr,sage,nage);
 	matrix seltotal(syr,proj_yr,sage,nage);
 	matrix yCatchtotalage(syr,proj_yr,sage,nage);

 	3darray selfisharea(syr,proj_yr,1,fisharea,sage-2,nage);
 	3darray selnation(syr,proj_yr,1,nations,sage-2,nage); 	
 		
	3darray Nage(1,ngroup,1,ntstp,sage-2,nage); 
 	3darray VulB(1,ngroup,1,ntstp,sage-2,nage);
	
 	3darray VBarea(1,ngroup,1,ntstp,sarea,narea);
 	3darray totB(1,ngroup,1,ntstp,sage,nage);
 	
 	3darray PosX(1,ngroup,1,ntstp,sage-1,nage);
 	3darray NAreaAge(1,ntstp,sarea,narea,sage,nage);
 	3darray CatchAreaAge(1,ntstp,sarea,narea,sage,nage);
 	3darray CatchNatAge(1,ntstp,1,fisharea,sage-2,nage);
 	//3darray EffNatAge(1,fisharea,1,ntstp,sage-2,nage);
 	3darray yCatchNatAge(syr,proj_yr,1,fisharea,sage-2,nage);
 	3darray yCatchStateAge(syr,proj_yr,1,nations,sage-2,nage);
 	

 	4darray propVBarea(1,ntstp,1,ngroup,sarea,narea,sage-3,nage);

 	matrix obsCatchNatAge(1,tot_pcat,sage-3,nage);

PRELIMINARY_CALCS_SECTION

	read_catlim();
	incidence_functions();
	initialization();
	move_grow_die();
	
	run_projections();
	
	output_true();

	save_OMrep();
	
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

	maxPos.initialize();
	calcmaxpos();
	varPos=maxPos*cvPos;
	varPosg=sqrt((varPos*varPos)/(ngroup*ngroup*4));

	//cout<<"Ok after incidence_functions"<<endl;  

FUNCTION void calc_numbers_at_age(const int& ii, const dvariable& expwt)

	dvar_vector propBarea(sarea,narea);
	
	for(int g=1;g<=ngroup;g++)
	{
			Nage(g)(ii)(sage-2)=ii;
			Nage(g)(ii)(sage-1)=g;	
			
			
			switch (indmonth(ii)){
            	
            	case 1:  
            		
            		Nage(g)(ii)(sage) =(So*SB(ii-nmon)/(1.+beta*SB(ii-nmon)))*mfexp(expwt*err)* prop_ng(g);//*1.0/dg;

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
            	
            		
            		yNage(indyr(ii))(sage,nage) += Nage(g)(ii)(sage,nage);
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
			VulBage(ii)(sage,nage) += VulB(g)(ii)(sage,nage);

			tnage(ii)(sage,nage) += Nage(g)(ii)(sage,nage);
			tB(ii) += Nage(g)(ii)(sage,nage)*wa;
			totB(g)(ii)(sage,nage) = elem_prod(Nage(g)(ii)(sage,nage),wa);
			


	}
	
	SB(ii) = elem_prod(tnage(ii),fa)*wa/2.0;


	//cout<<"Ok after calc_numbers_at_age"<<endl;
	
FUNCTION void calc_effarea(const int& ii,const int& ia,const dvar_vector& ctlim)

	dvar_vector tmp1(sarea,narea);
	dvar_vector tmp2(sarea,narea);
	tmp1.initialize();
	tmp2.initialize();

	for(int n=1; n<=nations;n++){
       totVBnation(ii,n) = sum(pow(tVBarea(ii)(ntmp1(n),ntmp1(n+1)-1.0)+0.00001,fbeta));		 
	}

	for(int r= sarea; r<=narea; r++)
	{
		if(sum(yCatchNatAge(indyr(ii))(indnatarea(r))(sage,nage))<ctlim(indnatarea(r))){
			tmp1(r)= (pow(tVBarea(ii)(r)+0.0000001,fbeta)/(totVBnation(ii)(indnatarea(r))+0.01)) * effPwr(r);
			tmp2(r) = tmp1(r)*TotEffyear(indfisharea(r))(indyr(ia));
			Effarea(ii)(r) = tmp2(r)*TotEffmonth(indfisharea(r))(indmonth(ii));
		}else{
			Effarea(ii)(r) = 0.0;
		}
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

			Catage(ii)(sage,nage) += CatchAreaAge(ii)(r)(sage,nage);
			
			yCatchNatAge(indyr(ii))(indfisharea(r))(sage-2) = indyr(ii);
			yCatchNatAge(indyr(ii))(indfisharea(r))(sage-1) = indfisharea(r);
			yCatchNatAge(indyr(ii))(indfisharea(r))(sage,nage) += CatchAreaAge(ii)(r)(sage,nage);			
			
			yCatchStateAge(indyr(ii))(indnatarea(r))(sage-2) = indyr(ii);
			yCatchStateAge(indyr(ii))(indnatarea(r))(sage-1) = indnatarea(r);
			yCatchStateAge(indyr(ii))(indnatarea(r))(sage,nage) += CatchAreaAge(ii)(r)(sage,nage);
			
			yCatchtotalage(indyr(ii))(sage,nage) += CatchAreaAge(ii)(r)(sage,nage);

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
		yNage(1)(sage,nage) += Nage(g)(1)(sage,nage);
	}

	
	SB(1) = elem_prod(tnage(1),fa)*wa/2.0;

	
	calc_position(1);
	
	dvar_vector catl(1,nations);
	catl.fill("{100000000,100000000}");
	calc_effarea(1,1,catl);
	
	calc_catage(1);

	//cout<<"Ok after initialization"<<endl;
 

FUNCTION move_grow_die

		//if varps and mu are time varying need to move these to inside the loop
		calc_InitPos_gtg();
		
		maxPos.initialize();		
		calcmaxpos();
		
		varPos=maxPos*cvPos;
		varPosg=sqrt((varPos*varPos)/(ngroup*ngroup*4));

	dvar_vector catl(1,nations);
	catl.fill("{100000000,100000000}");


	for(int ie=2;ie<=(rep_yr*tmon);ie++)
	{
		calc_numbers_at_age(ie, 0.0);		
			
		calc_position(ie);
		calc_effarea(ie,ie,catl);
		calc_catage(ie);
		//cout<<"i is "<<i<<endl;

	}

	int p;
 	p=1;

 	int svyr;
 	svyr = 1;

	for(int i=rep_yr*tmon+1;i<=ststp;i++){
		
		calc_numbers_at_age(i,wt(indyr(i)));
		

		//cout<<"maxPos "<<maxPos<<endl;
		calc_position(i);
	
		calc_effarea(i,i,catl);		

		calc_catage(i);

		clean_catage(i);

		if(indmonth(i)==nmon){
			//survey calculations

			if(surv_yrs(svyr)==indyr(i)){
				
				survey_data(svyr);
				svyr +=1;
			}

			catage_comm(indyr(i));	
		}
	}


	//cout<<"Ok after move_grow_die"<<endl;

FUNCTION void clean_catage(const int& ii)
	
	ivector tmppcat(1,fisharea);
	for(int n=1;n<=fisharea;n++)
	{	
		tmppcat(n)= sum(pmat(n)(1,ii-1));
	}
		
	int p;
			
	p =sum(tmppcat)+1;

	//for(int i=rep_yr*nmon+1;i<=ntstp;i++)
	//{
		for(int n=1;n<=fisharea;n++)
		{
							
			if((TotEffmonth(n)(indmonth(ii)))>0)
       		{
       			
       			dvector pa(nage,sage);
       			pa.initialize();

       			obsCatchNatAge(p)(sage-3) = ii;
       			obsCatchNatAge(p)(sage-2) = indmonth(ii);
				obsCatchNatAge(p)(sage-1) = n;
				
				pa = value((CatchNatAge(ii)(n)(sage,nage))/sum(CatchNatAge(ii)(n)(sage,nage)+0.01))+0.000001;
				obsCatchNatAge(p)(sage,nage) = rmvlogistic(pa,tau_c,seed+ii);
				p++;	
			}	
		}
	

	//cout<<"Ok after clean_catage"<<endl;
	
FUNCTION dvar_vector calcmaxpos()
	
		 					
		maxPos(sage,nage) = 1./(1.+mfexp(-(age-maxPos50)/maxPossd));
		maxPos(sage,nage) *= (narea-sarea);
		maxPos(sage,nage) += sarea;		
	
		return(maxPos);
		
	//cout<<"Ok after calcmaxpos"<<endl;



FUNCTION void catage_comm(const int& ii)

	
	//for(int i=rep_yr*nmon+1;i<=ntstp;i++)
	//{
	for(int i=ii*tmon-tmon+1;i<=ii*tmon;i++)
	{
		
		for(int n=1;n<=fisharea;n++)
		{
							
			if(TotEffmonth(n)(indmonth(i))>0)
       		{

       			tot_comm_obsCatage(indyr(i))(sage,nage) += CatchNatAge(i)(n)(sage,nage);
       		}
       		
		}
	}
			dvector pa(sage,nage);
       		pa.initialize();
      
    
		pa = value((tot_comm_obsCatage(ii)(sage,nage))/(sum(tot_comm_obsCatage(ii)(sage,nage))+0.01))+0.000001;			
		//pa = value((tot_comm_obsCatage(ii)(sage,nage))/(sum(tot_comm_obsCatage(ii)(sage,nage))));			
		
		comm_obsCatage(ii)(sage,nage) = rmvlogistic(pa,tau_c,seed+ii);	
		totcatch(ii) = 	sum(tot_comm_obsCatage(ii)(sage,nage));	
		Ut(ii) = totcatch(ii)/tB(ii*tmon-(nmon-smon));
		ytB(ii) = tB(ii*tmon-(nmon-smon));


FUNCTION  void survey_data(const int& ii)
		
		int ind_sv;
		ind_sv = surv_yrs(ii)*(tmon)-(nmon-surv_mon);
		
		//if(indmonth(ind_sv)==surv_mon)
       	//{	
       		survB(ii)=sum(VulBage(ind_sv)(sage,nage))* mfexp(epsilon(ii));
       		

       		dvector pp(sage,nage);
       		pp.initialize();
       		for(int n=1;n<=fisharea;n++)
			{    
				pp += value(CatchNatAge(ind_sv)(n)(sage,nage)); 			
				//pp += value((CatchNatAge(ind_sv)(n)(sage,nage))/(sum(CatchNatAge(ind_sv)(n)(sage,nage))+0.01))+0.000001;
			}	
			dvector ppp(nage,sage);
       		ppp.initialize();
       		ppp = ((pp)/(sum(pp)+0.01))+0.000001;
			surv_obsCatage(ii)(sage,nage) = rmvlogistic(ppp,tau_survey,seed+ii);
       	//}	
	
			
    cout<<"Ok after survey_data"<< endl;

FUNCTION calc_length_comps


    

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

		seltotal(i)(sage,nage) = elem_div(yCatchtotalage(i)(sage,nage),yNage(i)(sage,nage));
	}


FUNCTION calc_spr	

	int i, ii, a;

	
	dvector lz(sage,nage);
	
	
	for(ii=syr; ii<=proj_yr;ii++){


		yFatage(ii)(sage,nage) = elem_div(yCatchtotalage(ii)(sage,nage),yNage(ii)(sage,nage));

		lz.initialize();
		lz(sage) = 1.;
		for(a=sage+1; a<=nage;a++){
			lz(a)= value(lz(a-1)*mfexp(-m-yFatage(ii)(a-1)));
		}
		lz(nage) /= value(1.-mfexp(-m-yFatage(ii)(nage)));


		phie(ii)=elem_prod(lz,fa)*wa;
		spr(ii)=phie(ii)/phiE;

	
	}

	
	
	cout<<"Ok after calculating SPR"<<endl;


FUNCTION void run_stock_assessment(const int& ii,const int& svy)
		
	cout<<"running stock assessment"<<endl;

	switch (satype) {
            case 1:  
            //to be deprecated
            //output_pin_SA();
            output_datSA(ii,svy);
            //output_ctlSA(ii);

            #if defined __APPLE__ || defined __linux
            // cout<<m_est_fmsy<<endl;
            
            system("cd ../../stock_assessment/ && ./lagrangian_SA");

            #endif

            #if defined _WIN32 || defined _WIN64

            system("lagrangian_SS.exe");

            #endif

            calc_catlim();
			read_catlim();

			system("cd ../../../R/read_mse && make");

            break;

            case 2:

            write_iscam_data_file(ii,svy);
            write_iscam_ctl_file(ii,svy);

            

            #if defined __APPLE__ || defined __linux
            // cout<<m_est_fmsy<<endl;
            system("cd /Users/catarinawor/Documents/iSCAM/examples/hakelag/DATA && make clean && make");// && "exec make");
           
            #endif

            calc_catlim();
			read_catlim();
			

			system("cd ../../../R/read_mse && make readRdattwo");


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

	catlim = TAC * nationTACprop;
	cout<< "catlim" << endl << catlim<< endl;
	//cout<< "seltotal_SA" << endl << seltotal_SA<< endl;

	ofstream afs("catlim.txt");
	afs<<"catlim" << endl << catlim <<endl;





FUNCTION run_projections

	int svyr;
	svyr=1;

	
	for(int ib=syr;ib<=nyr;ib++){
		if(ib==surv_yrs(svyr)) svyr++;
	}

	catage_comm(nyr);
	run_stock_assessment(nyr,svyr-1);


	for(int ii=ststp+1;ii<=ntstp;ii++)
	{ 

		calc_numbers_at_age(ii,wt(indyr(ii)));
		//maxPos.initialize();		
		//calcmaxpos();

		//cout<<"maxPos "<<maxPos<<endl;
		calc_position(ii);
	
		calc_effarea(ii,ststp,catlim);

		calc_catage(ii);
		clean_catage(ii);

		

		if(indmonth(ii)==nmon){
			//survey calculations

			if(surv_yrs(svyr)==indyr(ii)){
				survey_data(svyr);
				svyr +=1;
			}
			
			catage_comm(indyr(ii));
			run_stock_assessment(indyr(ii),svyr-1);
		}
		

		
		
	}

	calc_spr();


FUNCTION save_OMrep

	system("cd ../../../R/read_mse && make readROM");
	// Terminal year of projection.
   

	
FUNCTION read_catlim

	cifstream ifs_clm("../TAC_input.dat");

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
	ofs<<"Bo " << endl << Bo <<endl;
	ofs<<"So " << endl << So <<endl;
	ofs<<"h " << endl << h <<endl;
	ofs<<"m " << endl << m <<endl;
	ofs<<"fisharea" << endl << fisharea <<endl;
	ofs<<"nations" << endl << nations <<endl;
	ofs<<"maxPos" << endl << maxPos <<endl;
	ofs<<"minPos" << endl << minPos <<endl;
	ofs<<"varPos" << endl << varPos <<endl;
	ofs<<"PosX" << endl << PosX <<endl;	
	ofs<<"SB" << endl << SB <<endl;
	ofs<<"tB" << endl << tB <<endl;
	ofs<<"Ut" << endl << Ut <<endl;
	ofs<<"ytB" << endl << ytB <<endl;
	ofs<<"survB" << endl << survB <<endl;
	ofs<<"VulB" << endl << VulB <<endl;
	ofs<<"Nage" << endl << Nage <<endl;
	ofs<<"VBarea" << endl << VBarea <<endl;
	ofs<<"Effarea"<< endl << Effarea <<endl;
	ofs<<"comm_obsCatage"<< endl << comm_obsCatage <<endl;
	ofs<<"surv_obsCatage"<< endl << surv_obsCatage <<endl;
	ofs<<"totVBnation" << endl << totVBnation <<endl;
	ofs<<"CatchNatAge"<< endl << CatchNatAge<<endl;
	ofs<<"CatchAreaAge"<< endl << CatchAreaAge<<endl;
	ofs<<"indyr"<< endl << indyr<<endl;
	ofs<<"indmonth"<< endl << indmonth<<endl;
	ofs<<"indnatarea"<< endl << indnatarea<<endl;
	ofs<<"propVBarea"<< endl << propVBarea <<endl;
	ofs<<"selfisharea"<< endl << selfisharea <<endl;
	ofs<<"selnation"<< endl << selnation <<endl;
	ofs<<"seltotal"<< endl << seltotal <<endl;
	ofs<<"yCatchtotalage"<< endl << yCatchtotalage <<endl;
	ofs<<"yCatchNatAge"<< endl << yCatchNatAge <<endl;
	ofs<<"yCatchStateAge"<< endl << yCatchStateAge <<endl;
	ofs<<"yNage"<< endl << yNage <<endl;
	ofs<<"Fatage"<< endl << Fatage <<endl;
	ofs<<"yFatage"<< endl << yFatage <<endl;
	ofs<<"spr"<< endl <<spr <<endl;
	ofs<<"phie"<< endl <<phie <<endl;
	ofs<<"phiE"<< endl <<phiE <<endl;



	
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



FUNCTION void output_datSA(const int& ii,const int& svy)

	int tmpts;
	tmpts = ii*tmon;

	int tot_tmppcat;
	ivector tmppcat(1,fisharea);

	for(int n=1;n<=fisharea;n++){
		 tmppcat(n)= sum(pmat(n)(1,tmpts));
	}
	 tot_tmppcat = sum(tmppcat);

	 
	ofstream mfs("../../stock_assessment/lagrangian_SA.dat");
	mfs<<"# syr " << endl << rep_yr+1 <<endl;
	mfs<<"# nyr " << endl << ii <<endl;
	mfs<<"# sage " << endl << sage <<endl;
	mfs<<"# nage " << endl << nage <<endl;
	mfs<<"# smon " << endl << smon <<endl;
	mfs<<"# nmon " << endl << nmon <<endl;
	mfs<<"# sarea " << endl << sarea <<endl;
	mfs<<"# narea " << endl << narea <<endl;
	mfs<<"# fisharea " << endl << fisharea <<endl;
	mfs<<"# fishbound " << endl << fishbound <<endl;
	mfs<<"# nations " << endl << nations <<endl;
	mfs<<"# border " << endl << border <<endl;
	mfs<<"# m " << endl << m <<endl;
	mfs<<"# fe " << endl << fe <<endl;
	mfs<<"# q " << endl << q <<endl;
	mfs<<"# fbeta " << endl << fbeta <<endl;
	mfs<<"# sigR " << endl << sigR <<endl;
	mfs<<"# weight at age " << endl << wa <<endl;
	mfs<<"# fecundity at age " << endl << fa <<endl;
	mfs<<"# vulnerability at age " << endl << va <<endl;
	mfs<<"# minPos "<< endl << minPos <<endl;
	mfs<<"# Fmult "<< endl << Fmult <<endl;
	mfs<<"# Total effort by country and year " << endl << trans(trans(TotEffyear_rep).sub(rep_yr+1,ii)) <<endl;
	mfs<<"# Total effort by country and month " << endl << TotEffmonth <<endl;
	mfs<<"# effPwr"<< endl << effPwr <<endl;
	mfs<<"# dMinP " << endl << 0.1e-15 <<endl;
	mfs<<"# tstp month area catage " << endl << (obsCatchNatAge).sub(1,tot_tmppcat)<<endl;	
	mfs<<"# nItNobs " << endl << svy <<endl;
	mfs<<"## iyr     index(it) gear area group sex log(SE) log(PE)   timing" << endl;
	for(int l=1;l<=svy;l++){
		//  	iyr     		index(it) 	log(SE)  	log(PE)   	timing
		mfs<<surv_yrs(l) <<"\t"<< survB(l) <<"\t" << 0.2<<"\t"<< 0.01<<"\t"<< surv_mon <<endl;
	}
	mfs<<"# eof " << endl << 999 <<endl;


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
    mfs<<"1" <<"\t"<<"4" <<"\t"<<"3.5" <<"\t"<<"0.45" <<"\t"<<"0" <<"\t"<<"7" <<"\t"<<"1" <<"\t"<<"2" <<"\t"<<"150.0" <<"\t"<<"1.0" <<"\t"<<"1.0" <<"\t"<<"71" <<"\t"<<"91"<<endl;
    mfs<<"1" <<"\t"<<"4" <<"\t"<<"3.5" <<"\t"<<"0.45" <<"\t"<<"0" <<"\t"<<"7" <<"\t"<<"1" <<"\t"<<"2" <<"\t"<<"150.0" <<"\t"<<"1.0" <<"\t"<<"1.0" <<"\t"<<"92" <<"\t"<<ii<<endl;
    mfs<<"2" <<"\t"<<"2" <<"\t"<<"1.5" <<"\t"<<"0.45" <<"\t"<<"0" <<"\t"<<"4" <<"\t"<<"5" <<"\t"<<"2" <<"\t"<<"200.0" <<"\t"<<"50.0"<<"\t"<<"50.0"<<"\t"<<"71" <<"\t"<<ii<<endl;
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
	mfs<<"71 "<< endl;
	mfs<<"71 "<< endl;
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
		ufs<<  i <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 0 <<"\t"<< 1 <<"\t"<<totcatch(i)<<"\t"<<0.1<<endl;
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
	ufs<< -99 <<endl;
	ufs<<"## year gear area stock sex |age columns (sage, nage) of weight at age data   ##" << endl;
	ufs<<"##1975 2     1     1    0  0.0550 0.1575 0.2987 0.3658 0.6143 0.6306 0.7873 0.8738 0.9678 0.9075 0.9700 1.6933 1.5000 1.9000 1.9555 2.7445 2.7445 2.7445 2.7445 2.7445 2.7445" << endl;
 	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<<"## MARKER FOR END OF DATA FILE (eof)                                         ##" << endl;
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<< 999 <<endl;                                
	

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

