

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
	init_int nations;

	init_vector border(1,nations-1);

	init_int ngroup;
	init_vector prop_ng(1,ngroup);

	//======================
	//model parameters -- fixed
	//======================

	init_number Ro;
	init_number h;
	init_number m;
	init_number fe;
	init_number q;
	init_number sigR;
	

	init_vector wa(sage,nage);
	init_vector fa(sage,nage);
	init_vector va(sage,nage);
	init_vector minPos(sage,nage);
	

	init_matrix TotEffyear(1,nations,syr,nyr);
	init_matrix TotEffmonth(1,nations,smon,nmon);

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


		ntstp = (nmon-smon+1) * (nyr-(syr-20)+1);


			age.fill_seqadd(sage,1);
			areas.fill_seqadd(sarea,1);

			nationareas.initialize();
			dvector natmp1(1,nations);
			
			natmp1(1)=sarea;


			for(int n=1; n<=nations-1; n++)
			{
				natmp1(n+1)=border(n);

				for(int a=sarea;a<=narea;a++)
				{
					if(areas(a)>=natmp1(n)&areas(a)<border(n))
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
		ivector pcat(1,nations);
		vector natmp(1,nations+1);

		int itsp;
		int tot_pcat;
		


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

       			natmp(1) = sarea;

       			for(int n=1;n<=nations;n++)
       			{
       				natmp(n+1)= natmp(n)+nationareas(n);
       				indnatarea(natmp(n),natmp(n+1)-1)=n;

       			}
       			natmp(nations+1) = narea;
       			indnatarea(narea) = nations;

       			

       			pcat.initialize();
       			for(int n=1;n<=nations;n++)
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

	init_int eof;
	
	
	LOC_CALCS

		
		if( eof != 999 )
		{
			cout<< "Error reading data.\n Fix it."<<endl;
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

	init_number log_mo;
	
	init_number log_cvPos;
	init_vector log_maxPos50(1,1);
	init_number log_maxPossd;

	//init_bounded_number log_mo(0,2.484907);
	
	//init_bounded_number log_cvPos(-2.995732,-1.203973);
	//init_bounded_vector log_maxPos50(1,ngroup,0,2.079442);
	//init_vector log_maxPos50(1,ngroup);
	//init_bounded_number log_maxPossd(-0.6931472,1.386294);
	
	init_vector wt(syr-20,nyr,-1);

	objective_function_value f;

	//=================================
	//derived quantities
	//=================================

	number kappa;
	number phie;
	number So;
	number Bo
	number beta;

	number m_tsp;
	
	number mo;
	vector maxPos50(1,ngroup);
	
	number maxPossd;
	number cvPos;
	//vector maxPossd(1,ngroup);
	//vector cvPos(1,ngroup);

	//vector wt(syr,nyr);

	vector lxo(sage,nage);
	vector za(sage,nage);
	vector SB(1,ntstp);
	vector tB(1,ntstp);

	matrix maxPos(1,1,sage,nage);
	matrix varPos(1,1,sage,nage);
	//matrix maxPos(1,ngroup,sage,nage);
	//matrix varPos(1,ngroup,sage,nage);
	
	
	vector nlvec(1,nations);
	
	//vector npvec(1,ngroup);
	
 	matrix Effarea(1,ntstp,sarea,narea);
 	matrix tnage(1,ntstp,sage,nage);
 	matrix tVBarea(1,ntstp,sarea,narea);


	3darray Nage(1,1,1,ntstp,sage-1,nage); 
 	3darray VulB(1,1,1,ntstp,sage-1,nage);
	3darray Effage(1,1,1,ntstp,sage,nage);
 	3darray VBarea(1,1,1,ntstp,sarea,narea);
 	//3darray Nage(1,ngroup,1,ntstp,sage-1,nage); 
 	//3darray VulB(1,ngroup,1,ntstp,sage-1,nage);
	//3darray Effage(1,ngroup,1,ntstp,sage,nage);
 	//3darray VBarea(1,ngroup,1,ntstp,sarea,narea);
	//3darray PosX(1,ngroup,1,ntstp,sage-1,nage);
	3darray PosX(1,1,1,ntstp,sage-1,nage);
 	
 	3darray NAreaAge(1,ntstp,sarea,narea,sage,nage);	
 	3darray CatchAreaAge(itsp,ntstp,sarea,narea,sage,nage);//
 	3darray CatchNatAge(itsp,ntstp,1,nations,sage,nage); //
 	3darray EffNatAge(1,nations,1,ntstp,sage-2,nage); //

 	//4darray propBarea(1,ntstp,1,ngroup,sarea,narea,sage-3,nage);
 	
 	matrix predCatchNatAge(1,tot_pcat,sage-3,nage);

PROCEDURE_SECTION

	incidence_functions();
	initialization();
	burn_in();
	
	move_grow_die();
	
	clean_catage();
	
	calc_obj_func();
	
	//cout<<"SB(1)"<<SB<<endl;
	//cout<<"(Nage(g)(1)(sage,nage)"<<Nage(1)(1)(sage,nage)<<endl;
	//cout<<"(sage,nage)"<<VulB(1)(1)(sage,nage)<<endl;
	//cout<<"fa"<<va<<endl;
	//cout<<"wa"<<wa<<endl;
	

	//exit(1);

FUNCTION dvar_vector cnorm(const double& x, const dvar_vector& mu, const dvar_vector& sd)

	dvar_vector rst(sage,nage);
	dvar_vector stx(sage,nage);


	for(int a= sage; a<= nage; a++)
	{
		stx(a) = (x-mu( a ))/sd( a );
		rst(a) = cumd_norm(stx( a ));
	}

	return(rst);


FUNCTION incidence_functions
	
	maxPos.initialize();

	lxo = mfexp(-m*age);
	lxo(nage) /= 1. - mfexp(-m); 

	kappa 	= 4*h/(1-h);
	phie	= lxo*fa;
	So 		= kappa/phie;
	Bo 		= kappa/So*Ro;
	beta 	= (kappa-1)/Bo;

	m_tsp 	= m/nmon;
	za 		= m_tsp+fe;

	maxPos50 = mfexp(log_maxPos50);
	maxPossd = mfexp(log_maxPossd);
	cvPos 	 = mfexp(log_cvPos);
	mo 	= mfexp(log_mo);


FUNCTION void calc_numbers_at_age(const int& ii, const dvariable& expwt )

	//for(int g=1;g<=ngroup;g++)
	//{
	//		Nage(g)(ii)(sage-1)=g;
	//		
	//		
	//		switch (indmonth(ii)){
    //        	
    //        	case 1:  
    //        		
    //        		Nage(g)(ii)(sage) =((So*SB(ii-nmon)/(1.+beta*SB(ii-nmon)))*expwt)* prop_ng(g);//*1.0/dg;
	//
    //        		for(int a = sage+1;a<=nage;a++)
    //        		{
    // 
    //        			Nage(g)(ii)(a) = Nage(g)(ii-1)(a-1)*mfexp(-(m_tsp+q*Effage(g)(ii-1)(a-1)*va(a-1)));
    //        			
    //        		}
    //        		Nage(g)(ii)(nage) /= (1.-mfexp(-(m_tsp+q*Effage(g)(ii-1)(nage)*va(nage))));
    //        	break;
    //        	
    //        	default: 
	//
    //        		Nage(g)(ii)(sage,nage) = elem_prod(Nage(g)(ii-1)(sage,nage),mfexp(-(m_tsp+q*elem_prod(Effage(g)(ii-1)(sage,nage),va))));
    //       
    //        	break;
    //    	}
            	

    //  	VulB(g)(ii)(sage-1) = g;
	//
	//		VulB(g)(ii)(sage,nage) = elem_prod(elem_prod(Nage(g)(ii)(sage,nage),va),wa);
	
	
	//		tnage(ii)(sage,nage) += Nage(g)(ii)(sage,nage);
	//		tB(ii) += Nage(g)(ii)(sage,nage)*wa;
	//}

	//SB(ii) = elem_prod(tnage(ii),fa)*wa/2;
	
	
			Nage(1)(ii)(sage-1)=g;
			
			
			switch (indmonth(ii)){
            	
            	case 1:  
            		
            		Nage(1)(ii)(sage) =((So*SB(ii-nmon)/(1.+beta*SB(ii-nmon)))*expwt)* prop_ng(1);//*1.0/dg;

            		for(int a = sage+1;a<=nage;a++)
            		{
     
            			Nage(1)(ii)(a) = Nage(1)(ii-1)(a-1)*mfexp(-(m_tsp+q*Effage(1)(ii-1)(a-1)*va(a-1)));
            			
            		}
            		Nage(1)(ii)(nage) /= (1.-mfexp(-(m_tsp+q*Effage(1)(ii-1)(nage)*va(nage))));
            	break;
            	
            	default: 

            		Nage(1)(ii)(sage,nage) = elem_prod(Nage(1)(ii-1)(sage,nage),mfexp(-(m_tsp+q*elem_prod(Effage(1)(ii-1)(sage,nage),va))));
            
            	break;
        	}
            	

        	VulB(1)(ii)(sage-1) = 1;
	
			VulB(1)(ii)(sage,nage) = elem_prod(elem_prod(Nage(1)(ii)(sage,nage),va),wa);
	
	
			tnage(ii)(sage,nage) += Nage(1)(ii)(sage,nage);
			tB(ii) += Nage(1)(ii)(sage,nage)*wa;
	}

	SB(ii) = elem_prod(tnage(ii),fa)*wa/2;
	
				



FUNCTION void calc_effarea(const int& ii,const int& ie)

	dvar_vector tmp1(sarea,narea);
	dvar_vector tmp2(sarea,narea);
	tmp1.initialize();
	tmp2.initialize();

	for(int r= sarea; r<=narea; r++)
	{

		tmp1(r)= (tVBarea(ii)(r)/ sum(tVBarea(ii)))* effPwr(r);
		tmp2(r) = tmp1(r)*TotEffyear(indnatarea(r))(indyr(ie));
		Effarea(ii)(r) = tmp2(r)*TotEffmonth(indnatarea(r))(indmonth(ie));

		//survey option
		//tmp1(rr)= 0.3;
		//tmp2(rr) = tmp1(rr)*TotEffyear(indnatarea(rr))(indyr(1));
		//Effarea(1)(rr) = tmp2(rr)*TotEffmonth(indnatarea(rr))(indmonth(1));
		
	}	

FUNCTION void calc_position(const int& ii)

	//int g, r, ig;

	//for(ig=1;ig<=n_rg ;ig++)
	//{
	//	
	//	
	//	g = n_group(ig);
	//	r = n_area(ig);
	//
	//	varPos(g) = maxPos(g)*cvPos;
	//
	//	PosX(g)(ii)(sage-1) = g;
	//	PosX(g)(ii)(sage,nage) = minPos + (maxPos(g) - minPos) * (0.5+0.5*sin(indmonth(ii)*PI/6 - mo*PI/6)); 
	//
	//	//VBarea(1,sarea) = VulB(1)* (cnorm(areas(sarea)+0.5,PosX(1),varPos));	
	//	
	//	VBarea(g)(ii)(r) = VulB(g)(ii)(sage,nage) * (cnorm(areas(r)+0.5,PosX(g)(ii),varPos(g))-cnorm(areas(r)-0.5,PosX(g)(ii),varPos(g)));
	//	NAreaAge(ii)(r) += elem_prod(Nage(g)(ii)(sage,nage),(cnorm(areas(r)+0.5,PosX(g)(ii),varPos(g))-cnorm(areas(r)-0.5,PosX(g)(ii),varPos(g))));
	//	tVBarea(ii)(r) += VBarea(g)(ii)(r);
	//}

	int r;

	for(int r= sarea; r<=narea; r++)
	{
		//g = n_group(ig);
		//r = n_area(ig);

		varPos(1) = maxPos(1)*cvPos;

		PosX(1)(ii)(sage-1) = 1;
		PosX(1)(ii)(sage,nage) = minPos + (maxPos(1) - minPos) * (0.5+0.5*sin(indmonth(ii)*PI/6 - mo*PI/6)); 
	
		//VBarea(1,sarea) = VulB(1)* (cnorm(areas(sarea)+0.5,PosX(1),varPos));	

		
		VBarea(g)(ii)(r) = VulB(1)(ii)(sage,nage) * (cnorm(areas(r)+0.5,PosX(g)(ii),varPos(1))-cnorm(areas(r)-0.5,PosX(g)(ii),varPos(1)));
		NAreaAge(ii)(r) += elem_prod(Nage(1)(ii)(sage,nage),(cnorm(areas(r)+0.5,PosX(g)(ii),varPos(1))-cnorm(areas(r)-0.5,PosX(g)(ii),varPos(1))));
		tVBarea(ii)(r) += VBarea(1)(ii)(r);
		
	
	}

FUNCTION void calc_effage(const int& ii)

		int g, a;
		
		for(a = sage; a<=nage;a++)
		{
			dvar_matrix propVBarea(1,1,sarea,narea);
			propVBarea.initialize();

			for(int r= sarea; r<=narea; r++)
			{
				
				dvariable temp1;
				dvariable temp2;

				temp1 = (VBarea(1)(ii)(r)+(0.00001* prop_ng(1)))/(tVBarea(ii)(r)+0.00001);
				temp2 = (cnorm(areas(r)+0.5,PosX(g)(ii),varPos(1))-cnorm(areas(r)-0.5,PosX(1)(ii),varPos(1)))(a-sage+1);

				propVBarea(g1)(r) = temp1*temp2;
				//propVBarea(g)(r) = temp2;

				//propBarea(i)(g)(r)(sage-3) = i;
				//propBarea(i)(g)(r)(sage-2) = g;
				//propBarea(i)(g)(r)(sage-1) = r;
				//propBarea(i)(g)(r)(a) =  elem_prod(totB(g)(i)(sage,nage), (cnorm(areas(r)+0.5,PosX(g)(i),varPos(g))-cnorm(areas(r)-0.5,PosX(g)(i),varPos(g))))(a-sage+1);
				
				
				EffNatAge(indnatarea(r))(ii)(sage-2) = ii;
				EffNatAge(indnatarea(r))(ii)(sage-1) = indnatarea(r);
				EffNatAge(indnatarea(r))(ii)(a) += Effarea(ii)(r)* propVBarea(1)(r);

				
				Effage(1)(ii)(a) = Effarea(ii)*propVBarea(1);	

			}

		}

		//for(int a = sage; a<=nage;a++)
		//{
		//	dvar_matrix propVBarea(1,ngroup,sarea,narea);
		//	propVBarea.initialize();
		//
		//	for(int ig=1;ig<=n_rg;ig++)
		//	{
		//		int g, r;
		//		g = n_group(ig);
		//		r = n_area(ig);	
		//		
		//		dvariable temp1;
		//		dvariable temp2;
		//
		//		temp1 = (VBarea(g)(ii)(r)+(0.00001* prop_ng(g)))/(tVBarea(ii)(r)+0.00001);
		//		temp2 = (cnorm(areas(r)+0.5,PosX(g)(ii),varPos(g))-cnorm(areas(r)-0.5,PosX(g)(ii),varPos(g)))(a-sage+1);
		//		
		//		propVBarea(g)(r) = temp1*temp2;
		//		//propVBarea(g)(r) = temp2;
		//
				//propBarea(i)(g)(r)(sage-3) = i;
				//propBarea(i)(g)(r)(sage-2) = g;
				//propBarea(i)(g)(r)(sage-1) = r;
				//propBarea(i)(g)(r)(a) =  elem_prod(totB(g)(i)(sage,nage), (cnorm(areas(r)+0.5,PosX(g)(i),varPos(g))-cnorm(areas(r)-0.5,PosX(g)(i),varPos(g))))(a-sage+1);
				
				
		//		EffNatAge(indnatarea(r))(ii)(sage-2) = ii;
		//		EffNatAge(indnatarea(r))(ii)(sage-1) = indnatarea(r);
		//		EffNatAge(indnatarea(r))(ii)(a) += Effarea(ii)(r)* propVBarea(g)(r);
		//
		//		
		//		Effage(g)(ii)(a) = Effarea(ii)*propVBarea(g);	
		//
		//	}

		//}
	


FUNCTION void calc_catage(const int& ii)

		int a,r;
		for(a = sage; a<=nage;a++)
		{
			for(int r=sarea;r<=narea;r++)
			{		
				CatchAreaAge(ii)(r)(a) = q*Effarea(ii)(r)*va(a)/(q*Effarea(ii)(r)*va(a)+m_tsp)*(1-mfexp(-(q*Effarea(ii)(r)*va(a)+m_tsp)))*NAreaAge(ii)(r)(a);
				CatchNatAge(ii)(indnatarea(r))(a) += CatchAreaAge(ii)(r)(a);

			}

		}
	

FUNCTION initialization
	

	NAreaAge.initialize();
 	CatchAreaAge.initialize();
 	CatchNatAge.initialize();
 	Nage.initialize();

 	//for(int g=1;g<=ngroup;g++)
	//{
	//	Nage(g)(1)(sage-1)=g;	
	//		
	//	Nage(g)(1)(sage) =(So*Bo/(1+beta*Bo))* prop_ng(g);
	//	
    //    for(int a=sage+1 ; a <= nage ; a++)
	//	{
	//		Nage(g)(1)(a) = Nage(g)(1)(a-1) * mfexp(-za(a-1));
	//	}
	//	Nage(g)(1)(nage) /= (1.-mfexp(-za(nage)));
    //        		
    //        	
	//	
    //    VulB(g)(1)(sage-1) = g;
	//
	//	VulB(g)(1)(sage,nage) = elem_prod(elem_prod(Nage(g)(1)(sage,nage),va),wa);
	//
	//
	//	tnage(1)(sage,nage) += Nage(g)(1)(sage,nage);
	//	tB(1) += Nage(g)(1)(sage,nage)*wa;
	//}

		Nage(1)(1)(sage-1)=1;	
			
		Nage(1)(1)(sage) =(So*Bo/(1+beta*Bo))* prop_ng(1);

        for(int a=sage+1 ; a <= nage ; a++)
		{
			Nage(1)(1)(a) = Nage(1)(1)(a-1) * mfexp(-za(a-1));
		}
		Nage(1)(1)(nage) /= (1.-mfexp(-za(nage)));
            		
            	

        VulB(1)(1)(sage-1) = 1;
	
		VulB(1)(1)(sage,nage) = elem_prod(elem_prod(Nage(1)(1)(sage,nage),va),wa);
	
	
		tnage(1)(sage,nage) += Nage(1)(1)(sage,nage);
		tB(1) += Nage(1)(1)(sage,nage)*wa;
	

	
	SB(1) = elem_prod(tnage(1),fa)*wa/2;

	maxPos.initialize();
	calcmaxpos(tB(1));
	
	calc_position(1);
	
	calc_effarea(1,itsp);
	
	calc_effage(1);

	
	

FUNCTION burn_in


	//dvariable dg(ngroup);

	for(int i=2;i<=itsp-1;i++)
	{
		
		
		calc_numbers_at_age(i, mfexp(wt(indyr(i))));		
		
		maxPos.initialize();		
		calcmaxpos(tB(i));	
		calc_position(i);
		calc_effarea(i,itsp);
		calc_effage(i);
	}		
		


FUNCTION move_grow_die

	
	for(int i=itsp;i<=ntstp;i++)
	{
		
		calc_numbers_at_age(i,mfexp(wt(indyr(i))));		
		
		maxPos.initialize();		
		calcmaxpos(tB(i));	
		
		calc_position(i);
		calc_effarea(i,i);
		calc_effage(i);
		calc_catage(i);
	}
        	 


FUNCTION clean_catage

	int p;
       
	p=1;
	for(int i=itsp;i<=ntstp;i++)
	{
		for(int n=1;n<=nations;n++)
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

	//cout<<"predCatchNatAge"<<predCatchNatAge<<endl;


FUNCTION calc_obj_func

	//double tau_c;

	nlvec.initialize();
	//npvec.initialize();

	
	

		for(int n = 1; n<=nations;n++)
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
				
					
				O(i) = (obsCatchNatAge(ii)(sage,nage)+0.1e-30)/sum(obsCatchNatAge(ii)(sage,nage)+0.1e-5);
				
				P(i) = (predCatchNatAge(ii)(sage,nage)+0.1e-30)/sum(predCatchNatAge(ii)(sage,nage)+0.1e-5);
				
				

			}			
			
			
			nlvec(n) =  dmvlogistic(O,P,nu,tau_c,dMinP);
			
		}

	
	
	f=sum(nlvec)/10000;
	//cout<<" nlvec "<<nlvec<<endl;

	
       			
	

FUNCTION dvar_matrix calcmaxpos(const dvariable& tb)

	//for(int g=1;g<=ngroup;g++)
	//{
	//	maxPos(g)(sage,nage) = 1./(1.+mfexp(-(age-maxPos50(g))/maxPossd));
	//	maxPos(g)(sage,nage) *= (narea-sarea);
	//	maxPos(g)(sage,nage) += sarea;			
	//}
	//sreturn(maxPos);


		maxPos(1)(sage,nage) = 1./(1.+mfexp(-(age-maxPos50(1))/maxPossd));
		maxPos(1)(sage,nage) *= (narea-sarea);
		maxPos(1)(sage,nage) += sarea;			
	
		
		return(maxPos);
			
	


REPORT_SECTION
	
	REPORT(mo);
	REPORT(maxPos50);
	REPORT(maxPossd);
	REPORT(cvPos);
	REPORT(syr);
	REPORT(nyr);
	REPORT(sage);
	REPORT(nage);
	REPORT(smon);
	REPORT(nmon);
	REPORT(sarea);
	REPORT(narea);
	REPORT(nations);
	REPORT(border);
	REPORT(ngroup);
	REPORT(maxPos);
	REPORT(minPos);
	REPORT(varPos);
	REPORT(SB);
	REPORT(VulB);
	REPORT(Nage);
	REPORT(VBarea);
	REPORT(Effage);
	REPORT(Effarea);
	REPORT(EffNatAge);
	REPORT(CatchNatAge);
	REPORT(CatchAreaAge);
	REPORT(predCatchNatAge);
	REPORT(indyr);
	REPORT(indmonth);
	REPORT(indnatarea);




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

