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

	

	//init_int ngroup;
	//init_vector prop_ng(1,ngroup);

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
	

	init_matrix TotEffyear(1,fisharea,syr,nyr);
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

	init_int eof;
	
	
	LOC_CALCS
		
		if( eof != 999 )
		{
			cout<< "minPos is: "<<minPos <<endl;
			cout<< "effPwr is: "<<effPwr <<endl;
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

	init_bounded_number log_mo(0,1.94591);
	//
	init_bounded_number log_cvPos(-3,-1.3);
	init_bounded_number log_maxPos50(1.00,2.1);
	init_bounded_number log_maxPossd(-0.7,1.609438);
	//init_number log_mo;
	//init_number log_cvPos7
	//init_number log_maxPos50;
	//init_number log_maxPossd;

	init_bounded_vector wt(syr,nyr,-3.0,3.0,-1);

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
	number maxPos50;
	number maxPossd;
	number cvPos;

	//vector wt(syr,nyr);

	vector lxo(sage,nage);
	vector za(sage,nage);
	vector SB(1,ntstp);
	vector tB(1,ntstp);

	vector maxPos(sage,nage);
	vector varPos(sage,nage);
	
	
	//vector maxPos(sage,nage);
	//vector varPos(sage,nage);
	
	vector nlvec(1,fisharea);
	vector npvec(1,1);

	matrix NationVulB(1,ntstp,1,nations);
	matrix Nage(1,ntstp,sage,nage);
 	matrix VulB(1,ntstp,sage,nage);
 	matrix PosX(1,ntstp,sage,nage);
 	matrix Effage(1,ntstp,sage,nage);
 	matrix VBarea(1,ntstp,sarea,narea);
 	
 	//matrix propVBarea(1,ntstp,sarea,narea);
 	matrix Effarea(1,ntstp,sarea,narea);
	
 	3darray NAreaAge(1,ntstp,sarea,narea,sage,nage);	
 	3darray CatchAreaAge(itsp,ntstp,sarea,narea,sage,nage);
 	3darray CatchNatAge(itsp,ntstp,1,fisharea,sage-2,nage);
 	3darray EffNatAge(1,fisharea,itsp,ntstp,sage-2,nage);
 	
 	matrix predCatchNatAge(1,tot_pcat,sage-3,nage);

PRELIMINARY_CALCS_SECTION
		

PROCEDURE_SECTION

	incidence_functions();

	initialization();
	burn_in();
	move_grow_die();
	clean_catage();
	
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
            	
            	Nage(ii)(sage) = (So*SB(ii-nmon)/(1.+beta*SB(ii-nmon)))*(mfexp(expwt));

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

		cout<<"OK after calc_numbers_at_age"<<endl;		

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

		}

		//cout<<"OK after calc_catage"<<endl;



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
	za 		= m_tsp+va*fe;

	maxPos50 = mfexp(log_maxPos50);
	maxPossd = mfexp(log_maxPossd);
	cvPos 	 = mfexp(log_cvPos);
	mo 	= mfexp(log_mo);
	
	

FUNCTION initialization
	
	dvariable tBo;

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

	VulB(1) = elem_prod(elem_prod(Nage(1),va),wa);
	SB(1) = elem_prod(Nage(1),fa)*wa/2;

	maxPos.initialize();
	tBo = Nage(1)*wa;
	
	calcmaxpos();
	
	calc_position(1);
	calc_effarea(1,itsp);	
	

FUNCTION burn_in

	dvariable tB;

	for(int i=2;i<=itsp-1;i++)
	{
		
		calc_numbers_at_age(i,0.0);	
		
		maxPos.initialize();
		dvariable tB;	
		tB = Nage(i)*wa;		
		calcmaxpos();

		calc_position(i);	
		calc_effarea(i,itsp);	
		
	}


FUNCTION move_grow_die

	dvariable tB;

	for(int i=itsp;i<=ntstp;i++)
	{
		
		calc_numbers_at_age(i,wt(indyr(i)));	
		
		maxPos.initialize();	
		tB = Nage(i)*wa;		
		calcmaxpos();	
		calc_position(i);
		calc_effarea(i,i);
		
		calc_catage(i);

		//cout<< "i is "<<i << endl;
		
	} 


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
				//cout<<"ii is "<<ii<<endl;
				//cout<<"ii"<<ii<<endl;
					
				//O(i) = (obsCatchNatAge(ii)(sage,nage)+0.1e-30)/sum(obsCatchNatAge(ii)(sage,nage)+0.1e-5);
				//P(i) = (predCatchNatAge(ii)(sage,nage)+0.1e-30)/sum(predCatchNatAge(ii)(sage,nage)+0.1e-5);
				O(i) = (obsCatchNatAge(ii)(sage,nage))/sum(obsCatchNatAge(ii)(sage,nage));
				P(i) = ((predCatchNatAge(ii)(sage,nage))/(sum(predCatchNatAge(ii)(sage,nage))+0.01))+0.000001;
				

				//cout<<"P("<<i<<")"<<P(i)<<endl;
				//cout<<"O("<<i<<")"<<O(i)<<endl;
				
				
			}			

			//cout<<"dmvlogistic(O,P,nu,tau_c(n),dMinP) is "<< dmvlogistic(O,P,nu,tau_c,dMinP)<<endl;
			//exit(1);									
			//nlvec(n) =  dmvlogistic(O,P,nu,tau_c(n),dMinP);
			cout<<"mo is"<<mo<<endl;
			cout<<"cvPos is"<<cvPos<<endl;
			cout<<"maxPos50 is"<<maxPos50<<endl;
			//cout<<"tau_c is"<<tau_c<<endl;
			//cout<<"dMinP is"<<dMinP<<endl;
			//cout<<"pcat(n) is"<<pcat(n)<<endl;
			cout<<"maxPossd is"<<maxPossd<<endl;
			//cout<<"dmvlogistic(O,P,nu,tau_c,dMinP) is "<<dmvlogistic(O,P,nu,tau_c,dMinP)<<endl;
			nlvec(n) =  dmvlogistic(O,P,nu,tau_c,dMinP);

		}
		
		cout<<"nlvec is"<<nlvec<<endl;
		//cout<<"maxPos50 is "<<maxPos50<<endl;
		//cout<<"maxPossd is "<<maxPossd<<endl;
		//cout<<"cvPos is "<<cvPos<<endl;
		//cout<<"mo is "<<mo<<endl;
	//exit(1);

	
	//f=sum(nlvec)+sum(npvec);
	f=sum(nlvec)/10000;
	


FUNCTION dvar_vector calcmaxpos()

	maxPos(sage,nage) = 1./(1.+mfexp(-(age-maxPos50)/maxPossd));
	maxPos(sage,nage) *= (narea-sarea);
	maxPos(sage,nage) += sarea;

			
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

