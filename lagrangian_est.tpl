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
   	

   	
   		LOC_CALCS		
			ntstp = (nmon-smon+1) * (nyr-syr+1);

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

		int tot_pcat;
		

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
       				for(int i=1;i<=ntstp;i++)
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
	init_number log_maxPos50;
	init_number log_maxPossd;
	init_vector wt(syr,nyr,-1);

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
	
	vector maxPos(sage,nage);
	vector varPos(sage,nage);
	
	vector nlvec(1,nations);
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
 	3darray CatchAreaAge(1,ntstp,sarea,narea,sage,nage);
 	3darray CatchNatAge(1,ntstp,1,nations,sage,nage);
 	3darray EffNatAge(1,nations,1,ntstp,sage-2,nage);
 	
 	matrix predCatchNatAge(1,tot_pcat,sage-3,nage);

PRELIMINARY_CALCS_SECTION
		

PROCEDURE_SECTION

	incidence_functions();
	initialization();
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

	Nage(1,1) = So*Bo/(1+beta*Bo);

	for(int i=sage+1 ; i <= nage ; i++)
	{
		Nage(1,i) = Nage(1,i-1) * mfexp(-za(i-1));
	}

	VulB(1) = elem_prod(elem_prod(Nage(1),va),wa);
	SB(1) = elem_prod(Nage(1),fa)*wa/2;

	maxPos.initialize();
	tBo = Nage(1)*wa;
	calcmaxpos(tBo);
	varPos = maxPos*cvPos;

	PosX(1) = minPos + (maxPos - minPos) * (0.5+0.5*sin(indmonth(1)*PI/6 - mo*PI/6)); 
	
	VBarea(1,sarea) = VulB(1)* (cnorm(areas(sarea)+0.5,PosX(1),varPos));	

	for(int r=sarea ; r <= narea ; r++)
	{
		VBarea(1,r) = VulB(1)* (cnorm(areas(r)+0.5,PosX(1),varPos)-cnorm(areas(r)-0.5,PosX(1),varPos));
		NAreaAge(1)(r) = elem_prod(Nage(1)(sage,nage),(cnorm(areas(r)+0.5,PosX(1),varPos)-cnorm(areas(r)-0.5,PosX(1),varPos)));
	}
	
	for(int n=1; n<=nations;n++){
       	NationVulB(1,n) = sum(VBarea(1)(natmp(n),natmp(n+1)));		 
	}
	
	
	dvar_vector tmp1(sarea,narea);
	dvar_vector tmp2(sarea,narea);
	dvar_vector tmp3(sarea,narea);

	
	for(int rr= sarea; rr<=narea; rr++)
	{
		//tmp1(rr)= pow((VBarea(1)(rr)/ (NationVulB(1)(indnatarea(rr)) + 0.0001))+1.0e-20,effPwr(rr));
		tmp1(rr)= (VBarea(1)(rr)/ (sum(VulB(1)) + 0.0001))* effPwr(rr);
		tmp2(rr) = tmp1(rr)*TotEffyear(indnatarea(rr))(indyr(1));
		Effarea(1)(rr) = tmp2(rr)*TotEffmonth(indnatarea(rr))(indmonth(1));
	}
	
	
	for(int a= sage; a<= nage;a++)
	{
		dvar_vector propVBarea(sarea,narea);
		for(int rr =sarea; rr<=narea; rr++)
		{
			propVBarea(rr) = (cnorm(areas(rr)+0.5,PosX(1),varPos)-cnorm(areas(rr)-0.5,PosX(1),varPos))(a-sage+1);
			CatchAreaAge(1)(rr)(a) = q*Effarea(1)(rr)*va(a)/(q*Effarea(1)(rr)*va(a)+m_tsp)*(1-mfexp(-(q*Effarea(1)(rr)*va(a)+m_tsp)))*NAreaAge(1)(rr)(a);
			CatchNatAge(1)(indnatarea(rr))(a) += CatchAreaAge(1)(rr)(a);

			EffNatAge(indnatarea(rr))(1)(sage-2) = 1;
			EffNatAge(indnatarea(rr))(1)(sage-1) = indnatarea(rr);
			EffNatAge(indnatarea(rr))(1)(a) += Effarea(1)(rr)*propVBarea(rr);

		}

		Effage(1)(a) = Effarea(1)* propVBarea;

	}

FUNCTION move_grow_die

	dvariable tB;

	for(int i=2;i<=ntstp;i++)
	{
		
		switch (indmonth(i)) {
            case 1:           	
            	
            	Nage(i)(sage) = (So*SB(i-nmon)/(1.+beta*SB(i-nmon)))*mfexp(wt(indyr(i)));


            	for(int a = sage+1;a<=nage;a++)
            	{
            		Nage(i)(a) = Nage(i-1)(a-1)*mfexp(-(m_tsp+q*Effage(i-1)(a-1)*va(a-1)));
            	}
            	
            	break;
            	
            default: 
            
            	Nage(i) = elem_prod(Nage(i-1),mfexp(-(m_tsp+q*elem_prod(Effage(i-1),va))));
            	break;
        }
		
		VulB(i) = elem_prod(elem_prod(Nage(i),va),wa);
		SB(i) = elem_prod(Nage(i),fa)*wa/2;
		
		maxPos.initialize();	
		tB = Nage(i)*wa;		
		calcmaxpos(tB);	
		varPos = maxPos*cvPos;

		PosX(i) = minPos + (maxPos - minPos) * (0.5+0.5*sin(indmonth(i)*PI/6 - mo*PI/6)); 

		VBarea(i,sarea) = VulB(i)* (cnorm(areas(sarea)+0.5,PosX(i),varPos));


		for(int r = sarea+1;r <= narea;r++)
		{
			VBarea(i)(r) = VulB(i)*(cnorm(areas(r)+0.5,PosX(i),varPos)-cnorm(areas(r)-0.5,PosX(i),varPos));
			NAreaAge(i)(r) = elem_prod(Nage(i)(sage,nage),(cnorm(areas(r)+0.5,PosX(i),varPos)-cnorm(areas(r)-0.5,PosX(i),varPos)));
		}	

		
		for(int n=1; n<=nations;n++){
       		NationVulB(i,n) = sum(VBarea(i)(natmp(n),natmp(n+1)));		 
		}

		dvar_vector tmp1(sarea,narea);
		dvar_vector tmp2(sarea,narea);

		for(int rr= sarea; rr<=narea; rr++)
		{
			//tmp1(rr)= pow((VBarea(i)(rr)/ (NationVulB(i)(indnatarea(rr)) + 0.0001))+ 1.0e-20,effPwr(rr));
			tmp1(rr)= (VBarea(i)(rr)/ (sum(VulB(i)) + 0.0001))* effPwr(rr);
				
			tmp2(rr) = tmp1(rr)*TotEffyear(indnatarea(rr))(indyr(i));
			Effarea(i)(rr) = tmp2(rr)*TotEffmonth(indnatarea(rr))(indmonth(i));
		}
		
		
		for(int a = sage; a<=nage;a++)
		{
			dvar_vector propVBarea(sarea,narea);
			for(int rr =sarea; rr<=narea; rr++)
			{
				propVBarea(rr) = (cnorm(areas(rr)+0.5,PosX(i),varPos)-cnorm(areas(rr)-0.5,PosX(i),varPos))(a-sage+1);
				
				EffNatAge(indnatarea(rr))(i)(sage-2) = i;
				EffNatAge(indnatarea(rr))(i)(sage-1) = indnatarea(rr);
				EffNatAge(indnatarea(rr))(i)(a) += Effarea(i)(rr)* propVBarea(rr);

			}
			//cout<<"propVBarea "<<propVBarea<<endl;
			//cout<<"Effarea(1) "<<Effarea(1)<<endl;
			Effage(i)(a) = Effarea(i)*propVBarea;
		}

		for(int r = sarea;r <= narea;r++)
		{
			for(int a = sage; a<=nage;a++)
			{
				CatchAreaAge(i)(r)(a) = q*Effarea(i)(r)*va(a)/(q*Effarea(i)(r)*va(a)+m_tsp)*(1-mfexp(-(q*Effarea(i)(r)*va(a)+m_tsp)))*NAreaAge(i)(r)(a);
				CatchNatAge(i)(indnatarea(r))(a)+= CatchAreaAge(i)(r)(a);
			}

		}
		
	}


FUNCTION clean_catage

	int p;
       
	p=1;
	for(int i=1;i<=ntstp;i++)
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



FUNCTION calc_obj_func

	//double tau_c;

	nlvec.initialize();
	npvec.initialize();
	
	
	
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
				//cout<<"ii"<<ii<<endl;
					
				//O(i) = (obsCatchNatAge(ii)(sage,nage)+0.1e-30)/sum(obsCatchNatAge(ii)(sage,nage)+0.1e-5);
				//P(i) = (predCatchNatAge(ii)(sage,nage)+0.1e-30)/sum(predCatchNatAge(ii)(sage,nage)+0.1e-5);
				O(i) = (obsCatchNatAge(ii)(sage,nage)+0.1e-30)/sum(obsCatchNatAge(ii)(sage,nage)+0.1e-5);
				P(i) = (predCatchNatAge(ii)(sage,nage)+0.1e-30)/sum(predCatchNatAge(ii)(sage,nage)+0.1e-5);
				
				//cout<<"O("<<i<<")"<<O(i)<<endl;
				//cout<<"P("<<i<<")"<<P(i)<<endl;
			}			

			//cout<<"dmvlogistic(O,P,nu,tau_c(n),dMinP) is "<< dmvlogistic(O,P,nu,tau_c,dMinP)<<endl;
			//exit(1);									
			//nlvec(n) =  dmvlogistic(O,P,nu,tau_c(n),dMinP);
			nlvec(n) =  dmvlogistic(O,P,nu,tau_c,dMinP);

		}
		
		cout<<"tau_c is"<<tau_c<<endl;
	//exit(1);

	
	//f=sum(nlvec)+sum(npvec);
	f=sum(nlvec)/100;

FUNCTION dvar_vector calcmaxpos(const dvariable& tb)

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

