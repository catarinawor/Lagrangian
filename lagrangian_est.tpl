//  ******************************************************************
//  Lagragian movement model
//  
//  Created by Catarina Wor on 2015-04-13.
//  Copyright (c) 2015. All rights reserved.
//  Comments:
//  ******************************************************************


DATA_SECTION


	// model dimensions
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


	//model parameters
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

	//control parameters
	init_number dMinP;
	


	// accessory quantities

	int ntstp;


   	vector age(sage,nage);
   	vector areas(sarea,narea);
   	ivector nationareas(1,nations);
   	vector wt(syr,nyr);

   	
   		LOC_CALCS		
			ntstp = (nmon-smon+1) * (nyr-syr+1);

			age.fill_seqadd(sage,1);
			areas.fill_seqadd(sarea,1);

			nationareas.initialize();

			dvector natmp(1,nations);
			

			natmp(1)=sarea;

			for(int n=1; n<=nations-1; n++)
			{
				natmp(n+1)=border(n);
				for(int a=sarea;a<=narea;a++)
				{
					if(areas(a)>=natmp(n)&areas(a)<border(n))
					{
						nationareas(n)++;
					}
				}
			}
		
			nationareas(nations)=narea-sarea+1 - sum(nationareas(1,nations-1));

	END_CALCS

	 	vector indyr(1,ntstp);
 		ivector indmonth(1,ntstp);
		vector indnatarea(sarea,narea);
		ivector pcat(1,nations);

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
      		
       			ivector natmp1(1,nations+1);

       			natmp1(1) = sarea;

       			for(int n=1;n<=nations;n++)
       			{
       				natmp1(n+1)= natmp1(n)+nationareas(n);
       				for(int b = natmp1(n); b <= natmp1(n+1)-1 ; b++)
       				{
       					indnatarea(b)=n;
       				}

       			}
       			indnatarea(narea)=nations;

       			pcat.initialize();
       			for(int n=1;n<=nations;n++)
       			{
       				for(int i=1;i<=ntstp;i++)
       				{

       					if(TotEffmonth(n)(indmonth(i))>0)
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
			cout<< "pcat is: "<<pcat<<endl;
			cout<< "eof is: "<<eof<<endl;
			ad_exit(1);
		}

	END_CALCS
	


PARAMETER_SECTION

	init_bounded_number mo(smon,nmon);
	init_number log_tau_c(-2);
	
	init_number log_cvPos;
	init_number log_maxPos50;
	init_number log_maxPossd;
	init_vector wt(syr,nyr,-1);
	//init_bounded_vector maxPos(sage,nage,sarea,narea);


	objective_function_value f;

	//derived quantities
	number kappa;
	number phie;
	number So;
	number Bo
	number beta;
	number maxPos50;
	number maxPossd;
	number cvPos;
	number tau_c;

	vector lxo(sage,nage);
	vector za(sage,nage);
	vector SB(1,ntstp);
	vector nlvec(1,nations);
	vector npvec(1,1);
	vector maxPos(sage,nage);
	vector varPos(sage,nage);


	matrix NationVulB(1,ntstp,1,nations);
	matrix Nage(1,ntstp,sage,nage);
 	matrix VulB(1,ntstp,sage,nage);
 	matrix PosX(1,ntstp,sage,nage);
 	matrix Effage(1,ntstp,sage,nage);
 	matrix VBarea(1,ntstp,sarea,narea);
 	matrix propVBarea(1,ntstp,sarea,narea);
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

	za 		= m+va*fe;

	maxPos50 = mfexp(log_maxPos50);
	maxPossd = mfexp(log_maxPossd);
	cvPos 	 = mfexp(log_cvPos);
	tau_c = mfexp(log_tau_c);

	maxPos(sage,nage) = 1./(1.+mfexp(-(age-maxPos50)/maxPossd));
	maxPos(sage,nage) *= (narea-sarea);
	maxPos(sage,nage) += sarea;

	varPos = maxPos *cvPos;


FUNCTION initialization

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

	for(int ii=1 ; ii <= ntstp ; ii++)
	{
		PosX(ii) = minPos + (maxPos - minPos) * (0.5+0.5*sin(indmonth(ii)*PI/6 - mo*PI/6)); 
	}

	for(int r=sarea ; r <= narea ; r++)
	{
		VBarea(1,r) = VulB(1)* (cnorm(areas(r)+0.5,PosX(1),varPos)-cnorm(areas(r)-0.5,PosX(1),varPos));
		NAreaAge(1)(r) = elem_prod(Nage(1)(sage,nage),(cnorm(areas(r)+0.5,PosX(1),varPos)-cnorm(areas(r)-0.5,PosX(1),varPos)));
	}
	
	NationVulB(1,1) = sum(VBarea(1)(sarea,sarea+nationareas(1)-1)); 
	NationVulB(1,2) = sum(VBarea(1)(sarea+nationareas(1),narea)); 
	

	dvar_vector tmp1(sarea,narea);
	dvar_vector tmp2(sarea,narea);
	dvar_vector tmp3(sarea,narea);

	
	for(int rr= sarea; rr<=narea; rr++)
	{
		tmp1(rr)= VBarea(1)(rr)/ (NationVulB(1)(indnatarea(rr)) + 0.0001);
		tmp2(rr) = tmp1(rr)*TotEffyear(indnatarea(rr))(indyr(1));
		Effarea(1)(rr) = tmp2(rr)*TotEffmonth(indnatarea(rr))(indmonth(1));
	}
	
	
	for(int a= sage; a<= nage;a++)
	{
		for(int rr =sarea; rr<=narea; rr++)
		{
			propVBarea(1)(rr) = (cnorm(areas(rr)+0.5,PosX(1),varPos)-cnorm(areas(rr)-0.5,PosX(1),varPos))(a-sage+1);
			CatchAreaAge(1)(rr)(a) = ((q*Effarea(1)(rr)*va(a))/12)/((q*Effarea(1)(rr)*va(a)+m)/12)*(1-mfexp(-(q*Effarea(1)(rr)*va(a)+m)/12))*NAreaAge(1)(rr)(a);
			CatchNatAge(1)(indnatarea(rr))(a) += CatchAreaAge(1)(rr)(a);

			EffNatAge(indnatarea(rr))(1)(sage-2) = 1;
			EffNatAge(indnatarea(rr))(1)(sage-1) = indnatarea(rr);
			EffNatAge(indnatarea(rr))(1)(a) += Effarea(1)(rr)*propVBarea(1)(rr);
		}
		Effage(1)(a) = Effarea(1)* propVBarea(1);
	}


FUNCTION move_grow_die

	for(int i=2;i<=ntstp;i++)
	{

		dvar_vector tmp1(sarea,narea);
		dvar_vector tmp2(sarea,narea);

		for(int rr= sarea; rr<=narea; rr++)
		{
			tmp1(rr)= VBarea(i-1)(rr)/ (NationVulB(i-1)(indnatarea(rr)) + 0.0001);
			
			tmp2(rr) = tmp1(rr)*TotEffyear(indnatarea(rr))(indyr(i));
			Effarea(i)(rr) = tmp2(rr)*TotEffmonth(indnatarea(rr))(indmonth(i));
		}
		

		for(int a = sage; a<=nage;a++)
		{
			for(int rr =sarea; rr<=narea; rr++)
			{
				propVBarea(i)(rr) = (cnorm(areas(rr)+0.5,PosX(i),varPos)-cnorm(areas(rr)-0.5,PosX(i),varPos))(a-sage+1);
				CatchAreaAge(i)(rr)(a) = ((q*Effarea(i)(rr)*va(a))/12)/((q*Effarea(i)(rr)*va(a)+m)/12)*(1-mfexp(-(q*Effarea(i)(rr)*va(a)+m)/12))*NAreaAge(i-1)(rr)(a);
				CatchNatAge(i)(indnatarea(rr))(a)+= CatchAreaAge(i)(rr)(a);


				EffNatAge(indnatarea(rr))(i)(sage-2) = i;
				EffNatAge(indnatarea(rr))(i)(sage-1) = indnatarea(rr);
				EffNatAge(indnatarea(rr))(i)(a) += Effarea(i)(rr)* propVBarea(i)(rr);
			}
			Effage(i)(a) = Effarea(i)*propVBarea(i);
		}

		//cout<<"cheguei aqui"<<endl;
		

		switch (indmonth(i)) {
            case 1:
            	
            	Nage(i,sage) = So*SB(i-nmon)/(1.+beta*SB(i-nmon))*mfexp(wt(indyr(i)));
            	for(int a = sage+1;a<=nage;a++)
            	{
            		Nage(i)(a) = Nage(i-1)(a-1)*mfexp(-(m+q*Effage(i)(a-1)*va(a-1))/12);
            	}
            	
            	break;

            default:
               Nage(i) = elem_prod(Nage(i-1),mfexp(-(m+q*elem_prod(Effage(i),va))/12));
               break;
        }
		
		VulB(i) = elem_prod(elem_prod(Nage(i),va),wa);
		SB(i) = elem_prod(Nage(i),fa)*wa/2;

		for(int r = sarea;r <= narea;r++)
		{
			VBarea(i)(r) = VulB(i)* (cnorm(areas(r)+0.5,PosX(i),varPos)-cnorm(areas(r)-0.5,PosX(i),varPos));
			NAreaAge(i)(r) = elem_prod(Nage(i)(sage,nage),(cnorm(areas(r)+0.5,PosX(i),varPos)-cnorm(areas(r)-0.5,PosX(i),varPos)));
		}


		NationVulB(i,1) = sum(VBarea(i)(sarea,sarea+nationareas(1)-1)); 
		NationVulB(i,2) = sum(VBarea(i)(sarea+nationareas(1),narea)); 
		
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
					
				O(i) = obsCatchNatAge(ii)(sage,nage)+0.1e-15/sum(obsCatchNatAge(ii)(sage,nage)+0.1e-15);
				P(i) = predCatchNatAge(ii)(sage,nage)+0.1e-15/sum(predCatchNatAge(ii)(sage,nage)+0.1e-15);
			}
				
								
			nlvec(n)=  dmvlogistic(O,P,nu,value(tau_c),dMinP);
		}
		
	//npvec(1) = dnorm(log_tau_c,0.1);

	
	f=sum(nlvec)+sum(npvec);




REPORT_SECTION
	
	REPORT(mo);
	REPORT(log_tau_c);
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

