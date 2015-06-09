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
		seed += 10; // add 10 to the seed
		ofstream ofs( "seed.txt" ); //put out to seed.txt
		ofs<<seed<<endl; //the new value of the seed
	END_CALCS

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
	init_number tau_c;
	init_number mo;

	init_vector wa(sage,nage);
	init_vector fa(sage,nage);
	init_vector va(sage,nage);
	init_vector minPos(sage,nage);
	//init_vector maxPos(sage,nage);
	init_number maxPos50;
	init_number maxPossd;
	init_number cvPos;
	

	init_matrix TotEffyear(1,nations,syr,nyr);
	init_matrix TotEffmonth(1,nations,smon,nmon);

	init_int eof;
	
	
	LOC_CALCS
		
		if( eof != 999 )
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			cout<< "eof is: "<<eof<<endl;
			ad_exit(1);
		}

	END_CALCS
	

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
			for(int n=1; n<=nations-1; n++)
			{
				for(int a=sarea;a<=narea;a++)
				{
					if(areas(a)<border(n))
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

	 	vector indyr(1,ntstp);
 		ivector indmonth(1,ntstp);
		vector indnatarea(sarea,narea);
		ivector pcat(sarea,narea);

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
      		
       			for(int b = sarea; b <= sarea+nationareas(1) ; b++)
       			{
       				indnatarea(b)=1;
       			}
     
       			for(int bb = 2;bb<=nations;bb++)
       			{
       				for(int dd = sarea+nationareas(bb-1)+1;dd<=sarea+sum(nationareas(1,bb))-1;dd++)
       				{
       					indnatarea(dd)=bb;
       				}	 	
       			}

       			
       			pcat.initialize();
       			for(int r=sarea;r<=narea;r++)
       			{
       				for(int i=1;i<=ntstp;i++)
       				{

       					if(TotEffmonth(indnatarea(r))(indmonth(i))>0)
       					{
       						pcat(r)++;
       					}
       							
       				}
       			}

       			tot_pcat=sum(pcat);
       		


       			
	END_CALCS


PARAMETER_SECTION

	objective_function_value no_f;

	//derived quantities
	number kappa;
	number phie;
	number So;
	number Bo
	number beta;

	vector lxo(sage,nage);
	vector za(sage,nage);

	vector SB(1,ntstp);
	vector varPos(sage,nage);
	vector maxPos(sage,nage);


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
 	3darray EffNatAge(1,nations,1,ntstp,sage-2,nage);

 	matrix obsCatchAreaAge(1,tot_pcat,sage-3,nage);


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
  

FUNCTION incidence_functions

	maxPos.initialize();
	
	lxo = exp(-m*age);
	lxo(nage) /= 1. - exp(-m); 

	kappa 	= 4*h/(1-h);
	phie	= lxo*fa;
	So 		= kappa/phie;
	Bo 		= kappa/So*Ro;
	beta 	= (kappa-1)/Bo;

	za 		= m+va*fe;

	maxPos(sage,nage) = 1./(1.+exp(-(age-maxPos50)/maxPossd));
	maxPos(sage,nage) *= (narea-sarea);
	maxPos(sage,nage) += sarea;

	varPos = maxPos*cvPos;

FUNCTION initialization
	
	NAreaAge.initialize();
 	CatchAreaAge.initialize();

	Nage(1,1) = So*Bo/(1+beta*Bo);

	for(int i=sage+1 ; i <= nage ; i++)
	{
		Nage(1,i) = Nage(1,i-1) * exp(-za(i-1));
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
			CatchAreaAge(1)(rr)(a) = q*Effarea(1)(rr)*va(a)/(q*Effarea(1)(rr)*va(a)+m)*(1-exp(-(q*Effarea(1)(rr)*va(a)+m)))*NAreaAge(1)(rr)(a);

			EffNatAge(indnatarea(rr))(1)(sage-2) = 1;
			EffNatAge(indnatarea(rr))(1)(sage-1) = indnatarea(rr);
			EffNatAge(indnatarea(rr))(1)(a) += Effarea(1)(rr)*propVBarea(1)(rr);
		}
		Effage(1)(a) = Effarea(1)* propVBarea(1);
	}


	//CatchAgeArea



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
				CatchAreaAge(i)(rr)(a) = q*Effarea(i)(rr)*va(a)/(q*Effarea(i)(rr)*va(a)+m)*(1-exp(-(q*Effarea(i)(rr)*va(a)+m)))*NAreaAge(i-1)(rr)(a);
				
				EffNatAge(indnatarea(rr))(i)(sage-2) = i;
				EffNatAge(indnatarea(rr))(i)(sage-1) = indnatarea(rr);
				EffNatAge(indnatarea(rr))(i)(a) += Effarea(i)(rr)* propVBarea(i)(rr);
			}
			Effage(i)(a) = Effarea(i)*propVBarea(i);
		}

		//cout<<"indmonth(i) is "<<indmonth(i)<<endl;
		

		switch (indmonth(i)) {
            case 1:
            	Nage(i) = elem_prod(Nage(i-1),exp(-(m+q*elem_prod(Effage(i),va))/12));
            	Nage(i,sage) = So*SB(i-nmon)/(1.+beta*SB(i-nmon))*exp(wt(indyr(i)));

            default:
               Nage(i) = elem_prod(Nage(i-1),exp(-(m+q*elem_prod(Effage(i),va))/12));

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
		for(int r=sarea;r<=narea;r++)
		{
							
			if(TotEffmonth(indnatarea(r))(indmonth(i))>0)
       		{
       			
       			dvector pa(nage,sage);
       			pa.initialize();

       			obsCatchAreaAge(p)(sage-3) = i;
       			obsCatchAreaAge(p)(sage-2) = indmonth(i);
				obsCatchAreaAge(p)(sage-1) = r;
       			pa = value((CatchAreaAge(i)(r)(sage,nage)+0.1e-30)/sum(CatchAreaAge(i)(r)(sage,nage)+0.1e-30));
				obsCatchAreaAge(p)(sage,nage) = rmvlogistic(pa,tau_c,seed+i);
				
				 p++;	
       		}	
		}
	}
			

FUNCTION output_true
	
	ofstream ofs("lagrangian_OM.rep");

	ofs<<"mo" << endl << mo <<endl;
	ofs<<"log_tau_c" << endl << log(tau_c) <<endl;
	ofs<<"maxPos50" << endl << maxPos50 <<endl;
	ofs<<"maxPossd" << endl << maxPossd <<endl;
	ofs<<"cvPos" << endl << cvPos <<endl;
	ofs<<"syr" << endl << syr <<endl;
	ofs<<"nyr" << endl << nyr <<endl;
	ofs<<"sage" << endl << sage <<endl;
	ofs<<"nage" << endl << nage <<endl;
	ofs<<"smon" << endl << smon <<endl;
	ofs<<"nmon" << endl << nmon <<endl;
	ofs<<"sarea" << endl << sarea <<endl;
	ofs<<"narea" << endl << narea <<endl;
	ofs<<"nations" << endl << nations <<endl;
	ofs<<"maxPos" << endl << maxPos <<endl;
	ofs<<"minPos" << endl << minPos <<endl;
	ofs<<"varPos" << endl << varPos <<endl;
	ofs<<"VBarea" << endl << VBarea <<endl;
	ofs<<"Effage" << endl << Effage <<endl;
	ofs<<"Effarea"<< endl << Effarea <<endl;
	ofs<<"EffNatAge"<< endl << EffNatAge<<endl;
	
	
FUNCTION output_pin
	
	ofstream ifs("lagrangian_est.pin");

	ifs<<"# mo " << endl << mo <<endl;
	ifs<<"# tau_c " << endl << log(tau_c) <<endl;
	ifs<<"# cvPos "<< endl << log(.1) <<endl;	
	//ifs<<"# maxPos "<< endl << minPos <<endl;
	ifs<<"# maxPos50 "<< endl << log(4) <<endl;
	ifs<<"# maxPossd "<< endl << log(.5) <<endl;


FUNCTION output_dat

	ofstream afs("lagrangian_est.dat");
	afs<<"# syr " << endl << syr <<endl;
	afs<<"# nyr " << endl << nyr <<endl;
	afs<<"# sage " << endl << sage <<endl;
	afs<<"# nage " << endl << nage <<endl;
	afs<<"# smon " << endl << smon <<endl;
	afs<<"# nmon " << endl << nmon <<endl;
	afs<<"# sarea " << endl << sarea <<endl;
	afs<<"# narea " << endl << narea <<endl;
	afs<<"# nations " << endl << nations <<endl;
	afs<<"# border " << endl << border <<endl;
	afs<<"# Ro " << endl << Ro <<endl;
	afs<<"# h " << endl << h <<endl;
	afs<<"# m " << endl << m <<endl;
	afs<<"# fe " << endl << fe <<endl;
	afs<<"# q " << endl << q <<endl;
	afs<<"# sigR " << endl << sigR <<endl;
	afs<<"# weight at age " << endl << wa <<endl;
	afs<<"# fecundity at age " << endl << fa <<endl;
	afs<<"# vulnerability at age " << endl << va <<endl;
	afs<<"# minPos "<< endl << minPos <<endl;
	afs<<"# Total effort by country and year " << endl << TotEffyear <<endl;
	afs<<"# Total effort by country and month " << endl << TotEffmonth <<endl;
	afs<<"# dMinP " << endl << 0.00001 <<endl;
	afs<<"#  tstp month area catage " << endl << obsCatchAreaAge <<endl;
	afs<<"# eof " << endl << 999 <<endl;
	


REPORT_SECTION

	REPORT(VBarea);
	REPORT(Effage);
	REPORT(Effarea);



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

