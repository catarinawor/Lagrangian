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
	init_number err;

	init_vector wa(sage,nage);
	init_vector fa(sage,nage);
	init_vector va(sage,nage);
	init_vector minPos(sage,nage);
	init_number maxPos501;
	init_number maxPos502;
	init_number maxPossd1;
	init_number maxPossd2;
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
		
			random_number_generator rng(seed);
			wt.fill_randn(rng);
			wt*=sigR;

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


PARAMETER_SECTION

	objective_function_value no_f;

	//derived quantities
	number kappa;
	number phie;
	number So;
	number Bo
	number beta;
	number tBo;

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
 	//matrix propVBarea(1,ntstp,sarea,narea);
 	matrix Effarea(1,ntstp,sarea,narea);
 
 	3darray NAreaAge(1,ntstp,sarea,narea,sage,nage);
 	3darray CatchAreaAge(1,ntstp,sarea,narea,sage,nage);
 	3darray CatchNatAge(1,ntstp,1,nations,sage,nage);
 	3darray EffNatAge(1,nations,1,ntstp,sage-2,nage);

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
	tBo = Nage(1)*wa;

	calcmaxpos(tBo);
	varPos = maxPos*cvPos;

	PosX(1) = minPos + (maxPos - minPos) * (0.5+0.5*sin(indmonth(1)*PI/6 - mo*PI/6)); 
		
	VBarea(1,sarea) = VulB(1)* (cnorm(areas(sarea)+0.5,PosX(1),varPos));
	
	for(int r=sarea+1 ; r <= narea-1 ; r++)
	{
		VBarea(1,r) = VulB(1)* (cnorm(areas(r)+0.5,PosX(1),varPos)-cnorm(areas(r)-0.5,PosX(1),varPos));
		NAreaAge(1)(r) = elem_prod(Nage(1)(sage,nage),(cnorm(areas(r)+0.5,PosX(1),varPos)-cnorm(areas(r)-0.5,PosX(1),varPos)));
	}

	//VBarea(1,narea) = VulB(1)* (1.0-cnorm(areas(narea)-0.5,PosX(1),varPos));

	
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
		dvar_vector propVBarea(sarea,narea);
		for(int rr =sarea; rr<=narea; rr++)
		{
			propVBarea(rr) = (cnorm(areas(rr)+0.5,PosX(1),varPos)-cnorm(areas(rr)-0.5,PosX(1),varPos))(a-sage+1);
			CatchAreaAge(1)(rr)(a) = q*Effarea(1)(rr)*va(a)/(q*Effarea(1)(rr)*va(a)+m)*(1-mfexp(-(q*Effarea(1)(rr)*va(a)+m)))*NAreaAge(1)(rr)(a);
			CatchNatAge(1)(indnatarea(rr))(a) += CatchAreaAge(1)(rr)(a);

			EffNatAge(indnatarea(rr))(1)(sage-2) = 1;
			EffNatAge(indnatarea(rr))(1)(sage-1) = indnatarea(rr);
			EffNatAge(indnatarea(rr))(1)(a) += Effarea(1)(rr)*propVBarea(rr);

		}

		//cout<<"propVBarea "<<propVBarea<<endl;
		//cout<<"Effarea(1) "<<Effarea(1)<<endl;
		Effage(1)(a) = Effarea(1)* propVBarea;

	}




FUNCTION move_grow_die

	dvariable tB;

	for(int i=2;i<=ntstp;i++)
	{
		//if(i>12)exit(1);

		
		switch (indmonth(i)) {
            case 1:           	
            	
            	Nage(i)(sage) = (So*SB(i-nmon)/(1.+beta*SB(i-nmon)))*mfexp(wt(indyr(i))*err);


            	for(int a = sage+1;a<=nage;a++)
            	{
            		Nage(i)(a) = Nage(i-1)(a-1)*mfexp(-(m+q*Effage(i-1)(a-1)*va(a-1))/12);
            	}
            	
            	break;
            	
            default: 
            
            	Nage(i) = elem_prod(Nage(i-1),mfexp(-(m+q*elem_prod(Effage(i-1),va))/12));
            	break;
        }
		
		VulB(i) = elem_prod(elem_prod(Nage(i),va),wa);
		SB(i) = elem_prod(Nage(i),fa)*wa/2;
		
		maxPos.initialize();
		tB = Nage(i)*wa;
		calcmaxpos(tB);
		//cout<<"maxPos "<<maxPos<<endl;
		varPos = maxPos*cvPos;

		PosX(i) = minPos + (maxPos - minPos) * (0.5+0.5*sin(indmonth(i)*PI/6 - mo*PI/6)); 

		VBarea(i,sarea) = VulB(i)* (cnorm(areas(sarea)+0.5,PosX(i),varPos));


		for(int r = sarea+1;r <= narea;r++)
		{
			VBarea(i)(r) = VulB(i)* (cnorm(areas(r)+0.5,PosX(i),varPos)-cnorm(areas(r)-0.5,PosX(i),varPos));
			NAreaAge(i)(r) = elem_prod(Nage(i)(sage,nage),(cnorm(areas(r)+0.5,PosX(i),varPos)-cnorm(areas(r)-0.5,PosX(i),varPos)));
		}	
			
		

		//VBarea(i,narea) = VulB(i)* (1.0-cnorm(areas(narea)-0.5,PosX(i),varPos));


		NationVulB(i,1) = sum(VBarea(i)(sarea,sarea+nationareas(1)-1)); 
		NationVulB(i,2) = sum(VBarea(i)(sarea+nationareas(1),narea)); 

		dvar_vector tmp1(sarea,narea);
		dvar_vector tmp2(sarea,narea);

		for(int rr= sarea; rr<=narea; rr++)
		{
			tmp1(rr)= VBarea(i)(rr)/ (NationVulB(i)(indnatarea(rr)) + 1);
			
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

		for(int r = sarea+1;r <= narea-1;r++)
		{
			for(int a = sage; a<=nage;a++)
			{
				CatchAreaAge(i)(r)(a) = q*Effarea(i)(r)*va(a)/(q*Effarea(i)(r)*va(a)+m)*(1-mfexp(-(q*Effarea(i)(r)*va(a)+m)))*NAreaAge(i)(r)(a);
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
       			
       			dvector pa(nage,sage);
       			pa.initialize();

       			obsCatchNatAge(p)(sage-3) = i;
       			obsCatchNatAge(p)(sage-2) = indmonth(i);
				obsCatchNatAge(p)(sage-1) = n;
       			pa = value((CatchNatAge(i)(n)(sage,nage)+0.1e-30)/sum(CatchNatAge(i)(n)(sage,nage)+0.1e-30));
				obsCatchNatAge(p)(sage,nage) = rmvlogistic(pa,tau_c,seed+i);
				
				p++;	
       		}	
		}
	}
	
FUNCTION dvar_vector calcmaxpos(const dvariable& tb)
	
	int pp = tb > 0.5*tBo?1:2;

	//cout<<"tb "<<tb<<endl;
	//cout<<"0.7tBo "<<0.7*tBo<<endl;
	//cout<<"pp "<<pp<<endl;

	switch (pp) {
        
        case 1:           	

			maxPos(sage,nage) = 1./(1.+mfexp(-(age-maxPos501)/maxPossd1));
			maxPos(sage,nage) *= (narea-sarea);
			maxPos(sage,nage) += sarea;

			cout<<"high B"<<endl;

    	break;
    
    	case 2: 

    		maxPos(sage,nage) = 1./(1.+mfexp(-(age-maxPos502)/maxPossd2));
			maxPos(sage,nage) *= (narea-sarea);
			maxPos(sage,nage) += sarea; 
			cout<<"low B"<<endl;

		break;
	}

	return(maxPos);


	


FUNCTION output_true
	
	ofstream ofs("lagrangian_OM.rep");

	ofs<<"mo" << endl << mo <<endl;
	ofs<<"log_tau_c" << endl << log(tau_c) <<endl;
	ofs<<"maxPos501" << endl << maxPos501 <<endl;
	ofs<<"maxPos502" << endl << maxPos502 <<endl;
	ofs<<"maxPossd1" << endl << maxPossd1 <<endl;
	ofs<<"maxPossd2" << endl << maxPossd2 <<endl;
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
	ofs<<"SB" << endl << SB <<endl;
	ofs<<"VulB" << endl << VulB <<endl;
	ofs<<"Nage" << endl << Nage <<endl;
	ofs<<"VBarea" << endl << VBarea <<endl;
	ofs<<"Effage" << endl << Effage <<endl;
	ofs<<"Effarea"<< endl << Effarea <<endl;
	ofs<<"EffNatAge"<< endl << EffNatAge<<endl;
	ofs<<"CatchNatAge"<< endl << CatchNatAge<<endl;

	
FUNCTION output_pin
	
	ofstream ifs("lagrangian_est.pin");

	ifs<<"# mo " << endl << 1 <<endl;
	ifs<<"# tau_c " << endl << log(.2) <<endl;
	ifs<<"# cvPos "<< endl << log(.1) <<endl;	
	//ifs<<"# maxPos "<< endl << minPos <<endl;
	ifs<<"# maxPos501 "<< endl << log(4) <<endl;
	ifs<<"# maxPos502 "<< endl << log(4) <<endl;
	ifs<<"# maxPossd1 "<< endl << log(.5) <<endl;
	ifs<<"# maxPossd2 "<< endl << log(4) <<endl;
	ifs<<"# wt "<< endl << wt <<endl;


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
	afs<<"#  tstp month area catage " << endl << obsCatchNatAge <<endl;
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
