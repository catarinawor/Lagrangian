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
		seed += 1; // add 10 to the seed
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
	init_number err;

	//=================================
	//parameters in estimation model 
	//
	//=================================
	init_number mo; 		
	init_number maxPos50;					
	init_number maxPossd;			
	init_number cvPos; 

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
			cout<<"TotEffmonth "<<TotEffmonth<<endl;
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

   	
   		LOC_CALCS		
			
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

       			//cout<<"Ok after indfisharea calcs"<< endl;

       			//calc for indnatarea
       			ntmp1(1) = sarea;
       			
       			for(int n=1;n<=nations;n++)
       			{
       				ntmp1(n+1)= ntmp1(n)+nationareas(n);	
       				indnatarea(ntmp1(n),ntmp1(n+1)-1)=n;
       			}

       			ntmp1(nations+1) = narea;
       			indnatarea(narea) = nations;

       			//cout<<"Ok after indnatarea calcs"<< endl;

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

	
	matrix Nage(1,ntstp,sage,nage);
 	matrix VulB(1,ntstp,sage,nage);
 	matrix PosX(1,ntstp,sage,nage);
 	//matrix Effage(1,ntstp,sage,nage);
 	matrix VBarea(1,ntstp,sarea,narea);
 	matrix totVBnation(1,ntstp,1,nations);
 	
 	//matrix propVBarea(1,ntstp,sarea,narea);
 	matrix Effarea(1,ntstp,sarea,narea);
 
 	3darray NAreaAge(1,ntstp,sarea,narea,sage,nage);
 	3darray CatchAreaAge(1,ntstp,sarea,narea,sage,nage);
 	3darray CatchNatAge(1,ntstp,1,fisharea,sage-2,nage);
 	//3darray EffNatAge(1,fisharea,1,ntstp,sage-2,nage);

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


	m_tsp = m/nmon;
	za 	= m_tsp+va*fe;
	//cout<<"Ok after incidence_functions"<< endl;
	
FUNCTION void calc_numbers_at_age(const int& ii)


		
		dvar_vector propBarea(sarea,narea);
		
		switch (indmonth(ii)) {
            case 1:           	
            	
            	Nage(ii)(sage) = (So*SB(ii-nmon)/(1.+beta*SB(ii-nmon)))*(mfexp(wt(indyr(ii))*err));

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
            					(Nage(ii-1)(nage-1)*(1-sum(propBarea))*mfexp(-m_tsp))/(1.-mfexp(-m_tsp));
            	
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
            					 Nage(ii-1)(a)*(1-sum(propBarea))*mfexp(-(m_tsp));
            	}

            	break;
        }
		
		VulB(ii) = elem_prod(elem_prod(Nage(ii),va),wa);
		SB(ii) = elem_prod(Nage(ii),fa)*wa/2;
		tB(ii) = Nage(ii)*wa;
		//cout<<"Ok after calc_numbers_at_age"<< endl;

FUNCTION void calc_effarea(const int& ii)

		dvar_vector tmp1(sarea,narea);
		dvar_vector tmp2(sarea,narea);

		//for(int r= sarea; r<=narea; r++)
		//{
		//	totVBnation(ii,indnatarea(r)) += VBarea(ii)(r);
		//}

		for(int n=1; n<=nations;n++){
       		totVBnation(ii,n) = sum(VBarea(ii)(ntmp1(n),ntmp1(n+1)-1));
       	}


		for(int rr= sarea; rr<=narea; rr++)
		{
			tmp1(rr)= pow((VBarea(ii)(rr)/(totVBnation(ii,indnatarea(rr))+0.01))+0.0000001,fbeta) * effPwr(rr);

			tmp2(rr) = tmp1(rr)*TotEffyear(indfisharea(rr))(indyr(ii));
			Effarea(ii)(rr) = tmp2(rr)*TotEffmonth(indfisharea(rr))(indmonth(ii));

		}
		

FUNCTION void calc_position(const int& ii)

		varPos = maxPos*cvPos;

		PosX(ii) = minPos + (maxPos - minPos) * (0.5+0.5*sin(indmonth(ii)*PI/6 - mo*PI/6)); 

		//VBarea(ii,sarea) = VulB(ii)* (cnorm(areas(sarea)+0.5,PosX(ii),varPos) - cnorm(areas(sarea)-0.5,PosX(ii),varPos));


		for(int r = sarea;r <= narea;r++)
		{
			VBarea(ii)(r) = VulB(ii)* (cnorm(areas(r)+0.5,PosX(ii),varPos)-cnorm(areas(r)-0.5,PosX(ii),varPos));
			NAreaAge(ii)(r) = elem_prod(Nage(ii)(sage,nage),(cnorm(areas(r)+0.5,PosX(ii),varPos)-cnorm(areas(r)-0.5,PosX(ii),varPos)));
		}	


		

FUNCTION void calc_catage(const int& ii)

	for(int r = sarea;r <= narea;r++)
		{
			CatchNatAge(ii)(indfisharea(r))(sage-2) = ii;
			CatchNatAge(ii)(indfisharea(r))(sage-1) = indfisharea(r);
			for(int a = sage; a<=nage;a++)
			{
				CatchAreaAge(ii)(r)(a) = q*Effarea(ii)(r)*va(a)/(q*Effarea(ii)(r)*va(a)+m_tsp)*(1-mfexp(-(q*Effarea(ii)(r)*va(a)+m_tsp)))*NAreaAge(ii)(r)(a);
				CatchNatAge(ii)(indfisharea(r))(a)+= CatchAreaAge(ii)(r)(a);
			}

		}
		//cout<<"Ok after calc_catage"<< endl;


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

	VulB(1) = elem_prod(elem_prod(Nage(1),va),wa);
	SB(1) = elem_prod(Nage(1),fa)*wa/2;

	maxPos.initialize();
	tB(1) = Nage(1)*wa;
	calcmaxpos(tB(1));
	calc_position(1);



	calc_effarea(1);


	

	calc_catage(1);


	
	

FUNCTION move_grow_die


	for(int i=2;i<=ntstp;i++)
	{
		
		calc_numbers_at_age(i);

		maxPos.initialize();		
		calcmaxpos(tB(i));

		//cout<<"maxPos "<<maxPos<<endl;
		calc_position(i);
	
		calc_effarea(i);
		

		calc_catage(i);
		
	}

FUNCTION clean_catage
	
	int p;
       
	p=1;
	for(int i=rep_yr*nmon+1;i<=ntstp;i++)
	{
		for(int n=1;n<=fisharea;n++)
		{
							
			if(TotEffmonth(n)(indmonth(i))>0)
       		{
       			
       			dvector pa(nage,sage);
       			pa.initialize();

       			obsCatchNatAge(p)(sage-3) = i;
       			obsCatchNatAge(p)(sage-2) = indmonth(i);
				obsCatchNatAge(p)(sage-1) = n;
       			
       			//non-survey option
				pa = value((CatchNatAge(i)(n)(sage,nage))/sum(CatchNatAge(i)(n)(sage,nage)));
				
				//survey option
				//pa = elem_prod(Nage(i)(sage,nage),va(sage,nage))
				
				obsCatchNatAge(p)(sage,nage) = rmvlogistic(pa,tau_c,seed+i);
				
				p++;	
       		}	
		}
	}
	
FUNCTION dvar_vector calcmaxpos(const dvariable& tb)
	
	maxPos(sage,nage) = 1./(1.+mfexp(-(age-maxPos50)/maxPossd));
	maxPos(sage,nage) *= (narea-sarea);
	maxPos(sage,nage) += sarea;			
	return(maxPos);


FUNCTION output_true
	
	ofstream ofs("lagrangian_OM.rep");

	ofs<<"mo" << endl << mo <<endl;
	ofs<<"tau_c" << endl << tau_c<<endl;
	ofs<<"maxPos50" << endl << maxPos50 <<endl;
	ofs<<"maxPossd" << endl << maxPossd <<endl;
	ofs<<"cvPos" << endl << cvPos <<endl;
	ofs<<"seed"<< endl << seed <<endl;
	ofs<<"syr" << endl << syr <<endl;
	ofs<<"nyr" << endl << nyr <<endl;
	ofs<<"sage" << endl << sage <<endl;
	ofs<<"nage" << endl << nage <<endl;
	ofs<<"smon" << endl << smon <<endl;
	ofs<<"nmon" << endl << nmon <<endl;
	ofs<<"sarea" << endl << sarea <<endl;
	ofs<<"narea" << endl << narea <<endl;
	ofs<<"fisharea" << endl << fisharea <<endl;
	ofs<<"nations" << endl << nations <<endl;
	ofs<<"maxPos" << endl << maxPos <<endl;
	ofs<<"minPos" << endl << minPos <<endl;
	ofs<<"varPos" << endl << varPos <<endl;
	ofs<<"PosX" << endl << PosX <<endl;	
	ofs<<"SB" << endl << SB <<endl;
	ofs<<"VulB" << endl << VulB <<endl;
	ofs<<"Nage" << endl << Nage <<endl;
	ofs<<"VBarea" << endl << VBarea <<endl;
	//ofs<<"Effage" << endl << Effage <<endl;
	ofs<<"Effarea"<< endl << Effarea <<endl;
	//ofs<<"EffNatAge"<< endl << EffNatAge<<endl;
	ofs<<"totVBnation" << endl << totVBnation <<endl;
	ofs<<"CatchNatAge"<< endl << CatchNatAge<<endl;
	ofs<<"CatchAreaAge"<< endl << CatchAreaAge<<endl;
	ofs<<"indyr"<< endl << indyr<<endl;
	ofs<<"indmonth"<< endl << indmonth<<endl;
	ofs<<"indnatarea"<< endl << indnatarea<<endl;



	
FUNCTION output_pin

	//Generate initial values at random

	//mo
	
	random_number_generator rngmo(seed);
	random_number_generator rngcvPos(seed);
	random_number_generator rngmaxPos50(seed);
	random_number_generator rngmaxPossd(seed);
	
	double tmp_mo;
	double tmp_cvPos;
	double tmp_maxPos50;
	double tmp_maxPossd;
	
	dvector guess_cvPos(1,6);
	dvector guess_maxPos50(1,10);
	dvector guess_maxPossd(1,8);

	
	guess_cvPos.fill_seqadd(0.05,0.05);
	guess_maxPos50.fill_seqadd(3,0.5);
	guess_maxPossd.fill_seqadd(0.5,0.5);


	

	tmp_mo 		= rand() % 6 + 1;
	tmp_cvPos	= rand() % 6 + 1;
	tmp_maxPos50= rand() % 10 + 1;
	tmp_maxPossd= rand() % 8 + 1;

	

	//cout<<tmp_mo<<endl;
	//cout<<guess_cvPos(tmp_cvPos)<<endl;
	//cout<<guess_maxPos50(tmp_maxPos50)<<endl;
	//cout<<guess_maxPossd(tmp_maxPossd)<<endl;


	
	ofstream ifs("lagrangian_est.pin");

	ifs<<"#log_mo \n "  << log(tmp_mo) <<endl;
	//ifs<<"#log_mo \n "  << log(mo) <<endl;
	ifs<<"#cvPos \n" << log(guess_cvPos(tmp_cvPos)) <<endl;	
	//ifs<<"#cvPos \n" << log(cvPos) <<endl;	
	//ifs<<"#maxPos "<< endl << minPos <<endl;
	ifs<<"# maxPos50 \n" << log(guess_maxPos50(tmp_maxPos50)) <<endl;
	//ifs<<"# maxPos50 \n" << log(maxPos50) <<endl;
	//ifs<<"#maxPos502 "<< endl << log(4) <<endl;
	ifs<<"# maxPossd \n"<< log(guess_maxPossd(tmp_maxPossd)) <<endl;
	//ifs<<"# maxPossd \n"<< log(maxPossd) <<endl;
	//ifs<<"#maxPossd2 "<< endl << log(4) <<endl;
	ifs<<"#wt \n" << wt(rep_yr+1,nyr)*err <<endl;

	


FUNCTION output_dat

	ofstream afs("lagrangian_est.dat");
	afs<<"# syr " << endl << rep_yr+1 <<endl;
	afs<<"# nyr " << endl << nyr <<endl;
	afs<<"# sage " << endl << sage <<endl;
	afs<<"# nage " << endl << nage <<endl;
	afs<<"# smon " << endl << smon <<endl;
	afs<<"# nmon " << endl << nmon <<endl;
	afs<<"# sarea " << endl << sarea <<endl;
	afs<<"# narea " << endl << narea <<endl;
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
	afs<<"# dMinP " << endl << 0.1e-15 <<endl;
	afs<<"# tstp month area catage " << endl << obsCatchNatAge <<endl;	
	afs<<"# eof " << endl << 999 <<endl;
	

REPORT_SECTION

	REPORT(VBarea);
	//REPORT(Effage);
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

