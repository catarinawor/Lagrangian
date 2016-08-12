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
	init_number tau_survey; 		
	init_number err;

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

	init_vector wa(sage,nage);
	init_vector fa(sage,nage);
	init_vector va(sage,nage);
	init_vector minPos(sage,nage);						

	init_matrix TotEffyear(1,fisharea,syr,nyr);
	init_matrix TotEffmonth(1,fisharea,smon,nmon);

	init_vector effPwr(sarea,narea);

	init_vector wt(rep_yr+1,nyr);

	init_int eof;
	
	
	LOC_CALCS
		
		if( eof != 999 )
		{
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
   	//vector wt(syr,nyr);
   	//vector wtsim(syr,rep_yr);

   	vector epsilon(1,surv_nobs)

   	
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
	vector itB(rep_yr+1,nyr);
	vector totcatch(rep_yr+1,nyr);
	vector Ut(rep_yr+1,nyr);
	vector spr(syr,nyr);
	vector phie(syr,nyr);


	

	
	matrix Nage(1,ntstp,sage,nage);
	matrix Catage(1,ntstp,sage,nage);
	matrix Fatage(1,ntstp,sage,nage);
 	matrix VulB(1,ntstp,sage,nage);
 	matrix PosX(1,ntstp,sage,nage);
 	//matrix Effage(1,ntstp,sage,nage);
 	matrix VBarea(1,ntstp,sarea,narea);
 	matrix totVBnation(1,ntstp,1,nations);	
 	matrix Effarea(1,ntstp,sarea,narea);
 	matrix tot_comm_obsCatage(rep_yr+1,nyr,sage,nage);
 	matrix comm_obsCatage(rep_yr+1,nyr,sage,nage);
 	matrix surv_obsCatage(1,surv_nobs,sage,nage);
 	matrix yNage(syr,nyr,sage,nage);
 	matrix yFatage(syr,nyr,sage,nage);
 	matrix seltotal(syr,nyr,sage,nage);
 	matrix yCatchtotalage(syr,nyr,sage,nage);

 	3darray selfisharea(syr,nyr,1,fisharea,sage-2,nage);
 	3darray selnation(syr,nyr,1,nations,sage-2,nage) 	
 	//3darray surv_obsCatage(1,surv_nobs,1,fisharea,sage,nage);
 	3darray propVBarea(1,ntstp,sarea,narea,sage-2,nage);
 	3darray NAreaAge(1,ntstp,sarea,narea,sage,nage);
 	3darray CatchAreaAge(1,ntstp,sarea,narea,sage,nage);
 	3darray CatchNatAge(1,ntstp,1,fisharea,sage-2,nage);
 	3darray yCatchNatAge(syr,nyr,1,fisharea,sage,nage);
 	3darray yCatchStateAge(syr,nyr,1,nations,sage,nage);
 	//3darray EffNatAge(1,fisharea,1,ntstp,sage-2,nage);

 	matrix obsCatchNatAge(1,tot_pcat,sage-3,nage);

PRELIMINARY_CALCS_SECTION

	incidence_functions();
	initialization();
	move_grow_die();
	clean_catage();

	survey_data();
	catage_comm();
	calc_selectivity();
	calc_spr();
	
	output_true();
	output_dat();
	write_iscam_data_file();
	output_pin();
	output_datSS();
	output_pin_SS();
	
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
	
	lxo(sage)=1;
	for(int a = sage+1; a<= nage; a++){
		lxo(a) = lxo(a-1)*mfexp(-m);
	}	
	lxo(nage) /= (1. - mfexp(-m)); 

	cout<<"lxo is "<< lxo<<endl;


	kappa 	= 4*h/(1-h);
	phiE	= elem_prod(lxo,fa)*wa;

	So 		= kappa/phiE;
	Bo 		= kappa/So*Ro;
	beta 	= (kappa-1)/Bo;


	m_tsp = m/nmon;
	za 	= m_tsp+va*fe;
	//cout<<"Ok after incidence_functions"<< endl;
	
FUNCTION void calc_numbers_at_age(const int& ii, const dvariable& expwt)


		
		dvar_vector propBarea(sarea,narea);
		
		switch (indmonth(ii)) {
            case 1:           	
            														//mfexp(wt(indyr(ii))
            	Nage(ii)(sage) = (So*SB(ii-nmon)/(1.+beta*SB(ii-nmon)))*(mfexp(expwt*err));

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

            	yNage(indyr(ii))(sage,nage)= Nage(ii)(sage,nage);
            	
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
		tB(ii) = Nage(ii)*wa;
		//cout<<"Ok after calc_numbers_at_age"<< endl;

FUNCTION void calc_effarea(const int& ii)

		dvar_vector tmp1(sarea,narea);
		dvar_vector tmp2(sarea,narea);
		tmp1.initialize();
		tmp2.initialize();

		//for(int r= sarea; r<=narea; r++)
		//{
		//	totVBnation(ii,indnatarea(r)) += pow(VBarea(ii)(r);
		//}

		for(int n=1; n<=nations;n++){
       		totVBnation(ii,n) = sum(pow(VBarea(ii)(ntmp1(n),ntmp1(n+1)-1.0) +0.00001,fbeta));
       	}


		for(int rr= sarea; rr<=narea; rr++)
		{
			//tmp1(rr)= (pow(VBarea(ii)(rr)+0.0000001,fbeta)/(totVBnation(ii,indnatarea(rr))+0.01)) * effPwr(rr);

			tmp1(rr)= (pow(VBarea(ii)(rr)+0.00001,fbeta)/(totVBnation(ii,indnatarea(rr))+0.01)) * effPwr(rr);

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
		

			propVBarea(ii)(r)(sage-2) = ii;
			propVBarea(ii)(r)(sage-1) = r;
			propVBarea(ii)(r)(sage,nage) =  elem_prod(VulB(ii)(sage,nage),(cnorm(areas(r)+0.5,PosX(ii),varPos)-cnorm(areas(r)-0.5,PosX(ii),varPos)));
			

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
			
			Catage(ii)(sage,nage) += CatchAreaAge(ii)(r)(sage,nage);
			yCatchNatAge(indyr(ii))(indfisharea(r))(sage,nage) += CatchAreaAge(ii)(r)(sage,nage);			
			yCatchStateAge(indyr(ii))(indnatarea(r))(sage,nage) += CatchAreaAge(ii)(r)(sage,nage);
			yCatchtotalage(indyr(ii))(sage,nage) += CatchAreaAge(ii)(r)(sage,nage);
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
	yNage(1)(sage,nage)= Nage(1)(sage,nage);

	maxPos.initialize();
	tB(1) = Nage(1)*wa;
	calcmaxpos();
	calc_position(1);



	calc_effarea(1);


	

	calc_catage(1);


	
	

FUNCTION move_grow_die

	for(int ie=2;ie<=((rep_yr)*(nmon-smon+1));ie++)
	{ 
		calc_numbers_at_age(ie,0.0);

		maxPos.initialize();		
		calcmaxpos();

		//cout<<"maxPos "<<maxPos<<endl;
		calc_position(ie);
	
		calc_effarea(ie);
		
		calc_catage(ie);
		
		
	}


	for(int i=(rep_yr)*(nmon-smon+1)+1;i<=ntstp;i++)
	{
		
		calc_numbers_at_age(i,wt(indyr(i)));
		//cout<<"wt(indyr(i)) "<<wt(indyr(i))<<endl;
		

		maxPos.initialize();		
		calcmaxpos();

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
				pa = value((CatchNatAge(i)(n)(sage,nage))/(sum(CatchNatAge(i)(n)(sage,nage))+0.01))+0.000001;
				
				//survey option
				//pa = elem_prod(Nage(i)(sage,nage),va(sage,nage))
				
				obsCatchNatAge(p)(sage,nage) = rmvlogistic(pa,tau_c,seed+i);
				//cout<<"obsCatchNatAge(p)(sage,nage) is "<<obsCatchNatAge(p)(sage,nage)<< endl;

				p++;	
       		}	
		}
	}
	
	//cout<<"Ok after clean_catage"<< endl;


FUNCTION catage_comm

	
	for(int i=rep_yr*nmon+1;i<=ntstp;i++)
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
      
    for(int ii=rep_yr+1;ii<=nyr;ii++)
	{	
		//pa = value((tot_comm_obsCatage(ii)(sage,nage))/(sum(tot_comm_obsCatage(ii)(sage,nage))+0.01))+0.000001;			
		pa = value((tot_comm_obsCatage(ii)(sage,nage))/(sum(tot_comm_obsCatage(ii)(sage,nage))));			
		
		comm_obsCatage(ii)(sage,nage) = rmvlogistic(pa,tau_c,seed+ii);	
		totcatch(ii) = 	sum(tot_comm_obsCatage(ii)(sage,nage));	
		Ut(ii) = totcatch(ii)/tB(ii*(nmon-smon+1)-(nmon-smon));
		itB(ii) = tB(ii*(nmon-smon+1)-(nmon-smon));
    }	
	
	//cout<<"Ok after catage_comm"<< endl;	





FUNCTION survey_data

	
	survB.initialize();
	surv_obsCatage.initialize();


	for(int i=1;i<=surv_nobs;i++)
	{	
		int ind_sv;
		ind_sv = surv_yrs(i)*(nmon-smon+1)-(nmon-surv_mon);

		
		if(indmonth(ind_sv)==surv_mon)
       	{	
       		survB(i)=sum(VulB(ind_sv)(sage,nage))* mfexp(epsilon(i));

       		dvector pp(sage,nage);
       		pp.initialize();
       		for(int n=1;n<=fisharea;n++)
			{    
				pp += value(CatchNatAge(ind_sv)(n)(sage,nage)); 			
				//pp += value((CatchNatAge(ind_sv)(n)(sage,nage))/(sum(CatchNatAge(ind_sv)(n)(sage,nage))+0.01))+0.000001;
			}	
			dvector ppp(nage,sage);
       		ppp.initialize();
       		ppp = (pp)/sum(pp);
			surv_obsCatage(i)(sage,nage) = rmvlogistic(ppp,tau_survey,seed+i);
       	}	
			
    }	
    //cout<<"Ok after survey_data"<< endl;

    



	
FUNCTION dvar_vector calcmaxpos()
	
	maxPos(sage,nage) = 1./(1.+mfexp(-(age-maxPos50)/maxPossd));
	maxPos(sage,nage) *= (narea-sarea);
	maxPos(sage,nage) += sarea;			
	//cout<<"age"<<age<<endl;
	//cout<<"maxPos"<<maxPos<<endl;
	
	return(maxPos);



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

	for(i=1; i<=ntstp;i++){

		Fatage(i)(sage,nage)=elem_div(Catage(i)(sage,nage),Nage(i)(sage,nage));
		yFatage(indyr(i))(sage,nage) += Fatage(i)(sage,nage);
	}
	
	
	for(ii=syr; ii<=nyr;ii++){
			
		lz.initialize();
		lz(1) = 1;


		for(a=sage+1; a<=nage;a++){
			lz(a)= value(lz(a-1)*mfexp(-m-yFatage(ii)(a-1)));
		}
		lz(nage) /=  value(1.-mfexp(-m-yFatage(ii)(nage)));

		
			//cout<<"Fatage(ii)(nage) is "<< Fatage(ii)(nage)<<endl;
			//cout<<"lz is "<< lz<<endl;
	


		phie(ii)=elem_prod(lz,fa)*wa;
		spr(ii)=phie(ii)/phiE;


	}




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
	ofs<<"Ro " << endl << Ro <<endl;
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
	ofs<<"itB" << endl << itB <<endl;
	ofs<<"survB" << endl << survB <<endl;
	ofs<<"VulB" << endl << VulB <<endl;
	ofs<<"Nage" << endl << Nage <<endl;
	ofs<<"VBarea" << endl << VBarea <<endl;
	ofs<<"propVBarea" << endl << propVBarea <<endl;
	ofs<<"Effarea"<< endl << Effarea <<endl;
	ofs<<"comm_obsCatage"<< endl << comm_obsCatage <<endl;
	ofs<<"surv_obsCatage"<< endl << surv_obsCatage <<endl;
	ofs<<"totVBnation" << endl << totVBnation <<endl;
	ofs<<"CatchNatAge"<< endl << CatchNatAge<<endl;
	ofs<<"CatchAreaAge"<< endl << CatchAreaAge<<endl;
	ofs<<"totcatch"<< endl << totcatch<<endl;
	ofs<<"indyr"<< endl << indyr<<endl;
	ofs<<"indmonth"<< endl << indmonth<<endl;
	ofs<<"indnatarea"<< endl << indnatarea<<endl;
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

	//Generate initial values at random

	//mo
	
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
	
	//dvector guess_wt(rep_yr+1,nyr);
	//guess_wt.initialize();
	
	ofstream ifs("lagrangian_est.pin");

	ifs<<"#log_mo \n "  << log(tmp_mo) <<endl;
	//ifs<<"#log_mo \n "  << log(mo) <<endl;
	ifs<<"#cvPos \n" << log(guess_cvPos(ceil(randu(rngcvPos)*5))) <<endl;	
	//ifs<<"#cvPos \n" << log(cvPos) <<endl;	
	ifs<<"# maxPos50 \n" << log(guess_maxPos50(ceil(randu(rngmaxPos50)*9))) <<endl;
	//ifs<<"# maxPos50 \n" << log(maxPos50) <<endl;
	ifs<<"# maxPossd \n"<< log(guess_maxPossd(ceil(randu(rngmaxPossd)*7))) <<endl;
	//ifs<<"# maxPossd \n"<< log(maxPossd) <<endl;
	ifs<<"#wt \n" << wt*err <<endl;


FUNCTION output_pin_SS

	//Generate initial values at random

	//mo
	
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
	
	//dvector guess_wt(rep_yr+1,nyr);
	//guess_wt.initialize();
	
	ofstream iifs("lagrangian_SS.pin");

	iifs<<"#log_mo \n "  << log(tmp_mo) <<endl;
	iifs<<"#cvPos \n" << log(guess_cvPos(ceil(randu(rngcvPos)*5))) <<endl;	
	iifs<<"# maxPos50 \n" << log(guess_maxPos50(ceil(randu(rngmaxPos50)*9))) <<endl;
	iifs<<"# maxPossd \n"<< log(guess_maxPossd(ceil(randu(rngmaxPossd)*7))) <<endl;
	iifs<<"#Ro \n" << Ro <<endl;
	iifs<<"#h \n" << h <<endl;
	iifs<<"#avg_rec \n" << log(5) <<endl;
	iifs<<"#sigma_r \n" << 3.0 <<endl;
	iifs<<"#wt \n" << wt*err <<endl;




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
	
	
FUNCTION output_datSS

	ofstream mfs("lagrangian_SS.dat");
	mfs<<"# syr " << endl << rep_yr+1 <<endl;
	mfs<<"# nyr " << endl << nyr <<endl;
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
	mfs<<"# Total effort by country and year " << endl << TotEffyear_rep <<endl;
	mfs<<"# Total effort by country and month " << endl << TotEffmonth <<endl;
	mfs<<"# effPwr"<< endl << effPwr <<endl;
	mfs<<"# dMinP " << endl << 0.1e-15 <<endl;
	mfs<<"# tstp month area catage " << endl << obsCatchNatAge <<endl;	
	mfs<<"# nItNobs " << endl << surv_nobs <<endl;
	mfs<<"## iyr     index(it) gear area group sex log(SE) log(PE)   timing" << endl;
	for(int l=1;l<=surv_nobs;l++){
		//  	iyr     		index(it) 	log(SE)  	log(PE)   	timing
		mfs<<surv_yrs(l) <<"\t"<< survB(l) <<"\t" << 0.2<<"\t"<< 0.01<<"\t"<< surv_mon <<endl;
	}
	mfs<<"# eof " << endl << 999 <<endl;

FUNCTION write_iscam_data_file

	ofstream ufs("/Users/catarinawor/Documents/iSCAM/examples/hakelag/DATA/hakelag.dat");
	ufs<<"# DATA FILE FOR iSCAM  " << endl;
	ufs<<"## ------------------------------------------------------------------------- ##"<< endl;
	ufs<<"## MODEL DIMENSIONS  " << endl;
	ufs<<"## ------------------------------------------------------------------------- ##"<< endl;
	ufs<< 1  	<< "\t" << "# -number of areas            (narea) "  <<endl;
	ufs<< 1 	<< "\t" << "# -number of groups or stocks (ngroup)"  <<endl;
	ufs<< 1  	<< "\t" << "# -number of sexes            (nsex)"  <<endl;
	ufs<< rep_yr+1 	<< "\t" << "# -first year of data         (syr)"  <<endl;
	ufs<< nyr 	<< "\t" << "# -last year of data          (nyr)"  <<endl;
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
	ufs<<  nyr-rep_yr <<endl;
	ufs<<"## Year gear area group sex type value" << endl;
	for(int i=rep_yr+1;i<=nyr;i++){
		//   year 		gear 		area 	  group       sex       type   		value             se
		ufs<<  i <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 0 <<"\t"<< 1 <<"\t"<<totcatch(i)<<"\t"<<0.1<<endl;
	}
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<<"## ABUNDANCE INDICES -A RAGGED ARRAY: (1,nit,1,nit_nobs,1,5)                 ##" << endl;
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<< 1 	<< "\t" << "# Number of abundance series         int(nit)"  <<endl;
	ufs<< surv_nobs 	<< "\t" << "# Number of observations in series   ivector(nit_nobs(1,nit))"  <<endl;
	ufs<< 2 	<< "\t" << "#Survey type (see key below)        ivector(survey_type(1,nit))"  <<endl;
	ufs<<"## 1 = survey is proportional to vulnerable numbers" << endl;
	ufs<<"## 2 = survey is proportional to vulnerable biomass" << endl;
	ufs<<"## 3 = survey is proportional to spawning biomass (e.g., a spawn survey)" << endl;
	ufs<<"## iyr     index(it) gear area group sex log(SE) log(PE)   timing" << endl;
	for(int l=1;l<=surv_nobs;l++){
		//  	iyr     		index(it) 		  gear       area 		group 		sex 	log(SE)  		log(PE)   	timing
		ufs<<surv_yrs(l) <<"\t"<< survB(l) <<"\t"<< 2 <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 0 <<"\t"<< 0.2<<"\t"<< 0.01<<"\t"<< 0.5 <<endl;
	}
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<<"## AGE COMPOSITION DATA (ROW YEAR, COL=AGE) Ragged object                    ##" << endl;
	ufs<<"## ------------------------------------------------------------------------- ##" << endl;
	ufs<< 2 	<< "\t" << "# Number of gears with age-comps int(na_gears)"  <<endl;
	ufs<< nyr-rep_yr << "\t" << surv_nobs	<< "\t" << "# Number of rows in the matrix   ivector(na_gears)"  <<endl;
	ufs<< sage << "\t" << sage	<< "\t" << "## ivector(na_gears) of youngest age-class"  <<endl;
	ufs<< nage << "\t" << nage	<< "\t" << "## ivector(na_gears) of oldest age-class + group"  <<endl;
	ufs<< 10 << "\t" << 10	<< "\t" << "## effective sample size for multinomial"  <<endl;
	ufs<< 1	<< "\t" << 1 << "\t" << "## Age composition flag (1==age comps, 0==length comps)"  <<endl;
	ufs<<"## year gear area group sex age_err | data columns (numbers or proportions)" << endl;
	for(int i=rep_yr+1;i<=nyr;i++){
		//  year 	  gear 		area 		group 		sex 	age_err | data columns (numbers or proportions)
		ufs<<i <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 1 <<"\t"<< 0 <<"\t"<< 1 <<"\t"<< comm_obsCatage(i)(sage,nage) <<endl;
	}
	for(int l=1;l<=surv_nobs;l++){	
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

