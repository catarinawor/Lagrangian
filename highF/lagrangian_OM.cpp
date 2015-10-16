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
	
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <lagrangian_OM.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
		ifstream ifs( "seed.txt" ); // if this file is available
		ifs>>seed; //read in the seed
		seed += 10; // add 10 to the seed
		ofstream ofs( "seed.txt" ); //put out to seed.txt
		ofs<<seed<<endl; //the new value of the seed
  syr.allocate("syr");
  nyr.allocate("nyr");
  sage.allocate("sage");
  nage.allocate("nage");
  smon.allocate("smon");
  nmon.allocate("nmon");
  sarea.allocate("sarea");
  narea.allocate("narea");
  nations.allocate("nations");
  border.allocate(1,nations-1,"border");
  Ro.allocate("Ro");
  h.allocate("h");
  m.allocate("m");
  fe.allocate("fe");
  q.allocate("q");
  sigR.allocate("sigR");
  tau_c.allocate("tau_c");
  mo.allocate("mo");
  err.allocate("err");
  wa.allocate(sage,nage,"wa");
  fa.allocate(sage,nage,"fa");
  va.allocate(sage,nage,"va");
  minPos.allocate(sage,nage,"minPos");
  maxPos501.allocate("maxPos501");
  maxPos502.allocate("maxPos502");
  maxPossd1.allocate("maxPossd1");
  maxPossd2.allocate("maxPossd2");
  cvPos.allocate("cvPos");
  TotEffyear.allocate(1,nations,syr,nyr,"TotEffyear");
  TotEffmonth.allocate(1,nations,smon,nmon,"TotEffmonth");
  effPwr.allocate(sarea,narea,"effPwr");
  eof.allocate("eof");
		
		if( eof != 999 )
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			cout<< "eof is: "<<eof<<endl;
			ad_exit(1);
		}
  age.allocate(sage,nage);
  areas.allocate(sarea,narea);
  nationareas.allocate(1,nations);
  wt.allocate(syr,nyr);
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
  indyr.allocate(1,ntstp);
  indmonth.allocate(1,ntstp);
  indnatarea.allocate(sarea,narea);
  pcat.allocate(1,nations);
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
    
       			
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  no_f.allocate("no_f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  kappa.allocate("kappa");
  #ifndef NO_AD_INITIALIZE
  kappa.initialize();
  #endif
  phie.allocate("phie");
  #ifndef NO_AD_INITIALIZE
  phie.initialize();
  #endif
  So.allocate("So");
  #ifndef NO_AD_INITIALIZE
  So.initialize();
  #endif
  Bo.allocate("Bo");
  #ifndef NO_AD_INITIALIZE
  Bo.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  tBo.allocate("tBo");
  #ifndef NO_AD_INITIALIZE
  tBo.initialize();
  #endif
  m_tsp.allocate("m_tsp");
  #ifndef NO_AD_INITIALIZE
  m_tsp.initialize();
  #endif
  lxo.allocate(sage,nage,"lxo");
  #ifndef NO_AD_INITIALIZE
    lxo.initialize();
  #endif
  za.allocate(sage,nage,"za");
  #ifndef NO_AD_INITIALIZE
    za.initialize();
  #endif
  SB.allocate(1,ntstp,"SB");
  #ifndef NO_AD_INITIALIZE
    SB.initialize();
  #endif
  varPos.allocate(sage,nage,"varPos");
  #ifndef NO_AD_INITIALIZE
    varPos.initialize();
  #endif
  maxPos.allocate(sage,nage,"maxPos");
  #ifndef NO_AD_INITIALIZE
    maxPos.initialize();
  #endif
  NationVulB.allocate(1,ntstp,1,nations,"NationVulB");
  #ifndef NO_AD_INITIALIZE
    NationVulB.initialize();
  #endif
  Nage.allocate(1,ntstp,sage,nage,"Nage");
  #ifndef NO_AD_INITIALIZE
    Nage.initialize();
  #endif
  VulB.allocate(1,ntstp,sage,nage,"VulB");
  #ifndef NO_AD_INITIALIZE
    VulB.initialize();
  #endif
  PosX.allocate(1,ntstp,sage,nage,"PosX");
  #ifndef NO_AD_INITIALIZE
    PosX.initialize();
  #endif
  Effage.allocate(1,ntstp,sage,nage,"Effage");
  #ifndef NO_AD_INITIALIZE
    Effage.initialize();
  #endif
  VBarea.allocate(1,ntstp,sarea,narea,"VBarea");
  #ifndef NO_AD_INITIALIZE
    VBarea.initialize();
  #endif
  Effarea.allocate(1,ntstp,sarea,narea,"Effarea");
  #ifndef NO_AD_INITIALIZE
    Effarea.initialize();
  #endif
  NAreaAge.allocate(1,ntstp,sarea,narea,sage,nage,"NAreaAge");
  #ifndef NO_AD_INITIALIZE
    NAreaAge.initialize();
  #endif
  CatchAreaAge.allocate(1,ntstp,sarea,narea,sage,nage,"CatchAreaAge");
  #ifndef NO_AD_INITIALIZE
    CatchAreaAge.initialize();
  #endif
  CatchNatAge.allocate(1,ntstp,1,nations,sage,nage,"CatchNatAge");
  #ifndef NO_AD_INITIALIZE
    CatchNatAge.initialize();
  #endif
  EffNatAge.allocate(1,nations,1,ntstp,sage-2,nage,"EffNatAge");
  #ifndef NO_AD_INITIALIZE
    EffNatAge.initialize();
  #endif
  obsCatchNatAge.allocate(1,tot_pcat,sage-3,nage,"obsCatchNatAge");
  #ifndef NO_AD_INITIALIZE
    obsCatchNatAge.initialize();
  #endif
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
	incidence_functions();
	initialization();
	move_grow_die();
	clean_catage();
	output_true();
	output_dat();
	output_pin();
	exit(1);
}

void model_parameters::userfunction(void)
{
  no_f =0.0;
}

dvar_vector model_parameters::cnorm(const double& x, const dvar_vector& mu, const dvar_vector& sd)
{
	dvar_vector rst(sage,nage);
	dvar_vector stx(sage,nage);
	for(int a= sage; a<= nage; a++)
	{
		stx(a) = (x-mu( a ))/sd( a );
		rst(a) = cumd_norm(stx( a ));
	}
	return(rst);
}

void model_parameters::incidence_functions(void)
{
	maxPos.initialize();
	lxo = mfexp(-m*age);
	lxo(nage) /= 1. - mfexp(-m); 
	kappa 	= 4*h/(1-h);
	phie	= lxo*fa;
	So 		= kappa/phie;
	Bo 		= kappa/So*Ro;
	beta 	= (kappa-1)/Bo;
	m_tsp 	= m/12;
	za 		= m+va*fe;
}

void model_parameters::initialization(void)
{
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
		//tmp1(rr)= pow((VBarea(1)(rr)/ (NationVulB(1)(indnatarea(rr)) + 0.0001))+ 1.0e-20,effPwr(rr));
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
		//cout<<"propVBarea "<<propVBarea<<endl;
		//cout<<"Effarea(1) "<<Effarea(1)<<endl;
		Effage(1)(a) = Effarea(1)* propVBarea;
	}
}

void model_parameters::move_grow_die(void)
{
	dvariable tB;
	for(int i=2;i<=ntstp;i++)
	{
		//if(i>12)exit(1);
		switch (indmonth(i)) {
            case 1:           	
            	Nage(i)(sage) = (So*SB(i-nmon)/(1.+beta*SB(i-nmon)))*mfexp(wt(indyr(i))*err);
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
			//tmp1(rr)= VBarea(i)(rr)/ (NationVulB(i)(indnatarea(rr)) + 1);
			//tmp1(rr)= pow((VBarea(1)(rr)/ (NationVulB(1)(indnatarea(rr)) + 0.0001))+ 1.0e-20,effPwr(rr));
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
		for(int r = sarea+1;r <= narea-1;r++)
		{
			for(int a = sage; a<=nage;a++)
			{
				CatchAreaAge(i)(r)(a) = q*Effarea(i)(r)*va(a)/(q*Effarea(i)(r)*va(a)+m_tsp)*(1-mfexp(-(q*Effarea(i)(r)*va(a)+m_tsp)))*NAreaAge(i)(r)(a);
				CatchNatAge(i)(indnatarea(r))(a)+= CatchAreaAge(i)(r)(a);
			}
		}
	}
}

void model_parameters::clean_catage(void)
{
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
}

dvar_vector model_parameters::calcmaxpos(const dvariable& tb)
{
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
}

void model_parameters::output_true(void)
{
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
}

void model_parameters::output_pin(void)
{
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
}

void model_parameters::output_dat(void)
{
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
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
	REPORT(VBarea);
	REPORT(Effage);
	REPORT(Effarea);
}

void model_parameters::final_calcs()
{
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
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
