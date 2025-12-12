#include "getPlanetvals.h"
#include "buildEvent.h"
#include "croin.h"

//void getPlanetvals(struct event* Event, struct obsfilekeywords World[], struct filekeywords *Paramfile, struct slcat* Sources, struct slcat* Lenses, vector<struct pcat>* Planets)
void getPlanetvals(struct event* Event, struct obsfilekeywords World[], struct filekeywords *Paramfile, struct slcat* Sources, struct slcat* Lenses, vector<vector<double> >* planet_data, vector<string>* planet_string)
{

  //extract and calculate the planet parameters

  int sdx = Event->id;
  int ln = Event->lens;
  //Event->params.resize(NPLANETINPUT+NPLANETDE);
  //Event->paramsHeader.resize(NPLANETINPUT+NPLANETDERIV);
  //Event->params.resize(planet_data[0].size() + 3);

  //Event->paramsHeader[PMASS] = string("mass");
  //Event->paramsHeader[AA] = string("semimajoraxis");
  //Event->paramsHeader[PHASE] = string("orbphase");
  //Event->paramsHeader[INC] = string("inclination");
  //Event->paramsHeader[QQ] = string("q");
  //Event->paramsHeader[SS] = string("s");
  //Event->paramsHeader[TT] = string("period");

  //For each planet there will be 6 parameters that define its orbit, etc.
  //1 mass
  //2 semimajor axis
  //3 eccentricity
  //4 inclination
  //5 longitude of perihelion
  //6 longitude of ascending node
  //Will also need dL, L0 for orbit

  //More planets/moons will come in additional columns
  //Each planet will have an extra column indicating what it orbits/can orbit
  //0 - no preference, orbits the chosen lens star or the binary depending on orbit size
  //1 - planet orbits the chosen lens star
  //2 - planet orbits the binary if its binary
  //3 - moon orbits the first planet
  //and not yet implemented
  //3 - no preference, orbits the chosen source star or the binary depending on orbit size
  //4 - orbits the chosen source star

  //so 7 columns per object

  if(int(planet_string.size())%7 != 0)
    {
      cerr << "Error: Bad number of columns in the planet file, should be 7 per object." << endl;
      exit(1);
    }
  vector<double> pd = (*planet_data)[Event->id];
  if(int(pd.size())%7 != 0)
    {
      cerr << "Error: Bad number of columns in the planet file on line " << Event->id << ", should be 7 per object." << endl;
      for(auto pdi : pd) cerr << pdi << " ";
      cerr << endl;
      exit(1);
    }

  int nplanets = int(pd.size()/7);
  int orbtype;
  Event->nplanets = nplanets;

  for(int i=0;i<nplanets;i++)
    {
      for(int j=0;j<7;j++)
	{
	  int col = 7*i+j;
	  if(planet_string[col].rfind("Mass",0)==0)
	    Event->p_mass.push_back(stod(pd[col]));
	  if(planet_string[col].rfind("SemimajorAxis",0)==0)
	    Event->p_a.push_back(stod(pd[col]));
	  if(planet_string[col].rfind("Eccentricity",0)==0)
	    Event->p_e.push_back(stod(pd[col]));
	  if(planet_string[col].rfind("Inclination",0)==0)
	    Event->p_I.push_back(stod(pd[col]));
	  if(planet_string[col].rfind("LongitudePerihelion",0)==0)
	    Event->p_w.push_back(stod(pd[col]));
	  if(planet_string[col].rfind("LongitudeAscNode",0)==0)
	    Event->p_O.push_back(stod(pd[col]));
	  if(planet_string[col].rfind("OrbitType",0)==0)
	    {
	      Event->p_orbtype.push_back(stod(pd[col]));
	      orbtype = stoi(pd[col]);
	    }
	}
      Event->p_L0.push_back(360.0*ran2(idum));

      if(orbcode==0)
	{
	  //Figure out the orbit type
	  if(Event->nlens==1)
	    {
	      orbcode=1;
	      Event->p_orbcode.last() = 1;
	    }
	  else
	    {
	      //binary star
	      if(Event->p_a.last()<Event->lcomp_a)
		{
		  orbcode=1;
		  Event->p_orbcode.last() = 1;
		}
	      else
		{
		  //cicumbinary
		  orbcode=2;
		  Event->p_orbcode.last() = 2;
		}
	    }
	}
      if(orbcode==3)
	{
	  //moon orbits a planet
	  double M1 = Event->p_mass[0];
	  double totmass = M1 + Event->p_mass.last();
	  double q = Event->p_mass.last()/M1;
	  Event->p_q.push_back(q);
	  Event->qsum += q*M1;
	  double acomb = Event->p_a.last() * (1+Event->p_q.last());
	  double period = DAYINYR * sqrt(cube(acomb)/totmass);
	  Event->p_period.push_back(period);
	  Event->p_dL.push_back(2*PI/period);
	}
      if(orbcode==2)
	{
	  //circumbinary planets
	  if(Event->lcompanions.size()==0)
	    {
	      //not a circumbinary with just one star
	      orbcode=1;
	      Event->p_orbcode.last()=1;
	    }
	  else
	    {
	      double M1 = Lenses->data[Event->lens][Lenses->datadict["Mass"]]
		+ Lenses->data[Event->lcompanions[0]][Lenses->datadict["Mass"]];
	      double totmass = M1 + Event->p_mass.last();
	      double q = Event->p_mass.last()/M1;
	      Event->p_q.push_back(q);
	      Event->qsum += q;
	      double acomb = Event->p_a.last() * (1+Event->p_q.last());
	      double period = DAYINYR * sqrt(cube(acomb)/totmass);
	      Event->p_period.push_back(period);
	      Event->p_dL.push_back(2*PI/period);
	    }
	}
      if(orbcode==1)
	{
	  //planet orbits star 1
	  double M1 = Lenses->data[Event->lens][Lenses->datadict["Mass"]];
	  double totmass = M1 + Event->p_mass.last();
	  double q = Event->p_mass.last()/M1;
	  Event->p_q.push_back(q);
	  Event->qsum += q;
	  double acomb = Event->p_a.last() * (1+Event->p_q.last());
	  double period = DAYINYR * sqrt(cube(acomb)/totmass);
	  Event->p_period.push_back(period);
	  Event->p_dL.push_back(2*PI/period);
	}

      if(orbcode>3)
	{
	  cerr << "Planet orbit codes above 3 aren't implemented." << endl;
	  exit(1);
	}
      
      //Event->lcomp_alpha.push_back(360.0*ran2(idum));
      //Event->lcomp_phase.push_back(360.0*ran2(idum));
    }

  Paramfile->parameterization=0;
  Paramfile->tref=Event->t0;  

  

}
