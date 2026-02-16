#include<string>
#include "getPlanetvals.h"
#include "buildEvent.h"
#include "croin.h"

using namespace std;

//void getPlanetvals(struct event* Event, struct obsfilekeywords World[], struct filekeywords *Paramfile, struct slcat* Sources, struct slcat* Lenses, vector<struct pcat>* Planets)
void getPlanetvals(struct event* Event, struct obsfilekeywords World[], struct filekeywords *Paramfile, struct slcat* Sources, struct slcat* Lenses, struct planetdata* Planets)
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

  Event->p_mass.clear();
  Event->p_a.clear();
  Event->p_e.clear();
  Event->p_I.clear();
  Event->p_L0.clear();
  Event->p_w.clear();
  Event->p_O.clear();
  Event->p_dL.clear();
  Event->p_orbtype.clear();
  Event->p_period.clear();
  Event->p_q.clear();
  Event->p_s0.clear();

  if(int(Planets->header.size())%7 != 0)
    {
      cerr << "Error: Bad number of columns in the planet file, should be 7 per object." << endl;
      exit(1);
    }
  vector<double> pd = Planets->data[Event->id];
  if(int(pd.size())%7 != 0)
    {
      cerr << "Error: Bad number of columns in the planet file on line " << Event->id << ", should be 7 per object." << endl;
      for(auto pdi : pd) cerr << pdi << " ";
      cerr << endl;
      exit(1);
    }

  int nplanets = int(pd.size()/7);
  Event->ncatalog_planets=nplanets;
  int orbtype;
  int skipped_planets=0;
  double first_inc=0;

  for(int i=0;i<nplanets;i++)
    {
      double mass;
      for(int j=0;j<7;j++)
	{
	  int col = 7*i+j;

	  if(Planets->header[col].rfind("Mass",0)==0) mass=pd[col];

	  if(mass>0.0)
	    {
	      if(Planets->header[col].rfind("Mass",0)==0)
		Event->p_mass.push_back(pd[col]);
	      if(Planets->header[col].rfind("SemimajorAxis",0)==0)
		Event->p_a.push_back(pd[col]);
	      if(Planets->header[col].rfind("Eccentricity",0)==0)
		Event->p_e.push_back(pd[col]);
	      if(Planets->header[col].rfind("Inclination",0)==0)
		{
		  double inc=pd[col];
		  if(inc>900) //If inclination is relative to the binary orbit, then it should have 1000 degrees added to it
		    {
		      if(Event->lcompanions.size()>0)
			{
			  inc = Event->lcomp_I[0] + (inc-1000.0);
			}
		      else
			{
			  if(i==0)
			    {
			      double rnd = ran2(Paramfile->seed);
			      inc = (180*(rnd<0.5?acos(2*rnd):-acos(2-2*rnd))/PI);
			      first_inc = inc;
			    }
			  else
			    {
			      inc=first_inc;
			    }
			}
		    }
		  Event->p_I.push_back(inc);
		}
	      if(Planets->header[col].rfind("LongitudePerihelion",0)==0)
		Event->p_w.push_back(pd[col]);
	      if(Planets->header[col].rfind("LongitudeAscNode",0)==0)
		Event->p_O.push_back(pd[col]);
	      if(Planets->header[col].rfind("OrbitType",0)==0)
		{
		  Event->p_orbtype.push_back(int(pd[col]));
		  orbtype = int(pd[col]);
		}
	    }
	   
	}
      if(mass>0.0)
	{
	  Event->p_L0.push_back(360.0*ran2(Paramfile->seed));
	  Event->p_s0.push_back(0.0);
	  

	  if(orbtype==0)
	    {
	      //Figure out the orbit type
	      if(Event->nlens==1)
		{
		  orbtype=1;
		  Event->p_orbtype.back() = 1;
		}
	      else
		{
		  //binary star
		  if(Event->p_a.back()<Event->lcomp_a.back())
		    {
		      orbtype=1;
		      Event->p_orbtype.back() = 1;
		    }
		  else
		    {
		      //cicumbinary
		      orbtype=2;
		      Event->p_orbtype.back() = 2;
		    }
		}
	    
	    }
	  if(orbtype==3)
	    {
	      //moon orbits a planet
	      double M1 = Lenses->data[Event->lens][Lenses->datadict["Mass"]];
	      double q = Event->p_mass.back()/M1;
	      Event->p_q.push_back(q);
	      Event->p_period.push_back(0.0);
	      Event->p_dL.push_back(360.0);
	      if(Paramfile->verbosity>=1) cout << "Planet orbtype=" << orbtype << " a=" << Event->p_a.back() << " mass=" << Event->p_mass.back() << endl;
	    }
	  if(orbtype==2)
	    {
	      //circumbinary planets
	      if(Event->lcompanions.size()==0)
		{
		  //not a circumbinary with just one star
		  orbtype=1;
		  Event->p_orbtype.back()=1;
		}
	      else
		{
		  double M1 = Lenses->data[Event->lens][Lenses->datadict["Mass"]];
		  double q = Event->p_mass.back()/M1;
		  Event->p_q.push_back(q);
		  Event->p_period.push_back(0.0);
		  Event->p_dL.push_back(360.0);
		  if(Paramfile->verbosity>=1) cout << "Planet orbtype=" << orbtype << " a=" << Event->p_a.back() << " mass=" << Event->p_mass.back() << endl;		  
		}
	    }
	  if(orbtype==1)
	    {
	      //planet orbits star 1
	      double M1 = Lenses->data[Event->lens][Lenses->datadict["Mass"]];
	      double q = Event->p_mass.back()/M1;
	      Event->p_q.push_back(q);
	      Event->p_period.push_back(0.0);
	      Event->p_dL.push_back(360.0);
	      if(Paramfile->verbosity>=1) cout << "Planet orbtype=" << orbtype << " a=" << Event->p_a.back() << " mass=" << Event->p_mass.back() << endl;		      
	    }
	  
	  if(orbtype>3)
	    {
	      cerr << "Planet orbit codes above 3 aren't implemented." << endl;
	      exit(1);
	    }

	} //end if mass >0
      else
	{
	  if(Paramfile->verbosity>=1) cout << "Skipping planet " << i << " with zero mass" << endl;
	  skipped_planets++;
	}
      
      //Event->lcomp_alpha.push_back(360.0*ran2(idum));
      //Event->lcomp_phase.push_back(360.0*ran2(idum));
    }

  Paramfile->parameterization=0;
  Event->tref=Event->t0;
  Event->nplanets = nplanets - skipped_planets;

  Event->nlens += Event->nplanets;
  if(Paramfile->verbosity>=2)
    {
      cout << "At end of planet setup, nlens = " << Event->nlens << ", nsrc = " << Event->nsrc << endl;
    }
  
}
