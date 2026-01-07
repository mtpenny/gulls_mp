#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include<cstdlib>

#include "readPlanets.h"
#include "split.h"

int readPlanets(struct filekeywords *Paramfile, struct planetdata *Planets)
{

  //There is expected to be one planet file per subrun and per field with the filename format
  //of <root>.<field>.<subrun>
  
  int ncols;
  int nlist=0;

  string line;
  string ignore = 
    string("ABCDFGHIJKLMNOPQRSTUVWXYZabcdfghijklmnopqrstuvwxyz*:;=(),[]{}#");
  vector<double> data;
  vector<string> header_data;
  string subrun = to_string(Paramfile->instance);

  ifstream pf;

  string fname = string(Paramfile->planetdir) + string(Paramfile->planetroot);

  if(Paramfile->choosefield>=0)
    {
      //char field[20];
      //sprintf(field,"%d",choosefield);
      fname += to_string(Paramfile->choosefield) + string(".") + subrun;
    }
  else
    {
      cerr << "Error: field not chosen. Use the -f flag" << endl;
      exit(1);
    }

  pf.open(fname.c_str());
  if(!pf)
    {
      cerr << "ERROR READING PLANETS FILE: " << fname << endl;
      exit(1);
      return 0;
    }

  int headerSet=0;

  //read in the planets
  while(!pf.eof())
    {
      getline(pf,line);
      //remove any comments

      if(line.substr(0,1).find("#")!=line.npos) continue;
      
      if(line.find_first_of(ignore)==line.npos)
	{
	  split(line,data);

	  if(data.size()>0)
	    {
	      Planets->data.push_back(data);
	      nlist++;
	    }

	}
      else
	{
	  if(headerSet==0)
	    {
	      split(line,header_data);
	      Planets->header = header_data;
	      headerSet=1;
	    }
	}
    }

  if(Planets->header.size()==0)
    {
      for(int i=0;i<Planets->data[0].size();i++)
	{
	  Planets->header.push_back("Planet_" + to_string(i));
	}
    }

  if(Planets->header.size()!=Planets->data[0].size())
    {
      cerr << "WARNING: Planet file header has different number of columns than the data." << endl;
      for(auto ph : Planets->header) cerr << ph << " ";
      cerr << endl;
      for(auto pd : Planets->data[0]) cerr << pd << " ";
      cerr << endl;
    }

  return nlist;
}

