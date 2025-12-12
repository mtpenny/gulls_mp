#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include<cstdlib>

#include "readPlanets.h"
#include "split.h"

int readPlanets(struct filekeywords *Paramfile, vector<vector<double> > *planets, string subrun, int choosefield, vector<string>* planet_header)
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

  ifstream pf;

  string fname = string(Paramfile->planetdir) + string(Paramfile->planetroot);

  if(choosefield>=0)
    {
      //char field[20];
      //sprintf(field,"%d",choosefield);
      fname += itos(choosefield) + string(".") + subrun;
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

  //read in the planets
  while(!pf.eof())
    {
      getline(pf,line);
      //remove any comments
      if(line.find_first_of(ignore)==line.npos)
	{
	  split(line,data);

	  planet_data->push_back(data);
	  nlist++;

	}
      else
	{
	  split(line,header_data);
	  (*planet_header) = header_data;
	  
	}
    }

  if(planet_header.size()==0)
    {
      for(int i=0;i<planet_data[0].size();i++)
	{
	  planet_header.push_back("Planet_" + to_string(i));
	}
    }

  if(planet_header.size()!=planet_data[0].size())
    {
      cerr << "WARNING: Planet file header has different number of columns than the data." << endl;
      for(auto i : planet_header.size()) cerr << planet_header[i] << " ";
      cerr << endl;
      for(auto i : planet_data[0].size()) cerr << planet_data[0][i] << " ";
      cerr << endl;
    }

  return nlist;
}

