#include<vector>
#include<iostream>
#include<fstream>
#include<cmath>
#include<fitsio.h>
#include<ctime>
#include<cstdlib>
#include<string>
#include<sstream>
// #include<algorithm>   //added initially for min/max elements in course vectors
// #include<string>   //added for double to string conversion in file naming

#include "psf.h"
#include "image.h"
#include "constants.h"
#include "split.h"
#include "integerPowers.h"

// Variables for use with course data
// Easier to define globally and read for each iteration of images, even if it's not
// the most efficient
vector<double> course_year, course_x, course_y, course_theta;
int bounds_xmin, bounds_xmax, bounds_ymin, bounds_ymax;
bool bounds_set = false;
double origin_year, origin_x, origin_y, origin_theta;


double fs(double Fin, double Fout)
{
  double Fsource = Fin-Fout;
  return Fsource/Fin;
}


double magnification(double t, double t0, double tE, double u0)
{
  double u = qAdd(u0,((t-t0)/tE));
  return (sqr(u)+2)/u/qAdd(u,2);
}

void load_course(string file_name)
{
  // clear course data
  course_year.clear();
  course_x.clear();
  course_y.clear();
  course_theta.clear();

  //load dithering course file
  ifstream coursein(file_name.c_str());
  if(!coursein)
  {
    cout << "Could not open dithering course file: " << file_name << endl;
    exit(1);
  }

  while(!coursein.eof())
  {
    string line;
    vector<string> data;

    getline(coursein,line);
    split(line,data);

    if(int(data.size())>=4)
    {
      course_year.push_back(atof(data[0].c_str())/365.25);
      course_x.push_back(atof(data[1].c_str()));
      course_y.push_back(atof(data[2].c_str()));
      course_theta.push_back(atof(data[3].c_str())*pi/180.0);
    }
  }

  coursein.close();

  //normalizes course_x and course_y so that first waypoint is treated as 0,0
  origin_year = course_year[0];
  origin_x = course_x[0];
  origin_y = course_y[0];
  origin_theta = course_theta[0];

  for(int i=0;i<int(course_x.size());i++)
  {
    course_year[i] -= origin_year;
    course_x[i] -= origin_x;
    course_y[i] -= origin_y;
    course_theta[i] -= origin_theta;
  }

  int course_xmin = int(floor(*min_element(course_x.begin(),course_x.end())));
  int course_xmax = int(ceil(*max_element(course_x.begin(),course_x.end())));
  int course_ymin = int(floor(*min_element(course_y.begin(),course_y.end())));
  int course_ymax = int(ceil(*max_element(course_y.begin(),course_y.end())));

  if(course_xmin<bounds_xmin)
  {
    bounds_xmin = course_xmin;
  }

  if(course_xmax>bounds_xmax)
  {
    bounds_xmax = course_xmax;
  }

  if(course_ymin<bounds_ymin)
  {
    bounds_ymin = course_ymin;
  }

  if(course_ymax>bounds_ymax)
  {
    bounds_ymax = course_ymax;
  }
}

int main(int argc, char* argv[])
{
//////////////////////////////////////////////////////////////////////////////
//                               USER OPTIONS
//

  //const int maximsize = 1024;            //maximum image size without 
                                         //receiving a warning requiring
                                         //user interaction
  //const bool output_image = true;        //output simulated images
  const double large_psf_mag = 15.0;     //magnitude where large psf is used
                                         //missed flux still gets distributed, 
                                         //but uniformly
  //const int subpix_method = 0;           //Method by which sub-pixel systematic
                                         //is computed (see below for options)
  //const bool zero_systematic = true;     //Do not include a systematic 
                                         //component in the raw error
                                         //can then add either const or 
                                         //calculated systematic

  //Subpixel computation:
  // 0 - the fractional rms fluctuation due the total flux within a slightly 
  //     expanded aperture. (USE THIS)
  // 1 - the fractional rms fluctuation due to the total flux of all stars 
  //     that lie within the same expanded aperture (DO NOT USE THIS)
  //
  //  i.e. method 1 gives a larger error if there are bright stars at the edge 
  //       of the aperture. Method 1 requires significantly more computation
  //       as for each photometered star it must look through all other stars
  //       to see if they are within the aperture. Method 0 just counts up the
  //       total flux in the expanded aperture
                                 
//  
//                               USER OPTIONS
/////////////////////////////////////////////////////////////////////////////

  if(argc<8 || (argc>=8 && argc%2==1))
  {
    cerr << "\nUsage:\n./transitfield <detectorlist> <mlevents> <nfilters> <nfiltersml> <field> <area> {<field> <area> {...} } <output filename root>\n" << endl;
    cerr << "argc = " << argc << endl;
    exit(1);
  }

  string detectorList = string(argv[1]);      //detector parameter file
  string outroot = string(argv[argc-1]);  //root of the output filename
  string mleventsfile = string(argv[2]);      //File containing microlensing events - follows same structure as gulls files
  int nfilters = atoi(argv[3]);      //number of filters in the galactic model files
    int nfiltersml = atoi(argv[4]);      //number of filters in the galactic model files
   
  //Output consists of:
  ////                   <outroot>.txt -the list of stars and photometry
  ////                   <outroot>_true.fits -true starfield image
  ////                   <outroot>_image.fits -single simulated starfield image
  ////                   <outroot>_stack.fits -stacked sim starfield image

  vector<string> detectorNames;
  vector<double> texp;
  vector<int> nstack;
  vector<int> xpix,ypix;
  vector<int> filter;
  vector<int> mlfilter;
  vector<string> dither;

  double meanl=0.0;
  double meanb=0.0;
  
  //load the detectors

  string line;
  vector<string> data;

  ifstream detin(detectorList.c_str());
  if(!detin)
  {
    cerr << "Could not open the list of detectors: " << detectorList << endl;
    exit(1);
  }

  //double largestx, largesty;
  int largestfilter = 0;

  while(!detin.eof())
  {
    getline(detin,line);
    split(line,data);

    if(int(data.size())>=7)
    {
      detectorNames.push_back(data[0]);
      texp.push_back(atof(data[1].c_str()));
      nstack.push_back(atoi(data[2].c_str()));
      xpix.push_back(atoi(data[3].c_str()));
      ypix.push_back(atoi(data[4].c_str()));
      filter.push_back(atoi(data[5].c_str()));
      mlfilter.push_back(atoi(data[6].c_str()));
      dither.push_back(data[7]);

      cout << "Detector parsing" << endl;
      cout << detectorNames.back() << endl;
      cout << texp.back() << endl;
      cout << nstack.back() << endl;
      cout << xpix.back() << " " << ypix.back() << endl;
      cout << filter.back() << endl;
      cout << mlfilter.back() << endl;
      cout << dither.back() << endl << endl;

      if(bounds_set==false)
      {
        //load dithering course file
        string file_name = dither[dither.size()-1];
        ifstream coursein(file_name.c_str());
        if(!coursein)
        {
          cout << "Could not open dithering course file: " << file_name << endl;
          exit(1);
        }

        string line;
        vector<string> data;

        getline(coursein,line);
        split(line,data);

        if(int(data.size())>=4)
        {
          bounds_xmin = atof(data[1].c_str());
          bounds_xmax = atof(data[1].c_str());
          bounds_ymin = atof(data[2].c_str());
          bounds_ymax = atof(data[2].c_str());
        }

        coursein.close();
        bounds_set = true;

      }
      load_course(dither[dither.size()-1]);   // For setting min/max bounds of dither area

      if(filter[filter.size()-1]>largestfilter) largestfilter = filter[filter.size()-1];
    }
  }

  string course_file, output_file;
  stringstream i_buff;
  for(int i=0; i<int(dither.size()); i++)
  {
    course_file = dither[i];
    ifstream src(course_file.c_str(), ios_base::binary);
    i_buff.clear();
    output_file = "";
    i_buff << "detector" << i << ".course";
    i_buff >> output_file;
    ofstream dst(output_file.c_str(), ios_base::binary);

    dst << src.rdbuf();
  }

  cout << detectorNames.size() << " observatories" << endl;

  vector<image> images(detectorNames.size());

  long seed;
  

  for(int im=0;im<int(detectorNames.size());im++)
  {
    if(images[im].load_detector(detectorNames[im])<0) exit(1);
    
    //free up the memory used by the psf loading process
    images[im].psf.flush_splines();

    //set the magnitude above which the small psf gets used
    images[im].set_largepsfmag(large_psf_mag);

    //set up the random seed
    images[im].pass_seed(&seed);
    if(im==0) images[im].reseed();

    images[im].set_image_properties(xpix[im],ypix[im]);

  }


  //Load in the microlensing events

  ifstream mlin(mleventsfile.c_str());
  if(!mlin)
  {
    cerr << "Could not open the list of microlensing events: " << mleventsfile << endl;
    exit(1);
  }

  vector<vector<double> > mlevents;
  
  while(!mlin.eof())
  {
    getline(mlin,line);
    split(line,data);

    if(data.size()<30)
      {
	continue;
      }

    vector<double> event;

    for(int i=0;i<int(data.size());i++)
      {
	//cerr << i << " " << data[i] << " " << atof(data[i].c_str()) << endl;
	event.push_back(atof(data[i].c_str()));
      }

    mlevents.push_back(event);
  }

  vector<int> mlids = vector<int> (mlevents.size());
  vector<int> mlidl = vector<int> (mlevents.size());
  

  vector<string> fields; //filename of all the fields

  vector<double> areas;  //areas of all the fields
  vector<int> nstars;    //number of stars in all the fields
  vector<int> sid;    //pointer between star and its data
  vector<int> sfn;    //pointer between star and its data

  vector<int> paddock; //TOO MANY FIELDS!!!! The field number of each 
                       //photometered star
  vector<int> ref;     //reference to the star in the input list 

  int nfields=0;       //number of starfields that will get read in

  int sourcefield=-1;
  int lensfield=-1;

  for(int i=5;i<argc-1;i+=2)
  {
    fields.push_back(string(argv[i]));
    areas.push_back(atof(argv[i+1]));
    nstars.push_back(0);
    nfields++;
  }

  //read in each starfield

  vector<vector<string> > stardata = vector<vector<string> >(nfields);

  //ignore any lines with these characters
  string ignore = 
    string("ABCDFGHIJKLMNOPQRSTUVWXYZabcdfghijklmnopqrstuvwxyz*:;=(),[]{}#");

  for(int i=0;i<nfields;i++)
  {
    ifstream in(fields[i].c_str());
    if(!in)
    {
      cerr << "Could not open field (" << fields[i] << ")." << endl;
      exit(1);
    }

    vector<double> data;
    string line;

    while(!in.eof())
    {
      getline(in,line);
      line = line.substr(0,line.find_first_of("#")); //remove comments

      if(line.find_first_of(ignore)==line.npos)
      {
        split(line,data);

        //keep track of the star's properties 
        //so we can do photometry with it later
        if(int(data.size())>=nfilters)
        {
          stardata[i].push_back(line);
          nstars[i]++;
        }
      }
    }

    cout << "Will use " << nstars[i] << " stars for photometry in field " << fields[i] << " over an area " << areas[i] << endl;

    in.close();
  }
  
  //start creating the image

  freestarlist stars;

  for(int k=0;k<int(detectorNames.size());k++) images[k].addbg();

  ofstream out((outroot + string("txt")).c_str());
  if(!out)
  {
    cerr << "Could not open output file (" << outroot + string(".txt") << ")" << endl;
    exit(1);
  }

  ofstream reg((outroot + string("reg")).c_str());
  if(!reg)
    {
      cerr << "Could not open region file (" << outroot + string(".reg") << ")" << endl;
      exit(1);
    }

  ofstream mlinfo((outroot + string("mlinfo")).c_str());
  if(!mlinfo)
    {
      cerr << "Could not open mlinfo file (" << outroot + string(".mlinfo") << ")" << endl;
      exit(1);
    }
		  

  double scaledif;
  vector<double> psfback; //for each image
  int xmin,xmax,ymin,ymax;
  double mag;
  int setdumb, ndraw, j;   //vars for poisson draw process
  

  //Now add stars to the first image
  for(int i=0;i<nfields;i++)
  {

    if(areas[i] == -1)
      {
	sourcefield=i;
      }
    if(areas[i] == -2)
      {
	lensfield=i;
      }
    areas[i]*=3600*3600; //convert to sq arcsec

    if(areas[i]<=0 && !(i==sourcefield || i==lensfield))
    {
      cerr << "Error: Nonsense solid angle inputed (" << areas[i] << "). No stars will be added to the image. If entering negative area for source or lens catalogs, use -1 for source -2 for lens." << endl;
      exit(1);
    }

    //calculate the dimensions of the star-field required

    double Afield;

    cout << "field " << i << endl;

    load_course(dither[0]);
    xmin = bounds_xmin;
    xmax = bounds_xmax;
    ymin = bounds_ymin;
    ymax = bounds_ymax;

    images[0].field_dimensions_roll(areas[i],xmin,xmax,ymin,ymax,Afield);
    
    cout << "stardata.size = " << stardata[i].size() << endl;

    //double x_cen = images[0].Xpix/2*images[0].psf.Nsub;
    //double y_cen = images[0].Ypix/2*images[0].psf.Nsub;
    double xsubmax = images[0].Xpix*images[0].psf.Nsub;
    double ysubmax = images[0].Ypix*images[0].psf.Nsub;

    //transit variables
    double period=7.13452;
    double tx1=0.5*(xmin+xmax);
    double ty1=0.5*(ymin+ymax);
    double tx2=tx1+50.0*9.0;
    double ty2=ty1;
    double tx3=tx2+0.5;
    double ty3=ty2+0.5;
    double m1=21.0;
    double delta1=0.01;
    double delta2=0.2;
    double m2 = m1 - 2.5*log10(delta1/delta2);
    

    //now we can generate the stars
    setdumb=0;
    if(areas[i]>0) ndraw = poisson(nstars[i]*Afield/areas[i],&seed);
    else ndraw = int(mlevents.size());
    for(int reps=0; reps < ndraw; reps++)
      {

	
	if(areas[i]>0)
	  {
	    //regular stars
	    j = randint(0,nstars[i]-1,&seed);
	    split(stardata[i][j],data);
	    mag = atof(data[filter[0]].c_str());
	    
	    images[0].freeaddstar(xmin,xmax,ymin,ymax,mag,&stars);
	  } //if areas>0
	else
	  { 
	    //microlenses

	    //Position of the source and lens at t=0
	    double xzero, yzero;
		
	    if(i==sourcefield)
	      {
		//Add the source first
		//09/18/20j = mlevents[reps][8];
		j = reps;
		split(stardata[i][j],data);
		cout << "Source mags " << j << " ";
		for(int ii=0; ii<5; ii++)
		  {
		    cout << data[54+ii].c_str() << " ";
		  }
		cout << "chosen " << atof(data[54+mlfilter[0]].c_str()) << endl;
		mag = atof(data[54+mlfilter[0]].c_str());
		//split(stardata[i][j],data);
		//mag = atof(data[filter[0]].c_str());

		//cerr << 54+mlfilter[0] << endl;

		//09/18/20 mag = mlevents[reps][54+mlfilter[0]];

		//cerr << "ML" << reps << " " << mag << " " << 54+mlfilter[0] << endl;

		//OK, need to solve for the position of the both lens and source at time 0 given t0, mul and mub of each
		
		//Choose the source's position at t=0 - inside the image
		xzero = xsubmax*ran2(&seed);
		yzero = ysubmax*ran2(&seed);
		

		reg << "MLs" << reps << " " << xzero << " " << yzero << " " << xzero/images[0].psf.Nsub << " " << yzero/images[0].psf.Nsub << " " << mag << endl << std::flush;
		cout << "MLs" << reps << " " << xzero << " " << yzero << " " << xzero/images[0].psf.Nsub << " " << yzero/images[0].psf.Nsub << " " << mag << endl << std::flush;
				
	      } //if i==sourcefield
	    else
	      {
		//Add the lens
		if(sourcefield==-1)
		  {
		    cerr << "Error: You must list the source field first, then the lens field" << endl;
		    exit(1);
		  }

		vector<double> sourcedata;
		//09/18/20int jsource = mlevents[reps][8];
		//09/18/20split(stardata[sourcefield][jsource],sourcedata);

		j = reps;
		split(stardata[i][j],data);
		cout << "Lens mags " << j << " ";
		for(int ii=0; ii<5; ii++)
		  {
		    cout << data[55+nfiltersml+ii].c_str() << " ";
		  }
		cout << "chosen " << atof(data[55+nfiltersml+mlfilter[0]].c_str()) << endl;
		
		mag = atof(data[55+nfiltersml+mlfilter[0]].c_str());		

		vector<double> lensdata;
		//09/18/20int jlens = mlevents[reps][17];
		//09/18/20split(stardata[lensfield][jlens],lensdata);

		/* If recovering from star lists
		double lmul = atof(lensdata[nfilters].c_str()) * scale;
		double lmub = atof(lensdata[nfilters+1].c_str()) * scale;

		double smul = atof(sourcedata[nfilters].c_str()) * scale;
		double smub = atof(sourcedata[nfilters+1].c_str()) * scale;
		*/

		/* If gathering from microlensing event */
		//arcsec to subpixels
		double mastosubpix = images[0].psf.Nsub/(images[0].psf.pixscale*1000.0);
		//Proper motions in subpixels per year
		/*09/18/20
		double lmul = mlevents[reps][20] * mastosubpix;
		double lmub = mlevents[reps][21] * mastosubpix;

		double smul = mlevents[reps][11] * mastosubpix;
		double smub = mlevents[reps][12] * mastosubpix;

		double t0 = mlevents[reps][32];
		double u0 = mlevents[reps][30];
		double thetaE = mlevents[reps][35];

	        mag = mlevents[reps][55+nfiltersml+mlfilter[0]];
		*/
		/* */

		double lmul = atof(data[20].c_str()) * mastosubpix;
		double lmub = atof(data[21].c_str()) * mastosubpix;

		double smul = atof(data[11].c_str()) * mastosubpix;
		double smub = atof(data[12].c_str()) * mastosubpix;

		double t0 = atof(data[32].c_str());
		double u0 = atof(data[30].c_str());
		double thetaE = atof(data[35].c_str());

	        //09/18/20mag = mlevents[reps][55+nfiltersml+mlfilter[0]];

		
		//Note a more natural formalisim is reversed l-s
		//double murel = qAdd(smul-lmul,smub-lmub);
		double amurel = atan2(smub-lmub,smul-lmul);

		//Source at t=t0
		double xs0 = stars.x[mlids[reps]] + (t0-0.0)/365.25 * smul;
		double ys0 = stars.y[mlids[reps]] + (t0-0.0)/365.25 * smub;

		//Lens at t=t0
		
		double xl0 = xs0 - u0*thetaE*mastosubpix*cos(amurel+pi/2.0);
		double yl0 = ys0 - u0*thetaE*mastosubpix*sin(amurel+pi/2.0);

		//Lens at t=0
		xzero = xl0 - (t0-0.0)/365.25 * lmul;
		yzero = yl0 - (t0-0.0)/365.25 * lmub;

		cout << "MLl" << reps << " " << xzero << " " << yzero << " " << xzero/images[0].psf.Nsub << " " << yzero/images[0].psf.Nsub << " " << mag << endl << std::flush;

		cout << "pred S-L @ t=t0 MLl" << reps << " " << xs0 << " " << ys0 << " " << xzero + (t0-0.0)/365.25*lmul << " " << yzero + (t0-0.0)/365.25*lmub << " " << xzero + (t0-0.0)/365.25*lmul-xs0 << " " << yzero + (t0-0.0)/365.25*lmub-ys0 << endl;

		mlinfo << reps << " " << t0 << " " << xs0 << " " << ys0 << " " << xl0 << " " << yl0 << " " << course_year[0]*365.25 << " " << stars.x[mlids[reps]] << " " << stars.y[mlids[reps]] << " " << xzero << " " << yzero << endl;

	      } //else if i==sourcefield

	    images[0].freeaddstar(xzero,yzero,mag,&stars);
	    if(i==sourcefield) mlids[reps] = int(stars.nstars-1);
	    else mlidl[reps] = int(stars.nstars-1);
	  } //end if areas>0 else

	//Store the star's location in lists
	sfn.push_back(i);
	sid.push_back(j);



	//Add each star to the other images	
	for(int im=1;im<int(detectorNames.size());im++)
	  {
	    //double x_cen = images[im].Xpix/2*images[im].psf.Nsub;
	    //double y_cen = images[im].Ypix/2*images[im].psf.Nsub;
	    //double xsubmax = images[im].Xpix*images[im].psf.Nsub;
	    //double ysubmax = images[im].Ypix*images[im].psf.Nsub;

	    
	    if(!setdumb)
	      {
		cout << images[im].psf.Nsub << " " << images[0].psf.Nsub << endl;
	      }
	    // load_course(dither[im]);
	    // xmin = bounds_xmin;
	    // xmax = bounds_xmax;
	    // ymin = bounds_ymin;
	    // ymax = bounds_ymax;
	    // images[im].field_dimensions(areas[i],xmin,xmax,ymin,ymax,Afield);
	    
	    scaledif=(images[im].psf.pixscale/images[0].psf.pixscale);
	    if(areas[i]>0)
	      {
		mag = atof(data[filter[im]].c_str());
	      }
	    else
	      {
		if(i==sourcefield) mag = mlevents[reps][54+mlfilter[im]];
		else mag = mlevents[reps][55+nfiltersml+mlfilter[im]];
	      }
	    
	    //why is this needed? it's just one star - MP: This is adding the same star to the different detectors
	    //add star at a random position
	    images[im].freeaddstar(stars.x[stars.nstars-1]/scaledif,stars.y[stars.nstars-1]/scaledif,mag);
	  }
	setdumb=1;
	
	out << stars.x[stars.nstars-1] << " " << stars.y[stars.nstars-1] << " " << stardata[i][j] << "\n";
      }
    
  } //end for reps


  //Add transits and a blended eclipsing binary
  images[0].freeaddstar(tx1,ty1,m1);
  images[0].freeaddstar(tx2,ty2,m1);
  images[0].freeaddstar(tx3,ty3,m2);
  cout << "Transit: " << tx1 << " " << tx2 << " " << tx1 << " " << tx2 << " " << tx1 << " " << tx2 << " " << m1 << " " << m2 << endl;;
  for(int im=1;im<int(detectorNames.size());im++)
    {
      scaledif=(images[im].psf.pixscale/images[0].psf.pixscale);
      images[im].freeaddstar(tx1/scaledif,ty1/scaledif,m1);
      images[im].freeaddstar(tx2/scaledif,ty2/scaledif,m1);
      images[im].freeaddstar(tx3/scaledif,ty3/scaledif,m2);
    }

  cout << "Added " << stars.nstars << " stars to the image" << endl;


  //Add in the microlensing events

  


  //add the background due to missed psf tails
  //Pretty sure this is currently wrong - should be the commented out line, but not sure why scaledif enters in
  for(int im=0;im<int(detectorNames.size());im++)
    {
      scaledif=(images[im].psf.pixscale/images[0].psf.pixscale);
      //psfback[im] = images[im].addpsfbg(xmin/scaledif,xmax/scaledif,ymin/scaledif,ymax/scaledif);
      images[im].addbg();
    }

  
  for(int im=0;im<int(detectorNames.size());im++)
  {
    char imno[12]; sprintf(imno,"%d",im);
    scaledif=(images[im].psf.pixscale/images[0].psf.pixscale);
    load_course(dither[im]);
    images[im].expose(texp[im],nstack[im]);
    images[im].set_time(course_year[0]*365.25);
    //images[im].set_crpix(images[im].Xpix/2,images[im].Ypix/2);

    //COMPUTE THE CRPIX center
    double cos_theta = cos(course_theta[0]);
    double sin_theta = sin(course_theta[0]);
    double x_cen = images[im].Xpix/2*images[im].psf.Nsub;
    double y_cen = images[im].Ypix/2*images[im].psf.Nsub;
    
    double crpx_t = images[0].Xpix*images[0].psf.Nsub/2/scaledif-course_x[0]*images[im].psf.Nsub;
    double crpy_t = images[0].Ypix*images[0].psf.Nsub/2/scaledif-course_y[0]*images[im].psf.Nsub;

    // 3. Rotate around image center
    double crpx = crpx_t*cos_theta+crpy_t*sin_theta-x_cen*cos_theta-y_cen*sin_theta+x_cen;
    double crpy = -crpx_t*sin_theta+crpy_t*cos_theta+x_cen*sin_theta-y_cen*cos_theta+y_cen;
	  
    //images[im].set_crpix(images[im].Xpix/2,images[im].Ypix/2);
    images[im].set_crpix(crpx/images[im].psf.Nsub,crpy/images[im].psf.Nsub);
    
    images[im].set_crval(meanl,meanb);
    images[im].set_wcs_rot(-course_theta[0]);
    images[im].write_fits(outroot+string("_detector")+string(imno)+string("_origin")+string("_stack.fits"),true);
    images[im].write_fits(outroot+string("_detector")+string(imno)+string("_origin")+string("_true.fits"),true);
  }



  //double x0,y0;
  double mul, mub;
  double scale;
  double star_x, star_y;

  // stringstream ss;
  string dither_index;
  int init_len;
  
  //not sure why these have to be initialized to a particular value
  double trans_x; //=30*9;
  double trans_y; //=30*9;
  double x; //=30*9;
  double y; //=30*9;

  vector<double> sx(mlevents.size()), sy(mlevents.size());
  
  for(int im=0;im<int(detectorNames.size());im++)
    {
      double x_cen = images[im].Xpix/2*images[im].psf.Nsub;
      double y_cen = images[im].Ypix/2*images[im].psf.Nsub;
      //double xsubmax = images[im].Xpix*images[im].psf.Nsub;
      //double ysubmax = images[im].Ypix*images[im].psf.Nsub;
    
      cout << "  Detector " << im << endl;
      
      load_course(dither[im]);

      for(int w=0;w<int(course_x.size());w++)
	{
	  cout << "    Waypoint " << w << endl;
	  
	  scale = images[0].psf.Nsub/(images[0].psf.pixscale*1000.0);
	  images[im].reset_image();
	  images[im].reset_detector();
	  scaledif=(images[im].psf.pixscale/images[0].psf.pixscale);
	  //y+=0.1372456;
	  //images[im].freeaddstar(x,y,mag);
	  
	  double cos_theta = cos(course_theta[w]);
	  double sin_theta = sin(course_theta[w]);
	  
	  for(int s=0;s<stars.nstars;s++)
	    {

	      if(s<mlids.front())
		{
	      //This data should be in an array by now, no need to redo the text processing.
		  split(stardata[sfn[s]][sid[s]],data);
		  mag = atof(data[filter[im]].c_str());
		  
		  // scale converts from milliarcseconds/year to arcseconds/year
		  // then we divide by 3600 to get degrees/year
		  mul = atof(data[nfilters].c_str()) * scale;
		  mub = atof(data[nfilters+1].c_str()) * scale;
		}


	      //These values are different for microlensing event stars


	      
	      //The sources
	      if(s>=mlids.front() && s<=mlids.back())
		{
		  int mlev = s-mlids.front();
		  double t0 = mlevents[mlev][32];
		  double u0 = mlevents[mlev][30];
		  double tE = mlevents[mlev][33];
		  mul = mlevents[mlev][11] * scale;
		  mub = mlevents[mlev][12] * scale;
		  double smag = mlevents[mlev][54+mlfilter[im]];
		  cout << im << " Source ";
		  for(int ii=0;ii<5;ii++)
		    {
		      cout << mlevents[mlev][54+ii] << " ";
		    }
		  cout << "actual " << smag << endl;
		  double mlmag = magnification(course_year[w]*365.25,t0,tE,u0);
		  mag = smag - 2.5*log10(mlmag);
		  //cerr << "ML event mag " << mag << endl;
		}

	      //The lenses
	      if(s>=mlidl.front() && s<=mlidl.back())
		{
		   int mlev = s-mlidl.front();
		   mul = mlevents[mlev][20] * scale;
		   mub = mlevents[mlev][21] * scale;
		   mag = mlevents[mlev][55+nfiltersml+mlfilter[im]];
		   cout << im << " Lens ";
		   for(int ii=0;ii<5;ii++)
		     {
		       cout << mlevents[mlev][55+nfiltersml+ii] << " ";
		     }
		   cout << "actual " << mag << endl;
		   //cerr << "ML lens mag " << mag << endl;
		}

	      
	      
	      // Updating star locations for a dither position has 3 steps:
	      
	      // 1. Calculate location after proper motion
	      star_x = stars.x[s] + course_year[w]*mul;
	      star_y = stars.y[s] + course_year[w]*mub;
	      
	      // 2. Adjust for dither x,y offset
	      trans_x = star_x/scaledif-course_x[w]*images[im].psf.Nsub;
	      trans_y = star_y/scaledif-course_y[w]*images[im].psf.Nsub;

	      // 3. Rotate around image center
	      x = trans_x*cos_theta+trans_y*sin_theta-x_cen*cos_theta-y_cen*sin_theta+x_cen;
	      y = trans_x*-1*sin_theta+trans_y*cos_theta+x_cen*sin_theta-y_cen*cos_theta+y_cen;

	      if(s>=mlids.front() && s<=mlids.back())
		{
		  int mlev = s-mlids.front();
		  sx[mlev] = x;
		  sy[mlev] = y;
		}
	      if(s>=mlidl.front() && s<=mlidl.back())
		{
		  int mlev = s-mlidl.front();
		  double dx = sx[mlev]-x;
		  double dy = sy[mlev]-y;

		  double t0 = mlevents[mlev][32];
		  double tE = mlevents[mlev][33];
		  double u0 = mlevents[mlev][30];
		  double thetaE = mlevents[mlev][35];

		  double t = course_year[w]*365.25;
		  
		  cout << "ML" << mlev << " " << t << " " << t/365.25 << " " << t-t0 << " " << dx << " " << dy << " " << qAdd(dx,dy)/images[im].psf.Nsub*images[im].psf.pixscale*1000.0 << " " << thetaE*qAdd(u0,((t-t0)/tE)) << " " << qAdd(dx,dy)/images[im].psf.Nsub*images[im].psf.pixscale*1000.0/thetaE << " " << (t-t0)/tE << endl << std::flush;
		}
	      
	      images[im].freeaddstar(x, y, mag);
	    }
	  
	  images[im].addbg();
	  char imno[12]; sprintf(imno,"%d",im);
	  
	  char dithno[15]; sprintf(dithno,"%d",w);
	  dither_index = string(dithno);
	  // ss.str(string());
	  // ss << w;
	  // dither_index = "";
	  // ss >> dither_index;
	  init_len = dither_index.length();
	  for(int sc=0; sc<6-init_len; sc++)
	    {
	      dither_index = "0" + dither_index;
	    }
	  
	  images[im].expose(texp[im],nstack[im]);
	  images[im].set_time(course_year[w]*365.25);

	  //COMPUTE THE CRPIX center
	  double crpx_t = images[0].Xpix*images[0].psf.Nsub/2/scaledif-course_x[w]*images[im].psf.Nsub;
	  double crpy_t = images[0].Ypix*images[0].psf.Nsub/2/scaledif-course_y[w]*images[im].psf.Nsub;

	  // 3. Rotate around image center
	  double crpx = crpx_t*cos_theta+crpy_t*sin_theta-x_cen*cos_theta-y_cen*sin_theta+x_cen;
	  double crpy = -crpx_t*sin_theta+crpy_t*cos_theta+x_cen*sin_theta-y_cen*cos_theta+y_cen;
	  
	  //images[im].set_crpix(images[im].Xpix/2,images[im].Ypix/2);
	  images[im].set_crpix(crpx/images[im].psf.Nsub,crpy/images[im].psf.Nsub);
	  images[im].set_crval(meanl,meanb);
	  images[im].set_wcs_rot(course_theta[w]);
	  images[im].write_fits(outroot+string("_detector")+string(imno)+string("_")+dither_index+string("_stack.fits"),true);
	  //images[im].freesubstar(x,y,mag);
    
	}
    }
  
  // cout << images[0].Xpix/2 << "   " << images[0].Ypix/2 << endl;
  
  cout << "Images written" << endl;
  
  // cout << "Course file had " << course_x.size() << " x entries and " << course_y.size() << " y entries" << endl;
  // for(int i=0;i<course_x.size();i++)
  // {
  //   cout << course_x[i] << " " << course_y[i] << endl;
  // }

}
