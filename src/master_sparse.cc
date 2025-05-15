#include <iostream>
#include <string>
#include <netcdf>
#include <cstdio>
//#include <omp.h>
#include "io.h"
#include "cmemt.h"
#include "atmosphere.h"
#include "math_tools.h"
//#include "sparse.h"
#include <mpi.h>
#include "comm.h"
#include "depthmodel.h"
#include "clte.h"
#include "crh.h"
#include "fpi.h"
#include <algorithm>
//
using namespace netCDF;
using namespace std;
//
void read_instruments(iput_t &iput)
{

  static const string no = "none";

  int const nReg = int(iput.regions.size());
  for(int ii=0; ii<nReg; ++ii){
    if(iput.regions[ii].inst == no) continue;

    io ifil(file_exists(iput.regions[ii].ifile), NcFile::read, false);
    ifil.read_Tstep<double>("iprof",iput.regions[ii].psf , 0, false);

    if(iput.regions[ii].inst == "specrebin"){
      mat<int> tmp;
      ifil.read_Tstep<int>("rebin", tmp, 0, false);

      iput.regions[ii].reb = tmp.d[0];

      fprintf(stderr,"read_instruments: nRebin=%d\n",iput.regions[ii].reb );
    }

    std::string siz = formatVect<int>(iput.regions[ii].psf.getdims());
    if(iput.verbose)
      fprintf(stderr,"read_instruments: region[%3d] <%s> %s\n", ii,iput.regions[ii].ifile.c_str(), siz.c_str() );

  }

}


void logSentReceived(int const dir, FILE* const io, int const ipix, int const nx)
{
  int const yy = ipix / nx;
  int const xx = ipix - yy*nx;

  if(dir == 0){
    fprintf(io, "0 %d %d\n", yy,xx);
  }else{
    fprintf(io, "1 %d %d\n", yy,xx);
  }
  fflush(io);
}

//
void slaveInversion(iput_t &iput, mdepthall_t &m, mat<double> &obs, mat<double> &x, mat<double> &chi2, mat<double> &dsyn){

  /* --- Init dimensions --- */
  int const nx = x.size(1);

  unsigned long ntot = (unsigned long)(x.size(0) * x.size(1));
  int ncom = (int)(std::floor(ntot / (double)iput.npack));
  if((unsigned long)(ncom * iput.npack) != ntot) ncom++;
  int nprocs = iput.nprocs;
  int iproc = 0;
  unsigned long ipix = 0;
  int tocom = ncom;
  FILE* io = fopen("sent_received_log.txt", "w");

  int compute_gradient = 0; // dummy parameter here
  chi2.set({x.size(0), x.size(1)});


  if(nprocs > 1){




    /* --- Init slaves with first package --- */
    for(int ss = 1; ss<=min(nprocs-1,ncom); ss++){
      logSentReceived(0,io,ipix,nx);
      comm_master_pack_data(iput, obs, x, ipix, ss, m, compute_gradient);
    }
    int per  = 0;
    int oper  = -1;
    float pno =  100.0 / double(std::max<int>(1, ntot - 1));
    unsigned long irec = 0;
    fprintf(stdout,"\rProcessed -> %d%s -> sent=%d, received=%d     ", per, "%", ipix, irec);
    fflush(stdout);


    /* --- manage packages as long as needed --- */
    while(irec < ntot){

      // Receive processed data from any slave (iproc)
      int rec_pix = comm_master_unpack_data(iproc, iput, obs, x, chi2, irec, dsyn, compute_gradient, m);
      logSentReceived(1,io,rec_pix,nx);

      per = irec * pno;

      // Send more data to that same slave (iproc)
      if(ipix < ntot){
	logSentReceived(0,io,ipix,nx);
	comm_master_pack_data(iput, obs, x, ipix, iproc, m, compute_gradient);
      }

      // Keep count of communications left
      tocom--;

      // Printout
      if(per > oper){
	oper = per;
	//cout << "\rProcessed -> "<<per<<"% -> sent="<<ipix<<", received="<<irec;
	fprintf(stdout,"\rProcessed -> %d%s -> sent=%d, received=%d", per, "%", ipix, irec);
	fflush(stdout);
      }
    }

    fprintf(stdout, "\n");
    fflush(stdout);
  }
  fclose(io);
}

void master_inverter(mdepthall_t &model, mat<double> &pars, mat<double> &obs, mat<double> &w, iput_t &input)
{

  int ndep = (int)model.ndep, nx = input.nx, ny = input.ny;
  mdepth_t m(ndep);
  atmos *atm;

  /* --- Init atmos --- */

  if(input.atmos_type == string("lte")) atm = new clte(input, 4.44);
  else if(input.atmos_type == string("rh")) atm = new crh(input, 4.44);
  else{
    cerr<<"master_invert: ERROR, atmos unknown -> "<<input.atmos_type<<endl;
    exit(0);
  }


  /* --- (TO-DO, change this!) --- */

  vector<instrument*> inst;
  int nreg = atm->input.regions.size();
  inst.resize(nreg);

  for(int kk = 0; kk<nreg; kk++){
    if(atm->input.regions[kk].inst == "spectral") inst[kk] = new spectral(atm->input.regions[kk], 1);
    else if(atm->input.regions[kk].inst == "fpi") inst[kk] = new     sfpi(atm->input.regions[kk], 1);
    else inst[kk] = new instrument();
  }
  atm->inst = &inst[0];


  /* --- Loop and invert --- */

  int per = 0, oper = -1, kk = 0;

  for(int yy = 0; yy<input.ny; yy++)
    for(int xx = 0; xx<input.nx; xx++){

      /* --- Copy data to single pixel model --- */

      memcpy(&m.cub.d[0], &model.cub(yy,xx,0,0), 12*ndep*sizeof(double));
      m.bound_val = model.boundary(yy,xx);

      // --- update instrumental info --- //

      for(int kk = 0; kk<nreg; kk++){
	int nd = input.regions[kk].psf.ndims();
	if(nd == 1)
	  inst[kk]->update(input.regions[kk].psf.d.size(), &input.regions[kk].psf.d[0]);
	else
	  inst[kk]->update(input.regions[kk].psf.size(nd-1), &input.regions[kk].psf(yy,xx,0));
      }

      /* --- invert --- */


      atm->fitModel2( m, input.npar, &pars(yy,xx,0),
		    (int)(input.nw_tot*input.ns), &obs(yy,xx,0,0), w);


      /* --- Copy inverted model back to model cube --- */

      memcpy( &model.cub(yy,xx,0,0), &m.cub.d[0], 12*ndep*sizeof(double));


      per = ++kk / std::max(1,nx*ny-1);
      if(per > oper){
	oper = per;
	fprintf(stderr,"\rmaster_invert: %d%s",per,"% ");
      }
    }

  fprintf(stderr,"master_invert: %d%s",per,"%");


  /* --- Clean-up --- */
  delete atm;
  for(auto &it: inst)
    delete it;
}



void do_master_sparse(int myrank, int nprocs,  char hostname[]){

  /* --- Printout number of processes --- */
  cerr << "STIC: Initialized with "<<nprocs <<" process(es)"<<endl;
  mat<double> model, obs, dobs, wav, w, syn, chi2;;
  mdepthall_t im;

  static const vector<string> vnames = {"temperature","vlos","vturb", "Blong", "Bhor", "azi","dens", "nne"};



  /* --- Read input files and fill in variables---*/
  iput_t input = read_input("input.cfg", true);
  read_lines("lines.cfg", input, input.verbose);
  wav.set(vector<int>{input.nw_tot});
  wav.d = fill_lambdas(input, false);


  input.nprocs = nprocs; input.myrank = myrank;
  input.myid = (string)"master, ";

  /* --- Set OpenMP threads in the master --- */

  // omp_set_num_threads(input.master_threads);


  /* --- Open input files with the models and profiles --- */

  io ipfile, omfile, imfile(file_exists(input.imodel), NcFile::read);

  bool inversion = false;
  if(input.mode == 1 || input.mode == 3) inversion = true;

  if(inversion)
    ipfile.initRead(file_exists(input.iprof), NcFile::read);



  /* --- Get dimensions and fill input struct --- */

  input.nt = 1;
  if(inversion){
    int off = 0;
    vector<int> odims = ipfile.dimSize(string("profiles"));
    if(odims.size() == 5) {
      off = 1;
      input.nt = odims[0];
    }
    input.ny = odims[0+off];
    input.nx = odims[1+off];
    input.ns = odims[3+off];
    if(input.nw_tot != odims[2+off]){
      cerr << input.myid <<"ERROR, number of wavelength points in regions != "<<input.iprof<<endl;
      exit(0);
    }
  }else{
    int off = 0;
    vector<int> odims = imfile.dimSize(string("temp"));
    if(odims.size() == 4) {
      off = 1;
      input.nt = odims[0];
    }
    input.ny = odims[0+off];
    input.nx = odims[1+off];
    input.ns = 4;

    obs.set(vector<int>{input.ny, input.nx, input.nw_tot, input.ns});
  }
  vector<int> dims = {input.ny, input.nx, input.nw_tot, input.ns};


  /* --- Read time-independent variables --- */

  if(inversion){
    ipfile.read_Tstep<double>(string("weights"), w);
    ipfile.read_Tstep<double>(string("wav"), wav);

    if(w.d.size() == 0) {
      w.set({input.nw_tot, input.ns});
      fill(w.d.begin(), w.d.end(), 3.0e-3); // Default noise?
    }
  }

  // Read model for tstep = 0
  input.boundary = im.read_model2(input, input.imodel, 0, true);
  input.ndep = im.ndep;

  /* ---
     Open output files and init vars to store results
     (dimension = 0 means unlimited)
     We decide here if the variables are going to be
     stored as floats or doubles, regardless of their
     type in memory
     --- */
  io opfile(input.oprof,  NcFile::replace);
  opfile.initDim({"time","ndep","vtype", "y", "x", "wav", "stokes"},{0, input.ndep, input.nresp, input.ny, input.nx, input.nw_tot, input.ns});
  opfile.initVar<double>(string("profiles"), {"time","y", "x", "wav", "stokes"});
  opfile.initVar<double>(string("wav"), {"wav"});

  opfile.write_Tstep<double>(string("wav"), wav);
  if(input.mode == 4){

    string vn;
    vector<string> vnv;
    mat<int> idx({input.nresp});
    int k = 0;

    for(int ii=0;ii<8;ii++)
      if(input.getResponse[ii]){
	vn += vnames[ii]+string(" ");
	vnv.push_back(vnames[ii]);
	idx(k++) = ii;
      }

    cerr<<"master: computing derivatives for ["<<vn<<"]"<<endl;
    opfile.initVar<int>(string("vidx"), {"vtype"});
    opfile.write_Tstep<int>(string("vidx"), idx);


    opfile.initVar<double>(string("derivatives") ,{"time","y", "x", "vtype", "ndep", "wav", "stokes"});
    opfile.varAttr("derivatives","units", vn);
    dobs.set({input.ny, input.nx, input.nresp, input.ndep, input.nw_tot, input.ns});
  }


  if(inversion){
    opfile.initVar<float>(string("weights"), {"wav", "stokes"});
    opfile.write_Tstep<double>(string("weights"), w);

    // omfile.initRead(input.omodel, NcFile::replace);

    vector<int> dims = ipfile.dimSize("profiles");
    if(dims.size() == 5) dims.erase(dims.begin(), dims.begin()+1); // remove time

    //
    // Propagate information to slaves via MPI
    //
    double mmax = im.cub(0,0,9,0), mmin = im.cub(0,0,9,0);
    for(int zz = 1; zz < input.ndep; zz++){
      if(im.cub(0,0,9,zz) > mmax) mmax = im.cub(0,0,9,zz);
      if(im.cub(0,0,9,zz) < mmin) mmin = im.cub(0,0,9,zz);
    }

    /* --- Place Nodes --- */
    //input.npar = set_nodes(input.nodes,  mmin, mmax, input.verbose);

    {
      vector<double> idep(input.ndep, 0.0);

      if(input.nodes.depth_t == 0)
	memcpy(&idep[0], &im.cub(0,0,9,0), input.ndep*sizeof(double));
      else
	memcpy(&idep[0], &im.cub(0,0,11,0), input.ndep*sizeof(double));

      input.npar = set_nodes(input.nodes, idep, input.dint, input.verbose);

    }

  } // if inversion

  read_instruments(input);


  /* --- MPI bussiness --- */

  if(input.npack < 0){
    input.npack = (int)((double)(input.ny*input.nx) / nprocs) + 1;
    if(input.verbose) cerr<<input.myid<<"Using NPACK="<<input.npack<<endl;
  }

  comm_get_buffer_size(input);
  MPI_Barrier(MPI_COMM_WORLD); // Wait until all processors reach this point
  comm_send_parameters(input);

  if(inversion) comm_send_weights(input, w); //




  /* --- Init sparse class --- */

  // sparse2d inv;
  //if(input.mode == 3)
  // inv.init(input,  dims, input.sparse_threshold,
  //		 string2wavelet(input.wavelet_type), input.wavelet_order,
  //		 spt_hard, input.master_threads);

  //
  // Main loop
  //
  for(int tt = 0; tt<input.nt; tt++){ // Loop in time

    /* --- Read Tstep data --- */

    mat<double> pweight;
    if(inversion){
      ipfile.read_Tstep(string("profiles"), obs, tt);
      if(ipfile.is_var_defined((string)"pixel_weights")){
	ipfile.read_Tstep<double>(string("pixel_weights"), pweight, tt);
      }
    }

    if(tt > 0) im.read_model2(input, input.imodel, tt, true);

    /* --- Check dimensions in inversion mode --- */

    if(inversion){
      if(obs.size(0) != im.cub.size(0) && obs.size(1) != im.cub.size(1)){
	cerr << input.myid <<"ERROR, the input model and the observations do not have the same dimensions in X,Y axes:"<<endl;
	cerr << "   -> "<<input.imodel<<" "<<formatVect<int>(im.temp.getdims())<<endl;
	cerr << "   -> "<<input.iprof<<" "<<formatVect<int>(obs.getdims())<<endl;

	comm_kill_slaves(input, nprocs);
	exit(0);
      }

      if(obs.size(2) != input.nw_tot){
	cerr << input.myid <<"ERROR, number of wavelenghts in input.cfg ["<<input.nw_tot<< "] does not match the observations ["<< obs.size(2) << "]"<<endl;
	comm_kill_slaves(input, nprocs);
	exit(0);

      }

      /* --- get free-parameters from input model --- */


    } // inversion mode

    im.model_parameters2(model, input.nodes);



    /* --- Invert data --- */

    if     (input.mode == 1){

      if(nprocs == 1)
	master_inverter(im, model, obs, w, input);
      else
	slaveInversion(input, im, obs, model, chi2, dobs); // implemented above!

    }else if(input.mode == 2) slaveInversion(input, im, obs, model, chi2, dobs); // it won't invert if mode == 2
    //else if(input.mode == 3) inv.SparseOptimization(obs, model, w, im, pweight);
    else if(input.mode == 4) slaveInversion(input, im, obs, model, chi2, dobs);

    if(inversion){

      /* --- Expand fitted parameters into depth-stratified atmos --- */

      // if(input.depth_model == 0)
      //im.expandAtmos(input.nodes, model, input.dint);


      /* --- init output for inverted model --- */

      if(tt == 0){
	vector<int> odim = model.getdims();
	odim.insert(odim.begin(), 0);
	//omfile.initDim({"time","y", "x", "par"},odim);
	//omfile.initVar<float>(string("model"), {"time","y", "x", "par"});
      }


      /* --- Write model parameters, profiles and depth-stratified atmos --- */

      //omfile.write_Tstep(string("model"), model, tt);
      im.write_model2(input, input.oatmos, tt);
    }

    opfile.write_Tstep(string("profiles"), obs, tt);
    if(input.mode == 4) opfile.write_Tstep(string("derivatives"), dobs, tt);

  }



  /* --- Tell slaves to exit while(1) loop --- */

  comm_kill_slaves(input, nprocs);

}

