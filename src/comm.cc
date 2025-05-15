#include <mpi.h>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <chrono>
#include <thread>
#include "input.h"
#include "comm.h"
#include "cmemt.h"
#include "depthmodel.h"
//
using namespace std;
//
int getNinstrumentData(std::vector<region_t> const &reg)
{
  int const nReg = int(reg.size());
  int ndata = 0, nRegwithdata = 0;
  
  for(int ii=0; ii<nReg; ++ii){
    if(reg[ii].inst == string("none")) continue;
    std::vector<int> dim = reg[ii].psf.getdims();
    int ndim = int(dim.size());
    ndata += dim[ndim-1];
    nRegwithdata+= 1;
  }
  return ndata+nRegwithdata+1;
}

std::vector<double> packInstrumentalData(iput_t &input, int yy, int xx)
{
  std::vector<double> res;
  int ndata = getNinstrumentData(input.regions);
  int offset = 0;
  if(ndata > 0){
    res.resize(ndata, 0.0);
    res[offset++] = double(ndata);
    int const nReg = int(input.regions.size());
    for(int ii=0; ii<nReg; ++ii){
      if(input.regions[ii].inst == string("none")) continue;
      int const ndim = input.regions[ii].psf.ndims();
      int const nd =  input.regions[ii].psf.size(ndim-1);
      res[offset++]= double(nd);

      if(ndim == 1)
	for(int ww=0; ww<nd; ++ww)
	  res[offset++] = input.regions[ii].psf.d[ww];
      else
	for(int ww=0; ww<nd; ++ww)
	  res[offset++] = input.regions[ii].psf(yy,xx,ww);

    }
  }
  return res;
}
void unpackInstrumentalData(iput_t &input, std::vector<double> const &buff)
{
  static const string no = "none";
  int const nReg = int(input.regions.size());
  int offset = 0;
  std::vector<int> idim(1,0);
  
  for(int ii=0; ii<nReg; ++ii){
    if(input.regions[ii].inst == no) continue;
    const long int &head = long(buff[offset++]+0.1);
    //cerr<<"detected "<<head<<" elements in region "<<ii<<endl;
    idim[0] = int(head);
    input.regions[ii].psf.set(idim);
    memcpy(&input.regions[ii].psf.d[0], &buff[offset], head*sizeof(double));
    offset += head;
  }
}
//
void comm_get_buffer_size(iput_t &input){
  int ninstrumentaldata = getNinstrumentData(input.regions);
  
  int maxlen = 0;
  string inam = "comm_get_buffer_size: ";
  unsigned long ntot, rest, tocom;
  //
  switch(input.mode)
    { 
    case 1: // Inversion
      input.buffer_size =
	(input.nw_tot*4*sizeof(double) +  // Profiles
	 input.npar*sizeof(double) + // Model
	 (12 * input.ndep +2+ 1)* sizeof(double) + //non-inverted quantities
	 2*sizeof(double)) * input.npack + // Chi2, boundary value
	6*sizeof(int) +      // xx, yy, iproc, pix, action, npacked
	ninstrumentaldata * sizeof(double);      
      
      input.buffer_size1 = input.buffer_size;
      break;
      //
    case 2: // Synthesis
      input.buffer_size =
	(13 * input.ndep * sizeof(double)) * input.npack + // depth-stratified quantities
	6*sizeof(int) + // xx, yy, iproc, pix, action, npacked
	ninstrumentaldata * sizeof(double);      
      input.buffer_size1 =
	(input.nw_tot*4*sizeof(double)) * input.npack + 6*sizeof(int);
      break;
      //
    case 3: // Synthesis + derivatives
      input.buffer_size =
	8*sizeof(int) + // xx, yy, iproc, pix, action, npacked, ndep, cgrad
	input.npar*sizeof(double) * input.npack + // Values of the nodes * npacked
	1*sizeof(double) + // perturbation to the parameter
	(12*input.ndep + 1+2)*sizeof(double)*input.npack; // atmospheric model
      
      input.buffer_size1 = (input.nw_tot*4*sizeof(double) +                    
			    input.nw_tot*4*input.npar*sizeof(double) + // Derivatives
			    1*sizeof(double)) * input.npack +// perturbation to the parameter
	(13*input.ndep+2)*input.npack*sizeof(double)+ // Send back the pressure scale
	                    6*sizeof(int); // xx, yy, ipix, npacked, iproc

    case 4:// Synthesis + derivatives at all heights
      input.buffer_size =
	(12 * input.ndep * sizeof(double)) * input.npack + // depth-stratified quantities
	6*sizeof(int)+ninstrumentaldata * sizeof(double);; // xx, yy, iproc, pix, action, npacked
      
      input.buffer_size1 = (input.nw_tot*4*sizeof(double) * (input.ndep*input.nresp+1)) + 6*sizeof(int);
      break;
    default:
      cout << input.myid<< inam <<"ERROR, work mode ("<<input.mode<<") not valid"<<endl;
      exit(0);
      break;
    }
  cout << input.myid << inam << "MPI_buffer_size = "<<input.buffer_size/(8.*1024.)<<" kbytes"<<endl;
  cout << input.myid << inam << "MPI_buffer_size1 = "<<input.buffer_size1/(8.*1024.)<<" kbytes"<<endl;
}

//
void comm_send_parameters(iput_t &input){


 // number of lines
  int nline = (int)input.lines.size();
  int status = 0;
  int nregions = (int)input.regions.size();
  //
  status = MPI_Bcast(&nline,     1,    MPI_INT, 0, MPI_COMM_WORLD);
  status = MPI_Bcast(&nregions,  1,    MPI_INT, 0, MPI_COMM_WORLD);  
  status = MPI_Bcast(&input.buffer_size,  2,    MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

  status = MPI_Bcast(&input.nt, 39,    MPI_INT, 0, MPI_COMM_WORLD); // We are sending 11 ints from the struct!
  status = MPI_Bcast(&input.nodes.regul_type, 9,    MPI_INT, 0, MPI_COMM_WORLD); // We are getting 15 ints from the struc
  status = MPI_Bcast(&input.nodes.rewe, 10,    MPI_DOUBLE, 0, MPI_COMM_WORLD); // We are getting 15 ints from the struct!
  // status = MPI_Bcast(&input.nodes.nregul,     1,    MPI_INT, 0, MPI_COMM_WORLD);


  status = MPI_Bcast(&input.mu,  8, MPI_DOUBLE, 0, MPI_COMM_WORLD); // We are sending 4 doubles from the struct!
  status = MPI_Bcast(&input.max_inv_iter,  1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    
  //int dummy = (int)input.verbose;
  //status = MPI_Bcast(&dummy,  1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  status = MPI_Bcast(&input.nodes.depth_t, 1,    MPI_INT, 0, MPI_COMM_WORLD); 


  
  // Instrument type
  //status = MPI_Bcast(const_cast<char *>(input.instrument.c_str()), input.inst_len+1,  MPI_CHAR, 0, MPI_COMM_WORLD);
  
  // Atmos type
  std::string tempo = input.atmos_type;
  //if(tempo.size() == 0) tempo = " ";
  status = MPI_Bcast(const_cast<char *>(tempo.c_str()), tempo.size()+1, MPI_CHAR, 0, MPI_COMM_WORLD);

  // Abund file
  tempo = input.abfile;
  //if(tempo.size() == 0) tempo = " ";
  status = MPI_Bcast(const_cast<char *>(tempo.c_str()), tempo.size()+1, MPI_CHAR, 0, MPI_COMM_WORLD);

  // Line structs
  if(nline > 0){
    for (int ll = 0;ll<nline;ll++) {
      status = MPI_Bcast(&input.lines[ll].elem[0], 8,   MPI_CHAR, 0, MPI_COMM_WORLD); // We are getting 23 chars
      status = MPI_Bcast(&input.lines[ll].Jup,     19, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      status = MPI_Bcast(&input.lines[ll].anum,     4,    MPI_INT, 0, MPI_COMM_WORLD); 
    }
  }

  // Regions
  for (int ll = 0;ll<nregions;ll++){
    status = MPI_Bcast(&input.regions[ll].w0, 3,   MPI_DOUBLE, 0, MPI_COMM_WORLD); 
     status = MPI_Bcast(&input.regions[ll].nw, 2,   MPI_INT, 0, MPI_COMM_WORLD);
    //
    int tmp = (int)input.regions[ll].ifile.size() + 1;
    status = MPI_Bcast(&tmp, 1,   MPI_INT, 0, MPI_COMM_WORLD);
    status = MPI_Bcast(const_cast<char *>(input.regions[ll].ifile.c_str()), tmp,   MPI_CHAR, 0, MPI_COMM_WORLD);
    //
    tmp = (int)input.regions[ll].inst.size() + 1;
    status = MPI_Bcast(&tmp, 1,   MPI_INT, 0, MPI_COMM_WORLD);
    status = MPI_Bcast(const_cast<char *>(input.regions[ll].inst.c_str()), tmp,   MPI_CHAR, 0, MPI_COMM_WORLD);
  }
  if(input.mode == 1 || input.mode == 3){
    // nodes in temp
    int nnodes = (int)input.nodes.temp.size();
    status = MPI_Bcast(&nnodes,     1,    MPI_INT, 0, MPI_COMM_WORLD);
    if(nnodes > 0) status = MPI_Bcast(&input.nodes.temp[0],     nnodes,    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // nodes in vturb
    nnodes = (int)input.nodes.v.size();
    status = MPI_Bcast(&nnodes,     1,    MPI_INT, 0, MPI_COMM_WORLD);
    if(nnodes > 0) status = MPI_Bcast(&input.nodes.v[0],     nnodes,    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // nodes in vlos
    nnodes = (int)input.nodes.vturb.size();
    status = MPI_Bcast(&nnodes,     1,    MPI_INT, 0, MPI_COMM_WORLD);
    if(nnodes > 0) status = MPI_Bcast(&input.nodes.vturb[0],     nnodes,    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // nodes in blo
    nnodes = (int)input.nodes.bl.size();
    status = MPI_Bcast(&nnodes,     1,    MPI_INT, 0, MPI_COMM_WORLD);
    if(nnodes > 0) status = MPI_Bcast(&input.nodes.bl[0],     nnodes,    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // nodes in bhor
    nnodes = (int)input.nodes.bh.size();
    status = MPI_Bcast(&nnodes,     1,    MPI_INT, 0, MPI_COMM_WORLD);
    if(nnodes > 0) status = MPI_Bcast(&input.nodes.bh[0],     nnodes,    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // nodes in azi
    nnodes = (int)input.nodes.azi.size();
    status = MPI_Bcast(&nnodes,     1,    MPI_INT, 0, MPI_COMM_WORLD);
    if(nnodes > 0) status = MPI_Bcast(&input.nodes.azi[0],     nnodes,    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
    // Total number of nodes
    status = MPI_Bcast(&input.nodes.nnodes,     13,    MPI_INT, 0, MPI_COMM_WORLD);
    
    
    // Variables to Invert
    status = MPI_Bcast(input.nodes.toinv,     8,    MPI_INT, 0, MPI_COMM_WORLD);

    
    
    
    // Send type of nodes
    status = MPI_Bcast(&input.nodes.ntype[0],   input.nodes.nnodes,    MPI_INT, 0, MPI_COMM_WORLD);
    

    input.nodes.bound = input.boundary;
    
  }

}

void comm_recv_parameters(iput_t &input){
  string inam = "comm_recv_parameters: ";
  //
  int status=0, nline=0, nregions=0;
  status = MPI_Bcast(&nline,     1,    MPI_INT, 0, MPI_COMM_WORLD);
  status = MPI_Bcast(&nregions,  1,    MPI_INT, 0, MPI_COMM_WORLD);
  status = MPI_Bcast(&input.buffer_size,  2,    MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  status = MPI_Bcast(&input.nt, 39,    MPI_INT, 0, MPI_COMM_WORLD); // We are getting 15 ints from the struct!

  status = MPI_Bcast(&input.nodes.regul_type, 9,    MPI_INT, 0, MPI_COMM_WORLD); // We are getting 15 ints from the struct!
  status = MPI_Bcast(&input.nodes.rewe, 10,    MPI_DOUBLE, 0, MPI_COMM_WORLD); // We are getting 15 ints from the struct!

  // status = MPI_Bcast(&input.nodes.nregul,     1,    MPI_INT, 0, MPI_COMM_WORLD);


  status = MPI_Bcast(&input.mu,  8, MPI_DOUBLE, 0, MPI_COMM_WORLD); // We are getting 4 doubles from the struct!
  status = MPI_Bcast(&input.max_inv_iter,  1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  //int dummy = 0;
  //status = MPI_Bcast(&dummy,  1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  //input.verbose = (bool)dummy;
  status = MPI_Bcast(&input.nodes.depth_t, 1,    MPI_INT, 0, MPI_COMM_WORLD); 


  
  //
  //  cerr << nline << endl;
  input.lines.resize(nline);
  input.regions.resize(nregions);
  
  
  { // Instrument type
    // input.instrument.clear();
    // char buf[input.inst_len+1];
    /// status = MPI_Bcast(&buf,     input.inst_len+1,    MPI_CHAR, 0, MPI_COMM_WORLD);
    //input.instrument = removeSpaces(string(buf));
  }
  { // Atmos type
    input.atmos_type.clear();
    std::vector<char> buf(input.atmos_len+2, 0);
    //char buf[input.atmos_len+1];
    status = MPI_Bcast(&buf[0],     input.atmos_len+1,    MPI_CHAR, 0, MPI_COMM_WORLD);
    input.atmos_type = removeSpaces(string(&buf[0]));
  }

  { // abund type
    input.abfile.clear();
    std::vector<char> buf(input.ab_len+2, 0);
    status = MPI_Bcast(&buf[0],     input.ab_len+1,    MPI_CHAR, 0, MPI_COMM_WORLD);
    input.abfile = removeSpaces(string(&buf[0]));
  }

  if(nline > 0){
    for (int ll = 0;ll<nline;ll++){
      status = MPI_Bcast(input.lines[ll].elem, 8,   MPI_CHAR, 0, MPI_COMM_WORLD); // We are getting 23 chars
      status = MPI_Bcast(&input.lines[ll].Jup,     19, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      status = MPI_Bcast(&input.lines[ll].anum,     4,    MPI_INT, 0, MPI_COMM_WORLD);
    }
  }
  // Regions
  for (int ll = 0;ll<nregions;ll++){
    status = MPI_Bcast(&input.regions[ll].w0, 3,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
    status = MPI_Bcast(&input.regions[ll].nw, 2,   MPI_INT, 0, MPI_COMM_WORLD);

    {
      int tmp = 0;
      status = MPI_Bcast(&tmp, 1,   MPI_INT, 0, MPI_COMM_WORLD);
      //char buf[tmp];
      std::vector<char> buf(tmp+2, 0);
      status = MPI_Bcast(&buf[0], tmp, MPI_CHAR, 0, MPI_COMM_WORLD);
      input.regions[ll].ifile =  removeSpaces(string(&buf[0]));
    }
    {
      int tmp = 0;
      status = MPI_Bcast(&tmp, 1,   MPI_INT, 0, MPI_COMM_WORLD);
      //char buf[tmp];
      std::vector<char> buf(tmp+2, 0);
      status = MPI_Bcast(&buf[0], tmp, MPI_CHAR, 0, MPI_COMM_WORLD);
      input.regions[ll].inst =  removeSpaces(string(&buf[0]));
    }

  } // ll

  if(input.mode == 1 || input.mode == 3){
    
    // nodes in temp
    int nnodes = 0;
    status = MPI_Bcast(&nnodes,     1,    MPI_INT, 0, MPI_COMM_WORLD);
    if(nnodes > 0){
      input.nodes.temp.resize(nnodes);
      status = MPI_Bcast(&input.nodes.temp[0],     nnodes,    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    
    status = MPI_Bcast(&nnodes,     1,    MPI_INT, 0, MPI_COMM_WORLD);
    if(nnodes > 0){
      input.nodes.v.resize(nnodes);
      status = MPI_Bcast(&input.nodes.v[0],     nnodes,    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    
    status = MPI_Bcast(&nnodes,     1,    MPI_INT, 0, MPI_COMM_WORLD);
    if(nnodes > 0){
      input.nodes.vturb.resize(nnodes);
      status = MPI_Bcast(&input.nodes.vturb[0],     nnodes,    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    
    status = MPI_Bcast(&nnodes,     1,    MPI_INT, 0, MPI_COMM_WORLD);
    if(nnodes > 0){
      input.nodes.bl.resize(nnodes);
      status = MPI_Bcast(&input.nodes.bl[0],     nnodes,    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    
    status = MPI_Bcast(&nnodes,     1,    MPI_INT, 0, MPI_COMM_WORLD);
    if(nnodes > 0){
      input.nodes.bh.resize(nnodes);
      status = MPI_Bcast(&input.nodes.bh[0],     nnodes,    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    
    status = MPI_Bcast(&nnodes,     1,    MPI_INT, 0, MPI_COMM_WORLD);
    if(nnodes > 0){
      input.nodes.azi.resize(nnodes);
      status = MPI_Bcast(&input.nodes.azi[0],     nnodes,    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    
    // total number of nodes
    status = MPI_Bcast(&input.nodes.nnodes,     13,    MPI_INT, 0, MPI_COMM_WORLD);
    
    // Variables to Invert
    status = MPI_Bcast(input.nodes.toinv,     8,    MPI_INT, 0, MPI_COMM_WORLD);
    
    
    // Send type of nodes
    input.nodes.ntype.resize(input.nodes.nnodes);
    status = MPI_Bcast(&input.nodes.ntype[0],   input.nodes.nnodes,    MPI_INT, 0, MPI_COMM_WORLD);

    input.nodes.bound = input.boundary;

  }

  
}
void comm_master_pack_data(iput_t &input, mat<double> &obs, mat<double> &model,
			   unsigned long &ipix, int proc, mdepthall_t &m, int cgrad,
			   int action){
  string inam = "comm_pack_data: ";
  double dum[m.ndep];
  memset(dum, 0, m.ndep * sizeof(double));

  // Buffer to pack the data
  //char buffer[input.buffer_size];
  char *buffer = new char [input.buffer_size]; 
  MPI_Request request = {};
  int status;
  int pos = 0;
  int xx, yy;
  unsigned long len = 0;
  //
  unsigned long ntot = input.nx * input.ny;

  status = MPI_Pack(&action, 1,     MPI_INT, &buffer[0], input.buffer_size, &pos, MPI_COMM_WORLD);
  //
  int init = ipix;
  int end  = min(ipix + input.npack-1, ntot-1); 
  int nPacked = end - init + 1;
  //
  status = MPI_Pack(&nPacked, 1,     MPI_INT, &buffer[0], input.buffer_size, &pos, MPI_COMM_WORLD);
  //
  switch(input.mode)
    {
    case 1:
      {
	/* ---  Pack all pixels --- */
	comm_get_xy(ipix, model.size(1), yy, xx);
	
	status = MPI_Pack(&ipix  , 1,     MPI_INT, &buffer[0], input.buffer_size,
			  &pos, MPI_COMM_WORLD);
	status = MPI_Pack(&xx    , 1,     MPI_INT, &buffer[0], input.buffer_size,
			  &pos, MPI_COMM_WORLD);
	status = MPI_Pack(&yy    , 1,     MPI_INT, &buffer[0], input.buffer_size,
			  &pos, MPI_COMM_WORLD);
	
	/* --- Pack model and obs --- */
	int len = nPacked * input.ns * input.nw_tot;
	status = MPI_Pack(&obs(yy,xx,0,0), len,  MPI_DOUBLE, &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);
	
	len = nPacked * input.npar;
	status = MPI_Pack(&model(yy,xx,0), len,  MPI_DOUBLE, &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);
	

	/* --- Pack full stratified atmos ---*/
	len = m.ndep * nPacked * 12;
	status = MPI_Pack(&m.cub(yy,xx,0,0), len, MPI_DOUBLE, &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);
	
	len = nPacked;
	status = MPI_Pack(&m.boundary(yy,xx), len, MPI_DOUBLE, &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);
	//
	ipix += nPacked; // Increase the pixel count

	// --- Instrumental profiles --- //
	{
	  
	  vector<double> ires = packInstrumentalData(input, yy, xx);
	  len = unsigned(ires.size());
	  status = MPI_Pack(&ires[0], len, MPI_DOUBLE, &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);
	}
	break;
      }
    case 2:
      { // Synthesis


	/* ---  Pack all pixels --- */
	comm_get_xy(ipix, model.size(1), yy, xx);
	status = MPI_Pack(&ipix  , 1,     MPI_INT, &buffer[0], input.buffer_size,
			  &pos, MPI_COMM_WORLD);
	status = MPI_Pack(&xx    , 1,     MPI_INT, &buffer[0], input.buffer_size,
			  &pos, MPI_COMM_WORLD);
	status = MPI_Pack(&yy    , 1,     MPI_INT, &buffer[0], input.buffer_size,
			  &pos, MPI_COMM_WORLD);

	
	/* --- Pack full stratified atmos ---*/
	len = m.ndep * nPacked * 12;
	status = MPI_Pack(&m.cub(yy,xx,0,0), len, MPI_DOUBLE, &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);

	

	/* --- Increment pixel count --- */
	ipix += nPacked;
	// --- Instrumental profiles --- //
	{
	  vector<double> ires = packInstrumentalData(input, yy, xx);
	  len = unsigned(ires.size());
	  status = MPI_Pack(&ires[0], len, MPI_DOUBLE, &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);
	}
	break;
      }
    case 3:
      { //Synthesis + Derivatives
	comm_get_xy(ipix, model.size(1), yy, xx);
	
	/* --- Pack data: ipix, xx, yy, compute_grad (?) --- */
	int pint[4] = {(int)ipix, (int)xx, (int)yy, (int)cgrad};
	status = MPI_Pack(&pint[0]       , 4,     MPI_INT,    &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);
	status = MPI_Pack(&input.dpar    , 1,     MPI_DOUBLE, &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);
	
	/* --- Pack node values --- */
	len = input.npar * nPacked;
	status = MPI_Pack(&model(yy,xx,0), len, MPI_DOUBLE, &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);

	/* --- Pack full stratified atmos ---*/
	len = m.ndep * nPacked * 12;
	status = MPI_Pack(&m.cub(yy,xx,0,0), len, MPI_DOUBLE, &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);

	len = nPacked;
	status = MPI_Pack(&m.boundary(yy,xx), len, MPI_DOUBLE, &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);
	
	/* --- Increase the count --- */
	ipix += nPacked; // Increase the pixel count
       
	break;
      }
    case 4:
      { // synthesis + derivatives at all heights
	comm_get_xy(ipix, model.size(1), yy, xx);
	int pint[3] = {(int)ipix, (int)xx, (int)yy};
	status = MPI_Pack(&pint[0]       , 3,     MPI_INT,    &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);

	/* --- pack ful model --- */
	len = input.ndep * 12*nPacked;
	status = MPI_Pack(&m.cub(yy,xx,0,0), len, MPI_DOUBLE, &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);
	ipix += nPacked; // Increase the pixel count

	// --- Instrumental profiles --- //
	{
	  vector<double> ires = packInstrumentalData(input, yy, xx);
	  len = unsigned(ires.size());
	  status = MPI_Pack(&ires[0], len, MPI_DOUBLE, &buffer[0],
			  input.buffer_size, &pos, MPI_COMM_WORLD);
	}
	
	break;
      }
    } // Switch case

  
  /* ---  Send data to slave --- */
  
  status = MPI_Send(&buffer[0], input.buffer_size, MPI_PACKED, proc, 1, MPI_COMM_WORLD);

  delete [] buffer;

}

void comm_slave_unpack_data(iput_t &input, int &action, mat<double> &obs, mat<double> &pars, vector<mdepth_t> &m, int &cgrad){
  string inam = "comm_slave_unpack_data: ";

  //char buffer[input.buffer_size];
  char *buffer = NULL;

  int  pos = 0, tag = 0;
  MPI_Status ierr = {}, stat = {};
  int nPacked, nbuff=0, status = 0;
  unsigned long len;
  

  /* --- Get DATA from master ---*/
  
  //int status = MPI_Recv(buffer, input.buffer_size, MPI_PACKED, 0, 1, MPI_COMM_WORLD, &ierr);

  while(1){
    int from_process = MPI_ANY_SOURCE, msg_tag = 1, flag = 0;
    MPI_Iprobe(from_process, msg_tag, MPI_COMM_WORLD, &flag, &stat);
    
    if(flag > 0){
      MPI_Get_count(&stat, MPI_CHAR, &nbuff);
      from_process = stat.MPI_SOURCE;
      tag = stat.MPI_TAG;
      buffer = new char [nbuff]();
      status = MPI_Recv(buffer, nbuff, MPI_PACKED, from_process, 1, MPI_COMM_WORLD, &ierr);
      break;
    }

    std::this_thread::sleep_for(std::chrono::microseconds(1000)); // Avoid polling all the time!

  }
    

  /* --- Unpack action --- */
  
  status = MPI_Unpack(buffer, input.buffer_size, &pos, &action, 1, MPI_INT, MPI_COMM_WORLD );
  
  if(action == 1){

    // Check mode and do whatever is required
    switch(input.mode)
      {
      case 1:
	{
	  // Number of packed pixels
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &nPacked, 1,
			      MPI_INT, MPI_COMM_WORLD );
	  

	  input.chi.resize(nPacked);
	  input.nPacked = nPacked;
	  
	  // Allocate arrays for data
	  obs.set({nPacked, input.nw_tot, input.ns});
	  pars.set({nPacked, input.npar});
	  m.resize(nPacked);
	  int ndep = input.ndep;
	  for(auto &it: m) {
	    it.setsize(ndep);
	    it.zero();
	  }
	  
	  // unpack ints
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &input.ipix, 1,
			      MPI_INT, MPI_COMM_WORLD );
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &input.xx,   1,
			      MPI_INT, MPI_COMM_WORLD );
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &input.yy,   1,
			      MPI_INT, MPI_COMM_WORLD );
	  
	  // Unpack arrays
	  len = nPacked * input.nw_tot * input.ns;
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &obs.d[0], len,
			      MPI_DOUBLE, MPI_COMM_WORLD );
	  
	  len = nPacked * input.npar;
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &pars.d[0],  len,
			      MPI_DOUBLE, MPI_COMM_WORLD );
	  
	  
	   /* --- Unpack model atmosphere --- */
	  len = 12 * input.ndep;
	  for(auto &it: m){
	    it.setsize(input.ndep);
	    status = MPI_Unpack(buffer, input.buffer_size, &pos, &it.cub.d[0], len,
				MPI_DOUBLE, MPI_COMM_WORLD );

	  }

	  for(auto &it: m)
	    status = MPI_Unpack(buffer, input.buffer_size, &pos, &it.bound_val, (int)1,
				MPI_DOUBLE, MPI_COMM_WORLD );
	  


	  {
	    double bla = 0.0;
	    status = MPI_Unpack(buffer, input.buffer_size, &pos, &bla, (int)1,
				MPI_DOUBLE, MPI_COMM_WORLD );
	    long ndata = long(bla+0.1);
	    
	    if(ndata > 0){
	      std::vector<double> dat(ndata-1, 0.0);
	      status = MPI_Unpack(buffer, input.buffer_size, &pos, &dat[0], (int)ndata-1,
				  MPI_DOUBLE, MPI_COMM_WORLD );
	      
	    
	      unpackInstrumentalData(input, dat);
	    }
	  }
	  
	  break;
	}
      case 2:
	  // Number of packed pixels
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &nPacked, 1,
			      MPI_INT, MPI_COMM_WORLD );
	  

	  input.nPacked = nPacked;
	  

	  
	  // unpack ints
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &input.ipix, 1,
			      MPI_INT, MPI_COMM_WORLD );
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &input.xx,   1,
			      MPI_INT, MPI_COMM_WORLD );
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &input.yy,   1,
			      MPI_INT, MPI_COMM_WORLD );
	  
	  
	  /* --- Unpack model atmosphere --- */
	  
	  len = 12 * input.ndep;
	  m.resize(nPacked);
	  obs.set({nPacked, input.nw_tot, input.ns});
	  
	  for(auto &it: m){
	    it.setsize(input.ndep);
	    status = MPI_Unpack(buffer, input.buffer_size, &pos, &it.cub.d[0], len,
				MPI_DOUBLE, MPI_COMM_WORLD );
	  }

	  {
	    double bla = 0.0;
	    status = MPI_Unpack(buffer, input.buffer_size, &pos, &bla, (int)1,
				MPI_DOUBLE, MPI_COMM_WORLD );
	    long ndata = long(bla+0.1);
	    
	    if(ndata > 0){
	      std::vector<double> dat(ndata-1, 0.0);
	      status = MPI_Unpack(buffer, input.buffer_size, &pos, &dat[0], (int)ndata-1,
				  MPI_DOUBLE, MPI_COMM_WORLD );
	      
	      unpackInstrumentalData(input, dat);
	    }
	  }
	break;
      case 3: // Synthesize + derivatives
	{
	  /* --- Number of packed pixels --- */
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &nPacked, 1, MPI_INT, MPI_COMM_WORLD );
	  input.nPacked = nPacked;

	  
	  /* --- Allocate arrays for data --- */
	  obs.set({nPacked, input.nw_tot, input.ns});
	  obs.zero();
	  pars.set({nPacked, input.npar});
	  m.resize(nPacked);
	  
	  
	  /* --- unpack scalars --- */
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &input.ipix, 1,
			      MPI_INT, MPI_COMM_WORLD );
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &input.xx,   1,
			      MPI_INT, MPI_COMM_WORLD );
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &input.yy,   1,
			      MPI_INT, MPI_COMM_WORLD );
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &cgrad,      1,
			      MPI_INT, MPI_COMM_WORLD );
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &input.dpar, 1,
			      MPI_DOUBLE, MPI_COMM_WORLD );
	  
	  /* ---  Model parameters --- */
	  len = input.npar*nPacked;
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &pars.d[0], len,
			      MPI_DOUBLE, MPI_COMM_WORLD );

	  
	  /* --- Unpack model atmosphere --- */
	  len = 12 * input.ndep;
	  for(auto &it: m){
	    it.setsize(input.ndep);
	    status = MPI_Unpack(buffer, input.buffer_size, &pos, &it.cub.d[0], len,
				MPI_DOUBLE, MPI_COMM_WORLD );
	  }

	  for(auto &it: m)
	    status = MPI_Unpack(buffer, input.buffer_size, &pos, &it.bound_val, (int)1,
				MPI_DOUBLE, MPI_COMM_WORLD );
	    
	  
	  break;
	}
      case 4:
	status = MPI_Unpack(buffer, input.buffer_size, &pos, &nPacked, 1, MPI_INT, MPI_COMM_WORLD );
	input.nPacked = nPacked;
	m.resize(nPacked);

	obs.set({nPacked, input.nw_tot, input.ns});
	status = MPI_Unpack(buffer, input.buffer_size, &pos, &input.ipix, 1,
			    MPI_INT, MPI_COMM_WORLD );
	status = MPI_Unpack(buffer, input.buffer_size, &pos, &input.xx,   1,
			    MPI_INT, MPI_COMM_WORLD );
	status = MPI_Unpack(buffer, input.buffer_size, &pos, &input.yy,   1,
			    MPI_INT, MPI_COMM_WORLD );

	len = 12 * input.ndep;
	for(auto &it: m){
	  it.setsize(input.ndep);
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &it.cub.d[0], len,
			      MPI_DOUBLE, MPI_COMM_WORLD );
	}
	{
	  double bla = 0.0;
	  status = MPI_Unpack(buffer, input.buffer_size, &pos, &bla, (int)1,
			      MPI_DOUBLE, MPI_COMM_WORLD );
	  long ndata = long(bla+0.1);
	  
	  if(ndata > 0){
	    std::vector<double> dat(ndata-1, 0.0);
	    status = MPI_Unpack(buffer, input.buffer_size, &pos, &dat[0], (int)ndata-1,
				MPI_DOUBLE, MPI_COMM_WORLD );
	    
	    unpackInstrumentalData(input, dat);
	  }
	}
	
	break;
      }
  } // action = 1
    
    
  delete [] buffer;
}


int comm_master_unpack_data(int &iproc, iput_t input, mat<double> &obs, mat<double> &pars, mat<double> &chi2, unsigned long &irec, mat<double> &dsyn, int cgrad, mdepthall_t &m){
  
  // char buffer[input.buffer_size];
  char *buffer;// = new char [input.buffer_size1];
  
  int  pos = 0;
  int nPacked = 0, nbuff=0, tag = 0;
  int status = 0;
  MPI_Status ierr = {}, stat = {};
  unsigned long len;
  int pix, xx, yy, pixout;

  // Get buffer from the slave
  //status = MPI_Recv(buffer, input.buffer_size1, MPI_PACKED, MPI_ANY_SOURCE,
  //		    MPI_ANY_TAG, MPI_COMM_WORLD, &ierr);

  while(1){
    int from_process = MPI_ANY_SOURCE, msg_tag = MPI_ANY_TAG, flag = 0;
    MPI_Iprobe(from_process, msg_tag, MPI_COMM_WORLD, &flag, &stat);
    
    if(flag > 0){
      MPI_Get_count(&stat, MPI_CHAR, &nbuff);
      from_process = stat.MPI_SOURCE;
      tag = stat.MPI_TAG;
      
      buffer = new char [nbuff]();
      status = MPI_Recv(buffer, nbuff, MPI_PACKED, from_process, tag, MPI_COMM_WORLD, &ierr);
      break;
    }

    std::this_thread::sleep_for(std::chrono::microseconds(40)); // Avoid polling all the time
    
  }


  
  // Get nPacked
  status = MPI_Unpack(buffer, input.buffer_size1, &pos, &nPacked, 1, MPI_INT,
		      MPI_COMM_WORLD );

  // Get iproc
  status = MPI_Unpack(buffer, input.buffer_size1, &pos, &iproc, 1, MPI_INT,
		      MPI_COMM_WORLD );

  switch(input.mode)
    {
    case 1:
      {
	// Unpack pixel data
	//for(int pp = 0; pp<nPacked;pp++){
	status = MPI_Unpack(&buffer[0], input.buffer_size1, &pos, &pix, 1,
			    MPI_INT, MPI_COMM_WORLD );
	pixout = pix;
	
	comm_get_xy(pix, input.nx, yy, xx);
	
	status = MPI_Unpack(&buffer[0], input.buffer_size1, &pos, &chi2(yy,xx),
			    nPacked, MPI_DOUBLE, MPI_COMM_WORLD );
	len = input.nw_tot * input.ns * nPacked;
	status = MPI_Unpack(&buffer[0], input.buffer_size1, &pos, &obs(yy,xx,0,0),
			    len, MPI_DOUBLE, MPI_COMM_WORLD );

  
	len = input.npar*nPacked;
	status = MPI_Unpack(&buffer[0], input.buffer_size1, &pos, &pars(yy,xx,0),
			    len, MPI_DOUBLE, MPI_COMM_WORLD );
	irec += nPacked;

	len = input.ndep*12;
	for(int ii=0;ii<nPacked;ii++){
	  comm_get_xy(pix++, input.nx, yy, xx);
	  status = MPI_Unpack(buffer, input.buffer_size1, &pos, &m.cub(yy,xx,0,0),
			      len, MPI_DOUBLE, MPI_COMM_WORLD);
	  
	  status = MPI_Unpack(buffer, input.buffer_size1, &pos, &m.tr_loc(yy,xx),
			      1, MPI_DOUBLE, MPI_COMM_WORLD);
	  status = MPI_Unpack(buffer, input.buffer_size1, &pos, &m.tr_amp(yy,xx),
			      1, MPI_DOUBLE, MPI_COMM_WORLD);
	  
	}
	break;
      }
    case 2:
      {
	
	status = MPI_Unpack(&buffer[0], input.buffer_size1, &pos, &pix, 1,
			    MPI_INT, MPI_COMM_WORLD );
	pixout = pix;
	comm_get_xy(pix, input.nx, yy, xx);
	
	len = input.nw_tot * input.ns * nPacked;
	status = MPI_Unpack(&buffer[0], input.buffer_size1, &pos, &obs(yy,xx,0,0),
			    len, MPI_DOUBLE, MPI_COMM_WORLD );

  
	irec += nPacked;
	
	break;
      }
    case 3:
      {
	// for(int pp = 0; pp<nPacked;pp++){
	int nx;
	nx = (int)obs.size(1);
	
	status = MPI_Unpack(buffer, input.buffer_size1, &pos, &pix, 1, MPI_INT,
			    MPI_COMM_WORLD );
	comm_get_xy(pix, nx, yy, xx);
	pixout = pix;
	len = input.nw_tot * input.ns * nPacked;
	status = MPI_Unpack(buffer, input.buffer_size1, &pos, &obs(yy,xx,0,0), len,
			    MPI_DOUBLE, MPI_COMM_WORLD );

	len = input.ndep*12+2;
	int mypix = pix, iyy=0, ixx=0;
	
	for(int ii = 0; ii< nPacked;ii++){
	  comm_get_xy(mypix++, nx, iyy, ixx);
	  status = MPI_Unpack(buffer, input.buffer_size1, &pos, &m.cub(yy,xx,0,0),
			      len, MPI_DOUBLE, MPI_COMM_WORLD);
	  status = MPI_Unpack(buffer, input.buffer_size1, &pos, &m.tr_loc(yy,xx),
			      1, MPI_DOUBLE, MPI_COMM_WORLD);
	  status = MPI_Unpack(buffer, input.buffer_size1, &pos, &m.tr_amp(yy,xx),
			      1, MPI_DOUBLE, MPI_COMM_WORLD);
	  //cerr<<mypix-1<<" "<<iyy<<" "<<ixx<<endl;
	}
	  
	if(cgrad > 0){
	  len = input.nw_tot * input.ns * input.npar * nPacked;
	  status = MPI_Unpack(buffer, input.buffer_size1, &pos, &dsyn(yy,xx,0,0,0), len,
			      MPI_DOUBLE, MPI_COMM_WORLD);
	}
	
	irec += nPacked;
	
	break;
      }
    case 4:
      {
	int nx;
	nx = (int)obs.size(1);
	status = MPI_Unpack(buffer, input.buffer_size1, &pos, &pix, 1, MPI_INT,
			    MPI_COMM_WORLD );
	comm_get_xy(pix, nx, yy, xx);
	pixout = pix;
	
	len = input.nw_tot * input.ns * nPacked;
	status = MPI_Unpack(buffer, input.buffer_size1, &pos, &obs(yy,xx,0,0), len,
			    MPI_DOUBLE, MPI_COMM_WORLD );
	
	len = input.nw_tot * input.ns * input.ndep * input.nresp * nPacked;
	status = MPI_Unpack(buffer, input.buffer_size1, &pos, &dsyn(yy,xx,0,0,0,0), len,
			    MPI_DOUBLE, MPI_COMM_WORLD);
	
	irec += nPacked;

	break;
      }
    default:
      break;
    }
  

  delete [] buffer;
  
  return pixout;
}


void comm_slave_pack_data(iput_t &input, mat<double> &obs, mat<double> &pars, mat<double> &dobs, int cgrad, vector<mdepth_t> &m){

  char *buffer = new char [input.buffer_size1];

  int  pos = 0;
  int nPacked = input.nPacked;
  unsigned long len;
  //
  int status = MPI_Pack(&nPacked, 1,     MPI_INT, &buffer[0], input.buffer_size1,
			&pos, MPI_COMM_WORLD);
  status = MPI_Pack(&input.myrank, 1,     MPI_INT, &buffer[0], input.buffer_size1,
		    &pos, MPI_COMM_WORLD);

  switch(input.mode)
    {
    case 1:
      status = MPI_Pack(&input.ipix, 1,   MPI_INT, &buffer[0], input.buffer_size1,
			&pos, MPI_COMM_WORLD);
      
      status = MPI_Pack(&input.chi[0], nPacked, MPI_DOUBLE, &buffer[0], input.buffer_size1,
			&pos, MPI_COMM_WORLD);
      //   for(int pp = 0; pp<nPacked; pp++){
	
      len = input.nw_tot*input.ns*nPacked;
      status = MPI_Pack(&obs.d[0], len,  MPI_DOUBLE, &buffer[0], input.buffer_size1,
			&pos, MPI_COMM_WORLD);

      
      len = input.npar * nPacked;
      status = MPI_Pack(&pars.d[0], len,  MPI_DOUBLE, &buffer[0], input.buffer_size1,
			&pos, MPI_COMM_WORLD);
      
      // model
      len = input.ndep*12+2;
      for(int pp=0; pp<nPacked; pp++){
	vector<double> vec = m[pp].model2vector();
	status = MPI_Pack(&vec[0], len,  MPI_DOUBLE, &buffer[0], input.buffer_size1,
			  &pos, MPI_COMM_WORLD);
      }
      break;
      
    case 2:

      status = MPI_Pack(&input.ipix, 1,   MPI_INT, &buffer[0], input.buffer_size1,
			&pos, MPI_COMM_WORLD);
      	
      len = input.nw_tot*input.ns*nPacked;
      status = MPI_Pack(&obs.d[0], len,  MPI_DOUBLE, &buffer[0], input.buffer_size1,
			&pos, MPI_COMM_WORLD);

      
      break;
      
    case 3:
      //  for(int pp = 0; pp<nPacked; pp++){
      status = MPI_Pack(&input.ipix, 1,   MPI_INT, &buffer[0], input.buffer_size1, &pos, MPI_COMM_WORLD);
      
      // Synthetic profiles
      len = input.nw_tot*input.ns*nPacked;
      status = MPI_Pack(&obs.d[0], len,  MPI_DOUBLE, &buffer[0], input.buffer_size1, &pos, MPI_COMM_WORLD);

      // gas pressure
      len = input.ndep*12+2;
      for(int pp=0; pp<nPacked; pp++){
	vector<double> vec = m[pp].model2vector();
	status = MPI_Pack(&vec[0], len,  MPI_DOUBLE, &buffer[0], input.buffer_size1,
			  &pos, MPI_COMM_WORLD);
      }
      
      // Derivatives
      if(cgrad > 0){
	len = input.nw_tot*input.ns*input.npar*nPacked;
	status = MPI_Pack(&dobs.d[0], len,  MPI_DOUBLE, &buffer[0], input.buffer_size1, &pos, MPI_COMM_WORLD);
      }
      break;

    case 4:
      {
	status = MPI_Pack(&input.ipix, 1,   MPI_INT, &buffer[0], input.buffer_size1, &pos, MPI_COMM_WORLD);
	// Synthetic profiles
	len = input.nw_tot*input.ns*nPacked;
	status = MPI_Pack(&obs.d[0], len,  MPI_DOUBLE, &buffer[0], input.buffer_size1, &pos, MPI_COMM_WORLD);
	len = input.nw_tot*input.ns*input.ndep*input.nresp*nPacked;
	status = MPI_Pack(&dobs.d[0], len,  MPI_DOUBLE, &buffer[0], input.buffer_size1, &pos, MPI_COMM_WORLD);

	
	break;
      
      }
    default:
      break;
    } // Switch case

  
  /* --- Send data to master --- */
  
  status = MPI_Send(&buffer[0], input.buffer_size1, MPI_PACKED, 0, 3, MPI_COMM_WORLD);

  delete [] buffer;

}



void comm_kill_slaves(iput_t &input, int nprocs){
  char buffer[input.buffer_size];
  int  pos = 0;
  int action = 0;
  string inam = "comm_kill_slaves: ";
  
  // Pack kill command
  int status = MPI_Pack(&action, 1,     MPI_INT, &buffer[0], input.buffer_size, &pos, MPI_COMM_WORLD);
  
  // Send kill command to slaves
  for(int ss = 1; ss<nprocs; ss++) status = MPI_Send(&buffer[0], input.buffer_size, MPI_PACKED, ss, 1, MPI_COMM_WORLD);

  cout << " "<<endl;
  cout << input.myid << inam << "Killing slaves" << endl;
}


void comm_send_weights(iput_t &input, mat<double> &w){

  if(input.myrank > 0) w.set({input.nw_tot, input.ns});  
  int len = (int)w.d.size();

  MPI_Bcast(&w.d[0],     len,    MPI_DOUBLE, 0, MPI_COMM_WORLD);

  
}
