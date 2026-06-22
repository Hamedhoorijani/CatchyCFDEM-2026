/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    UPDATE:
    25-07-2024 - CHECK THE CONDITIONS AND ACTIVATE THE Conservation law function

    Copyright 2024-     LCT - UGENT, Ghent
    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include "fix_heat_gran_electric.h"
#include "domain.h"
#include "comm.h"

#include "atom.h"
#include "compute_pair_gran_local.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "force.h"
#include "math_extra.h"
#include "math_const.h"
#include "properties.h"
#include "modify.h"
#include "neigh_list.h"
#include "pair_gran.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <iomanip>
#include "update.h"
#include <thread>
#include <fstream>
#include <cstdlib>
#include "memory.h"
#include "mpi.h"
#include <queue>
#include "random_mars.h"
#include "neighbor.h"
#include <cassert>
#include <ctime>
#include <omp.h>
#include <utility>

#include <unordered_map>
#include <unordered_set>
#include <numeric>

#include <limits>  // For numeric limits


using namespace LAMMPS_NS;
using namespace FixConst;

using MathConst::MY_4PI;
using namespace MathExtra;

const double SMALL = 1.0E-30;

// contact area modes
enum{ CONDUCTION_CONTACT_AREA_OVERLAP,
      CONDUCTION_CONTACT_AREA_CONSTANT,
      CONDUCTION_CONTACT_AREA_PROJECTION};

/* ---------------------------------------------------------------------- */

FixHeatGranElectric::FixHeatGranElectric(class LAMMPS *lmp, int narg, char **arg) :
  FixHeatGran(lmp, narg, arg),
  fix_eConductivity_(0),
  eConductivity_(0),
  store_contact_data_(false),
  fix_conduction_contact_area_(0),
  fix_n_conduction_contacts_(0),
  fix_wall_heattransfer_coeff_(0),
  fix_wall_temperature_(0),
  conduction_contact_area_(0),
  n_conduction_contacts_(0),
  wall_heattransfer_coeff_(0),
  wall_temp_(0),
  fix_wall_electricPotential_(0),       //electric part
  wall_ep_(0),
  fix_eP(NULL),
  fix_eP_old(NULL),
  fix_tot_I(NULL),
  eP(0),
  eP_old(0),
  tot_I(0),
  fix_total_C(NULL),
  fix_countCheck(NULL),
  eP0(0), //initial_electricPotential
  eP_needs_init_(false),
  qj(0),   // fix to save joule heating
  qr(0),   // fix to save radiation
  hasEPArgs(true),
  electricMode_(false),               //electric flag
  area_calculation_mode_(CONDUCTION_CONTACT_AREA_OVERLAP),
  fixed_contact_area_(0.),
  area_correction_flag_(0),
  joule_heat_startTime(0),
  temp_control_flag_(false),		//flag temperature control
  temp_control_(0),					//control temperature
  aspot_fraction_global_(1.0),
  fix_aspotFrac_(NULL),
  aspotFrac_(NULL),
  Cmin_branch_(0.0),     // 0 → no clipping; you can set e.g. 1e-20 S if desired
  sigma0_(NULL),
  alpha_(NULL),
  Tref_(293.15),        // 20°C default
  sigma_floor_(1e-12),  // safe tiny floor
  TB(0), //Background Temp
  Sigma(5.670373E-8),
  electric_log_stride_(1),	//initialized to export electric results every step
  cutGhost(0.10),  //hardcoded cutoff distance for radiation
  cutGhostsq(0.01),
  radiation_flag_(false),
  compID(NULL), // ----
  deltan_ratio_(0)
{
  iarg_ = 5;
  bool hasargs = true;
  if(strcmp(style,"heat/gran/electric") == 0) electricMode_ = true;
  int rank;
  MPI_Comm_rank(world, &rank);

  while(iarg_ < narg && hasargs)
  {
    hasargs = false;
    if(strcmp(arg[iarg_],"contact_area") == 0) {

      if(strcmp(arg[iarg_+1],"overlap") == 0)
        area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_OVERLAP;
      else if(strcmp(arg[iarg_+1],"projection") == 0)
        area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_PROJECTION;
      else if(strcmp(arg[iarg_+1],"constant") == 0)
      {
        if (iarg_+3 > narg)
            error->fix_error(FLERR,this,"not enough arguments for keyword 'contact_area constant'");
        area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_CONSTANT;
        fixed_contact_area_ = force->numeric(FLERR,arg[iarg_+2]);
        if (fixed_contact_area_ <= 0.)
            error->fix_error(FLERR,this,"'contact_area constant' value must be > 0");
        iarg_++;
      }
      else error->fix_error(FLERR,this,"expecting 'overlap', 'projection' or 'constant' after 'contact_area'");
      iarg_ += 2;
      hasargs = true;
    } else if(strcmp(arg[iarg_],"area_correction") == 0) {
      if (iarg_+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'area_correction'");
      if(strcmp(arg[iarg_+1],"yes") == 0)
        area_correction_flag_ = 1;
      else if(strcmp(arg[iarg_+1],"no") == 0)
        area_correction_flag_ = 0;
      else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'area_correction'");
      iarg_ += 2;
      hasargs = true;
    } else if(strcmp(arg[iarg_],"store_contact_data") == 0) {
      if (iarg_+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'store_contact_data'");
      if(strcmp(arg[iarg_+1],"yes") == 0)
        store_contact_data_ = true;
      else if(strcmp(arg[iarg_+1],"no") == 0)
        store_contact_data_ = false;
      else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'store_contact_data'");
      iarg_ += 2;
      hasargs = true;
    } else if(strcmp(arg[iarg_],"initial_electricPotential") == 0 && electricMode_) {
        if (iarg_+1 > narg){
          error->fix_error(FLERR,this,"not enough arguments for keyword 'initial_electricPotential'");
        }
        eP0 = atof(arg[iarg_+1]); //initial electricPotential 
        iarg_ += 2;
        hasargs = true;
    } else if(strcmp(arg[iarg_],"start_time") == 0 && electricMode_) {
        if (iarg_+2 > narg)      error->fix_error(FLERR,this,"not enough arguments for keyword 'start_time'");
        joule_heat_startTime = atof(arg[iarg_+1]);
        if (joule_heat_startTime < 0.)  error->fix_error(FLERR,this,"'start_time' value must be >= 0");
        iarg_ += 2;
        hasargs = true;
    } else if(strcmp(arg[iarg_],"temp_control") == 0 && electricMode_) {
        if (iarg_+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'temp_control'");
        if(strcmp(arg[iarg_+1],"yes") == 0){
          if (iarg_+3 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'temp_control yes'");
          temp_control_flag_ = true;
          temp_control_ = atof(arg[iarg_+2]);
          iarg_ += 3;
		} else if(strcmp(arg[iarg_+1],"no") == 0){
          temp_control_flag_ = false;
          iarg_ += 2;
		} else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'radiation'");
        hasargs = true;
	} else if (strcmp(arg[iarg_],"aspot_fraction")==0) {
		if (iarg_+2 > narg) error->fix_error(FLERR,this,"aspot_fraction needs a scalar in (0,1]");
		aspot_fraction_global_ = atof(arg[iarg_+1]);
		if (aspot_fraction_global_ <= 0.0 || aspot_fraction_global_ > 1.0)
			error->fix_error(FLERR,this,"aspot_fraction must be in (0,1]");
		iarg_ += 2; 
		hasargs = true;
	} else if (strcmp(arg[iarg_],"Gmin")==0) {
		if (iarg_+2 > narg) error->fix_error(FLERR,this,"Gmin needs a scalar ≥ 0");
		Cmin_branch_ = std::max(0.0, atof(arg[iarg_+1]));
		iarg_ += 2; 
		hasargs = true;
	} else if(strcmp(arg[iarg_],"export_every") == 0) {
		if (iarg_+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'export_every'");
		const int stride = atoi(arg[iarg_+1]);
		if (stride <= 0) error->fix_error(FLERR,this,"'export_every' must be a positive integer");
		electric_log_stride_ = stride;
		iarg_ += 2;
		hasargs = true;
	} else if(strcmp(arg[iarg_],"radiation") == 0 && electricMode_) {
        if (iarg_+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'radiation'");
        if(strcmp(arg[iarg_+1],"yes") == 0){
          if (iarg_+3 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'radiation yes'");
          radiation_flag_ = true;
          TB = atof(arg[iarg_+2]);
          iarg_ += 3;
        } else if(strcmp(arg[iarg_+1],"no") == 0){
          radiation_flag_ = false;
          iarg_ += 2;
        } else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'radiation'");
        hasargs = true;
    } else if(strcmp(style,"heat/gran/electric") == 0)   error->fix_error(FLERR,this,"unknown keyword");
  }

  if(CONDUCTION_CONTACT_AREA_OVERLAP != area_calculation_mode_ && 1 == area_correction_flag_)
    error->fix_error(FLERR,this,"can use 'area_correction' only for 'contact_area = overlap'");

  // Deterministic seed for RanMars
  std::string seed_full = std::to_string(86028157);
  RGen = new RanMars(lmp, seed_full.c_str());

  // for optimization of trace() preallocate these
  raypoint = new double[3];
}

/* ---------------------------------------------------------------------- */

FixHeatGranElectric::~FixHeatGranElectric()
{
  if (eConductivity_)    delete []eConductivity_;
  if (raypoint) delete[] raypoint;
  if (stencilLength) delete[] stencilLength;
  if (binStencildx) delete[] binStencildx;
}

/* ---------------------------------------------------------------------- */

void FixHeatGranElectric::post_create()
{
  FixHeatGran::post_create();

  // register contact storage
  fix_conduction_contact_area_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("contactAreaConduction","property/atom","scalar",0,0,this->style,false));
  if(!fix_conduction_contact_area_ && store_contact_data_)
  {
    const char* fixarg[10];
    fixarg[0]="contactAreaConduction";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="contactAreaConduction";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fix_conduction_contact_area_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_n_conduction_contacts_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("nContactsConduction","property/atom","scalar",0,0,this->style,false));
  if(!fix_n_conduction_contacts_ && store_contact_data_)
  {
    const char* fixarg[10];
    fixarg[0]="nContactsConduction";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="nContactsConduction";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fix_n_conduction_contacts_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_qj = static_cast<FixPropertyAtom*>(modify->find_fix_property("qj","property/atom","scalar",0,0,this->style,false));
  if(!fix_qj)
  {
    const char* fixarg[10];
    fixarg[0]="qj";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="qj";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_qj = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
  
  fix_tot_I = static_cast<FixPropertyAtom*>(modify->find_fix_property("tot_I","property/atom","scalar",0,0,this->style,false));
  if(!fix_tot_I)
  {
    const char* fixarg[10];
    fixarg[0]="tot_I";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="tot_I";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_tot_I = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_qj_pp = static_cast<FixPropertyAtom*>(modify->find_fix_property("qj_pp","property/atom","scalar",0,0,this->style,false));
  if(!fix_qj_pp)
  {
    const char* fixarg[10];
    fixarg[0]="qj_pp";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="qj_pp";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_qj_pp = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_qj_pw = static_cast<FixPropertyAtom*>(modify->find_fix_property("qj_pw","property/atom","scalar",0,0,this->style,false));
  if(!fix_qj_pw)
  {
    const char* fixarg[10];
    fixarg[0]="qj_pw";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="qj_pw";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_qj_pw = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_qr = static_cast<FixPropertyAtom*>(modify->find_fix_property("qr","property/atom","scalar",0,0,this->style,false));
  if(!fix_qr)
  {
    const char* fixarg[10];
    fixarg[0]="qr";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="qr";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_qr = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_wall_heattransfer_coeff_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("wallHeattransferCoeff","property/atom","scalar",0,0,this->style,false));
  if(!fix_wall_heattransfer_coeff_ && store_contact_data_)
  {
    const char* fixarg[10];
    fixarg[0]="wallHeattransferCoeff";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="wallHeattransferCoeff";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fix_wall_heattransfer_coeff_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
  fix_wall_temperature_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("wallTemp","property/atom","scalar",0,0,this->style,false));
  if(!fix_wall_temperature_ && store_contact_data_)
  {
    const char* fixarg[10];
    fixarg[0]="wallTemp";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="wallTemp";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fix_wall_temperature_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // register electricalPotential
  fix_eP = static_cast<FixPropertyAtom*>(modify->find_fix_property("eP","property/atom","scalar",1,0,this->style,false));
  if(!fix_eP)
  {
    const char* fixarg[11];
    fixarg[0]="eP";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="eP";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_eP = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    eP_needs_init_ = true; //fresh run 
  } else
  {
  	eP_needs_init_ = false;
  }

  // register deltaE (optional diagnostics)
  fix_deltaE  = static_cast<FixPropertyAtom*>(modify->find_fix_property("deltaE","property/atom","scalar",1,0,this->style,false));
  if(!fix_deltaE )
  {
    const char* fixarg[11];
    fixarg[0]="deltaE";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="deltaE";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_deltaE  = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // register electricalPotential (old)
  fix_eP_old = static_cast<FixPropertyAtom*>(modify->find_fix_property("eP_old","property/atom","scalar",1,0,this->style,false));
  if(!fix_eP_old)
  {
    const char* fixarg[11];
    fixarg[0]="eP_old";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="eP_old";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_eP_old = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_current = static_cast<FixPropertyAtom*>(modify->find_fix_property("current","property/atom","vector",0,0,this->style,false));
  if(!fix_current)
  {
    const char* fixarg[11];
    fixarg[0]="current";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="current";
    fixarg[4]="vector";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fixarg[9]="0.";
    fixarg[10]="0.";
    fix_current = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
  }

  fix_compID = static_cast<FixPropertyAtom*>(modify->find_fix_property("compID","property/atom","scalar",0,0,this->style,false));
  if(!fix_compID)
  {
    const char* fixarg[10];
    fixarg[0]="compID";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="compID";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_compID = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_rhsG = static_cast<FixPropertyAtom*>(modify->find_fix_property("rhsG","property/atom","scalar",0,0,this->style,false));
  if(!fix_rhsG)
  {
    const char* fixarg[10];
    fixarg[0]="rhsG";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="rhsG";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_rhsG = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_diagFlag = static_cast<FixPropertyAtom*>(modify->find_fix_property("diagFlag","property/atom","scalar",0,0,this->style,false));
  if(!fix_diagFlag)
  {
    const char* fixarg[10];
    fixarg[0]="diagFlag";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="diagFlag";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_diagFlag = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

 fix_total_C = static_cast<FixPropertyAtom*>(modify->find_fix_property("total_C","property/atom","scalar",0,0,this->style,false));
  if(!fix_total_C)
  {
    const char* fixarg[10];
    fixarg[0]="total_C";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="total_C";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_total_C = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

 fix_countCheck= static_cast<FixPropertyAtom*>(modify->find_fix_property("countCheck","property/atom","scalar",0,0,this->style,false));
  if(!fix_countCheck)
  {
    const char* fixarg[10];
    fixarg[0]="countCheck";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="countCheck";
    fixarg[4]="scalar";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fix_countCheck = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  if(store_contact_data_ && (!fix_conduction_contact_area_ || !fix_n_conduction_contacts_ ||
                             !fix_wall_heattransfer_coeff_ || !fix_eP || !fix_eP_old ||
                             !fix_wall_temperature_ ||!fix_wall_electricPotential_))
    error->one(FLERR,"internal error");
}

/* ---------------------------------------------------------------------- */

void FixHeatGranElectric::pre_delete(bool unfixflag)
{
  if(cpl && unfixflag) cpl->reference_deleted();
}

/* ---------------------------------------------------------------------- */

int FixHeatGranElectric::setmask()
{
  int mask = FixHeatGran::setmask();
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeatGranElectric::updatePtrs()
{
  FixHeatGran::updatePtrs();
  eP = fix_eP->vector_atom;
  eP_old = fix_eP_old->vector_atom;

  if(store_contact_data_)
  {
    conduction_contact_area_ = fix_conduction_contact_area_->vector_atom;
    n_conduction_contacts_ = fix_n_conduction_contacts_->vector_atom;
    wall_heattransfer_coeff_ = fix_wall_heattransfer_coeff_->vector_atom;
    wall_temp_ = fix_wall_temperature_->vector_atom;
    if(electricMode_) wall_ep_ = fix_wall_electricPotential_->vector_atom;
  }
}

/* ---------------------------------------------------------------------- */

void FixHeatGranElectric::init()
{
  // initialize base class
  FixHeatGran::init();
  const double *Y, *nu, *Y_orig;
  double expo, Yeff_ij, Yeff_orig_ij, ratio;
  int max_type = atom->get_properties()->max_type();

  if (eConductivity_) delete []eConductivity_;
  eConductivity_ = new double[max_type];
  fix_eConductivity_ =
    static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalConductivity","property/global","peratomtype",max_type,0,style));

  if(electricMode_){
    e_cond = static_cast<FixPropertyGlobal*>(modify->find_fix_property("electricalConductivity","property/global","peratomtype",max_type,0,style))->get_values();
    fix_eP = static_cast<FixPropertyAtom*>(modify->find_fix_property("eP","property/atom","scalar",1,0,style));
    fix_eP_old = static_cast<FixPropertyAtom*>(modify->find_fix_property("eP_old","property/atom","scalar",1,0,style));
    fix_qj = static_cast<FixPropertyAtom*>(modify->find_fix_property("qj","property/atom","scalar",1,0,style));
    fix_qr = static_cast<FixPropertyAtom*>(modify->find_fix_property("qr","property/atom","scalar",1,0,style));
    fix_cMatrix = static_cast<FixPropertyAtom*>(modify->find_fix_property("cMatrix","property/atom","scalar",1,0,style));
    fix_pwER = static_cast<FixPropertyAtom*>(modify->find_fix_property("pwER","property/atom","scalar",1,0,style));
    fix_pwEP = static_cast<FixPropertyAtom*>(modify->find_fix_property("pwEP","property/atom","scalar",1,0,style));

    fix_current=static_cast<FixPropertyAtom*>(modify->find_fix_property("current","property/atom","vector",1,0,style));
    fix_compID = static_cast<FixPropertyAtom*>(modify->find_fix_property("compID","property/atom","scalar",1,0,style));
    fix_rhsG = static_cast<FixPropertyAtom*>(modify->find_fix_property("rhsG","property/atom","scalar",1,0,style));
    fix_diagFlag = static_cast<FixPropertyAtom*>(modify->find_fix_property("diagFlag","property/atom","scalar",1,0,style));
	fix_total_C= static_cast<FixPropertyAtom*>(modify->find_fix_property("total_C","property/atom","scalar",1,0,style));
    fix_countCheck= static_cast<FixPropertyAtom*>(modify->find_fix_property("countCheck","property/atom","scalar",1,0,style));

	fix_aspotFrac_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("electricalAspotFraction","property/global","peratomtype", max_type, 0, style));
	aspotFrac_ = fix_aspotFrac_ ? fix_aspotFrac_->get_values() : NULL;

	// Prefer new knobs; fall back to legacy electricalConductivity (e_cond) if σ0 not provided
	if (auto *fxS0 = static_cast<FixPropertyGlobal*>(modify->find_fix_property("electricalSigma0","property/global","peratomtype",max_type,0,style,false))) {
	  sigma0_ = fxS0->get_values();
	} else {
	  sigma0_ = e_cond; // fallback to your existing per-type constant conductivity
	}

	if (auto *fxAl = static_cast<FixPropertyGlobal*>(modify->find_fix_property("electricalAlpha","property/global","peratomtype",max_type,0,style,false))) {
	  alpha_ = fxAl->get_values();
	} else {
	  alpha_ = nullptr; // α defaults to 0
	}
	
	
	// Optional: scalar Tref
	if (auto *fix_Tref = static_cast<FixPropertyGlobal*>(modify->find_fix_property("electricalTref","property/global","scalar",0,0,style,false))) {
	  Tref_ = fix_Tref->compute_scalar();
	}

	// Optional: floor as scalar
	if (auto *fix_floor = static_cast<FixPropertyGlobal*>(modify->find_fix_property("electricalSigmaFloor","property/global","scalar",0,0,style,false))) {
	  sigma_floor_ = std::max(0.0, fix_floor->compute_scalar());
	}

    if(radiation_flag_){
      emissivity_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalEmissivity","property/global","peratomtype",max_type,0,style))->get_values();
    }

    wallfix = NULL;
    int nfix = modify->nfix;
    for (int i = 0; i < nfix; i++) {
        if (strcmp(modify->fix[i]->style, "wall/gran") == 0) {
            wallfix = (FixWallGran *)modify->fix[i];
            break;
        }
    }
    if (wallfix == NULL) {
        error->all(FLERR, "FixHeatGranElectric requires fix wall/gran");
    }
    if (force->newton_pair)
      error->all(FLERR, "FixHeatGranElectric assumes 'newton pair off'");

    // --- initialize eP once from eP0 (warm start) ---
    //fix_eP->set_all(eP0);
    if (eP_needs_init_) {
    	fix_eP->set_all(eP0);
    	eP_needs_init_ = false;   // only once, even across multiple run commands
    }
  }

  // pre-calculate conductivity combinations
  for(int i=1;i< max_type+1; i++)
      for(int j=1;j<max_type+1;j++)
      {
          eConductivity_[i-1] = fix_eConductivity_->compute_vector(i-1);
          if(eConductivity_[i-1] < 0.)
            error->all(FLERR,"Fix heat/gran/electric: Thermal conductivity must not be < 0");
      }

  // contact area correction table
  if(area_correction_flag_)
  {
    if(!force->pair_match("gran",0))
        error->fix_error(FLERR,this,"area correction only works with using granular pair styles");

    expo = 1./pair_gran->stressStrainExponent();

    Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
    nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
    Y_orig = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->get_values();

    static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->new_array(max_type,max_type);

    for(int i = 1; i < max_type+1; i++)
    {
      for(int j = 1; j < max_type+1; j++)
      {
        Yeff_ij      = 1./((1.-pow(nu[i-1],2.))/Y[i-1]     +(1.-pow(nu[j-1],2.))/Y[j-1]);
        Yeff_orig_ij = 1./((1.-pow(nu[i-1],2.))/Y_orig[i-1]+(1.-pow(nu[j-1],2.))/Y_orig[j-1]);
        ratio = pow(Yeff_ij/Yeff_orig_ij,expo);

        static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->array_modify(i-1,j-1,ratio);
      }
    }
    deltan_ratio_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->get_array_modified();
  }
  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixHeatGranElectric::pre_force(int vflag)
{
  // IMPORTANT: zero heatFlux/heatSource/etc. every step
  FixHeatGran::pre_force(vflag);

  if(store_contact_data_)
  {
    fix_wall_heattransfer_coeff_->set_all(0.);
    fix_wall_temperature_->set_all(0.);
    fix_wall_electricPotential_->set_all(0.);
  }

  fix_total_C->set_all(0.);

  fix_qr->set_all(0.);
  fix_rhsG->set_all(0.);
  fix_diagFlag->set_all(0.);
  

  const int nloc = atom->nlocal;
  double *V    = fix_eP->vector_atom;
  double *Vold = fix_eP_old->vector_atom;
  for (int i = 0; i < nloc; ++i) Vold[i] = V[i];
  fix_eP_old->do_forward_comm();   // ghosts consistent while iterating
  
  // clear wall ↔ particle electrical accumulators once per step
  
  
  if (fix_pwER)    fix_pwER->set_all(0.);
  if (fix_cMatrix) fix_cMatrix->set_all(0.);
  if (fix_pwEP)    fix_pwEP->set_all(0.);
  
}

/* ---------------------------------------------------------------------- */

void FixHeatGranElectric::post_force(int vflag)
{
  if(history_flag == 0 && CONDUCTION_CONTACT_AREA_OVERLAP == area_calculation_mode_)
    post_force_eval<0,CONDUCTION_CONTACT_AREA_OVERLAP>(vflag,0);

  if(history_flag == 1 && CONDUCTION_CONTACT_AREA_OVERLAP == area_calculation_mode_)
    post_force_eval<1,CONDUCTION_CONTACT_AREA_OVERLAP>(vflag,0);

  if(history_flag == 0 && CONDUCTION_CONTACT_AREA_CONSTANT == area_calculation_mode_)
    post_force_eval<0,CONDUCTION_CONTACT_AREA_CONSTANT>(vflag,0);

  if(history_flag == 1 && CONDUCTION_CONTACT_AREA_CONSTANT == area_calculation_mode_)
    post_force_eval<1,CONDUCTION_CONTACT_AREA_CONSTANT>(vflag,0);

  if(history_flag == 0 && CONDUCTION_CONTACT_AREA_PROJECTION == area_calculation_mode_)
    post_force_eval<0,CONDUCTION_CONTACT_AREA_PROJECTION>(vflag,0);

  if(history_flag == 1 && CONDUCTION_CONTACT_AREA_PROJECTION == area_calculation_mode_)
    post_force_eval<1,CONDUCTION_CONTACT_AREA_PROJECTION>(vflag,0);
}

/* ---------------------------------------------------------------------- */

void FixHeatGranElectric::cpl_evaluate(ComputePairGranLocal *caller)
{
  if(history_flag == 0 && CONDUCTION_CONTACT_AREA_OVERLAP == area_calculation_mode_)
    post_force_eval<0,CONDUCTION_CONTACT_AREA_OVERLAP>(0,1);
  if(history_flag == 1 && CONDUCTION_CONTACT_AREA_OVERLAP == area_calculation_mode_)
    post_force_eval<1,CONDUCTION_CONTACT_AREA_OVERLAP>(0,1);

  if(history_flag == 0 && CONDUCTION_CONTACT_AREA_CONSTANT == area_calculation_mode_)
    post_force_eval<0,CONDUCTION_CONTACT_AREA_CONSTANT>(0,1);
  if(history_flag == 1 && CONDUCTION_CONTACT_AREA_CONSTANT == area_calculation_mode_)
    post_force_eval<1,CONDUCTION_CONTACT_AREA_CONSTANT>(0,1);

  if(history_flag == 0 && CONDUCTION_CONTACT_AREA_PROJECTION == area_calculation_mode_)
    post_force_eval<0,CONDUCTION_CONTACT_AREA_PROJECTION>(0,1);
  if(history_flag == 1 && CONDUCTION_CONTACT_AREA_PROJECTION == area_calculation_mode_)
    post_force_eval<1,CONDUCTION_CONTACT_AREA_PROJECTION>(0,1);
}

/* ---------------------------------------------------------------------- */

template <int HISTFLAG,int CONTACTAREA>
void FixHeatGranElectric::post_force_eval(int vflag,int cpl_flag)
{
  int myrank;
  MPI_Comm_rank(world,&myrank);

  double erij, ecoi,ecoj;   //electric

  double hc,contactArea,delta_n,flux,dirFlux[3];
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r,tcoi,tcoj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *contact_flag,**first_contact_flag;

  int newton_pair = force->newton_pair;

  if (strcmp(force->pair_style,"hybrid")==0)
    error->warning(FLERR,"Fix heat/gran/electric implementation may not be valid for pair style hybrid");
  if (strcmp(force->pair_style,"hybrid/overlay")==0)
    error->warning(FLERR,"Fix heat/gran/electric implementation may not be valid for pair style hybrid/overlay");

  inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
  firstneigh = pair_gran->list->firstneigh;
  if(HISTFLAG) first_contact_flag = pair_gran->listgranhistory->firstneigh;

  double *radius = atom->radius;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *tag = atom->tag;
  std::vector<int> spanningClusters;

  double localSum = 0;
  for(int i=0; i<nlocal ; i++) localSum+=Temp[i]; 
  double globalSum = 0.0;
  MPI_Allreduce(&localSum,&globalSum,1,MPI_DOUBLE,MPI_SUM,world);
  int natoms = atom->natoms;
  double aveT = 1.0* globalSum/natoms;

  if(radiation_flag_) {
     calculateRadiation();
  }
    
  updatePtrs();

  if(store_contact_data_)
  {
    fix_conduction_contact_area_->set_all(0.);
    fix_n_conduction_contacts_->set_all(0.);
  }

  // loop over neighbors of my atoms (thermal conduction bookkeeping)
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    if(HISTFLAG) contact_flag = first_contact_flag[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

      if(!HISTFLAG)
      {
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radj = radius[j];
        radsum = radi + radj;
      }

      if ((HISTFLAG && contact_flag[jj]) || (!HISTFLAG && (rsq < radsum*radsum))) {  //contact

        if(HISTFLAG)
        {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          radj = radius[j];
          radsum = radi + radj;
          if(rsq >= radsum*radsum) continue;
        }

        r = sqrt(rsq);

        if(CONTACTAREA == CONDUCTION_CONTACT_AREA_OVERLAP)
        {
          if(area_correction_flag_)
          {
            delta_n = radsum - r;
            delta_n *= deltan_ratio_[type[i]-1][type[j]-1];
            r = radsum - delta_n;
          }

          if (r < fmax(radi, radj)) // one sphere is inside the other
          {
              contactArea = fmin(radi,radj);
              contactArea *= contactArea * M_PI;
          }
          else
              contactArea = - M_PI/4.0 * ( (r-radi-radj)*(r+radi-radj)*(r-radi+radj)*(r+radi+radj) )/(r*r);
        } else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_CONSTANT)
              contactArea = fixed_contact_area_;
        else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_PROJECTION)
        {
            double rmax = std::max(radi,radj);
            contactArea = M_PI*rmax*rmax;
        }

        tcoi = eConductivity_[type[i]-1];
        tcoj = eConductivity_[type[j]-1];
        if (tcoi < SMALL_FIX_HEAT_GRAN || tcoj < SMALL_FIX_HEAT_GRAN) hc = 0.;
        else hc = 4.*tcoi*tcoj/(tcoi+tcoj)*sqrt(contactArea);

        flux = (Temp[j]-Temp[i])*hc;

        dirFlux[0] = flux*delx;
        dirFlux[1] = flux*dely;
        dirFlux[2] = flux*delz;
        if(!cpl_flag)
        {
          heatFlux[i] += flux;
          directionalHeatFlux[i][0] += 0.50 * dirFlux[0];
          directionalHeatFlux[i][1] += 0.50 * dirFlux[1];
          directionalHeatFlux[i][2] += 0.50 * dirFlux[2];

          if(store_contact_data_)
          {
              conduction_contact_area_[i] += contactArea;
              n_conduction_contacts_[i] += 1.;
          }
          if (newton_pair || j < nlocal)
          {
            heatFlux[j] -= flux;
            directionalHeatFlux[j][0] += 0.50 * dirFlux[0];
            directionalHeatFlux[j][1] += 0.50 * dirFlux[1];
            directionalHeatFlux[j][2] += 0.50 * dirFlux[2];

            if(store_contact_data_)
            {
                conduction_contact_area_[j] += contactArea;
                n_conduction_contacts_[j] += 1.;
            }
          }
        }

        if(cpl_flag && cpl) cpl->add_heat(i,j,flux);
      }
    }
  }

  if(newton_pair)
  {
    fix_heatFlux->do_reverse_comm();
    fix_directionalHeatFlux->do_reverse_comm();
    fix_conduction_contact_area_->do_reverse_comm();
    fix_n_conduction_contacts_->do_reverse_comm();
  }

  if(!cpl_flag && store_contact_data_)
  for(int i = 0; i < nlocal; i++)
  {
     if(n_conduction_contacts_[i] > 0.5)
        conduction_contact_area_[i] /= n_conduction_contacts_[i];
  }
}
/* ----------------------------------------------------------------------
   register and unregister callback to compute
------------------------------------------------------------------------- */
void FixHeatGranElectric::end_of_step()
{
  if (!electricMode_) return;
  if (update->atime < joule_heat_startTime) return;

	// Gate on the MAX temperature among *electrically active* particles in this group.
  int *mask = atom->mask;
  int *type = atom->type;

  double localSum = 0.0;
  long long localCnt = 0;

  for (int i = 0; i < atom->nlocal; ++i) {
    //if ((mask[i] & groupbit) && e_cond[type[i]-1] > 1e-30) {		//this is just to calculate the average temp of conductive particles 
      localSum += Temp[i];
      localCnt += 1;
   // }
  }

  double gSum = 0.0;
  long long gCnt = 0;
  MPI_Allreduce(&localSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&localCnt, &gCnt, 1, MPI_LONG_LONG, MPI_SUM, world);
  const double aveActive = gCnt ? (gSum / (double)gCnt) : 0.0; // optional diag

  if (temp_control_flag_ && aveActive >= temp_control_) {
    int me; MPI_Comm_rank(world,&me);
    if (me==0) {
      printf("[ELECTRIC] Temp gate: Tave(active)=%.3f K >= %.3f K  → skip electric + Joule this step\n",
             aveActive, temp_control_);
    }
    return;
  }
  
  // make sure wall accumulators / ghosts are visible here
  if (fix_pwER)    fix_pwER->do_forward_comm();
  if (fix_cMatrix) fix_cMatrix->do_forward_comm();
  if (fix_eP)      fix_eP->do_forward_comm();
  if (fix_eP_old)  fix_eP_old->do_forward_comm();

  // dispatch to the right templated path
  if (history_flag == 0 && area_calculation_mode_ == CONDUCTION_CONTACT_AREA_OVERLAP)
    electric_end_step_impl<0, CONDUCTION_CONTACT_AREA_OVERLAP>();
  else if (history_flag == 1 && area_calculation_mode_ == CONDUCTION_CONTACT_AREA_OVERLAP)
    electric_end_step_impl<1, CONDUCTION_CONTACT_AREA_OVERLAP>();
  else if (history_flag == 0 && area_calculation_mode_ == CONDUCTION_CONTACT_AREA_CONSTANT)
    electric_end_step_impl<0, CONDUCTION_CONTACT_AREA_CONSTANT>();
  else if (history_flag == 1 && area_calculation_mode_ == CONDUCTION_CONTACT_AREA_CONSTANT)
    electric_end_step_impl<1, CONDUCTION_CONTACT_AREA_CONSTANT>();
  else if (history_flag == 0 && area_calculation_mode_ == CONDUCTION_CONTACT_AREA_PROJECTION)
    electric_end_step_impl<0, CONDUCTION_CONTACT_AREA_PROJECTION>();
  else if (history_flag == 1 && area_calculation_mode_ == CONDUCTION_CONTACT_AREA_PROJECTION)
    electric_end_step_impl<1, CONDUCTION_CONTACT_AREA_PROJECTION>();
}

template<int HISTFLAG, int CONTACTAREA>
inline void FixHeatGranElectric::electric_end_step_impl()
{
  // everything here happens after all post_force() calls of this step
  detectedComponentsParallel<HISTFLAG, CONTACTAREA>();
  iterativeDomainDecompSolve<HISTFLAG, CONTACTAREA>();
  calculateJouleHeat<HISTFLAG, CONTACTAREA>();
}

void FixHeatGranElectric::register_compute_pair_local(ComputePairGranLocal *ptr)
{
  if(cpl != NULL)
    error->all(FLERR,"Fix heat/gran/electric allows only one compute of type pair/local");
  cpl = ptr;
}

void FixHeatGranElectric::unregister_compute_pair_local(ComputePairGranLocal *ptr)
{
  if(cpl != ptr)
      error->all(FLERR,"Illegal situation in FixHeatGranElectric::unregister_compute_pair_local");
  cpl = NULL;
}


inline bool FixHeatGranElectric::type_is_conductive(int itype) const {
  const double s0 = (sigma0_ ? sigma0_[itype] : e_cond[itype]);
  return std::isfinite(s0) && (s0 > 1e-30);
}

inline double FixHeatGranElectric::sigma_at_typeT(int itype, double T) const {
  // If the material is nominally insulating, it's insulating. Period.
  if (!type_is_conductive(itype)) return 0.0;

  const double s0 = (sigma0_ ? sigma0_[itype] : e_cond[itype]);
  const double a  = (alpha_  ? alpha_[itype]  : 0.0);

  double denom = 1.0 + a*(T - Tref_);
  if (denom < 1e-6) denom = 1e-6;

  double s = s0 / denom;

  // Floor ONLY applies to already-conductive materials, to avoid zero divisions.
  if (!std::isfinite(s) || s < sigma_floor_) s = sigma_floor_;
  return s;
}

template<int HISTFLAG, int CONTACTAREA>
void FixHeatGranElectric::calculateJouleHeat()
{
  // ---------- per-step bookkeeping (for heatSource accumulation)
  static long long last_ts = -1;
  static int calls = 0;
  const long long ts = (long long)update->ntimestep;
  if (ts != last_ts) { last_ts = ts; calls = 0; }
  const bool first_call = (++calls == 1);

  // ---------- bring fields local
  fix_eP->do_forward_comm();
  if (fix_pwER)    fix_pwER->do_forward_comm();
  if (fix_cMatrix) fix_cMatrix->do_forward_comm();
  if (fix_pwEP)    fix_pwEP->do_forward_comm();
  if (fix_total_C) fix_total_C->do_forward_comm();

  // ---------- handy aliases
  const int nlocal = atom->nlocal;
  double *V   = fix_eP->vector_atom;
  double *GwA = fix_pwER->vector_atom;     // Σ_w G_iw
  double *cM  = fix_cMatrix->vector_atom;  // Σ_w G_iw Vw
  double *VwT = fix_pwEP->vector_atom;     // diagnostic (last wall V seen)

  double *qj     = fix_qj->vector_atom;
  double *qj_pp  = fix_qj_pp->vector_atom;
  double *qj_pw  = fix_qj_pw->vector_atom;
  double *tot_I  = fix_tot_I->vector_atom;
  double *totalC = fix_total_C ? fix_total_C->vector_atom : nullptr;
  double **Iviz  = fix_current->array_atom;

  // ---------- numeric guards
  const double epsV  = 1e-12;
  const double epsC  = 1e-40;
  const double epsR2 = 1e-30;

  // ---------- zero viz and per-step scratch
  for (int i=0;i<nlocal+atom->nghost;++i) { Iviz[i][0]=Iviz[i][1]=Iviz[i][2]=0.0; }

  // use two scratch property/atom slots already allocated in the fix
  fix_rhsG->set_all(0.0);      // PP half power bucket for "j"
  fix_diagFlag->set_all(0.0);  // PP diagonal accumulator
  double *ppHalfToJ = fix_rhsG->vector_atom;
  double *diagPP    = fix_diagFlag->vector_atom;

  // ---------- find GLOBAL terminal voltages and GLOBAL wall coupling sum
  double vmin_loc =  std::numeric_limits<double>::infinity();
  double vmax_loc = -std::numeric_limits<double>::infinity();
  double Gw_sum_loc = 0.0;

  for (int i=0;i<nlocal;++i) {
    if (GwA[i] > epsC && std::isfinite(VwT[i])) {
      vmin_loc = std::min(vmin_loc, VwT[i]);
      vmax_loc = std::max(vmax_loc, VwT[i]);
      Gw_sum_loc += GwA[i];
    }
  }

  double V0=vmin_loc, V1=vmax_loc, Gw_sum=Gw_sum_loc;
  MPI_Allreduce(&vmin_loc,&V0,   1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&vmax_loc,&V1,   1,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&Gw_sum_loc,&Gw_sum,1,MPI_DOUBLE,MPI_SUM,world);

  const bool any_wall = (Gw_sum > epsC);
  if (!std::isfinite(V0)) V0 = 0.0;
  if (!std::isfinite(V1)) V1 = V0;

  const double dVw = V1 - V0;
  const bool two_terminal = (std::fabs(dVw) > epsV);

  // ---------- EASY OUT: no wall contacts anywhere ⇒ floating network
  if (!any_wall) {
    double Ppp_loc = 0.0; // but we will *not* compute PP power on a floating network
    for (int i=0;i<nlocal;++i) {
      if (totalC) totalC[i] = 0.0;
      qj_pp[i]=qj_pw[i]=0.0;
      if (first_call) { if (heatSource) heatSource[i] += (0.0 - (std::isfinite(qj[i])?qj[i]:0.0)); }
      qj[i]=0.0; tot_I[i]=0.0;
    }
    double Ppp=0.0;
    MPI_Allreduce(&Ppp_loc,&Ppp,1,MPI_DOUBLE,MPI_SUM,world);
    /*if (comm->me==0 && screen) {
      fprintf(screen,"[ELECTRIC] floating network (no wall contacts): P=0, I=0\n");
      fflush(screen);
    }*/
    return;
  }

  // ---------- PARTICLE-PARTICLE (unique edge, split 50/50)
  double Ppp_loc = 0.0;
  std::vector<double> qpp_new(nlocal,0.0), qpw_new(nlocal,0.0);

  if (two_terminal) {
    // Only compute PP dissipation when there is a driving EMF;
    // otherwise rounding in V may create spurious heat.
    int inum = pair_gran->list->inum;
    int *ilist = pair_gran->list->ilist;
    int *numneigh = pair_gran->list->numneigh;
    int **firstneigh = pair_gran->list->firstneigh;
    int **firstHist  = (HISTFLAG && pair_gran->listgranhistory)
                     ? pair_gran->listgranhistory->firstneigh : nullptr;
    double **x = atom->x;
    int *type = atom->type;
    tagint *tags = atom->tag_enable ? atom->tag : nullptr;

    for (int ii=0; ii<inum; ++ii) {
      const int i = ilist[ii];
      if (i>=nlocal) continue;
      const int ti = type[i]-1;
      if (ti<0) continue;

      int *jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (int jj=0; jj<jnum; ++jj) {
        // unique edge filter
        bool unique = true;
        if (HISTFLAG && firstHist) unique = (firstHist[i][jj] != 0);
        else if (tags)             unique = (tags[i] < atom->tag[jlist[jj] & NEIGHMASK]);
        else                       unique = (i < (jlist[jj] & NEIGHMASK));
        if (!unique) continue;

        const int j = jlist[jj] & NEIGHMASK; //add uniqueness of the contact by tag i < tag j
        const int tj = type[j]-1;
        if (tj<0) continue;

        const double Gij = conductanceLocalIJ<HISTFLAG,CONTACTAREA>(i,j);
        if (!(Gij > epsC) || !std::isfinite(Gij)) continue;

        diagPP[i] += Gij;
        diagPP[j] += Gij;

        const double dV = V[i] - V[j];
        if (!std::isfinite(dV)) continue;

        const double Qij = Gij * dV * dV;
        if (!(Qij >= 0.0) || !std::isfinite(Qij)) continue;

        Ppp_loc += Qij;
        const double halfQ = 0.5*Qij;
        qpp_new[i]   += halfQ;
        ppHalfToJ[j] += halfQ;

        // current visualization (i→j)
        const double dx=x[i][0]-x[j][0], dy=x[i][1]-x[j][1], dz=x[i][2]-x[j][2];
        const double r2 = dx*dx+dy*dy+dz*dz;
        if (r2 > epsR2) {
          const double rinv = 1.0/std::sqrt(r2);
          const double Iij  = Gij * dV;
          if (std::isfinite(Iij)) {
            Iviz[i][0] += -Iij*dx*rinv;
            Iviz[i][1] += -Iij*dy*rinv;
            Iviz[i][2] += -Iij*dz*rinv;
            
            // add to j (may be ghost on this rank)
	    Iviz[j][0] +=  Iij * dx*rinv;
	    Iviz[j][1] +=  Iij * dy*rinv;
	    Iviz[j][2] +=  Iij * dz*rinv;
          }
        }
      }
    }
  } // else: single-terminal ⇒ do not accumulate PP Joule

  // send halves/diagonals back to owners
  fix_rhsG->do_reverse_comm();
  //fix_current->do_reverse_comm();	// bring j’s ghost contributions home 
  //fix_current->do_forward_comm();	//Suggested by GPT
  fix_diagFlag->do_reverse_comm();
  for (int i=0;i<nlocal;++i) qpp_new[i] += ppHalfToJ[i];

  // ---------- PARTICLE-WALL using recovered (G0,G1)
  double I0_loc=0.0, I1_loc=0.0, Ppw_loc=0.0;

  for (int i=0;i<nlocal;++i) {
    const double Gw = GwA[i];
    if (totalC) totalC[i] = diagPP[i] + (std::isfinite(Gw) ? std::max(0.0,Gw) : 0.0);

    if (!(Gw > epsC) || !std::isfinite(Gw)) continue;
    double cMi = std::isfinite(cM[i]) ? cM[i] : 0.0;

    const double Vi = V[i];
    if (!std::isfinite(Vi)) continue;

    double G0=0.0, G1=0.0;
    if (two_terminal) {
      // clamp cM into [Gw*min, Gw*max] to avoid pathological splits
      const double cLo = Gw*std::min(V0,V1);
      const double cHi = Gw*std::max(V0,V1);
      if (cMi < cLo) cMi = cLo;
      if (cMi > cHi) cMi = cHi;

      const double denom = (std::fabs(dVw) > epsV) ? dVw : (dVw>=0?epsV:-epsV);
      G1 = (cMi - Gw*V0)/denom;
      G0 = Gw - G1;

      // numeric clamping
      if (G1 < 0.0) G1 = 0.0;
      if (G1 > Gw ) G1 = Gw;
      if (G0 < 0.0) G0 = 0.0;
      if (G0 > Gw ) G0 = Gw;
    } else {
      // one terminal globally: anchor exists but no EMF → no device power
      G0 = Gw; G1 = 0.0;
    }

    const double I0 = G0*(V0 - Vi);
    const double I1 = G1*(V1 - Vi);
    if (std::isfinite(I0)) I0_loc += I0;
    if (std::isfinite(I1)) I1_loc += I1;

    // PW Joule at node
    double Qpw = 0.0;
    if (two_terminal) {
      Qpw = G0*(Vi - V0)*(Vi - V0) + G1*(Vi - V1)*(Vi - V1);
    } else {
      // one terminal, no EMF ⇒ ideally Vi≈V0 ⇒ tiny; but keep it safe
      Qpw = Gw*(Vi - V0)*(Vi - V0);
    }
    if (std::isfinite(Qpw) && Qpw >= 0.0) {
      qpw_new[i] += Qpw;
      Ppw_loc    += Qpw;
    }
  }

  // ---------- finalize per-atom power & heatSource accumulation
  for (int i=0;i<nlocal;++i) {
    const double q_old = std::isfinite(qj[i]) ? qj[i] : 0.0;
    const double qpp = qpp_new[i];
    const double qpw = qpw_new[i];
    const double qtot = qpp + qpw;

    qj_pp[i] = qpp;
    qj_pw[i] = qpw;
    qj[i]    = qtot;

    if (first_call && heatSource) {
      const double dq = qtot - q_old;
      if (std::isfinite(dq)) heatSource[i] += dq;

    }
  }

  // ---------- global reductions & logging
  double Ppp=0.0, Ppw=0.0, I0=0.0, I1=0.0;
  MPI_Allreduce(&Ppp_loc,&Ppp,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&Ppw_loc,&Ppw,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&I0_loc, &I0, 1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&I1_loc, &I1, 1,MPI_DOUBLE,MPI_SUM,world);

  const double Idev = std::fabs(I1);     // two-terminal case; single-terminal ⇒ ~0
  for (int i=0;i<nlocal;++i) tot_I[i] = Idev;

  const double Pterm = V0*I0 + V1*I1;
  const double Psum  = Ppp + Ppw;
  const double denom = std::max(1e-12, std::fabs(Psum) + std::fabs(Pterm));
  const double relMis = std::fabs(Psum - Pterm)/denom;

  /* 
  //  ---- Terminal print ------ //
  if (comm->me==0 && screen) {
    fprintf(screen,
      "[ELECTRIC] V0=%.6g V, V1=%.6g V | I=%.6g A | "
      "Ppp=%.3f W, Ppw=%.3f W, Psum=%.3f W | Pterm=%.3f W | relMis=%.3e%s\n",
      V0, V1, Idev, Ppp, Ppw, Psum, Pterm, relMis,
      first_call ? "" : "  [extra call; no accumulation]");
    fflush(screen);
  }
  */

  // ---------- rank-0 CSV logging (once per timestep, minimal overhead)   (COMMENTING THIS SECTION TO REDUCE THE COMPUTATIONAL LOAD)
  
  /*if (comm->me == 0 && first_call && (electric_log_stride_ >= 1)) {
    const long long step = (long long)update->ntimestep;
    if ((step % (long long)electric_log_stride_) == 0) {

    static FILE* fp = nullptr;
    static bool wrote_header = false;

    if (!fp) {
      char fname[256];
      snprintf(fname, sizeof(fname), "electric_device_%s.csv", id);
      fp = fopen(fname, "w");
      if (fp) setvbuf(fp, nullptr, _IOFBF, 1<<20); // 1MB buffered I/O
    }
    if (fp) {
      if (!wrote_header) {
	// step,time_s,V[V],I[A],R[Ohm],P[W]
	fprintf(fp, "step,time_s,V[V],I[A],R[Ohm],P[W]\n");
	wrote_header = true;
      }
      const double Vabs = std::fabs(V1 - V0);
      const double Iabs = std::fabs(I1);
      const double Rtot = (Iabs > 1e-20) ? (Vabs / Iabs)
	                                 : std::numeric_limits<double>::infinity();
      const double Ptot = Pterm; // terminal power (equals dissipated power ideally)
      const double t    = update->ntimestep * update->dt;

      fprintf(fp, "%lld,%.17g,%.9g,%.9g,%.9g,%.9g\n",
	      (long long)update->ntimestep, t, Vabs, Iabs, Rtot, Ptot);

      if ((update->ntimestep & 0x3FF) == 0) fflush(fp); // flush every 1024 steps
      }
    }
  }*/
}


// NEW: a helper to fetch the a-spot fraction for a type pair
inline double FixHeatGranElectric::aspot_fraction_for_types(int ti, int tj) const {
  double f = aspot_fraction_global_;
  if (aspotFrac_) {
    // per-atom-type values are scalars; take the mean for a pair
    f = 0.5*(std::max(0.0, aspotFrac_[ti]) + std::max(0.0, aspotFrac_[tj]));
  }
  if (f <= 0.0) return 0.0;
  if (f > 1.0)  return 1.0;
  return f;
}

template<int HISTFLAG, int CONTACTAREA>
void FixHeatGranElectric::detectedComponentsParallel()
{
  const int me     = comm->me;
  const int nprocs = comm->nprocs;

  const int nlocal = atom->nlocal;
  tagint   *tag    = atom->tag;
  int      *type   = atom->type;
  int      *mask   = atom->mask;
  const int gb     = groupbit;

  int  *numneigh    = pair_gran->list->numneigh;
  int **firstneigh  = pair_gran->list->firstneigh;
  int **firstHist   = (HISTFLAG && pair_gran->listgranhistory)
                      ? pair_gran->listgranhistory->firstneigh : NULL;

  double *pwER   = fix_pwER   ? fix_pwER->vector_atom   : NULL;   // Σ G_iw
  double *compID = fix_compID ? fix_compID->vector_atom : NULL;   // output

  if (!compID) error->all(FLERR,"detectedComponentsParallel: compID missing");
  if (!pwER)   error->all(FLERR,"detectedComponentsParallel: pwER missing");

  // Clear compID on local atoms
  for (int i = 0; i < nlocal; ++i) compID[i] = 0.0;

  const double tiny_sigma = 1.0e-30;

  struct Edge { unsigned long long a, b; };

  std::vector<unsigned long long> nodes_local;
  nodes_local.reserve(nlocal);

  std::vector<Edge> edges_local;
  edges_local.reserve(8 * std::max(1,nlocal));

  // ---- nodes: conductive atoms in group
  for (int i = 0; i < nlocal; i++) {
    const int itype = type[i] - 1;
    if (!(mask[i] & gb)) continue;
    if (itype >= 0 && e_cond[itype] > tiny_sigma) {
      nodes_local.push_back((unsigned long long) tag[i]);
    }
  }

  const int inum = pair_gran->list->inum;
  int *ilist     = pair_gran->list->ilist;

  // ---- edges: unique-edge rule (history owner when available; else tag ordering)
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    const int itype = type[i] - 1;
    if (!(mask[i] & gb)) continue;
    if (itype < 0 || e_cond[itype] <= tiny_sigma) continue;

    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & NEIGHMASK;
      if (j == i) continue;

      const int jtype = type[j] - 1;
      if (!(mask[j] & gb)) continue;
      if (jtype < 0 || e_cond[jtype] <= tiny_sigma) continue;

      // Unique owner for the edge:
      bool unique = true;
      if (HISTFLAG && firstHist) {
        unique = (firstHist[i][jj] != 0);
      } else {
        // fallback: tag ordering
        const unsigned long long ti = (unsigned long long) tag[i];
        const unsigned long long tj = (unsigned long long) tag[j];
        unique = (ti < tj);
      }
      if (!unique) continue;

      const double Cij = conductanceLocalIJ<HISTFLAG,CONTACTAREA>(i,j);
      if (!(Cij > 0.0) || !std::isfinite(Cij)) continue;

      const unsigned long long a = (unsigned long long) tag[i];
      const unsigned long long b = (unsigned long long) tag[j];
      if (a == b) continue;

      // With tag ordering fallback, (a<b) is guaranteed; with HISTFLAG it may not,
      // so enforce ordering for safety.
      if (a < b) edges_local.push_back({a,b});
      else       edges_local.push_back({b,a});
    }
  }

  // ---- Allgather helper arrays
  std::vector<int> cnt(nprocs), dsp(nprocs);

  // ---- gather nodes, dedup
  std::vector<unsigned long long> nodes_all;
  {
    int nloc_nodes = (int) nodes_local.size();
    MPI_Allgather(&nloc_nodes,1,MPI_INT,cnt.data(),1,MPI_INT,world);
    int tot = 0;
    for (int p = 0; p < nprocs; ++p) { dsp[p] = tot; tot += cnt[p]; }
    nodes_all.resize(tot);

    MPI_Allgatherv((void*) nodes_local.data(), nloc_nodes, MPI_UNSIGNED_LONG_LONG,
                   (void*) nodes_all.data(), cnt.data(), dsp.data(),
                   MPI_UNSIGNED_LONG_LONG, world);

    std::sort(nodes_all.begin(), nodes_all.end());
    nodes_all.erase(std::unique(nodes_all.begin(), nodes_all.end()), nodes_all.end());
  }
  const int N = (int) nodes_all.size();

  // ---- gather edges, dedup
  std::vector<Edge> E;
  {
    std::vector<unsigned long long> edges_flat;
    edges_flat.reserve(edges_local.size()*2);
    for (size_t k=0;k<edges_local.size();++k){
      edges_flat.push_back(edges_local[k].a);
      edges_flat.push_back(edges_local[k].b);
    }

    int nloc_edges = (int) edges_flat.size();
    MPI_Allgather(&nloc_edges,1,MPI_INT,cnt.data(),1,MPI_INT,world);
    int tot = 0;
    for (int p=0;p<nprocs;++p){ dsp[p]=tot; tot += cnt[p]; }

    std::vector<unsigned long long> edges_all(tot);
    MPI_Allgatherv((void*) edges_flat.data(), nloc_edges, MPI_UNSIGNED_LONG_LONG,
                   (void*) edges_all.data(), cnt.data(), dsp.data(),
                   MPI_UNSIGNED_LONG_LONG, world);

    E.reserve(edges_all.size()/2);
    for (size_t k=0;k+1<edges_all.size();k+=2){
      unsigned long long a = edges_all[k];
      unsigned long long b = edges_all[k+1];
      if (a==b) continue;
      if (a>b) std::swap(a,b);
      E.push_back({a,b});
    }

    std::sort(E.begin(),E.end(),[](const Edge& x,const Edge& y){
      if (x.a<y.a) return true;
      if (x.a>y.a) return false;
      return x.b<y.b;
    });
    E.erase(std::unique(E.begin(),E.end(),[](const Edge& x,const Edge& y){
      return x.a==y.a && x.b==y.b;
    }), E.end());
  }

  // ---- map tag -> compact index
  std::unordered_map<unsigned long long,int> idx;
  idx.reserve((size_t)N*2);
  for (int i=0;i<N;++i) idx[nodes_all[i]] = i;

  // ---- DSU
  struct DSU {
    std::vector<int> p,r;
    explicit DSU(int n): p(n), r(n,0) { for(int i=0;i<n;++i) p[i]=i; }
    int find(int a){ while(p[a]!=a){ p[a]=p[p[a]]; a=p[a]; } return a; }
    void unite(int a,int b){
      a=find(a); b=find(b); if(a==b) return;
      if(r[a]<r[b]) std::swap(a,b);
      p[b]=a; if(r[a]==r[b]) r[a]++;
    }
  } uf(N);

  for (size_t e=0;e<E.size();++e){
    const int ia = idx[E[e].a];
    const int ib = idx[E[e].b];
    uf.unite(ia,ib);
  }

  std::vector<unsigned long long> rep_min(N, std::numeric_limits<unsigned long long>::max());
  std::vector<int> comp_sz(N,0);
  for (int i=0;i<N;++i){
    const int r = uf.find(i);
    rep_min[r] = std::min(rep_min[r], nodes_all[i]);
    comp_sz[r]++;
  }

  std::unordered_map<unsigned long long,unsigned long long> tag2rep;
  tag2rep.reserve((size_t)N*2);
  for (int i=0;i<N;++i){
    tag2rep[nodes_all[i]] = rep_min[ uf.find(i) ];
  }

  // ---- reps that touch a wall (Σ G_iw > 0 anywhere)
  std::unordered_set<unsigned long long> wall_touch;
  {
    std::unordered_set<unsigned long long> local;
    local.reserve((size_t)nlocal/2+8);

    for (int i=0;i<nlocal;++i){
      const int itype = type[i]-1;
      if (!(mask[i] & gb)) continue;
      if (itype<0 || e_cond[itype]<=tiny_sigma) continue;
      if (!(pwER[i] > 0.0)) continue;

      const unsigned long long t = (unsigned long long) tag[i];
      auto it = tag2rep.find(t);
      if (it != tag2rep.end()) local.insert(it->second);
    }

    std::vector<unsigned long long> send(local.begin(), local.end());
    int nloc_w = (int)send.size();
    MPI_Allgather(&nloc_w,1,MPI_INT,cnt.data(),1,MPI_INT,world);
    int tot=0; for(int p=0;p<nprocs;++p){ dsp[p]=tot; tot+=cnt[p]; }
    std::vector<unsigned long long> recv(tot);

    MPI_Allgatherv((void*)send.data(), nloc_w, MPI_UNSIGNED_LONG_LONG,
                   (void*)recv.data(), cnt.data(), dsp.data(),
                   MPI_UNSIGNED_LONG_LONG, world);

    std::sort(recv.begin(), recv.end());
    recv.erase(std::unique(recv.begin(), recv.end()), recv.end());
    wall_touch.insert(recv.begin(), recv.end());
  }

  // ---- keep interesting components (size>=2 or touches wall)
  std::vector<unsigned long long> kept;
  kept.reserve(N);
  for (int i=0;i<N;++i){
    if (rep_min[i] == std::numeric_limits<unsigned long long>::max()) continue;
    if (comp_sz[i] >= 2 || wall_touch.count(rep_min[i]))
      kept.push_back(rep_min[i]);
  }
  std::sort(kept.begin(), kept.end());

  std::unordered_map<unsigned long long,int> rep2dense;
  rep2dense.reserve(kept.size()*2);
  for (int k=0;k<(int)kept.size();++k) rep2dense[kept[k]] = k+1;

  // ---- write compID on local atoms
  for (int i=0;i<nlocal;++i){
    const int itype = type[i]-1;
    if (!(mask[i] & gb) || itype<0 || e_cond[itype]<=tiny_sigma) { compID[i]=0.0; continue; }

    const unsigned long long t = (unsigned long long) tag[i];
    auto it_rep = tag2rep.find(t);
    if (it_rep == tag2rep.end()) { compID[i]=0.0; continue; }

    auto it_dense = rep2dense.find(it_rep->second);
    compID[i] = (it_dense == rep2dense.end()) ? 0.0 : (double) it_dense->second;
  }

  fix_compID->do_forward_comm();
}


template<int HISTFLAG, int CONTACTAREA>
void FixHeatGranElectric::iterativeDomainDecompSolve()
{
  // --- sync inputs ---
  if (fix_pwER)    fix_pwER->do_forward_comm();
  if (fix_cMatrix) fix_cMatrix->do_forward_comm();
  if (fix_compID)  fix_compID->do_forward_comm();
  if (fix_pwEP)    fix_pwEP->do_forward_comm();
  fix_eP->do_forward_comm();

  const int nloc  = atom->nlocal;
  const int nall  = atom->nlocal + atom->nghost;

  const int inum       = pair_gran->list->inum;
  int *     ilist      = pair_gran->list->ilist;
  int *     numneigh   = pair_gran->list->numneigh;
  int **    firstneigh = pair_gran->list->firstneigh;
  int **    firstHist  = (HISTFLAG && pair_gran->listgranhistory)
                         ? pair_gran->listgranhistory->firstneigh : NULL;

  int   *mask = atom->mask;
  int   *type = atom->type;
  tagint*tag  = atom->tag;
  const int gb = groupbit;

  // Primary fields
  double *V    = fix_eP->vector_atom;      // potentials
  double *Gw   = fix_pwER->vector_atom;    // Σ_w G_iw
  double *cM   = fix_cMatrix->vector_atom; // Σ_w (G_iw * Vw)
  double *cid  = fix_compID->vector_atom;  // component IDs
  double *Vtag = fix_pwEP ? fix_pwEP->vector_atom : NULL; // last wall potential

  // Scratch (property/atom arrays are sized for nall internally)
  double *sumPP_acc  = fix_rhsG->vector_atom;     // Σ_j G_ij V_j
  double *diagPP_acc = fix_diagFlag->vector_atom; // Σ_j G_ij

  fix_rhsG->set_all(0.0);
  fix_diagFlag->set_all(0.0);

  // Numerics
  const double OMEGA    = 0.75;
  const double RTOL_V   = 1e-10;
  const double ATOL_V   = 1e-14;
  const double ITOL_REL = 1e-10;
  const double ITOL_ABS = 1e-12;
  const int    KMAX     = 500;
  const double tinyC    = 1e-40;
  const double tinyD    = 1e-300;
  const double tinySig  = 1e-30;
  const double Veps     = 1e-12;

  auto gmax = [&](double l){ double g=l; MPI_Allreduce(&l,&g,1,MPI_DOUBLE,MPI_MAX,world); return g; };
  auto gsum = [&](double l){ double g=l; MPI_Allreduce(&l,&g,1,MPI_DOUBLE,MPI_SUM,world); return g; };

  // -------------------- number of components --------------------
  int Nc_loc = 0;
  for (int i=0;i<nloc;++i){
    if (!(mask[i] & gb)) continue;
    const int c = (int) llround(cid[i]);
    Nc_loc = std::max(Nc_loc, c);
  }
  int Nc = 0;
  MPI_Allreduce(&Nc_loc,&Nc,1,MPI_INT,MPI_MAX,world);

  if (Nc == 0) {
    for (int i=0;i<nloc;++i){ cid[i]=0.0; V[i]=eP0; }
    fix_compID->do_forward_comm();
    fix_eP->do_forward_comm();
    return;
  }

  // -------------------- per-component terminal info --------------------
  std::vector<int>    touch_loc(Nc+1, 0);
  std::vector<int>    hasEdge_loc(Nc+1,0);
  std::vector<double> Vmin_loc(Nc+1,  std::numeric_limits<double>::infinity());
  std::vector<double> Vmax_loc(Nc+1, -std::numeric_limits<double>::infinity());

  // wall touch + wall potential extrema
  for (int i=0;i<nloc;++i){
    if (!(mask[i] & gb)) continue;
    const int c = (int) llround(cid[i]);
    if (c <= 0) continue;

    const int ti = type[i]-1;
    if (ti<0 || e_cond[ti] <= tinySig) continue;

    const double Gwi = std::isfinite(Gw[i]) ? Gw[i] : 0.0;
    if (Gwi > tinyC) {
      touch_loc[c] = 1;
      if (Vtag) {
        const double Vw = Vtag[i];
        if (std::isfinite(Vw)) {
          Vmin_loc[c] = std::min(Vmin_loc[c], Vw);
          Vmax_loc[c] = std::max(Vmax_loc[c], Vw);
        }
      }
    }
  }

  // any PP edge in component (unique owner if history available; else tag ordering)
  for (int ii=0; ii<inum; ++ii) {
    const int i = ilist[ii];
    if (!(mask[i] & gb)) continue;

    const int ti = type[i]-1;
    if (ti<0 || e_cond[ti] <= tinySig) continue;

    const int c = (int) llround(cid[i]);
    if (c <= 0) continue;

    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    for (int jj=0; jj<jnum; ++jj) {
      const int j = jlist[jj] & NEIGHMASK;
      if (!(mask[j] & gb)) continue;

      const int tj = type[j]-1;
      if (tj<0 || e_cond[tj] <= tinySig) continue;

      bool unique = true;
      if (HISTFLAG && firstHist) unique = (firstHist[i][jj] != 0);
      else                       unique = ( (unsigned long long)tag[i] < (unsigned long long)tag[j] );
      if (!unique) continue;

      const double Cij = conductanceLocalIJ<HISTFLAG,CONTACTAREA>(i,j);
      if (Cij > tinyC && std::isfinite(Cij)) { hasEdge_loc[c] = 1; }
    }
  }

  // global reduce
  std::vector<int>    touch(Nc+1,0), hasEdge(Nc+1,0);
  std::vector<double> Vmin(Nc+1,  std::numeric_limits<double>::infinity());
  std::vector<double> Vmax(Nc+1, -std::numeric_limits<double>::infinity());

  MPI_Allreduce(touch_loc.data(),   touch.data(),   Nc+1, MPI_INT,    MPI_MAX, world);
  MPI_Allreduce(hasEdge_loc.data(), hasEdge.data(), Nc+1, MPI_INT,    MPI_MAX, world);
  MPI_Allreduce(Vmin_loc.data(),    Vmin.data(),    Nc+1, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(Vmax_loc.data(),    Vmax.data(),    Nc+1, MPI_DOUBLE, MPI_MAX, world);

  // classify components
  std::vector<char>   compSolve(Nc+1,0);
  std::vector<double> VminC(Nc+1, -std::numeric_limits<double>::infinity());
  std::vector<double> VmaxC(Nc+1,  std::numeric_limits<double>::infinity());

  for (int c=1;c<=Nc;++c){
    if (!hasEdge[c]) continue;
    if (!touch[c])   continue;
    if (!std::isfinite(Vmin[c]) || !std::isfinite(Vmax[c])) continue;

    if ((Vmax[c] - Vmin[c]) > Veps) {
      compSolve[c] = 1;
      VminC[c] = Vmin[c];
      VmaxC[c] = Vmax[c];
    } else {
      compSolve[c] = 0; // single-terminal
    }
  }

  // -------------------- build active set, sanitize compID / eP --------------------
  std::vector<char> active(nall,0);
  int nSolveLoc = 0;

  for (int i=0;i<nloc;++i){
    const int ti = type[i]-1;
    const int c  = (int) llround(cid[i]);

    if (!(mask[i] & gb) || ti<0 || e_cond[ti] <= tinySig || c<=0) {
      cid[i]=0.0; V[i]=eP0; active[i]=0; continue;
    }

    if (!compSolve[c]) {
      if (touch[c]) {
        // NOTE: keep your current behavior; consider pin-to-Vmin for neutrality
        double Vfix = std::numeric_limits<double>::quiet_NaN();
        if (std::isfinite(Vmin[c])) Vfix = Vmin[c];
        if (std::isfinite(Vmax[c])) Vfix = Vmax[c];
        V[i] = std::isfinite(Vfix) ? Vfix : eP0;
      } else {
        V[i] = eP0;
      }
      cid[i]=0.0;
      active[i]=0;
      continue;
    }

    // two-terminal solve
    active[i]=1;
    ++nSolveLoc;

    if (std::isfinite(VminC[c])) {
      if (V[i] < VminC[c]) V[i] = VminC[c];
      if (V[i] > VmaxC[c]) V[i] = VmaxC[c];
    }
  }

  int nSolveGlb=0;
  MPI_Allreduce(&nSolveLoc,&nSolveGlb,1,MPI_INT,MPI_SUM,world);

  // Make ghosts consistent
  fix_compID->do_forward_comm();
  fix_eP->do_forward_comm();

  // Fill active flags for ghosts safely (so active[j] is valid)
  for (int i=nloc;i<nall;++i){
    const int ti = type[i]-1;
    const int c  = (int) llround(cid[i]);
    active[i] = ((mask[i] & gb) && ti>=0 && e_cond[ti] > tinySig && c>0) ? 1 : 0;
  }

  if (nSolveGlb == 0) return;

  // -------------------- assemble PP diagonal once --------------------
  // We compute diagPP_acc (Σ Gij) once; it doesn't change during Jacobi.
  for (int ii=0; ii<inum; ++ii) {
    const int i = ilist[ii];
    if (!active[i]) continue;

    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj=0; jj<jnum; ++jj) {
      const int j = jlist[jj] & NEIGHMASK;
      if (!active[j]) continue;

      bool unique = true;
      if (HISTFLAG && firstHist) unique = (firstHist[i][jj] != 0);
      else                       unique = ( (unsigned long long)tag[i] < (unsigned long long)tag[j] );
      if (!unique) continue;

      const double Cij = conductanceLocalIJ<HISTFLAG,CONTACTAREA>(i,j);
      if (!(Cij > tinyC) || !std::isfinite(Cij)) continue;

      diagPP_acc[i] += Cij;
      diagPP_acc[j] += Cij;
    }
  }
  fix_diagFlag->do_reverse_comm();

  // diagonals = PP + PW (local only)
  std::vector<double> diag(nloc,0.0);
  for (int i=0;i<nloc;++i) if (active[i]) {
    const double Gwi = std::isfinite(Gw[i]) ? std::max(0.0,Gw[i]) : 0.0;
    const double dpp = std::isfinite(diagPP_acc[i]) ? diagPP_acc[i] : 0.0;
    diag[i] = dpp + Gwi;
    if (diag[i] < tinyD) diag[i] = tinyD;
  }

  std::vector<double> Vnew(nloc,0.0);

  // -------------------- damped Jacobi iterations --------------------
  for (int k=0; k<KMAX; ++k) {

    fix_eP->do_forward_comm();

    // Build sumPP_acc = Σ Gij * Vj (depends on V, so rebuild each sweep)
    fix_rhsG->set_all(0.0);

    for (int ii=0; ii<inum; ++ii) {
      const int i = ilist[ii];
      if (!active[i]) continue;

      int *jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (int jj=0; jj<jnum; ++jj) {
        const int j = jlist[jj] & NEIGHMASK;
        if (!active[j]) continue;

        bool unique = true;
        if (HISTFLAG && firstHist) unique = (firstHist[i][jj] != 0);
        else                       unique = ( (unsigned long long)tag[i] < (unsigned long long)tag[j] );
        if (!unique) continue;

        const double Cij = conductanceLocalIJ<HISTFLAG,CONTACTAREA>(i,j);
        if (!(Cij > tinyC) || !std::isfinite(Cij)) continue;

        sumPP_acc[i] += Cij * V[j];
        sumPP_acc[j] += Cij * V[i];
      }
    }

    fix_rhsG->do_reverse_comm();

    double du_max_loc = 0.0;
    double v_max_loc  = 0.0;

    for (int i=0;i<nloc;++i) if (active[i] && cid[i] > 0.0) {
      const double si  = std::isfinite(sumPP_acc[i]) ? sumPP_acc[i] : 0.0;
      const double bi  = std::isfinite(cM[i])        ? cM[i]        : 0.0;
      const double rhs = si + bi;

      const double Vi = V[i];
      const double Vc = rhs / diag[i];
      const double Vr = Vi + OMEGA * (Vc - Vi);

      const int c = (int) llround(cid[i]);
      double Vcl = Vr;
      if (c > 0 && std::isfinite(VminC[c])) {
        if (Vcl < VminC[c]) Vcl = VminC[c];
        if (Vcl > VmaxC[c]) Vcl = VmaxC[c];
      }

      Vnew[i] = Vcl;
      du_max_loc = std::max(du_max_loc, std::fabs(Vcl - Vi));
      v_max_loc  = std::max(v_max_loc,  std::fabs(Vcl));
    }

    for (int i=0;i<nloc;++i) if (active[i] && cid[i] > 0.0) V[i] = Vnew[i];

    const double g_du = gmax(du_max_loc);
    const double g_vm = gmax(v_max_loc);

    // residual r_i = d_i V_i − sumPP_i − cM_i
    double rmax_loc=0.0, rsum_loc=0.0;
    for (int i=0;i<nloc;++i) if (active[i] && cid[i] > 0.0) {
      const double si = std::isfinite(sumPP_acc[i]) ? sumPP_acc[i] : 0.0;
      const double bi = std::isfinite(cM[i])        ? cM[i]        : 0.0;
      const double ri = diag[i]*V[i] - si - bi;
      const double ar = std::fabs(ri);
      rmax_loc = std::max(rmax_loc, ar);
      rsum_loc += ar;
    }

    const double rmax = gmax(rmax_loc);
    const double rsum = gsum(rsum_loc);
    const double rref = std::max(1e-30, rsum);

    const bool convV = (g_du <= std::max(ATOL_V, RTOL_V*std::max(1.0, g_vm)));
    const bool convI = (rmax <= ITOL_ABS) || ((rmax / rref) <= ITOL_REL);
    if (convV && convI) break;
  }

  // final clamp for safety
  for (int i=0;i<nloc;++i) if (active[i] && cid[i] > 0.0) {
    const int c = (int) llround(cid[i]);
    if (c > 0 && std::isfinite(VminC[c])) {
      if (V[i] < VminC[c]) V[i] = VminC[c];
      if (V[i] > VmaxC[c]) V[i] = VmaxC[c];
    }
  }

  fix_eP->do_forward_comm();
}



/*
template<int HISTFLAG, int CONTACTAREA>
void FixHeatGranElectric::detectedComponentsParallel()
{
  const int me     = comm->me;
  const int nprocs = comm->nprocs;

  const int nlocal = atom->nlocal;
  tagint   *tag    = atom->tag;
  int      *type   = atom->type;
  int      *mask   = atom->mask;
  const int gb     = groupbit;

  // neighbor lists
  int  *numneigh    = pair_gran->list->numneigh;
  int **firstneigh  = pair_gran->list->firstneigh;
  int **firstHist   = (HISTFLAG && pair_gran->listgranhistory)
                      ? pair_gran->listgranhistory->firstneigh : NULL;

  // fields
  double *pwER   = fix_pwER   ? fix_pwER->vector_atom   : NULL;   // Σ G_iw
  double *compID = fix_compID ? fix_compID->vector_atom : NULL;   // output

  if (!compID)
    error->all(FLERR,"detectedComponentsParallel: compID missing");
  if (!pwER)
    error->all(FLERR,"detectedComponentsParallel: pwER missing");

  // clear compID on all local atoms; we'll rewrite for conductive ones in group
  for (int i = 0; i < nlocal; ++i)
    compID[i] = 0.0;

  const double tiny_sigma = 1.0e-30;	//non-conductive

  struct Edge {
    unsigned long long a, b;
  };

  std::vector<unsigned long long> nodes_local;
  nodes_local.reserve(nlocal);

  std::vector<Edge> edges_local;
  edges_local.reserve(8 * std::max(1,nlocal));

  // --- nodes: only conductive atoms in fix group ---
  for (int i = 0; i < nlocal; i++) {
    const int itype = type[i] - 1;
    if (!(mask[i] & gb)) continue;
    if (itype >= 0 && e_cond[itype] > tiny_sigma) {
      nodes_local.push_back((unsigned long long) tag[i]);
    }
  }

  const int inum = pair_gran->list->inum;
  int *ilist     = pair_gran->list->ilist;

  // --- edges: only between conductive atoms in fix group ---
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    const int itype = type[i] - 1;
    if (!(mask[i] & gb)) continue;
    if (itype < 0 || e_cond[itype] <= tiny_sigma) continue;

    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & NEIGHMASK;
      if (j == i) continue;

      const int jtype = type[j] - 1;
      if (!(mask[j] & gb)) continue;
      if (jtype < 0 || e_cond[jtype] <= tiny_sigma) continue;

      //if (HISTFLAG && firstHist && firstHist[i][jj] == 0) continue;

      const double Cij = conductanceLocalIJ<HISTFLAG,CONTACTAREA>(i,j);
      if (Cij <= 0.0 || !std::isfinite(Cij)) continue;

      unsigned long long a = (unsigned long long) tag[i];
      unsigned long long b = (unsigned long long) tag[j];
      if (a == b) continue;
      if (a > b) std::swap(a,b);
      edges_local.push_back({a,b});
    }
  }

  // --- Allgather helpers ---
  std::vector<int> cnt(nprocs), dsp(nprocs);

  // ----- gather nodes and dedup -----
  std::vector<unsigned long long> nodes_all;
  {
    int nloc_nodes = (int) nodes_local.size();
    MPI_Allgather(&nloc_nodes,1,MPI_INT,cnt.data(),1,MPI_INT,world);
    int tot = 0;
    for (int p = 0; p < nprocs; ++p) {
      dsp[p] = tot;
      tot += cnt[p];
    }
    nodes_all.resize(tot);
    MPI_Allgatherv((void*) nodes_local.data(), nloc_nodes, MPI_UNSIGNED_LONG_LONG,
                   (void*) nodes_all.data(), cnt.data(), dsp.data(),
                   MPI_UNSIGNED_LONG_LONG, world);

    std::sort(nodes_all.begin(), nodes_all.end());
    nodes_all.erase(std::unique(nodes_all.begin(), nodes_all.end()), nodes_all.end());
  }
  const int N = (int) nodes_all.size();

  // ----- gather edges and dedup -----
  std::vector<Edge> E;
  {
    std::vector<unsigned long long> edges_flat;
    edges_flat.reserve(edges_local.size()*2);
    for (size_t k = 0; k < edges_local.size(); ++k) {
      edges_flat.push_back(edges_local[k].a);
      edges_flat.push_back(edges_local[k].b);
    }

    int nloc_edges = (int) edges_flat.size();
    MPI_Allgather(&nloc_edges,1,MPI_INT,cnt.data(),1,MPI_INT,world);
    int tot = 0;
    for (int p = 0; p < nprocs; ++p) {
      dsp[p] = tot;
      tot += cnt[p];
    }

    std::vector<unsigned long long> edges_all;
    edges_all.resize(tot);
    MPI_Allgatherv((void*) edges_flat.data(), nloc_edges, MPI_UNSIGNED_LONG_LONG,
                   (void*) edges_all.data(), cnt.data(), dsp.data(),
                   MPI_UNSIGNED_LONG_LONG, world);

    E.reserve(edges_all.size()/2);
    for (size_t k = 0; k+1 < edges_all.size(); k += 2) {
      unsigned long long a = edges_all[k];
      unsigned long long b = edges_all[k+1];
      if (a == b) continue;
      if (a > b) std::swap(a,b);
      E.push_back({a,b});
    }

    std::sort(E.begin(), E.end(),
      [](const Edge& x, const Edge& y) {
        if (x.a < y.a) return true;
        if (x.a > y.a) return false;
        return x.b < y.b;
      });
    E.erase(std::unique(E.begin(),E.end(),
      [](const Edge& x, const Edge& y){ return x.a==y.a && x.b==y.b; }), E.end());
  }

  // ----- map tag → compact index -----
  std::unordered_map<unsigned long long,int> idx;
  idx.reserve(N*2);
  for (int i = 0; i < N; ++i)
    idx[nodes_all[i]] = i;

  // ----- DSU -----
  struct DSU {
    std::vector<int> p,r;
    explicit DSU(int n): p(n), r(n,0) {
      for (int i=0; i<n; ++i) p[i]=i;
    }
    int find(int a){
      while(p[a]!=a){ p[a]=p[p[a]]; a=p[a]; }
      return a;
    }
    void unite(int a,int b){
      a=find(a); b=find(b); if(a==b) return;
      if(r[a]<r[b]) std::swap(a,b);
      p[b]=a;
      if(r[a]==r[b]) r[a]++;
    }
  } uf(N);

  for (size_t e = 0; e < E.size(); ++e) {
    const Edge &ed = E[e];
    const int ia = idx[ed.a];
    const int ib = idx[ed.b];
    uf.unite(ia,ib);
  }

  std::vector<unsigned long long> rep_min(N, std::numeric_limits<unsigned long long>::max());
  std::vector<int> comp_sz(N,0);
  for (int i = 0; i < N; ++i) {
    const int r = uf.find(i);
    if (nodes_all[i] < rep_min[r]) rep_min[r] = nodes_all[i];
    comp_sz[r]++;
  }

  std::unordered_map<unsigned long long,unsigned long long> tag2rep;
  tag2rep.reserve(N*2);
  for (int i = 0; i < N; ++i) {
    const int r = uf.find(i);
    tag2rep[nodes_all[i]] = rep_min[r];
  }

  // ----- reps that touch a wall (Σ G_iw > 0 anywhere in group) -----
  std::unordered_set<unsigned long long> wall_touch;
  {
    std::unordered_set<unsigned long long> local;
    local.reserve(nlocal/2+8);
    for (int i = 0; i < nlocal; ++i) {
      const int itype = type[i]-1;
      if (!(mask[i] & gb)) continue;
      if (itype < 0 || e_cond[itype] <= tiny_sigma) continue;
      if (pwER[i] <= 0.0) continue;
      const unsigned long long t = (unsigned long long) tag[i];
      auto it = tag2rep.find(t);
      if (it != tag2rep.end())
        local.insert(it->second);
    }

    std::vector<unsigned long long> send(local.begin(), local.end());
    std::vector<unsigned long long> recv;
    int nloc_w = (int) send.size();
    MPI_Allgather(&nloc_w,1,MPI_INT,cnt.data(),1,MPI_INT,world);
    int tot = 0;
    for (int p=0; p<nprocs; ++p) { dsp[p]=tot; tot += cnt[p]; }
    recv.resize(tot);
    MPI_Allgatherv((void*) send.data(), nloc_w, MPI_UNSIGNED_LONG_LONG,
                   (void*) recv.data(), cnt.data(), dsp.data(),
                   MPI_UNSIGNED_LONG_LONG, world);
    std::sort(recv.begin(), recv.end());
    recv.erase(std::unique(recv.begin(), recv.end()), recv.end());
    wall_touch.insert(recv.begin(), recv.end());
  }

  // ----- keep components of interest -----
  std::vector<unsigned long long> kept;
  kept.reserve(N);
  for (int i = 0; i < N; ++i) {
    const unsigned long long rep = rep_min[i];
    if (rep == std::numeric_limits<unsigned long long>::max()) continue;
    if (comp_sz[i] >= 2 || wall_touch.count(rep))
      kept.push_back(rep);
  }
  std::sort(kept.begin(), kept.end());

  std::unordered_map<unsigned long long,int> rep2dense;
  rep2dense.reserve(kept.size()*2);
  for (int k = 0; k < (int)kept.size(); ++k)
    rep2dense[kept[k]] = k+1;   // IDs 1..Nc

  // ----- write compID on local atoms (only for nodes in fix group & conductive) -----
  for (int i = 0; i < nlocal; ++i) {
    const int itype = type[i]-1;
    if (!(mask[i] & gb) || itype < 0 || e_cond[itype] <= tiny_sigma) {
      compID[i] = 0.0;
      continue;
    }
    const unsigned long long t = (unsigned long long) tag[i];
    auto it_rep = tag2rep.find(t);
    if (it_rep == tag2rep.end()) {
      compID[i] = 0.0;
      continue;
    }
    auto it_dense = rep2dense.find(it_rep->second);
    compID[i] = (it_dense == rep2dense.end()) ? 0.0 : (double) it_dense->second;
  }

  // make ghosts consistent for later stages
  fix_compID->do_forward_comm();
}
*/


/*
template<int HISTFLAG, int CONTACTAREA>
void FixHeatGranElectric::iterativeDomainDecompSolve()
{
  // --- sync inputs ---
  if (fix_pwER)    fix_pwER->do_forward_comm();
  if (fix_cMatrix) fix_cMatrix->do_forward_comm();
  if (fix_compID)  fix_compID->do_forward_comm();
  if (fix_pwEP)    fix_pwEP->do_forward_comm();
  fix_eP->do_forward_comm();

  const int nloc       = atom->nlocal;
  const int inum       = pair_gran->list->inum;
  int *     ilist      = pair_gran->list->ilist;
  int *     numneigh   = pair_gran->list->numneigh;
  int **    firstneigh = pair_gran->list->firstneigh;
  int **    firstHist  = (HISTFLAG && pair_gran->listgranhistory)
                         ? pair_gran->listgranhistory->firstneigh : NULL;

  int   *mask = atom->mask;
  int   *type = atom->type;
  const int gb = groupbit;

  // Primary fields
  double *V    = fix_eP->vector_atom;      // potentials (physical V)
  double *Gw   = fix_pwER->vector_atom;    // Σ_w G_iw
  double *cM   = fix_cMatrix->vector_atom; // Σ_w (G_iw * Vw)
  double *cid  = fix_compID->vector_atom;  // component IDs (from detectedComponentsParallel)
  double *Vtag = fix_pwEP ? fix_pwEP->vector_atom : NULL; // last wall potential (diagnostic)

  // Scratch (reused arrays)
  double *sumPP_acc  = fix_rhsG->vector_atom;   // Σ_j G_ij V_j
  double *diagPP_acc = fix_diagFlag->vector_atom; // Σ_j G_ij
  fix_rhsG->set_all(0.0);
  fix_diagFlag->set_all(0.0);

  // Numerics
  const double OMEGA    = 0.75;
  const double RTOL_V   = 1e-10;
  const double ATOL_V   = 1e-14;
  const double ITOL_REL = 1e-10;
  const double ITOL_ABS = 1e-12;
  const int    KMAX     = 500;
  const double tinyC    = 1e-40;
  const double tinyD    = 1e-300;
  const double tinySig  = 1e-30;
  const double Veps     = 1e-12;

  auto gmax = [&](double l){ double g=l; MPI_Allreduce(&l,&g,1,MPI_DOUBLE,MPI_MAX,world); return g; };
  auto gsum = [&](double l){ double g=l; MPI_Allreduce(&l,&g,1,MPI_DOUBLE,MPI_SUM,world); return g; };

  // -------------------- number of components --------------------
  int Nc_loc = 0;
  for (int i = 0; i < nloc; ++i) {
    if (!(mask[i] & gb)) continue;
    const int c = (int) llround(cid[i]);
    if (c > Nc_loc) Nc_loc = c;
  }
  int Nc = 0;
  MPI_Allreduce(&Nc_loc,&Nc,1,MPI_INT,MPI_MAX,world);
  if (Nc == 0) {
    // no conductive components in this group; make sure eP is benign
    for (int i=0; i<nloc; ++i) {
      cid[i] = 0.0;
      V[i]   = eP0;
    }
    fix_compID->do_forward_comm();
    fix_eP->do_forward_comm();
    return;
  }

  // -------------------- per-component terminal info --------------------
  std::vector<int>    touch_loc(Nc+1, 0);  // touches *some* wall
  std::vector<int>    hasEdge_loc(Nc+1,0); // has at least one PP edge
  std::vector<double> Vmin_loc(Nc+1,  std::numeric_limits<double>::infinity());
  std::vector<double> Vmax_loc(Nc+1, -std::numeric_limits<double>::infinity());

  // wall touch + wall potential extrema
  for (int i = 0; i < nloc; ++i) {
    if (!(mask[i] & gb)) continue;
    const int c = (int) llround(cid[i]);
    if (c <= 0) continue;

    const int ti = type[i]-1;
    if (ti < 0 || e_cond[ti] <= tinySig) continue;

    const double Gwi = std::isfinite(Gw[i]) ? Gw[i] : 0.0;
    if (Gwi > tinyC) {
      touch_loc[c] = 1;
      if (Vtag) {
        const double Vw = Vtag[i];
        if (std::isfinite(Vw)) {
          if (Vw < Vmin_loc[c]) Vmin_loc[c] = Vw;
          if (Vw > Vmax_loc[c]) Vmax_loc[c] = Vw;
        }
      }
    }
  }

  // any PP edge in component (unique owner when history available)
  for (int ii = 0; ii < inum; ++ii) {
    const int i = ilist[ii];
    if (!(mask[i] & gb)) continue;

    const int ti = type[i]-1;
    if (ti < 0 || e_cond[ti] <= tinySig) continue;

    const int c = (int) llround(cid[i]);
    if (c <= 0) continue;

    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; ++jj) {
      if (HISTFLAG && firstHist && firstHist[i][jj] == 0) continue;
      const int j = jlist[jj] & NEIGHMASK;
      if (!(mask[j] & gb)) continue;

      const int tj = type[j]-1;
      if (tj < 0 || e_cond[tj] <= tinySig) continue;

      const double Cij = conductanceLocalIJ<HISTFLAG,CONTACTAREA>(i,j);
      if (Cij > tinyC && std::isfinite(Cij)) {
        hasEdge_loc[c] = 1;
      }
    }
  }

  // global reduce
  std::vector<int>    touch(Nc+1,0), hasEdge(Nc+1,0);
  std::vector<double> Vmin(Nc+1,  std::numeric_limits<double>::infinity());
  std::vector<double> Vmax(Nc+1, -std::numeric_limits<double>::infinity());

  MPI_Allreduce(touch_loc.data(),   touch.data(),   Nc+1, MPI_INT,    MPI_MAX, world);
  MPI_Allreduce(hasEdge_loc.data(), hasEdge.data(), Nc+1, MPI_INT,    MPI_MAX, world);
  MPI_Allreduce(Vmin_loc.data(),    Vmin.data(),    Nc+1, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(Vmax_loc.data(),    Vmax.data(),    Nc+1, MPI_DOUBLE, MPI_MAX, world);

  // classify components
  std::vector<char>   compSolve(Nc+1,0);   // two-terminal, solved
  std::vector<double> VminC(Nc+1, -std::numeric_limits<double>::infinity());
  std::vector<double> VmaxC(Nc+1,  std::numeric_limits<double>::infinity());

  for (int c = 1; c <= Nc; ++c) {
    if (!hasEdge[c]) continue;      // must have some PP connectivity
    if (!touch[c]) continue;        // no wall contact ⇒ floating network
    if (!std::isfinite(Vmin[c]) || !std::isfinite(Vmax[c])) continue;

    if ((Vmax[c] - Vmin[c]) > Veps) {
      // two distinct electrode potentials present ⇒ solve this component
      compSolve[c] = 1;
      VminC[c] = Vmin[c];
      VmaxC[c] = Vmax[c];
    } else {
      // single-terminal component, we will pin it and then drop compID
      compSolve[c] = 0;
    }
  }

  // -------------------- build active set, sanitize compID / eP --------------------
  const int nall = atom->nlocal + atom->nghost;
  std::vector<char> active(nall,0); 	//active(nloc, 0);
  int nSolveLoc = 0;

  for (int i = 0; i < nloc; ++i) {
    const int ti = type[i]-1;
    const int c  = (int) llround(cid[i]);

    // outside group / non-conductive / no component: flush
    if (!(mask[i] & gb) || ti < 0 || e_cond[ti] <= tinySig || c <= 0) {
      cid[i] = 0.0;
      V[i]   = eP0;
      active[i] = 0;
      continue;
    }

    // here: in group, conductive, c>0 from detection
    if (!compSolve[c]) {
      // not two-terminal: either single-terminal or floating
      if (touch[c]) {
        // single-terminal: pin to that wall potential (prefer any finite value)
        double Vfix = std::numeric_limits<double>::quiet_NaN();
        if (std::isfinite(Vmin[c])) Vfix = Vmin[c];
        if (std::isfinite(Vmax[c])) Vfix = Vmax[c];
        V[i] = std::isfinite(Vfix) ? Vfix : eP0;
      } else {
        // floating: reset to reference
        V[i] = eP0;
      }
      cid[i] = 0.0;   // from here on compID encodes *only* bridging components
      active[i] = 0;
      continue;
    }

    // two-terminal, to be solved
    active[i] = 1;
    ++nSolveLoc;

    // clamp initial guess into terminal range for stability
    if (std::isfinite(VminC[c])) {
      if (V[i] < VminC[c]) V[i] = VminC[c];
      if (V[i] > VmaxC[c]) V[i] = VmaxC[c];
    }
  }

  int nSolveGlb = 0;
  MPI_Allreduce(&nSolveLoc,&nSolveGlb,1,MPI_INT,MPI_SUM,world);
  fix_compID->do_forward_comm();   // we changed cid

  if (nSolveGlb == 0) {
    // nothing to solve electrically this step; potentials are now sanitized
    fix_eP->do_forward_comm();
    return;
  }

  // -------------------- assemble symmetric PP contributions --------------------
  for (int ii = 0; ii < inum; ++ii) {
    const int i = ilist[ii];
    if (!active[i]) continue;

    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; ++jj) {
      if (HISTFLAG && firstHist && firstHist[i][jj] == 0) continue;
      const int j = jlist[jj] & NEIGHMASK;
      if (!active[j]) continue;

      const double Cij = conductanceLocalIJ<HISTFLAG,CONTACTAREA>(i,j);
      if (!(Cij > tinyC) || !std::isfinite(Cij)) continue;

      sumPP_acc[i]  += Cij * V[j];
      diagPP_acc[i] += Cij;

      sumPP_acc[j]  += Cij * V[i];
      diagPP_acc[j] += Cij;
    }
  }

  fix_rhsG->do_reverse_comm();
  fix_diagFlag->do_reverse_comm();

  // diagonals = PP + PW
  std::vector<double> diag(nloc,0.0);
  for (int i = 0; i < nloc; ++i) if (active[i]) {
    const double Gwi = std::isfinite(Gw[i]) ? Gw[i] : 0.0;
    const double dpp = std::isfinite(diagPP_acc[i]) ? diagPP_acc[i] : 0.0;
    diag[i] = dpp + ((Gwi > 0.0) ? Gwi : 0.0);
  }

  std::vector<double> Vnew(nloc,0.0);

  // -------------------- damped Jacobi iterations --------------------
  for (int k = 0; k < KMAX; ++k) {

    fix_eP->do_forward_comm(); // refresh ghost V

    double du_max_loc = 0.0;
    double v_max_loc  = 0.0;

    for (int i = 0; i < nloc; ++i) if (active[i] && cid[i]) { //added cid[i] condition
      const double d   = (diag[i] > tinyD) ? diag[i] : tinyD;
      const double si  = std::isfinite(sumPP_acc[i]) ? sumPP_acc[i] : 0.0;
      const double bi  = std::isfinite(cM[i])        ? cM[i]        : 0.0;
      const double rhs = si + bi;

      const double Vi  = V[i];
      const double Vc  = rhs / d;
      const double Vr  = Vi + OMEGA * (Vc - Vi);

      const int c = (int) llround(cid[i]);
      double Vcl = Vr;
      if (c > 0 && std::isfinite(VminC[c])) {
        if (Vcl < VminC[c]) Vcl = VminC[c];
        if (Vcl > VmaxC[c]) Vcl = VmaxC[c];
      }

      Vnew[i] = Vcl;
      du_max_loc = std::max(du_max_loc, std::fabs(Vcl - Vi));
      v_max_loc  = std::max(v_max_loc, std::fabs(Vcl));
    }

    for (int i = 0; i < nloc; ++i) if (active[i] && cid[i]) V[i] = Vnew[i]; //added cid[i] condition

    const double g_du = gmax(du_max_loc);
    const double g_vm = gmax(v_max_loc);

    // residual r_i = d_i V_i − sumPP_i − cM_i
    double rmax_loc = 0.0, rsum_loc = 0.0;
    for (int i = 0; i < nloc; ++i) if (active[i] && cid[i]) {  //added cid[i] condition
      const double d  = (diag[i] > tinyD) ? diag[i] : tinyD;
      const double si = std::isfinite(sumPP_acc[i]) ? sumPP_acc[i] : 0.0;
      const double bi = std::isfinite(cM[i])        ? cM[i]        : 0.0;
      const double ri = d*V[i] - si - bi;
      const double ar = std::fabs(ri);
      if (ar > rmax_loc) rmax_loc = ar;
      rsum_loc += ar;
    }
    const double rmax = gmax(rmax_loc);
    const double rsum = gsum(rsum_loc);
    const double rref = std::max(1e-30, rsum);

    const bool convV = (g_du <= std::max(ATOL_V, RTOL_V*std::max(1.0, g_vm)));
    const bool convI = (rmax <= ITOL_ABS) || ((rmax / rref) <= ITOL_REL);
    if (convV && convI) break;

    // rebuild PP sums for next sweep with updated V
    fix_rhsG->set_all(0.0);
    fix_diagFlag->set_all(0.0);

    for (int ii = 0; ii < inum; ++ii) {
      const int i = ilist[ii];
      if (!active[i]) continue; 

      int *jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (int jj = 0; jj < jnum; ++jj) {
        if (HISTFLAG && firstHist && firstHist[i][jj] == 0) continue;
        const int j = jlist[jj] & NEIGHMASK;
        if (!active[j]) continue; 

        const double Cij = conductanceLocalIJ<HISTFLAG,CONTACTAREA>(i,j);
        if (!(Cij > tinyC) || !std::isfinite(Cij)) continue;

        sumPP_acc[i]  += Cij * V[j];
        diagPP_acc[i] += Cij;
        sumPP_acc[j]  += Cij * V[i];
        diagPP_acc[j] += Cij;
      }
    }

    fix_rhsG->do_reverse_comm();
    fix_diagFlag->do_reverse_comm();
  }

  // final clamp for safety
  for (int i = 0; i < nloc; ++i) if (active[i]&& cid[i]) {
    const int c = (int) llround(cid[i]);
    if (c > 0 && std::isfinite(VminC[c])) {
      if (V[i] < VminC[c]) V[i] = VminC[c];
      if (V[i] > VmaxC[c]) V[i] = VmaxC[c];
    }
  }

  fix_eP->do_forward_comm();
}
*/

template<int HISTFLAG, int CONTACTAREA>
double FixHeatGranElectric::conductanceLocalIJ(const int i, const int j)
{
  // ---- 1) Conductivities (S/m), temperature-dependent ----
  const int ti = atom->type[i] - 1;     // 0-based
  const int tj = atom->type[j] - 1;     // 0-based

  const double Ti   = Temp[i];
  const double Tj   = Temp[j];

  if (!type_is_conductive(ti) || !type_is_conductive(tj)) return 0.0;
  const double sigI = sigma_at_typeT(ti, Temp[i]);
  const double sigJ = sigma_at_typeT(tj, Temp[j]);
  if (!(sigI > 0.0) || !(sigJ > 0.0)) return 0.0;

  // ---- 2) Geometry / contact check (cheap) ----
  const double *xi = atom->x[i], *xj = atom->x[j];
  const double ri = atom->radius[i], rj = atom->radius[j];

  const double dx = xj[0]-xi[0], dy = xj[1]-xi[1], dz = xj[2]-xi[2];
  const double r  = std::sqrt(dx*dx + dy*dy + dz*dz);
  double delta    = (ri + rj) - r;                        // normal overlap
  if (!(delta > 0.0)) return 0.0;

  // optional overlap correction (same table as thermal)
  if (CONTACTAREA == CONDUCTION_CONTACT_AREA_OVERLAP && deltan_ratio_) {
    // ti,tj are already 0-based; do NOT subtract 1 again
    delta *= deltan_ratio_[ti][tj];
    if (!(delta > 0.0)) return 0.0;
  }

  // ---- 3) Nominal contact radius a_nom ----
  double a_nom = 0.0;
  if (CONTACTAREA == CONDUCTION_CONTACT_AREA_OVERLAP) {
    // Hertz-like: a = sqrt(Reff * delta), Reff = (ri rj)/(ri + rj)
    const double denom = ri + rj;
    if (!(denom > 0.0)) return 0.0;
    const double Reff = (ri * rj) / denom;
    if (!(Reff > 0.0) || !std::isfinite(Reff)) return 0.0;

    a_nom = std::sqrt(delta * Reff);

    // conservative geometric cap (match p–p logic)
    const double a_cap = 0.35 * std::min(ri, rj);
    if (a_nom > a_cap) a_nom = a_cap;

  } else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_CONSTANT) {
    const double A = fixed_contact_area_;
    if (!(A > 0.0)) return 0.0;
    a_nom = std::sqrt(A / M_PI);

  } else { // CONDUCTION_CONTACT_AREA_PROJECTION
    a_nom = std::min(ri, rj);
  }

  if (!(a_nom > 0.0) || !std::isfinite(a_nom)) return 0.0;

  // ---- 4) Metallic a-spot shrinkage ----
  // a_real = sqrt(fA) * a_nom
  const double fA = aspot_fraction_for_types(ti, tj);
  if (!(fA > 0.0) || !std::isfinite(fA)) return 0.0;

  const double a = a_nom * std::sqrt(fA);
  if (!(a > 0.0) || !std::isfinite(a)) return 0.0;

  // ---- 5) Holm constriction (no film) ----
  // Ri = 1/(2 a σi), Rj = 1/(2 a σj)  =>  G = 1 / (Ri + Rj)
  const double Ri = 1.0 / (2.0 * a * sigI);
  const double Rj = 1.0 / (2.0 * a * sigJ);
  const double R  = Ri + Rj;

  if (!(R > 0.0) || !std::isfinite(R)) return 0.0;

  const double G = 1.0 / R;

  // Optional: quick debug
  // std::cout << "pp G="<<G<<" a="<<a<<" a_nom="<<a_nom
  //           <<" fA="<<fA<<" sigI="<<sigI<<" sigJ="<<sigJ<<"\n";

  if (!std::isfinite(G) || G < Cmin_branch_) return 0.0;
  return G;
}
// --------------------------------------- RADIATIVE HEAT TRANSFER FUNCTIONS ---------------------------------------- //
void FixHeatGranElectric::calculateRadiation(){

  stencilLength = NULL;
  binStencildx = NULL;
  binStencilmdx = NULL;
  binStencildy = NULL;
  binStencilmdy = NULL;
  binStencildz = NULL;
  binStencilmdz = NULL;

  createStencils();

  int i;
  int ibin;
  int maxBounces = 100;

  double **x;
  double *radius;
  int *type;
  int nlocal, nghost;
  double *heatFlux;

  double hitEmis;
  int hitId;
  int hitBin;
  double hitp[3];
  double nextNormal[3];

  double *ci;
  double flux;
  double buffer3[3];
  double d[3];
  double o[3];

  double areai;
  double emisi;
  double radi;
  double tempi;

  nlocal = atom->nlocal;
  nghost = atom->nghost;
  radius = atom->radius;
  type   = atom->type;
  x      = atom->x;
  heatFlux = fix_heatFlux->vector_atom;
  qr = fix_qr->vector_atom;

  updatePtrs();

  for (i = 0; i < nlocal + nghost; i++)
  {
    radi  = radius[i];
    ci    = x[i];
    emisi = emissivity_[type[i]-1];
    tempi = Temp[i];
    const double cg_i  = force->cg(atom->type[i]);
    
    ibin = neighbor->public_coord2bin(ci);

    areai = MY_4PI * radi * radi;
    flux  = areai * emisi * Sigma * tempi * tempi * tempi * tempi * cg_i;

    qr[i] += -flux;
    heatFlux[i] -= flux;

    randOnSphere(ci, radi, o, buffer3);
    randDir(buffer3, d);

    hitId = trace(i, ibin, o, d, buffer3, hitp);

    if (hitId != -1){
      const double& sendflux = flux;
      hitEmis = emissivity_[type[hitId]-1];
      heatFlux[hitId] += hitEmis * sendflux;
      qr[hitId] += hitEmis * sendflux;

      hitBin   = neighbor->public_coord2bin(x[hitId]);

      sub3(hitp, x[hitId], buffer3);
      normalize3(buffer3, nextNormal);

      reflect(i, hitId, hitBin, hitp, nextNormal, sendflux, 1.0-hitEmis, maxBounces, buffer3);

    } else {
      heatFlux[i] += cg_i * (areai * emisi * Sigma * TB * TB * TB * TB);
    }
  }
}

void FixHeatGranElectric::randOnSphere(const double *c, double r, double *ansP, double *ansD)
{
  double s;
  do
  {
    ansD[0] = 2.0*RGen->uniform() - 1.0;
    ansD[1] = 2.0*RGen->uniform() - 1.0;
    s = ansD[0]*ansD[0] + ansD[1]*ansD[1];

  } while (s > 1.0);
  ansD[2] = 2.0*s - 1.0;
  const double a = 2.0 * sqrt(1.0 - s);
  ansD[0] *= a;
  ansD[1] *= a;

  addscaled3(c,ansD,r,ansP);
}

void FixHeatGranElectric::randDir(const double *n, double *d){
  double s;

  do
  {
    d[0] = 2.0*RGen->uniform() - 1.0;
    d[1] = 2.0*RGen->uniform() - 1.0;
    s = d[0]*d[0] + d[1]*d[1];

  } while (s > 1.0);

  d[2] = 2.0*s - 1.0;
  const double a = 2.0 * sqrt(1.0 - s);

  d[0] *= a;
  d[1] *= a;

  const double side = dot3(d, n);
  if (side < 0.0){
    d[0] = -d[0];
    d[1] = -d[1];
    d[2] = -d[2];
  }
}

int FixHeatGranElectric::trace(int orig_id, int ibin, double *o, double *d, double *buffer3, double *hitp)
{
  int *binhead = neighbor->binhead;
  int *bins = neighbor->bins;
  int stencilbin, stbX, stbY, stbZ;
  int currX, currY, currZ;
  int sx = neighbor->sx;
  int sy = neighbor->sy;
  int sz = neighbor->sz;
  int mbins = neighbor->mbins;

  double distsq, distx, disty, distz;
  bool hit;
  double t;

  bool hitFlag = false;
  double hitT  = 1.0e300;
  int hitId    = -1;

  double **x     = atom->x;
  double *radius = atom->radius;

  int i;
  int currentBin = neighbor->public_coord2bin(o);
  int dx = 0;
  int dy = 0;
  int dz = 0;

  int var_nstencil;
  int *stencil = pair_gran->list->stencil;
  int nstencil = pair_gran->list->nstencil;
  int nstencil2D = 0;

  int *currentStencil = stencil;

  bool check_boundary_only = false;

  while ((currentBin != -1) && (hitFlag == false)){

    var_nstencil = check_boundary_only ? nstencil2D : nstencil;

    for (int k = 0; k < var_nstencil; k++){

      stencilbin = currentBin + currentStencil[k];
      neighbor->public_bin2XYZ(stencilbin, stbX, stbY, stbZ);
      neighbor->public_bin2XYZ(currentBin, currX, currY, currZ);
      if (stencilbin < 0 || stencilbin >= mbins || abs(stbX - currX) > sx || abs(stbY - currY) > sy || abs(stbZ - currZ) > sz)
        continue;

      for (i = binhead[stencilbin]; i >= 0; i = bins[i]){

        if (i == orig_id){
          continue;
        }

        hit = intersectRaySphere(o, d, x[i], radius[i], t, buffer3);
        if (hit){
          hitFlag = true;
          if (t < hitT){
            hitT = t;
            hitId = i;
          }
        }
      }
    }

    if (hitFlag){
      addscaled3(o,d,hitT,hitp);
      return hitId;
    }

    ibin       = currentBin;
    currentBin = nextBin(ibin, o, d, raypoint, dx, dy, dz);

    check_boundary_only = true;
    if (dx == 1){
      currentStencil = binStencildx;
      nstencil2D = stencilLength[0];
    }
    else if(dx == -1){
      currentStencil = binStencilmdx;
      nstencil2D = stencilLength[3];
    }
    else if(dy == 1){
      currentStencil = binStencildy;
      nstencil2D = stencilLength[1];
    }
    else if(dy == -1){
      currentStencil = binStencilmdy;
      nstencil2D = stencilLength[4];
    }
    else if(dz == 1){
      currentStencil = binStencildz;
      nstencil2D = stencilLength[2];
    }
    else{
      currentStencil = binStencilmdz;
      nstencil2D = stencilLength[5];
    }

    distx = raypoint[0] - o[0];
    disty = raypoint[1] - o[1];
    distz = raypoint[2] - o[2];
    distsq = distx*distx + disty*disty + distz*distz;
    if (distsq >= cutGhostsq){
      return -1;
    }
  }

  return -1;
}

void FixHeatGranElectric::reflect(int radID, int orig_id, int ibin, double *o, double *d,double flux, double accum_eps, int n, double *buffer3)
{
  const double influx = flux * accum_eps;
  if (n == 0){
    heatFlux[radID] += influx;
    return;
  }
  if (accum_eps < 0.001){
    heatFlux[radID] += influx;
    return;
  }
  double hitp[3];
  int hitId;
  {
    double dd[3];
    randDir(d, dd);

    hitId = trace(orig_id, ibin, o, dd, buffer3, hitp);
  }
  int *type = atom->type;
  if (hitId != -1){
    double **x = atom->x;
    const double& sendflux = influx;
    const double hitEmis = emissivity_[type[hitId]-1];

    heatFlux[hitId] += hitEmis * sendflux;

    int hitBin = neighbor->public_coord2bin(x[hitId]);
    double nextNormal[3];

    sub3(hitp, x[hitId], nextNormal);
    norm3(nextNormal);

    reflect(radID, hitId, hitBin, hitp, nextNormal, flux, (1.0-hitEmis) * accum_eps, n-1, buffer3);
  }
  else {
    const double *radius = atom->radius;
    const double radRad  = radius[radID];
    const double cg_i  = force->cg(atom->type[radID]);
    const double radArea = MY_4PI * radRad * radRad;
    const double radEmis = emissivity_[type[radID]-1];
    heatFlux[radID] += (cg_i * radArea * radEmis * accum_eps * Sigma * TB * TB * TB * TB);
  }
}

int FixHeatGranElectric::nextBin(int ibin, const double *o, const double *d, double *p, int &dx, int &dy, int &dz)
{
  double s;
  double smax = 0.0;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  int nextBinId;

  dx = dy = dz = 0;
  neighbor->public_binBorders(ibin, xlo, xhi, ylo, yhi, zlo, zhi);
  if (d[0] != 0.0){
    s = (xlo - o[0]) / d[0];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];
    if ((s > smax) && (ylo <= p[1] && p[1] <= yhi) && (zlo <= p[2] && p[2] <= zhi)){
      smax = s;
      dy = dz = 0;
      dx = -1;
    }
    s = (xhi - o[0]) / d[0];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];
    if ((s > smax) && (ylo <= p[1] && p[1] <= yhi) && (zlo <= p[2] && p[2] <= zhi)){
      smax = s;
      dy = dz = 0;
      dx = 1;
    }
  }
  if (d[1] != 0.0){
    s = (ylo - o[1]) / d[1];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];
    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (zlo <= p[2] && p[2] <= zhi)){

      smax = s;
      dx = dz = 0;
      dy = -1;
    }
    s = (yhi - o[1]) / d[1];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];
    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (zlo <= p[2] && p[2] <= zhi)){
      smax = s;
      dx = dz = 0;
      dy = 1;
    }
  }
  if (d[2] != 0.0){
    s = (zlo - o[2]) / d[2];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];
    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (ylo <= p[1] && p[1] <= yhi)){
      smax = s;
      dx = dy = 0;
      dz = -1;
    }
    s = (zhi - o[2]) / d[2];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];
    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (ylo <= p[1] && p[1] <= yhi)){
      smax = s;
      dx = dy = 0;
      dz = 1;
    }
  }
  p[0] = o[0] + smax*d[0];
  p[1] = o[1] + smax*d[1];
  p[2] = o[2] + smax*d[2];
  nextBinId = neighbor->public_binHop(ibin,dx,dy,dz);
  if (nextBinId != ibin){
    return nextBinId;
  }
  else {
    error->one(FLERR, "FixHeatGranElectric::nextBin() did not find a suitable next bin.");
    return -1;
  }
}

bool FixHeatGranElectric::intersectRaySphere(const double *o, const double *d, const double *center, double radius, double &t, double *buffer3)
{
  double A, B, C;
  double discr, q;
  double t0, t1, ttemp;
  t = 0.0;
  A = lensq3(d);
  sub3(o, center, buffer3);
  B = 2.0 * dot3(buffer3, d);
  C = lensq3(buffer3) - radius*radius;
  discr = B*B - 4.0*A*C;

  if (discr < 0.0){
    t = 0.0;
    return false;
  }

  q = 0.0;
  if (B < 0.0){
    q = (-B - sqrt(discr)) * 0.5;
  } else {
    q = (-B + sqrt(discr)) * 0.5;
  }

  t0 = q / A;
  t1 = C / q;

  if (t0 > t1){
    ttemp = t0;
    t0 = t1;
    t1 = ttemp;
  }

  if (t1 <= 0.0){
    t = 0.0;
    return false;
  }
  t = t0;
  return true;
}

void FixHeatGranElectric::createStencils()
{
  if(stencilLength == NULL)
    stencilLength = new int[6];

  int mbinx = neighbor->mbinx;
  int mbiny = neighbor->mbiny;
  int sx = neighbor->sx;
  int sy = neighbor->sy;
  int sz = neighbor->sz;
  double cutneighmaxsq = neighbor->cutneighmaxsq;

  int n;
  int i,j,k;
  int nstencil;

  bool hit;

  // binStencildx
  n = (2*sy+1)*(2*sz+1);
  delete [] binStencildx;
  binStencildx = new int[n];
  nstencil = 0;
  for (k = -sz; k <= sz; k++){
    for (j = -sy; j <= sy; j++){
      hit = false;
      i = sx;
      while (!hit && i >= -sz){
        if (neighbor->public_bin_distance(i,j,k) < cutneighmaxsq){
          binStencildx[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        }
        else
          i--;
      }
    }
  }
  stencilLength[0] = nstencil;

  // binStencildy
  n = (2*sx+1)*(2*sz+1);
  delete [] binStencildy;
  binStencildy = new int[n];
  nstencil = 0;
  for (k = -sz; k <= sz; k++){
    for (i = -sx; i <= sx; i++){
      hit = false;
      j = sy;
      while (!hit && j >= -sy){
        if (neighbor->public_bin_distance(i,j,k) < cutneighmaxsq){
          binStencildy[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        }
        else
          j--;
      }
    }
  }
  stencilLength[1] = nstencil;

  // binStencildz
  n = (2*sx+1)*(2*sy+1);
  delete [] binStencildz;
  binStencildz = new int[n];
  nstencil = 0;
  for (j = -sy; j <= sy; j++){
    for (i = -sx; i <= sx; i++){
      hit = false;
      k = sz;
      while (!hit && k >= -sz){
        if (neighbor->public_bin_distance(i,j,k) < cutneighmaxsq){
          binStencildz[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        }
        else
          k--;
      }
    }
  }
  stencilLength[2] = nstencil;

  // binStencilmdx
  n = (2*sy+1)*(2*sz+1);
  delete [] binStencilmdx;
  binStencilmdx = new int[n];
  nstencil = 0;
  for (k = -sz; k <= sz; k++){
    for (j = -sy; j <= sy; j++){
      hit = false;
      i = -sz;
      while (!hit && i <= sz){
        if (neighbor->public_bin_distance(i,j,k) < cutneighmaxsq){
          binStencilmdx[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        }
        else
          i++;
      }
    }
  }
  stencilLength[3] = nstencil;

  // binStencilmdy
  n = (2*sx+1)*(2*sz+1);
  delete [] binStencilmdy;
  binStencilmdy = new int[n];
  nstencil = 0;
  for (k = -sz; k <= sz; k++){
    for (i = -sx; i <= sx; i++){
      hit = false;
      j = -sy;
      while (!hit && j <= sy){
        if (neighbor->public_bin_distance(i,j,k) < cutneighmaxsq){
          binStencilmdy[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        }
        else
          j++;
      }
    }
  }
  stencilLength[4] = nstencil;

  // binStencilmdz
  n = (2*sx+1)*(2*sy+1);
  delete [] binStencilmdz;
  binStencilmdz = new int[n];
  nstencil = 0;
  for (j = -sy; j <= sy; j++){
    for (i = -sx; i <= sx; i++){
      hit = false;
      k = -sz;
      while (!hit && k <= sz){
        if (neighbor->public_bin_distance(i,j,k) < cutneighmaxsq){
          binStencilmdz[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        }
        else
          k++;
      }
    }
  }
  stencilLength[5] = nstencil;
}

void FixHeatGranElectric::addscaled3(const double *v1, const double *v2, double s, double *ans)
{
    ans[0] = v1[0] + s * v2[0];
    ans[1] = v1[1] + s * v2[1];
    ans[2] = v1[2] + s * v2[2];
}
