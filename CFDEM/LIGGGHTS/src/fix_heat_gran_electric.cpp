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

    Copyright 2024-     LCT - UGENT, Ghent
    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include "fix_heat_gran_electric.h"
#include "atom.h"
#include "comm.h"
#include "compute_pair_gran_local.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair_gran.h"
#include "properties.h"
#include "random_mars.h"
#include "update.h"
#include "mpi.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;

using MathConst::MY_4PI;
using namespace MathExtra;

enum{ CONDUCTION_CONTACT_AREA_OVERLAP,
      CONDUCTION_CONTACT_AREA_CONSTANT,
      CONDUCTION_CONTACT_AREA_PROJECTION};

static void zero_with_ghosts(FixPropertyAtom *f, int nall)
{
  if (!f || !f->vector_atom) return;
  double *v = f->vector_atom;
  for (int i = 0; i < nall; ++i) v[i] = 0.0;
}

static FixPropertyAtom *find_or_add_property(Modify *modify, const char *caller, const char *name,
                                             bool isVector, int len1, bool restartAndReverse, bool create)
{
  const char *kind = isVector ? "vector" : "scalar";
  FixPropertyAtom *fix = static_cast<FixPropertyAtom*>(modify->find_fix_property(name,"property/atom",kind,len1,0,caller,false));
  if (fix || !create) return fix;

  const char *flag = restartAndReverse ? "yes" : "no";
  const char *fixarg[11];
  fixarg[0] = name;
  fixarg[1] = "all";
  fixarg[2] = "property/atom";
  fixarg[3] = name;
  fixarg[4] = kind;
  fixarg[5] = flag;
  fixarg[6] = "yes";
  fixarg[7] = flag;
  fixarg[8] = "0.";
  int nfixarg = 9;
  if (isVector) {
    fixarg[9]  = "0.";
    fixarg[10] = "0.";
    nfixarg = 11;
  }
  return modify->add_fix_property_atom(nfixarg,const_cast<char**>(fixarg),caller);
}

FixHeatGranElectric::FixHeatGranElectric(class LAMMPS *lmp, int narg, char **arg) :
  FixHeatGran(lmp, narg, arg)
{
  elec_contact_factor_   = 1.0;
  elec_min_overlap_frac_ = 0.0;
  elec_min_overlap_set_  = false;

  iarg_ = 5;
  bool hasargs = true;
  if(strcmp(style,"heat/gran/electric") == 0) electricMode_ = true;

  while(iarg_ < narg && hasargs)
  {
    hasargs = false;
    if(strcmp(arg[iarg_],"contact_area") == 0) {
      if(strcmp(arg[iarg_+1],"overlap") == 0)
        area_calculation_mode_ = CONDUCTION_CONTACT_AREA_OVERLAP;
      else if(strcmp(arg[iarg_+1],"projection") == 0)
        area_calculation_mode_ = CONDUCTION_CONTACT_AREA_PROJECTION;
      else if(strcmp(arg[iarg_+1],"constant") == 0)
      {
        if (iarg_+3 > narg)
          error->fix_error(FLERR,this,"not enough arguments for keyword 'contact_area constant'");
        area_calculation_mode_ = CONDUCTION_CONTACT_AREA_CONSTANT;
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
      if (iarg_+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'initial_electricPotential'");
      eP0 = atof(arg[iarg_+1]);
      iarg_ += 2;
      hasargs = true;
    } else if(strcmp(arg[iarg_],"start_time") == 0 && electricMode_) {
      if (iarg_+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'start_time'");
      joule_heat_startTime = atof(arg[iarg_+1]);
      if (joule_heat_startTime < 0.) error->fix_error(FLERR,this,"'start_time' value must be >= 0");
      iarg_ += 2;
      hasargs = true;
    } else if(strcmp(arg[iarg_],"temp_control") == 0 && electricMode_) {
      if (iarg_+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'temp_control'");
      if(strcmp(arg[iarg_+1],"yes") == 0) {
        if (iarg_+3 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'temp_control yes'");
        temp_control_flag_ = true;
        temp_control_ = atof(arg[iarg_+2]);
        iarg_ += 3;
      } else if(strcmp(arg[iarg_+1],"no") == 0) {
        temp_control_flag_ = false;
        iarg_ += 2;
      } else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'temp_control'");
      hasargs = true;
    } else if(strcmp(arg[iarg_],"aspot_fraction") == 0) {
      if (iarg_+2 > narg) error->fix_error(FLERR,this,"aspot_fraction needs a scalar in (0,1]");
      aspot_fraction_global_ = atof(arg[iarg_+1]);
      if (aspot_fraction_global_ <= 0.0 || aspot_fraction_global_ > 1.0)
        error->fix_error(FLERR,this,"aspot_fraction must be in (0,1]");
      iarg_ += 2;
      hasargs = true;
    } else if(strcmp(arg[iarg_],"Gmin") == 0) {
      if (iarg_+2 > narg) error->fix_error(FLERR,this,"Gmin needs a scalar >= 0");
      Cmin_branch_ = std::max(0.0, atof(arg[iarg_+1]));
      iarg_ += 2;
      hasargs = true;
    } else if(strcmp(arg[iarg_],"neighbor_dia") == 0) {
      if (iarg_+2 > narg) error->fix_error(FLERR,this,"neighbor_dia needs a scalar >= 1");
      elec_contact_factor_ = atof(arg[iarg_+1]);
      if (elec_contact_factor_ < 1.0) error->fix_error(FLERR,this,"neighbor_dia must be >= 1.0");
      iarg_ += 2;
      hasargs = true;
    } else if(strcmp(arg[iarg_],"min_overlap") == 0) {
      if (iarg_+2 > narg) error->fix_error(FLERR,this,"min_overlap needs a scalar in [0,1)");
      elec_min_overlap_frac_ = atof(arg[iarg_+1]);
      if (elec_min_overlap_frac_ < 0.0 || elec_min_overlap_frac_ >= 1.0)
        error->fix_error(FLERR,this,"min_overlap must be in [0,1)");
      elec_min_overlap_set_ = true;
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
      if(strcmp(arg[iarg_+1],"yes") == 0) {
        if (iarg_+3 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'radiation yes'");
        radiation_flag_ = true;
        TB = atof(arg[iarg_+2]);
        iarg_ += 3;
      } else if(strcmp(arg[iarg_+1],"no") == 0) {
        radiation_flag_ = false;
        iarg_ += 2;
      } else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'radiation'");
      hasargs = true;
    } else if(strcmp(style,"heat/gran/electric") == 0) error->fix_error(FLERR,this,"unknown keyword");
  }

  if(CONDUCTION_CONTACT_AREA_OVERLAP != area_calculation_mode_ && 1 == area_correction_flag_)
    error->fix_error(FLERR,this,"can use 'area_correction' only for 'contact_area = overlap'");

  std::string seed = std::to_string(86028157 + comm->me);
  RGen = new RanMars(lmp, seed.c_str());

  raypoint = new double[3];
}

FixHeatGranElectric::~FixHeatGranElectric()
{
  if (eConductivity_) delete [] eConductivity_;
  if (raypoint)       delete [] raypoint;
  if (stencilLength)  delete [] stencilLength;
  if (binStencildx)   delete [] binStencildx;
  if (binStencilmdx)  delete [] binStencilmdx;
  if (binStencildy)   delete [] binStencildy;
  if (binStencilmdy)  delete [] binStencilmdy;
  if (binStencildz)   delete [] binStencildz;
  if (binStencilmdz)  delete [] binStencilmdz;
  if (RGen)           delete RGen;
}

void FixHeatGranElectric::post_create()
{
  FixHeatGran::post_create();

  fix_conduction_contact_area_ = find_or_add_property(modify,style,"contactAreaConduction",false,0,false,store_contact_data_);
  fix_n_conduction_contacts_   = find_or_add_property(modify,style,"nContactsConduction",  false,0,false,store_contact_data_);
  fix_wall_heattransfer_coeff_ = find_or_add_property(modify,style,"wallHeattransferCoeff",false,0,false,store_contact_data_);
  fix_wall_temperature_        = find_or_add_property(modify,style,"wallTemp",             false,0,false,store_contact_data_);

  fix_qj      = find_or_add_property(modify,style,"qj",     false,0,true,true);
  fix_qj_pp   = find_or_add_property(modify,style,"qj_pp",  false,0,true,true);
  fix_qj_pw   = find_or_add_property(modify,style,"qj_pw",  false,0,true,true);
  fix_qr      = find_or_add_property(modify,style,"qr",     false,0,true,true);
  fix_tot_I   = find_or_add_property(modify,style,"tot_I",  false,0,true,true);
  fix_current = find_or_add_property(modify,style,"current",true, 0,true,true);

  fix_eP = find_or_add_property(modify,style,"eP",false,1,true,false);
  eP_needs_init_ = (fix_eP == NULL);
  if (!fix_eP) fix_eP = find_or_add_property(modify,style,"eP",false,1,true,true);

  fix_eP_old = find_or_add_property(modify,style,"eP_old",false,1,true,true);
  fix_deltaE = find_or_add_property(modify,style,"deltaE",false,1,true,true);

  fix_compID     = find_or_add_property(modify,style,"compID",    false,0,true,true);
  fix_rhsG       = find_or_add_property(modify,style,"rhsG",      false,0,true,true);
  fix_diagFlag   = find_or_add_property(modify,style,"diagFlag",  false,0,true,true);
  fix_total_C    = find_or_add_property(modify,style,"total_C",   false,0,true,true);
  fix_countCheck = find_or_add_property(modify,style,"countCheck",false,0,true,true);

  if(store_contact_data_ && (!fix_conduction_contact_area_ || !fix_n_conduction_contacts_ ||
                             !fix_wall_heattransfer_coeff_ || !fix_eP || !fix_eP_old ||
                             !fix_wall_temperature_))
    error->one(FLERR,"internal error");
}

void FixHeatGranElectric::pre_delete(bool unfixflag)
{
  if(cpl && unfixflag) cpl->reference_deleted();
}

int FixHeatGranElectric::setmask()
{
  int mask = FixHeatGran::setmask();
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

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
    if(electricMode_ && fix_wall_electricPotential_)
      wall_ep_ = fix_wall_electricPotential_->vector_atom;
  }
}

void FixHeatGranElectric::init()
{
  FixHeatGran::init();
  const double *Y, *nu, *Y_orig;
  double expo, Yeff_ij, Yeff_orig_ij, ratio;
  int max_type = atom->get_properties()->max_type();

  if (eConductivity_) delete [] eConductivity_;
  eConductivity_ = new double[max_type];
  fix_eConductivity_ =
    static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalConductivity","property/global","peratomtype",max_type,0,style));

  if(electricMode_)
  {
    FixPropertyGlobal *fix_econd = static_cast<FixPropertyGlobal*>(modify->find_fix_property("electricalConductivity","property/global","peratomtype",max_type,0,style,false));
    e_cond = fix_econd ? fix_econd->get_values() : NULL;

    fix_eP       = static_cast<FixPropertyAtom*>(modify->find_fix_property("eP","property/atom","scalar",1,0,style));
    fix_eP_old   = static_cast<FixPropertyAtom*>(modify->find_fix_property("eP_old","property/atom","scalar",1,0,style));
    fix_qj       = static_cast<FixPropertyAtom*>(modify->find_fix_property("qj","property/atom","scalar",1,0,style));
    fix_qr       = static_cast<FixPropertyAtom*>(modify->find_fix_property("qr","property/atom","scalar",1,0,style));
    fix_cMatrix  = static_cast<FixPropertyAtom*>(modify->find_fix_property("cMatrix","property/atom","scalar",1,0,style));
    fix_pwER     = static_cast<FixPropertyAtom*>(modify->find_fix_property("pwER","property/atom","scalar",1,0,style));
    fix_pwEP     = static_cast<FixPropertyAtom*>(modify->find_fix_property("pwEP","property/atom","scalar",1,0,style));
    fix_current  = static_cast<FixPropertyAtom*>(modify->find_fix_property("current","property/atom","vector",1,0,style));
    fix_compID   = static_cast<FixPropertyAtom*>(modify->find_fix_property("compID","property/atom","scalar",1,0,style));
    fix_rhsG     = static_cast<FixPropertyAtom*>(modify->find_fix_property("rhsG","property/atom","scalar",1,0,style));
    fix_diagFlag = static_cast<FixPropertyAtom*>(modify->find_fix_property("diagFlag","property/atom","scalar",1,0,style));
    fix_total_C  = static_cast<FixPropertyAtom*>(modify->find_fix_property("total_C","property/atom","scalar",1,0,style));
    fix_countCheck = static_cast<FixPropertyAtom*>(modify->find_fix_property("countCheck","property/atom","scalar",1,0,style));

    fix_aspotFrac_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("electricalAspotFraction","property/global","peratomtype",max_type,0,style,false));
    aspotFrac_ = fix_aspotFrac_ ? fix_aspotFrac_->get_values() : NULL;

    FixPropertyGlobal *fix_sigma0 = static_cast<FixPropertyGlobal*>(modify->find_fix_property("electricalSigma0","property/global","peratomtype",max_type,0,style,false));
    sigma0_ = fix_sigma0 ? fix_sigma0->get_values() : e_cond;
    if (!sigma0_)
      error->all(FLERR,"Fix heat/gran/electric needs fix property/global 'electricalSigma0' or 'electricalConductivity' (peratomtype)");

    FixPropertyGlobal *fix_alpha = static_cast<FixPropertyGlobal*>(modify->find_fix_property("electricalAlpha","property/global","peratomtype",max_type,0,style,false));
    alpha_ = fix_alpha ? fix_alpha->get_values() : NULL;

    FixPropertyGlobal *fix_Tref = static_cast<FixPropertyGlobal*>(modify->find_fix_property("electricalTref","property/global","scalar",0,0,style,false));
    if (fix_Tref) Tref_ = fix_Tref->compute_scalar();

    FixPropertyGlobal *fix_floor = static_cast<FixPropertyGlobal*>(modify->find_fix_property("electricalSigmaFloor","property/global","scalar",0,0,style,false));
    if (fix_floor) sigma_floor_ = std::max(0.0, fix_floor->compute_scalar());

    if(radiation_flag_)
      emissivity_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalEmissivity","property/global","peratomtype",max_type,0,style))->get_values();

    wallfix = NULL;
    for (int i = 0; i < modify->nfix; i++) {
      if (strcmp(modify->fix[i]->style, "wall/gran") == 0) {
        wallfix = (FixWallGran *)modify->fix[i];
        break;
      }
    }
    if (wallfix == NULL)
      error->all(FLERR,"FixHeatGranElectric requires fix wall/gran");
    if (force->newton_pair)
      error->all(FLERR,"FixHeatGranElectric assumes 'newton pair off'");

    if (eP_needs_init_) {
      fix_eP->set_all(eP0);
      eP_needs_init_ = false;
    }

    if (elec_contact_factor_ > 1.0 && !elec_min_overlap_set_) {
      elec_min_overlap_frac_ = 1.0e-4;
      if (comm->me == 0 && screen)
        fprintf(screen,"[ELECTRIC] neighbor_dia > 1 with no min_overlap: using min_overlap %.3g\n", elec_min_overlap_frac_);
    }

    double rmax_loc = 0.0;
    for (int i = 0; i < atom->nlocal; ++i) rmax_loc = std::max(rmax_loc, atom->radius[i]);
    double rmax_glb = 0.0;
    MPI_Allreduce(&rmax_loc, &rmax_glb, 1, MPI_DOUBLE, MPI_MAX, world);

    if (rmax_glb > 0.0 && comm->me == 0 && screen) {
      const double f_max = 1.0 + neighbor->skin / (2.0 * rmax_glb);
      if (elec_contact_factor_ > f_max)
        fprintf(screen,"[ELECTRIC] WARNING: neighbor_dia=%.4g is larger than the neighbor list allows "
                       "(max %.4g for skin=%.4g, r_max=%.4g). Increase the neighbor skin.\n",
                elec_contact_factor_, f_max, neighbor->skin, rmax_glb);
      fprintf(screen,"[ELECTRIC] contact model: neighbor_dia=%.4g  min_overlap=%.4g  skin=%.4g  r_max=%.4g\n",
              elec_contact_factor_, elec_min_overlap_frac_, neighbor->skin, rmax_glb);
      fflush(screen);
    }
  }

  for(int i = 1; i < max_type+1; i++)
  {
    eConductivity_[i-1] = fix_eConductivity_->compute_vector(i-1);
    if(eConductivity_[i-1] < 0.)
      error->all(FLERR,"Fix heat/gran/electric: Thermal conductivity must not be < 0");
  }

  if(area_correction_flag_)
  {
    if(!force->pair_match("gran",0))
      error->fix_error(FLERR,this,"area correction only works with using granular pair styles");

    expo = 1./pair_gran->stressStrainExponent();

    Y      = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
    nu     = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
    FixPropertyGlobal *fix_Yorig = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style));
    Y_orig = fix_Yorig->get_values();

    fix_Yorig->new_array(max_type,max_type);
    for(int i = 1; i < max_type+1; i++)
      for(int j = 1; j < max_type+1; j++)
      {
        Yeff_ij      = 1./((1.-pow(nu[i-1],2.))/Y[i-1]     +(1.-pow(nu[j-1],2.))/Y[j-1]);
        Yeff_orig_ij = 1./((1.-pow(nu[i-1],2.))/Y_orig[i-1]+(1.-pow(nu[j-1],2.))/Y_orig[j-1]);
        ratio = pow(Yeff_ij/Yeff_orig_ij,expo);
        fix_Yorig->array_modify(i-1,j-1,ratio);
      }
    deltan_ratio_ = fix_Yorig->get_array_modified();
  }
  updatePtrs();
}

void FixHeatGranElectric::pre_force(int vflag)
{
  FixHeatGran::pre_force(vflag);

  if(store_contact_data_)
  {
    fix_wall_heattransfer_coeff_->set_all(0.);
    fix_wall_temperature_->set_all(0.);
    if (fix_wall_electricPotential_) fix_wall_electricPotential_->set_all(0.);
  }

  fix_qr->set_all(0.);
  fix_rhsG->set_all(0.);
  fix_diagFlag->set_all(0.);

  const int nloc = atom->nlocal;
  double *V    = fix_eP->vector_atom;
  double *Vold = fix_eP_old->vector_atom;
  for (int i = 0; i < nloc; ++i) Vold[i] = V[i];
  fix_eP_old->do_forward_comm();

  if (fix_pwER)    fix_pwER->set_all(0.);
  if (fix_cMatrix) fix_cMatrix->set_all(0.);
  if (fix_pwEP)    fix_pwEP->set_all(0.);
}

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

template <int HISTFLAG,int CONTACTAREA>
void FixHeatGranElectric::post_force_eval(int vflag,int cpl_flag)
{
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

  if(radiation_flag_ && !cpl_flag) calculateRadiation();

  updatePtrs();

  if(store_contact_data_)
  {
    fix_conduction_contact_area_->set_all(0.);
    fix_n_conduction_contacts_->set_all(0.);
  }

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

      if ((HISTFLAG && contact_flag[jj]) || (!HISTFLAG && (rsq < radsum*radsum))) {

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

          if (r < fmax(radi, radj))
          {
            contactArea = fmin(radi,radj);
            contactArea *= contactArea * M_PI;
          }
          else
            contactArea = - M_PI/4.0 * ( (r-radi-radj)*(r+radi-radj)*(r-radi+radj)*(r+radi+radj) )/(r*r);
        }
        else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_CONSTANT)
          contactArea = fixed_contact_area_;
        else
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
    for(int k = 0; k < nlocal; k++)
      if(n_conduction_contacts_[k] > 0.5)
        conduction_contact_area_[k] /= n_conduction_contacts_[k];
}

void FixHeatGranElectric::end_of_step()
{
  if (!electricMode_) return;

  auto removeJouleHeat = [&]() {
    updatePtrs();
    const int nlocal = atom->nlocal;
    const int nall   = nlocal + atom->nghost;
    double *q = fix_qj->vector_atom;
    for (int i = 0; i < nlocal; ++i) {
      const double qold = std::isfinite(q[i]) ? q[i] : 0.0;
      if (heatSource && qold != 0.0) heatSource[i] -= qold;
      q[i] = 0.0;
    }
    if (fix_qj_pp) fix_qj_pp->set_all(0.0);
    if (fix_qj_pw) fix_qj_pw->set_all(0.0);
    if (fix_tot_I) fix_tot_I->set_all(0.0);
    if (fix_current) {
      double **I = fix_current->array_atom;
      for (int i = 0; i < nall; ++i) I[i][0] = I[i][1] = I[i][2] = 0.0;
    }
    edges_.clear();
    esolve_.clear();
  };

  if (update->atime < joule_heat_startTime) {
    removeJouleHeat();
    return;
  }

  double localSum = 0.0;
  long long localCnt = 0;
  for (int i = 0; i < atom->nlocal; ++i) {
    if (!(atom->mask[i] & groupbit)) continue;
    localSum += Temp[i];
    localCnt += 1;
  }

  double gSum = 0.0;
  long long gCnt = 0;
  MPI_Allreduce(&localSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&localCnt, &gCnt, 1, MPI_LONG_LONG, MPI_SUM, world);
  const double Tave = gCnt ? (gSum / (double)gCnt) : 0.0;

  if (temp_control_flag_ && Tave >= temp_control_) {
    if (comm->me == 0 && screen) {
      fprintf(screen,"[ELECTRIC] Temp gate: Tave=%.3f K >= %.3f K, skipping electric and Joule heating this step\n",
              Tave, temp_control_);
      fflush(screen);
    }
    removeJouleHeat();
    return;
  }

  if (fix_pwER)    fix_pwER->do_forward_comm();
  if (fix_cMatrix) fix_cMatrix->do_forward_comm();
  if (fix_eP)      fix_eP->do_forward_comm();
  if (fix_eP_old)  fix_eP_old->do_forward_comm();

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
  buildEdgeCache<HISTFLAG, CONTACTAREA>();

  const int nlocal = atom->nlocal;
  const int nall   = nlocal + atom->nghost;
  int *type = atom->type;
  int *mask = atom->mask;
  double *Rtot = fix_total_C->vector_atom;
  const double *Gw = fix_pwER ? fix_pwER->vector_atom : NULL;

  for (int i = 0; i < nall; ++i) Rtot[i] = 0.0;
  for (size_t e = 0; e < edges_.size(); ++e) {
    Rtot[edges_[e].i] += edges_[e].G;
    Rtot[edges_[e].j] += edges_[e].G;
  }
  fix_total_C->do_reverse_comm();

  for (int i = 0; i < nlocal; ++i) {
    double Gsum = (std::isfinite(Rtot[i]) && Rtot[i] > 0.0) ? Rtot[i] : 0.0;
    if (Gw && (mask[i] & groupbit) && type_is_conductive(type[i]-1) &&
        std::isfinite(Gw[i]) && Gw[i] > 1.0e-40)
      Gsum += Gw[i];
    Rtot[i] = (Gsum > 0.0) ? 1.0/Gsum : 0.0;
  }
  fix_total_C->do_forward_comm();

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

inline bool FixHeatGranElectric::type_is_conductive(int itype) const
{
  const double s0 = (sigma0_ ? sigma0_[itype] : e_cond[itype]);
  return std::isfinite(s0) && (s0 > 1e-30);
}

inline double FixHeatGranElectric::sigma_at_typeT(int itype, double T) const
{
  if (!type_is_conductive(itype)) return 0.0;

  const double s0 = (sigma0_ ? sigma0_[itype] : e_cond[itype]);
  const double a  = (alpha_  ? alpha_[itype]  : 0.0);

  double denom = 1.0 + a*(T - Tref_);
  if (denom < 1e-6) denom = 1e-6;

  double s = s0 / denom;
  if (!std::isfinite(s) || s < sigma_floor_) s = sigma_floor_;
  return s;
}

inline double FixHeatGranElectric::aspot_fraction_for_types(int ti, int tj) const
{
  double f = aspot_fraction_global_;
  if (aspotFrac_)
    f = 0.5*(std::max(0.0, aspotFrac_[ti]) + std::max(0.0, aspotFrac_[tj]));
  if (f <= 0.0) return 0.0;
  if (f > 1.0)  return 1.0;
  return f;
}

void FixHeatGranElectric::recoverG0G1(double Gw, double cM, double V0, double V1,
                                      double& G0, double& G1)
{
  G0 = 0.0;
  G1 = 0.0;
  if (!(Gw > 1.0e-40) || !std::isfinite(cM)) return;
  const double dV = V1 - V0;
  if (std::fabs(dV) < 1.0e-12) { G0 = Gw; return; }

  double g1 = (cM - V0 * Gw) / dV;
  if (g1 < 0.0) g1 = 0.0;
  if (g1 > Gw ) g1 = Gw;
  G1 = g1;
  G0 = Gw - g1;
}

template<int HISTFLAG, int CONTACTAREA>
double FixHeatGranElectric::conductanceLocalIJ(const int i, const int j, double &fi)
{
  fi = 0.5;
  const int ti = atom->type[i] - 1;
  const int tj = atom->type[j] - 1;

  if (!type_is_conductive(ti) || !type_is_conductive(tj)) return 0.0;
  const double sigI = sigma_at_typeT(ti, Temp[i]);
  const double sigJ = sigma_at_typeT(tj, Temp[j]);
  if (!(sigI > 0.0) || !(sigJ > 0.0)) return 0.0;

  const double *xi = atom->x[i], *xj = atom->x[j];
  const double ri = atom->radius[i], rj = atom->radius[j];
  const double rsum = ri + rj;
  if (!(rsum > 0.0)) return 0.0;

  const double dx = xj[0]-xi[0], dy = xj[1]-xi[1], dz = xj[2]-xi[2];
  const double r  = std::sqrt(dx*dx + dy*dy + dz*dz);
  double delta    = rsum - r;

  if (!(delta > 0.0)) {
    if (!(r < elec_contact_factor_ * rsum)) return 0.0;
    delta = 0.0;
  } else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_OVERLAP && deltan_ratio_) {
    delta *= deltan_ratio_[ti][tj];
  }

  const double Reff = (ri * rj) / rsum;
  if (!(Reff > 0.0) || !std::isfinite(Reff)) return 0.0;

  const double delta_floor = elec_min_overlap_frac_ * Reff;
  if (delta < delta_floor) delta = delta_floor;
  if (!(delta > 0.0)) return 0.0;

  double a_nom = 0.0;
  if (CONTACTAREA == CONDUCTION_CONTACT_AREA_OVERLAP) {
    a_nom = std::sqrt(delta * Reff);
    const double a_cap = 0.35 * std::min(ri, rj);
    if (a_nom > a_cap) a_nom = a_cap;
  } else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_CONSTANT) {
    const double A = fixed_contact_area_;
    if (!(A > 0.0)) return 0.0;
    a_nom = std::sqrt(A / M_PI);
  } else {
    a_nom = std::min(ri, rj);
  }
  if (!(a_nom > 0.0) || !std::isfinite(a_nom)) return 0.0;

  const double fA = aspot_fraction_for_types(ti, tj);
  if (!(fA > 0.0) || !std::isfinite(fA)) return 0.0;

  const double a = a_nom * std::sqrt(fA);
  if (!(a > 0.0) || !std::isfinite(a)) return 0.0;

  const double Ri = 1.0 / (4.0 * a * sigI);
  const double Rj = 1.0 / (4.0 * a * sigJ);
  const double R  = Ri + Rj;
  if (!(R > 0.0) || !std::isfinite(R)) return 0.0;

  const double G = 1.0 / R;
  fi = Ri / R;

  if (!std::isfinite(G) || G < Cmin_branch_) return 0.0;
  return G;
}

template<int HISTFLAG, int CONTACTAREA>
void FixHeatGranElectric::buildEdgeCache()
{
  const int nloc = atom->nlocal;
  int    *mask = atom->mask;
  int    *type = atom->type;
  tagint *tag  = atom->tag;

  const int inum   = pair_gran->list->inum;
  int *ilist       = pair_gran->list->ilist;
  int *numneigh    = pair_gran->list->numneigh;
  int **firstneigh = pair_gran->list->firstneigh;

  edges_.clear();
  esolve_.clear();
  if ((int)edges_.capacity() < 8*nloc) edges_.reserve(8*nloc);

  for (int ii = 0; ii < inum; ++ii) {
    const int i = ilist[ii];
    if (i >= nloc) continue;
    if (!(mask[i] & groupbit)) continue;
    if (!type_is_conductive(type[i]-1)) continue;

    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; ++jj) {
      const int j = jlist[jj] & NEIGHMASK;
      if (!(mask[j] & groupbit)) continue;
      if (!type_is_conductive(type[j]-1)) continue;

      bool owner;
      if (j < nloc) owner = true;
      else          owner = (tag[i] < tag[j]);
      if (!owner) continue;

      double fi = 0.5;
      const double G = conductanceLocalIJ<HISTFLAG,CONTACTAREA>(i, j, fi);
      if (!(G > 0.0) || !std::isfinite(G)) continue;

      edges_.push_back(EdgeCG{i, j, G, fi});
    }
  }
}

int FixHeatGranElectric::uf_find(std::vector<int>& p, int x)
{
  int r = x;
  while (p[r] != r) r = p[r];
  while (p[x] != r) {
    const int next = p[x];
    p[x] = r;
    x = next;
  }
  return r;
}

void FixHeatGranElectric::uf_union(std::vector<int>& p, std::vector<int>& sz, int a, int b)
{
  a = uf_find(p,a);
  b = uf_find(p,b);
  if (a == b) return;
  if (sz[a] < sz[b]) std::swap(a,b);
  p[b] = a;
  sz[a] += sz[b];
}

template<int HISTFLAG, int CONTACTAREA>
void FixHeatGranElectric::detectedComponentsParallel()
{
  if (fix_pwER)    fix_pwER->do_forward_comm();
  if (fix_pwEP)    fix_pwEP->do_forward_comm();
  if (fix_cMatrix) fix_cMatrix->do_forward_comm();

  const int nlocal = atom->nlocal;
  const int nall   = nlocal + atom->nghost;
  int *type = atom->type;
  int *mask = atom->mask;
  const bool logStep = (update->ntimestep % std::max(1, electric_log_stride_) == 0);

  double *Gw  = fix_pwER->vector_atom;
  double *Vw  = fix_pwEP->vector_atom;
  double *cM  = fix_cMatrix->vector_atom;
  double *cid = fix_compID->vector_atom;

  for (int i = 0; i < nall; ++i) cid[i] = 0.0;

  double vmin_loc =  std::numeric_limits<double>::infinity();
  double vmax_loc = -std::numeric_limits<double>::infinity();
  int nwall_loc = 0;

  for (int i = 0; i < nlocal; ++i) {
    if (!(mask[i] & groupbit))          continue;
    if (!type_is_conductive(type[i]-1)) continue;
    if (!(Gw[i] > 1.0e-40))             continue;
    if (!std::isfinite(Vw[i]))          continue;
    vmin_loc = std::min(vmin_loc, Vw[i]);
    vmax_loc = std::max(vmax_loc, Vw[i]);
    ++nwall_loc;
  }

  double v0, v1;
  int nwall_glb;
  MPI_Allreduce(&vmin_loc,  &v0,        1, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(&vmax_loc,  &v1,        1, MPI_DOUBLE, MPI_MAX, world);
  MPI_Allreduce(&nwall_loc, &nwall_glb, 1, MPI_INT,    MPI_SUM, world);

  if (!std::isfinite(v0)) { v0 = 0.0; v1 = 0.0; }
  if (!std::isfinite(v1)) v1 = v0;

  V0_global_    = v0;
  V1_global_    = v1;
  two_terminal_ = (std::fabs(v1 - v0) > 1.0e-12);

  if (comm->me == 0 && screen && (update->ntimestep % 1000 == 0)) {
    fprintf(screen,"[ELECTRIC] step=%lld  wall_contacts=%d  V0=%.4g  V1=%.4g  two_terminal=%d\n",
            (long long)update->ntimestep, nwall_glb, v0, v1, (int)two_terminal_);
    fflush(screen);
  }

  if (!two_terminal_) {
    fix_compID->do_forward_comm();
    return;
  }

  fix_rhsG->set_all(0.0);
  fix_diagFlag->set_all(0.0);
  double *reachAnode   = fix_rhsG->vector_atom;
  double *reachCathode = fix_diagFlag->vector_atom;

  for (int i = 0; i < nlocal; ++i) {
    if (!(mask[i] & groupbit))          continue;
    if (!type_is_conductive(type[i]-1)) continue;
    if (!(Gw[i] > 1.0e-40))             continue;
    double G0, G1;
    recoverG0G1(Gw[i], cM[i], v0, v1, G0, G1);
    if (G0 > 0.0) reachCathode[i] = 1.0;
    if (G1 > 0.0) reachAnode[i]   = 1.0;
  }
  fix_rhsG->do_forward_comm();
  fix_diagFlag->do_forward_comm();

  uf_p_.resize(nall);
  uf_sz_.assign(nall, 1);
  for (int i = 0; i < nall; ++i) uf_p_[i] = i;
  for (size_t e = 0; e < edges_.size(); ++e)
    uf_union(uf_p_, uf_sz_, edges_[e].i, edges_[e].j);

  const int RMAX = 1000;
  int rounds = 0;
  long long prev_total = -1;

  while (rounds < RMAX) {
    uf_rootA_.assign(nall, 0.0);
    uf_rootK_.assign(nall, 0.0);
    for (int m = 0; m < nall; ++m) {
      const int r = uf_find(uf_p_, m);
      if (reachAnode[m]   > 0.5) uf_rootA_[r] = 1.0;
      if (reachCathode[m] > 0.5) uf_rootK_[r] = 1.0;
    }
    for (int m = 0; m < nall; ++m) {
      const int r = uf_find(uf_p_, m);
      if (uf_rootA_[r] > 0.5) reachAnode[m]   = 1.0;
      if (uf_rootK_[r] > 0.5) reachCathode[m] = 1.0;
    }

    fix_rhsG->do_reverse_comm();
    fix_diagFlag->do_reverse_comm();
    for (int i = 0; i < nlocal; ++i) {
      reachAnode[i]   = (reachAnode[i]   > 0.5) ? 1.0 : 0.0;
      reachCathode[i] = (reachCathode[i] > 0.5) ? 1.0 : 0.0;
    }
    fix_rhsG->do_forward_comm();
    fix_diagFlag->do_forward_comm();

    long long cnt_loc = 0;
    for (int i = 0; i < nlocal; ++i) {
      if (reachAnode[i]   > 0.5) ++cnt_loc;
      if (reachCathode[i] > 0.5) ++cnt_loc;
    }
    long long cnt_glb = 0;
    MPI_Allreduce(&cnt_loc, &cnt_glb, 1, MPI_LONG_LONG, MPI_SUM, world);

    ++rounds;
    if (cnt_glb == prev_total) break;
    prev_total = cnt_glb;
  }

  if (rounds >= RMAX && comm->me == 0 && screen)
    fprintf(screen,"[ELECTRIC] warning: component merge hit cap %d\n", RMAX);

  std::vector<int> degree(nall, 0);
  for (size_t e = 0; e < edges_.size(); ++e) {
    ++degree[edges_[e].i];
    ++degree[edges_[e].j];
  }

  double *connClass = fix_deltaE ? fix_deltaE->vector_atom : NULL;

  int n_active = 0, n_wall_active = 0, n_isolated = 0, n_floating = 0, n_oneSide = 0;
  long long degree_sum = 0;

  for (int i = 0; i < nlocal; ++i) {
    const bool a = (reachAnode[i]   > 0.5);
    const bool k = (reachCathode[i] > 0.5);

    int cls;
    if (a && k)            cls = 3;
    else if (a || k)       cls = 2;
    else if (degree[i] > 0) cls = 1;
    else                   cls = 0;

    if (connClass) connClass[i] = (double)cls;
    degree_sum += degree[i];

    if (cls == 3) {
      cid[i] = 1.0;
      ++n_active;
      if (Gw[i] > 1.0e-40) ++n_wall_active;
    }
    else if (cls == 2) ++n_oneSide;
    else if (cls == 1) ++n_floating;
    else               ++n_isolated;
  }
  fix_compID->do_forward_comm();
  if (fix_deltaE) fix_deltaE->do_forward_comm();

  int counts_loc[5] = { n_active, n_wall_active, n_oneSide, n_floating, n_isolated };
  int counts_glb[5];
  MPI_Allreduce(counts_loc, counts_glb, 5, MPI_INT, MPI_SUM, world);
  long long degree_sum_glb = 0;
  MPI_Allreduce(&degree_sum, &degree_sum_glb, 1, MPI_LONG_LONG, MPI_SUM, world);

  if (comm->me == 0 && screen && logStep) {
    const long long ntot = (long long)atom->natoms;
    const double zbar = ntot ? (double)degree_sum_glb / (double)ntot : 0.0;
    fprintf(screen,"[ELECTRIC-DET] active=%d  wall_active=%d  rounds=%d | stub=%d  floating=%d  isolated=%d  <z_elec>=%.2f\n",
            counts_glb[0], counts_glb[1], rounds, counts_glb[2], counts_glb[3], counts_glb[4], 2.0*zbar);
    fflush(screen);
  }
}

void FixHeatGranElectric::buildLocalIC0()
{
  const int nloc = atom->nlocal;

  ic_map_.assign(nloc, -1);
  ic_inv_.clear();
  for (int i = 0; i < nloc; ++i)
    if (active_[i]) {
      ic_map_[i] = (int)ic_inv_.size();
      ic_inv_.push_back(i);
    }

  const int n = (int)ic_inv_.size();
  ic_n_  = n;
  ic_ok_ = false;
  if (n == 0) return;

  std::vector<int> rowCount(n, 0);
  std::vector<int> rowOf, colOf;
  std::vector<double> gOf;
  rowOf.reserve(edges_.size());
  colOf.reserve(edges_.size());
  gOf.reserve(edges_.size());

  const size_t nedge = edges_.size();
  for (size_t t = 0; t < esolve_.size(); ++t) {
    const int eidx = esolve_[t];
    if (eidx < 0 || (size_t)eidx >= nedge) continue;
    const EdgeCG& ed = edges_[eidx];
    if (ed.j >= nloc) continue;
    int a = ic_map_[ed.i], b = ic_map_[ed.j];
    if (a < 0 || b < 0 || a == b) continue;
    if (a < b) std::swap(a,b);
    rowOf.push_back(a);
    colOf.push_back(b);
    gOf.push_back(ed.G);
    ++rowCount[a];
  }

  icLp_.assign(n+1, 0);
  for (int r = 0; r < n; ++r) icLp_[r+1] = icLp_[r] + rowCount[r];
  const int nnz = icLp_[n];
  icLj_ .assign(nnz, 0);
  icLx0_.assign(nnz, 0.0);
  icLx_ .assign(nnz, 0.0);
  {
    std::vector<int> pos(icLp_.begin(), icLp_.begin()+n);
    for (size_t t = 0; t < rowOf.size(); ++t) {
      const int p = pos[rowOf[t]]++;
      icLj_[p]  = colOf[t];
      icLx0_[p] = -gOf[t];
    }
  }

  for (int r = 0; r < n; ++r) {
    for (int p = icLp_[r]+1; p < icLp_[r+1]; ++p) {
      const int    cj = icLj_[p];
      const double cv = icLx0_[p];
      int q = p-1;
      while (q >= icLp_[r] && icLj_[q] > cj) {
        icLj_[q+1]  = icLj_[q];
        icLx0_[q+1] = icLx0_[q];
        --q;
      }
      icLj_[q+1]  = cj;
      icLx0_[q+1] = cv;
    }
  }

  icCp_.assign(n+1, 0);
  for (int p = 0; p < nnz; ++p) ++icCp_[icLj_[p]+1];
  for (int k = 0; k < n; ++k) icCp_[k+1] += icCp_[k];
  icCi_.assign(nnz, 0);
  icCq_.assign(nnz, 0);
  {
    std::vector<int> pos(icCp_.begin(), icCp_.begin()+n);
    for (int r = 0; r < n; ++r)
      for (int p = icLp_[r]; p < icLp_[r+1]; ++p) {
        const int k = icLj_[p];
        const int q = pos[k]++;
        icCi_[q] = r;
        icCq_[q] = p;
      }
  }

  icAd_.assign(n, 0.0);
  for (int c = 0; c < n; ++c) icAd_[c] = cg_diag_[ic_inv_[c]];

  icLd_.assign(n, 0.0);
  ic_w_.assign(n, 0.0);
  ic_y_.assign(n, 0.0);
  ic_mark_.assign(n, -1);

  double shift = 0.0;
  for (int attempt = 0; attempt < 6 && !ic_ok_; ++attempt) {
    icLx_ = icLx0_;
    std::fill(ic_mark_.begin(), ic_mark_.end(), -1);
    bool ok = true;

    for (int i = 0; i < n && ok; ++i) {
      for (int p = icLp_[i]; p < icLp_[i+1]; ++p) {
        const int k = icLj_[p];
        ic_w_[k] = icLx0_[p];
        ic_mark_[k] = i;
      }
      double d = icAd_[i] * (1.0 + shift);

      for (int p = icLp_[i]; p < icLp_[i+1]; ++p) {
        const int k = icLj_[p];
        const double lik = ic_w_[k] / icLd_[k];
        ic_w_[k] = lik;
        d -= lik*lik;
        for (int q = icCp_[k]; q < icCp_[k+1]; ++q) {
          const int j = icCi_[q];
          if (j >= i) break;
          if (ic_mark_[j] != i) continue;
          ic_w_[j] -= lik * icLx_[icCq_[q]];
        }
      }

      if (!(d > 1.0e-300) || !std::isfinite(d)) { ok = false; break; }
      icLd_[i] = std::sqrt(d);
      for (int p = icLp_[i]; p < icLp_[i+1]; ++p) icLx_[p] = ic_w_[icLj_[p]];
    }

    ic_ok_ = ok;
    if (!ok) shift = (shift == 0.0) ? 1.0e-3 : shift*10.0;
  }
}

void FixHeatGranElectric::applyPrecond(const std::vector<double>& r, std::vector<double>& z)
{
  const int nloc = atom->nlocal;
  std::fill(z.begin(), z.begin()+nloc, 0.0);

  if (!ic_ok_ || ic_n_ == 0) {
    for (int i = 0; i < nloc; ++i)
      if (active_[i] && cg_diag_[i] > 0.0) z[i] = r[i] / cg_diag_[i];
    return;
  }

  const int n = ic_n_;
  for (int c = 0; c < n; ++c) ic_y_[c] = r[ic_inv_[c]];

  for (int i = 0; i < n; ++i) {
    double s = ic_y_[i];
    for (int p = icLp_[i]; p < icLp_[i+1]; ++p) s -= icLx_[p] * ic_y_[icLj_[p]];
    ic_y_[i] = s / icLd_[i];
  }
  for (int i = n-1; i >= 0; --i) {
    const double s = ic_y_[i] / icLd_[i];
    ic_y_[i] = s;
    for (int p = icLp_[i]; p < icLp_[i+1]; ++p) ic_y_[icLj_[p]] -= icLx_[p] * s;
  }

  for (int c = 0; c < n; ++c) z[ic_inv_[c]] = ic_y_[c];
}

template<int HISTFLAG, int CONTACTAREA>
void FixHeatGranElectric::iterativeDomainDecompSolve()
{
  esolve_.clear();
  if (!two_terminal_) return;

  fix_pwER->do_forward_comm();
  fix_cMatrix->do_forward_comm();
  fix_compID->do_forward_comm();
  fix_eP->do_forward_comm();

  const int nloc = atom->nlocal;
  const int nall = nloc + atom->nghost;
  const bool logStep = (update->ntimestep % std::max(1, electric_log_stride_) == 0);

  double *V   = fix_eP->vector_atom;
  double *Gw  = fix_pwER->vector_atom;
  double *cM  = fix_cMatrix->vector_atom;
  double *cid = fix_compID->vector_atom;

  double *sumPP  = fix_rhsG->vector_atom;
  double *diagPP = fix_diagFlag->vector_atom;
  zero_with_ghosts(fix_rhsG,     nall);
  zero_with_ghosts(fix_diagFlag, nall);

  const double V0    = V0_global_, V1 = V1_global_;
  const double Vmid  = 0.5*(V0 + V1);
  const double Vlo   = std::min(V0, V1);
  const double Vhi   = std::max(V0, V1);
  const double tinyC = 1.0e-40;

  active_.assign(nall, 0);
  for (int i = 0; i < nall; ++i) active_[i] = (cid[i] > 0.5) ? 1 : 0;

  int nsolve_loc = 0;
  for (int i = 0; i < nloc; ++i) nsolve_loc += active_[i];
  int nsolve_glb = 0;
  MPI_Allreduce(&nsolve_loc, &nsolve_glb, 1, MPI_INT, MPI_SUM, world);
  if (nsolve_glb == 0) {
    ic_n_  = 0;
    ic_ok_ = false;
    if (comm->me == 0 && screen && logStep) {
      fprintf(screen,"[ELECTRIC-SOL] no percolating path this step, skipping solve\n");
      fflush(screen);
    }
    return;
  }

  for (int i = 0; i < nloc; ++i) {
    if (!active_[i]) {
      if (Gw[i] > tinyC) {
        double G0, G1;
        recoverG0G1(Gw[i], cM[i], V0, V1, G0, G1);
        if      (G0 > 0.0 && G1 == 0.0) V[i] = V0;
        else if (G1 > 0.0 && G0 == 0.0) V[i] = V1;
      } else V[i] = eP0;
      continue;
    }
    if (!std::isfinite(V[i]) || V[i] < Vlo || V[i] > Vhi) {
      if (Gw[i] > tinyC) {
        double G0, G1;
        recoverG0G1(Gw[i], cM[i], V0, V1, G0, G1);
        const double Gs = G0 + G1;
        V[i] = (Gs > 0.0) ? (0.7*(G0*V0 + G1*V1)/Gs + 0.3*Vmid) : Vmid;
      } else V[i] = Vmid;
    }
  }
  fix_eP->do_forward_comm();

  esolve_.reserve(edges_.size());
  for (size_t e = 0; e < edges_.size(); ++e)
    if (active_[edges_[e].i] && active_[edges_[e].j]) esolve_.push_back((int)e);

  for (size_t t = 0; t < esolve_.size(); ++t) {
    const EdgeCG& ed = edges_[esolve_[t]];
    diagPP[ed.i] += ed.G;
    diagPP[ed.j] += ed.G;
  }
  fix_diagFlag->do_reverse_comm();

  cg_diag_.assign(nloc, 0.0);
  int n_zero_diag = 0;
  for (int i = 0; i < nloc; ++i) {
    if (!active_[i]) continue;
    const double dpp = std::isfinite(diagPP[i]) ? std::max(0.0, diagPP[i]) : 0.0;
    const double Gwi = std::isfinite(Gw[i])     ? std::max(0.0, Gw[i])     : 0.0;
    cg_diag_[i] = dpp + Gwi;
    if (!(cg_diag_[i] > 0.0)) ++n_zero_diag;
  }
  int n_zero_diag_glb = 0;
  MPI_Allreduce(&n_zero_diag, &n_zero_diag_glb, 1, MPI_INT, MPI_SUM, world);
  if (n_zero_diag_glb > 0 && comm->me == 0 && screen) {
    fprintf(screen,"[ELECTRIC-SOL] WARNING: %d active particles have diag=0\n", n_zero_diag_glb);
    fflush(screen);
  }

  buildLocalIC0();

  const int    KMAX = 2000;
  const double RTOL = 1.0e-8;

  cg_x_ .assign(nloc, 0.0);
  cg_r_ .assign(nloc, 0.0);
  cg_z_ .assign(nloc, 0.0);
  cg_p_ .assign(nloc, 0.0);
  cg_Ap_.assign(nloc, 0.0);
  cg_Vbk_.assign(nloc, 0.0);

  for (int i = 0; i < nloc; ++i) if (active_[i]) cg_x_[i] = V[i];

  auto matvec = [&](const std::vector<double>& x_in, std::vector<double>& y_out)
  {
    for (int i = 0; i < nloc; ++i) cg_Vbk_[i] = V[i];
    for (int i = 0; i < nloc; ++i) if (active_[i]) V[i] = x_in[i];
    fix_eP->do_forward_comm();

    for (int i = 0; i < nall; ++i) sumPP[i] = 0.0;
    for (size_t t = 0; t < esolve_.size(); ++t) {
      const EdgeCG& ed = edges_[esolve_[t]];
      sumPP[ed.i] += ed.G * V[ed.j];
      sumPP[ed.j] += ed.G * V[ed.i];
    }
    fix_rhsG->do_reverse_comm();

    for (int i = 0; i < nloc; ++i) {
      if (!active_[i]) { y_out[i] = 0.0; continue; }
      const double sp = std::isfinite(sumPP[i]) ? sumPP[i] : 0.0;
      y_out[i] = cg_diag_[i] * x_in[i] - sp;
    }
    for (int i = 0; i < nloc; ++i) V[i] = cg_Vbk_[i];
  };

  auto dotprod = [&](const std::vector<double>& a, const std::vector<double>& b) -> double
  {
    double s_loc = 0.0;
    for (int i = 0; i < nloc; ++i) if (active_[i]) s_loc += a[i]*b[i];
    double s_glb = 0.0;
    MPI_Allreduce(&s_loc, &s_glb, 1, MPI_DOUBLE, MPI_SUM, world);
    return s_glb;
  };

  matvec(cg_x_, cg_Ap_);
  for (int i = 0; i < nloc; ++i) {
    if (!active_[i]) continue;
    const double bi = std::isfinite(cM[i]) ? cM[i] : 0.0;
    cg_r_[i] = bi - cg_Ap_[i];
  }
  applyPrecond(cg_r_, cg_z_);
  for (int i = 0; i < nloc; ++i) if (active_[i]) cg_p_[i] = cg_z_[i];

  double rz_old = dotprod(cg_r_, cg_z_);
  const double rz0 = std::max(rz_old, 1.0e-30);

  int iter = 0;
  for (iter = 0; iter < KMAX; ++iter) {
    matvec(cg_p_, cg_Ap_);
    const double pAp = dotprod(cg_p_, cg_Ap_);
    if (!(pAp > 0.0)) break;

    const double alpha = rz_old / pAp;
    for (int i = 0; i < nloc; ++i) {
      if (!active_[i]) continue;
      cg_x_[i] += alpha * cg_p_[i];
      cg_r_[i] -= alpha * cg_Ap_[i];
    }

    applyPrecond(cg_r_, cg_z_);
    const double rz_new = dotprod(cg_r_, cg_z_);
    if (rz_new < RTOL * RTOL * rz0) { ++iter; break; }

    const double beta = rz_new / rz_old;
    for (int i = 0; i < nloc; ++i)
      if (active_[i]) cg_p_[i] = cg_z_[i] + beta * cg_p_[i];
    rz_old = rz_new;
  }

  for (int i = 0; i < nloc; ++i) if (active_[i]) V[i] = cg_x_[i];
  fix_eP->do_forward_comm();

  matvec(cg_x_, cg_Ap_);
  double rmax_loc = 0.0;
  for (int i = 0; i < nloc; ++i)
    if (active_[i]) rmax_loc = std::max(rmax_loc, std::fabs(cM[i] - cg_Ap_[i]));
  double rmax_glb = 0.0;
  MPI_Allreduce(&rmax_loc, &rmax_glb, 1, MPI_DOUBLE, MPI_MAX, world);

  if (comm->me == 0 && screen && (logStep || iter >= KMAX)) {
    fprintf(screen,"[ELECTRIC] CG solve: iters=%d  edges=%zu  nact=%d  max|r|=%.3e  precond=%s%s\n",
            iter, esolve_.size(), nsolve_glb, rmax_glb,
            ic_ok_ ? "IC0" : "Jacobi", (iter >= KMAX) ? "  (NOT CONVERGED)" : "");
    fflush(screen);
  }
}

template<int HISTFLAG, int CONTACTAREA>
void FixHeatGranElectric::calculateJouleHeat()
{
  static long long last_ts = -1;
  static int calls = 0;
  const long long ts = (long long)update->ntimestep;
  if (ts != last_ts) { last_ts = ts; calls = 0; }
  const bool first_call = (++calls == 1);

  updatePtrs();

  const int nlocal = atom->nlocal;
  const int nall   = nlocal + atom->nghost;
  int *type = atom->type;
  int *mask = atom->mask;

  double  *V    = fix_eP->vector_atom;
  double  *Gw   = fix_pwER->vector_atom;
  double  *cM   = fix_cMatrix->vector_atom;
  double  *cid  = fix_compID->vector_atom;
  double  *qj   = fix_qj->vector_atom;
  double **Iviz = fix_current->array_atom;

  zero_with_ghosts(fix_qj_pp, nall);
  zero_with_ghosts(fix_qj_pw, nall);
  for (int i = 0; i < nall; ++i) Iviz[i][0] = Iviz[i][1] = Iviz[i][2] = 0.0;
  double *qpp = fix_qj_pp->vector_atom;
  double *qpw = fix_qj_pw->vector_atom;

  if (!two_terminal_) {
    for (int i = 0; i < nlocal; ++i) {
      const double qold = std::isfinite(qj[i]) ? qj[i] : 0.0;
      qj[i] = 0.0;
      if (first_call && heatSource) heatSource[i] -= qold;
    }
    fix_tot_I->set_all(0.0);
    return;
  }

  fix_eP->do_forward_comm();
  fix_compID->do_forward_comm();

  const double V0 = V0_global_, V1 = V1_global_;
  double **x = atom->x;

  const size_t nedge = edges_.size();
  double Ppp_loc = 0.0;
  for (size_t t = 0; t < esolve_.size(); ++t) {
    const int eidx = esolve_[t];
    if (eidx < 0 || (size_t)eidx >= nedge) continue;
    const EdgeCG& ed = edges_[eidx];
    const int i = ed.i, j = ed.j;
    if (i < 0 || i >= nall || j < 0 || j >= nall) continue;

    const double dV = V[i] - V[j];
    if (!std::isfinite(dV)) continue;
    const double Q = ed.G * dV * dV;
    if (!(Q >= 0.0)) continue;

    Ppp_loc += Q;
    const double fi = (ed.fi >= 0.0 && ed.fi <= 1.0) ? ed.fi : 0.5;
    qpp[i] += fi * Q;
    qpp[j] += (1.0 - fi) * Q;

    const double dx = x[i][0]-x[j][0], dy = x[i][1]-x[j][1], dz = x[i][2]-x[j][2];
    const double r2 = dx*dx + dy*dy + dz*dz;
    if (r2 > 1.0e-30) {
      const double rinv = 1.0/std::sqrt(r2);
      const double Iij  = ed.G * dV;
      Iviz[i][0] -= Iij*dx*rinv; Iviz[i][1] -= Iij*dy*rinv; Iviz[i][2] -= Iij*dz*rinv;
      Iviz[j][0] -= Iij*dx*rinv; Iviz[j][1] -= Iij*dy*rinv; Iviz[j][2] -= Iij*dz*rinv;
    }
  }
  fix_qj_pp->do_reverse_comm();
  fix_current->do_reverse_comm();
  qpp = fix_qj_pp->vector_atom;

  double I0_loc = 0.0, I1_loc = 0.0, Ppw_loc = 0.0;
  for (int i = 0; i < nlocal; ++i) {
    if (!(mask[i] & groupbit))          continue;
    if (cid[i] <= 0.5)                  continue;
    if (!type_is_conductive(type[i]-1)) continue;
    if (!(Gw[i] > 1.0e-40))             continue;

    double G0, G1;
    recoverG0G1(Gw[i], cM[i], V0, V1, G0, G1);

    const double dV0 = V0 - V[i];
    const double dV1 = V1 - V[i];
    if (!std::isfinite(dV0) || !std::isfinite(dV1)) continue;

    I0_loc += G0 * dV0;
    I1_loc += G1 * dV1;

    const double Q = G0*dV0*dV0 + G1*dV1*dV1;
    if (std::isfinite(Q) && Q >= 0.0) {
      qpw[i] = Q;
      Ppw_loc += Q;
    }
  }

  for (int i = 0; i < nlocal; ++i) {
    const double qpp_i = std::isfinite(qpp[i]) ? qpp[i] : 0.0;
    const double qpw_i = std::isfinite(qpw[i]) ? qpw[i] : 0.0;
    const double qnew  = qpp_i + qpw_i;
    const double qold  = std::isfinite(qj[i]) ? qj[i] : 0.0;
    qj[i] = qnew;
    if (first_call && heatSource) {
      const double dq = qnew - qold;
      if (std::isfinite(dq)) heatSource[i] += dq;
    }
  }

  double sums_loc[4] = { Ppp_loc, Ppw_loc, I0_loc, I1_loc };
  double sums_glb[4];
  MPI_Allreduce(sums_loc, sums_glb, 4, MPI_DOUBLE, MPI_SUM, world);
  const double Ppp = sums_glb[0], Ppw = sums_glb[1];
  const double I0  = sums_glb[2], I1  = sums_glb[3];

  const double Idev = 0.5*(std::fabs(I0) + std::fabs(I1));
  double *totI = fix_tot_I->vector_atom;
  for (int i = 0; i < nlocal; ++i) totI[i] = Idev;

  if (comm->me == 0 && first_call && screen &&
      (update->ntimestep % std::max(1, electric_log_stride_) == 0)) {
    const double Pterm  = V0*I0 + V1*I1;
    const double Pdiss  = Ppp + Ppw;
    const double relMis = std::fabs(Pdiss - Pterm)/std::max(1.0e-20, std::fabs(Pdiss));
    const double kcl    = std::fabs(I0 + I1)/std::max(1.0e-20, Idev);
    fprintf(screen,"[ELECTRIC] V0=%.4g V1=%.4g I=%.4g (I0=%+.4g I1=%+.4g)  Ppp=%.4g Ppw=%.4g  "
                   "Pdiss=%.4g Pterm=%.4g  relMis=%.2e  KCL=%.2e\n",
            V0, V1, Idev, I0, I1, Ppp, Ppw, Pdiss, Pterm, relMis, kcl);
    fflush(screen);
  }
}

void FixHeatGranElectric::calculateRadiation()
{
  createStencils();

  int maxBounces = 100;

  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  double *radius = atom->radius;
  int *type = atom->type;
  double **x = atom->x;
  double *hf = fix_heatFlux->vector_atom;
  qr = fix_qr->vector_atom;

  int hitId, hitBin, ibin;
  double hitEmis, flux, areai, emisi, radi, tempi;
  double hitp[3], nextNormal[3], buffer3[3], d[3], o[3];
  double *ci;

  updatePtrs();

  for (int i = 0; i < nlocal + nghost; i++)
  {
    radi  = radius[i];
    ci    = x[i];
    emisi = emissivity_[type[i]-1];
    tempi = Temp[i];
    const double cg_i = force->cg(atom->type[i]);

    ibin = neighbor->public_coord2bin(ci);

    areai = MY_4PI * radi * radi;
    flux  = areai * emisi * Sigma * tempi * tempi * tempi * tempi * cg_i;

    qr[i] += -flux;
    hf[i] -= flux;

    randOnSphere(ci, radi, o, buffer3);
    randDir(buffer3, d);

    hitId = trace(i, ibin, o, d, buffer3, hitp);

    if (hitId != -1) {
      const double& sendflux = flux;
      hitEmis = emissivity_[type[hitId]-1];
      hf[hitId] += hitEmis * sendflux;
      qr[hitId] += hitEmis * sendflux;

      hitBin = neighbor->public_coord2bin(x[hitId]);

      sub3(hitp, x[hitId], buffer3);
      normalize3(buffer3, nextNormal);

      reflect(i, hitId, hitBin, hitp, nextNormal, sendflux, 1.0-hitEmis, maxBounces, buffer3);
    } else {
      hf[i] += cg_i * (areai * emisi * Sigma * TB * TB * TB * TB);
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

void FixHeatGranElectric::randDir(const double *n, double *d)
{
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

  if (dot3(d, n) < 0.0) {
    d[0] = -d[0];
    d[1] = -d[1];
    d[2] = -d[2];
  }
}

int FixHeatGranElectric::trace(int orig_id, int ibin, double *o, double *d, double *buffer3, double *hitp)
{
  int *binhead = neighbor->binhead;
  int *bins = neighbor->bins;
  int sx = neighbor->sx;
  int sy = neighbor->sy;
  int sz = neighbor->sz;
  int mbins = neighbor->mbins;

  int stencilbin, stbX, stbY, stbZ;
  int currX, currY, currZ;
  double distsq, distx, disty, distz, t;
  bool hit;

  bool hitFlag = false;
  double hitT  = 1.0e300;
  int hitId    = -1;

  double **x     = atom->x;
  double *radius = atom->radius;

  int *stencil = pair_gran->list->stencil;
  int nstencil = pair_gran->list->nstencil;
  int nstencil2D = 0;
  int var_nstencil;

  const int nall_atoms = atom->nlocal + atom->nghost;
  if (!binhead || !bins || mbins <= 0 || nstencil <= 0 || !stencil) {
    if (comm->me == 0 && screen) {
      fprintf(screen,"[RAD] bin data unusable (mbins=%d nstencil=%d), radiation ray dropped\n", mbins, nstencil);
      fflush(screen);
    }
    return -1;
  }
  if (!stencilLength || !binStencildx || !binStencilmdx ||
      !binStencildy  || !binStencilmdy || !binStencildz || !binStencilmdz)
    return -1;

  int currentBin = neighbor->public_coord2bin(o);
  if (currentBin < 0 || currentBin >= mbins) return -1;

  int dx = 0, dy = 0, dz = 0;
  int *currentStencil = stencil;
  bool check_boundary_only = false;

  while ((currentBin != -1) && (hitFlag == false)) {

    var_nstencil = check_boundary_only ? nstencil2D : nstencil;

    for (int k = 0; k < var_nstencil; k++) {
      stencilbin = currentBin + currentStencil[k];
      if (stencilbin < 0 || stencilbin >= mbins) continue;
      neighbor->public_bin2XYZ(stencilbin, stbX, stbY, stbZ);
      neighbor->public_bin2XYZ(currentBin, currX, currY, currZ);
      if (abs(stbX - currX) > sx || abs(stbY - currY) > sy || abs(stbZ - currZ) > sz)
        continue;

      for (int i = binhead[stencilbin]; i >= 0; i = bins[i]) {
        if (i >= nall_atoms) break;
        if (i == orig_id) continue;

        hit = intersectRaySphere(o, d, x[i], radius[i], t, buffer3);
        if (hit) {
          hitFlag = true;
          if (t < hitT) {
            hitT = t;
            hitId = i;
          }
        }
      }
    }

    if (hitFlag) {
      addscaled3(o,d,hitT,hitp);
      return hitId;
    }

    ibin       = currentBin;
    currentBin = nextBin(ibin, o, d, raypoint, dx, dy, dz);
    if (currentBin < 0 || currentBin >= mbins) return -1;

    check_boundary_only = true;
    if (dx == 1)       { currentStencil = binStencildx;  nstencil2D = stencilLength[0]; }
    else if (dx == -1) { currentStencil = binStencilmdx; nstencil2D = stencilLength[3]; }
    else if (dy == 1)  { currentStencil = binStencildy;  nstencil2D = stencilLength[1]; }
    else if (dy == -1) { currentStencil = binStencilmdy; nstencil2D = stencilLength[4]; }
    else if (dz == 1)  { currentStencil = binStencildz;  nstencil2D = stencilLength[2]; }
    else               { currentStencil = binStencilmdz; nstencil2D = stencilLength[5]; }

    distx = raypoint[0] - o[0];
    disty = raypoint[1] - o[1];
    distz = raypoint[2] - o[2];
    distsq = distx*distx + disty*disty + distz*distz;
    if (distsq >= cutGhostsq) return -1;
  }

  return -1;
}

void FixHeatGranElectric::reflect(int radID, int orig_id, int ibin, double *o, double *d, double flux,
                                  double accum_eps, int n, double *buffer3)
{
  const double influx = flux * accum_eps;
  if (n == 0 || accum_eps < 0.001) {
    heatFlux[radID] += influx;
    return;
  }

  double hitp[3];
  double dd[3];
  randDir(d, dd);
  int hitId = trace(orig_id, ibin, o, dd, buffer3, hitp);

  int *type = atom->type;
  if (hitId != -1) {
    double **x = atom->x;
    const double hitEmis = emissivity_[type[hitId]-1];
    heatFlux[hitId] += hitEmis * influx;

    int hitBin = neighbor->public_coord2bin(x[hitId]);
    double nextNormal[3];
    sub3(hitp, x[hitId], nextNormal);
    norm3(nextNormal);

    reflect(radID, hitId, hitBin, hitp, nextNormal, flux, (1.0-hitEmis) * accum_eps, n-1, buffer3);
  }
  else {
    const double radRad  = atom->radius[radID];
    const double cg_i    = force->cg(atom->type[radID]);
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

  dx = dy = dz = 0;
  neighbor->public_binBorders(ibin, xlo, xhi, ylo, yhi, zlo, zhi);

  if (d[0] != 0.0) {
    s = (xlo - o[0]) / d[0];
    p[0] = o[0] + s*d[0]; p[1] = o[1] + s*d[1]; p[2] = o[2] + s*d[2];
    if ((s > smax) && (ylo <= p[1] && p[1] <= yhi) && (zlo <= p[2] && p[2] <= zhi)) {
      smax = s; dy = dz = 0; dx = -1;
    }
    s = (xhi - o[0]) / d[0];
    p[0] = o[0] + s*d[0]; p[1] = o[1] + s*d[1]; p[2] = o[2] + s*d[2];
    if ((s > smax) && (ylo <= p[1] && p[1] <= yhi) && (zlo <= p[2] && p[2] <= zhi)) {
      smax = s; dy = dz = 0; dx = 1;
    }
  }
  if (d[1] != 0.0) {
    s = (ylo - o[1]) / d[1];
    p[0] = o[0] + s*d[0]; p[1] = o[1] + s*d[1]; p[2] = o[2] + s*d[2];
    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (zlo <= p[2] && p[2] <= zhi)) {
      smax = s; dx = dz = 0; dy = -1;
    }
    s = (yhi - o[1]) / d[1];
    p[0] = o[0] + s*d[0]; p[1] = o[1] + s*d[1]; p[2] = o[2] + s*d[2];
    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (zlo <= p[2] && p[2] <= zhi)) {
      smax = s; dx = dz = 0; dy = 1;
    }
  }
  if (d[2] != 0.0) {
    s = (zlo - o[2]) / d[2];
    p[0] = o[0] + s*d[0]; p[1] = o[1] + s*d[1]; p[2] = o[2] + s*d[2];
    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (ylo <= p[1] && p[1] <= yhi)) {
      smax = s; dx = dy = 0; dz = -1;
    }
    s = (zhi - o[2]) / d[2];
    p[0] = o[0] + s*d[0]; p[1] = o[1] + s*d[1]; p[2] = o[2] + s*d[2];
    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (ylo <= p[1] && p[1] <= yhi)) {
      smax = s; dx = dy = 0; dz = 1;
    }
  }
  p[0] = o[0] + smax*d[0];
  p[1] = o[1] + smax*d[1];
  p[2] = o[2] + smax*d[2];

  int nextBinId = neighbor->public_binHop(ibin,dx,dy,dz);
  if (nextBinId != ibin) return nextBinId;
  return -1;
}

bool FixHeatGranElectric::intersectRaySphere(const double *o, const double *d, const double *center, double radius,
                                             double &t, double *buffer3)
{
  double A, B, C, discr, q, t0, t1;
  t = 0.0;
  A = lensq3(d);
  sub3(o, center, buffer3);
  B = 2.0 * dot3(buffer3, d);
  C = lensq3(buffer3) - radius*radius;
  discr = B*B - 4.0*A*C;

  if (discr < 0.0) return false;

  if (B < 0.0) q = (-B - sqrt(discr)) * 0.5;
  else         q = (-B + sqrt(discr)) * 0.5;

  t0 = q / A;
  t1 = C / q;
  if (t0 > t1) std::swap(t0, t1);

  if (t1 <= 0.0) {
    t = 0.0;
    return false;
  }
  t = t0;
  return true;
}

void FixHeatGranElectric::createStencils()
{
  if (stencilLength == NULL) stencilLength = new int[6];

  int mbinx = neighbor->mbinx;
  int mbiny = neighbor->mbiny;
  int sx = neighbor->sx;
  int sy = neighbor->sy;
  int sz = neighbor->sz;
  double cutneighmaxsq = neighbor->cutneighmaxsq;

  int i, j, k, nstencil;
  bool hit;

  delete [] binStencildx;
  binStencildx = new int[(2*sy+1)*(2*sz+1)];
  nstencil = 0;
  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++) {
      hit = false;
      i = sx;
      while (!hit && i >= -sx) {
        if (neighbor->public_bin_distance(i,j,k) < cutneighmaxsq) {
          binStencildx[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        } else i--;
      }
    }
  stencilLength[0] = nstencil;

  delete [] binStencildy;
  binStencildy = new int[(2*sx+1)*(2*sz+1)];
  nstencil = 0;
  for (k = -sz; k <= sz; k++)
    for (i = -sx; i <= sx; i++) {
      hit = false;
      j = sy;
      while (!hit && j >= -sy) {
        if (neighbor->public_bin_distance(i,j,k) < cutneighmaxsq) {
          binStencildy[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        } else j--;
      }
    }
  stencilLength[1] = nstencil;

  delete [] binStencildz;
  binStencildz = new int[(2*sx+1)*(2*sy+1)];
  nstencil = 0;
  for (j = -sy; j <= sy; j++)
    for (i = -sx; i <= sx; i++) {
      hit = false;
      k = sz;
      while (!hit && k >= -sz) {
        if (neighbor->public_bin_distance(i,j,k) < cutneighmaxsq) {
          binStencildz[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        } else k--;
      }
    }
  stencilLength[2] = nstencil;

  delete [] binStencilmdx;
  binStencilmdx = new int[(2*sy+1)*(2*sz+1)];
  nstencil = 0;
  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++) {
      hit = false;
      i = -sx;
      while (!hit && i <= sx) {
        if (neighbor->public_bin_distance(i,j,k) < cutneighmaxsq) {
          binStencilmdx[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        } else i++;
      }
    }
  stencilLength[3] = nstencil;

  delete [] binStencilmdy;
  binStencilmdy = new int[(2*sx+1)*(2*sz+1)];
  nstencil = 0;
  for (k = -sz; k <= sz; k++)
    for (i = -sx; i <= sx; i++) {
      hit = false;
      j = -sy;
      while (!hit && j <= sy) {
        if (neighbor->public_bin_distance(i,j,k) < cutneighmaxsq) {
          binStencilmdy[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        } else j++;
      }
    }
  stencilLength[4] = nstencil;

  delete [] binStencilmdz;
  binStencilmdz = new int[(2*sx+1)*(2*sy+1)];
  nstencil = 0;
  for (j = -sy; j <= sy; j++)
    for (i = -sx; i <= sx; i++) {
      hit = false;
      k = -sz;
      while (!hit && k <= sz) {
        if (neighbor->public_bin_distance(i,j,k) < cutneighmaxsq) {
          binStencilmdz[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        } else k++;
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
