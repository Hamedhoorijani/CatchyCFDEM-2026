/* ----------------------------------------------------------------------
   LIGGGHTS / LAMMPS — fix_heat_gran_electric6.h
   Electric + thermal coupling for granular materials
------------------------------------------------------------------------- */

#include <vector>
#include <unordered_map>

#include <map>
#include <set>

#ifdef FIX_CLASS
FixStyle(heat/gran/electric,FixHeatGranElectric)
FixStyle(heat/gran,FixHeatGranElectric)
#else

#ifndef LMP_FIX_HEATGRAN_ELECTRIC_H
#define LMP_FIX_HEATGRAN_ELECTRIC_H

#include "fix_heat_gran.h"

namespace LAMMPS_NS {

class ComputePairGranLocal;
class FixWallGran;
class FixPropertyAtom;
class FixPropertyGlobal;
class RanMars;

class FixHeatGranElectric : public FixHeatGran {
public:
  FixHeatGranElectric(class LAMMPS *, int, char **);
  ~FixHeatGranElectric() override;

  void post_create() override;
  void pre_delete(bool) override;

  int  setmask() override;
  void init() override;
  void pre_force(int vflag) override;
  void post_force(int vflag) override;

  void cpl_evaluate(ComputePairGranLocal *caller);
  void register_compute_pair_local(ComputePairGranLocal *);
  void unregister_compute_pair_local(ComputePairGranLocal *);
  
  void updatePtrs();

  void end_of_step();  // new
  // Electric core (kept intact in your .cpp)
  template<int HISTFLAG,int CONTACTAREA>
  double conductanceLocalIJ(int i, int j);

  template<int HISTFLAG,int CONTACTAREA>
  void detectedComponentsParallel();

  template<int HISTFLAG,int CONTACTAREA>
  void iterativeDomainDecompSolve();

  template<int HISTFLAG,int CONTACTAREA>
  void calculateJouleHeat();

  // Radiation
  void   calculateRadiation();
  void   addscaled3(const double *, const double*, double, double*);
  bool   intersectRaySphere(const double *, const double *, const double *, double, double &, double *);
  int    nextBin(int, const double *, const double *, double *, int &, int &, int &);
  int    trace(int, int, double *, double *, double *, double *);
  void   randDir(const double *, double *);
  void   randOnSphere(const double *, double, double *, double *);
  void   reflect(int, int, int, double *, double *, double, double, int, double *);
  void   createStencils();

  // in class FixHeatGranElectric:
  double aspot_fraction_global_;                // default 1.0 (no reduction)
  FixPropertyGlobal *fix_aspotFrac_;            // optional per-atom-type values
  const double *aspotFrac_;                           // pointer returned by fix_aspotFrac_
  double Cmin_branch_;                          // optional tiny cutoff for G (stability)

  FixPropertyAtom *fix_total_C  = nullptr;
  FixPropertyAtom *fix_countCheck = nullptr;
  
  const double *sigma0_;       // peratomtype σ at Tref  [S/m]
  const double *alpha_;        // peratomtype temp coeff [1/K]
  double        Tref_;         // scalar reference temp  [K]
  double        sigma_floor_;  // floor, e.g. 1e-12 S/m

  inline bool type_is_conductive(int) const; 
  inline double sigma_at_typeT(int, double) const;

  // logging stride for CSV export (>=1)
  int electric_log_stride_;

private:

  template<int HISTFLAG, int CONTACTAREA>
  inline void electric_end_step_impl();  // new helper
  
  inline double aspot_fraction_for_types(int, int) const;

  // -------- order matches constructor init-list in your .cpp ------------
  // global/per-type thermal
  FixPropertyGlobal *fix_eConductivity_ = nullptr;
  double *eConductivity_ = nullptr;

  bool    store_contact_data_ = false;
  FixPropertyAtom *fix_conduction_contact_area_ = nullptr;
  FixPropertyAtom *fix_n_conduction_contacts_   = nullptr;
  FixPropertyAtom *fix_wall_heattransfer_coeff_ = nullptr;
  FixPropertyAtom *fix_wall_temperature_        = nullptr;
  double *conduction_contact_area_ = nullptr;
  double *n_conduction_contacts_   = nullptr;
  double *wall_heattransfer_coeff_ = nullptr;
  double *wall_temp_               = nullptr;

  // electric / wall-electrics
  FixPropertyAtom *fix_wall_electricPotential_ = nullptr;
   
  FixWallGran     *wallfix = nullptr;
  double *wall_ep_ = nullptr;

  // particle potentials (new/old) and flags
  FixPropertyAtom *fix_eP     = nullptr;
  FixPropertyAtom *fix_eP_old = nullptr;
  double *eP     = nullptr;
  double *eP_old = nullptr;
  double  eP0    = 0.0;
  bool eP_needs_init_;	//flag for restart and initialization of Temp

  // Joule / radiation accumulators
  FixPropertyAtom *fix_qj     = nullptr;
  FixPropertyAtom *fix_qj_pp  = nullptr;
  FixPropertyAtom *fix_qj_pw  = nullptr;
  FixPropertyAtom *fix_qr     = nullptr;
  FixPropertyAtom *fix_tot_I     = nullptr;
  double *qj = nullptr;
  double *qr = nullptr;
  double *tot_I = nullptr;

  // misc flags / diagnostics
  bool hasEPArgs = true;
  bool electricMode_ = false;

  // area modes for thermal contact
  enum { CONDUCTION_CONTACT_AREA_OVERLAP=0,
         CONDUCTION_CONTACT_AREA_CONSTANT,
         CONDUCTION_CONTACT_AREA_PROJECTION };
  int    area_calculation_mode_ = CONDUCTION_CONTACT_AREA_OVERLAP;
  double fixed_contact_area_    = 0.0;
  int    area_correction_flag_  = 0;

  double joule_heat_startTime = 0.0;

  bool   temp_control_flag_ = false;
  double temp_control_      = 0.0;

  // radiation params
  double TB          = 0.0;
  double Sigma       = 5.670373e-8;
  double cutGhost    = 0.10;
  double cutGhostsq  = 0.01;
  bool   radiation_flag_ = false;

  // 🔧 MISSING BEFORE → used in init(), calculateRadiation(), reflect()
  FixPropertyGlobal *fix_emissivity_ = nullptr; // (optional handle)
  const double *emissivity_ = nullptr;          // per-type emissivity

  // component labeling helpers
  int *compID = nullptr;

  // area-correction look-up
  const double * const * deltan_ratio_ = nullptr;

  // parser index
  int iarg_ = 0;

  // per-type electrical conductivity
  FixPropertyGlobal *fix_econductivity_ = nullptr;
  const double *e_cond = nullptr;

  // per-atom props referenced in .cpp
  FixPropertyAtom *fix_current  = nullptr;  // vector
  FixPropertyAtom *fix_compID   = nullptr;  // scalar
  FixPropertyAtom *fix_rhsG     = nullptr;  // scalar
  FixPropertyAtom *fix_diagFlag = nullptr;  // scalar
  FixPropertyAtom *fix_pwER     = nullptr;  // Σ Cw
  FixPropertyAtom *fix_pwEP     = nullptr;  // eff Vw (diag)
  FixPropertyAtom *fix_pwEC     = nullptr;  // conductive flag
  FixPropertyAtom *fix_cMatrix  = nullptr;  // Σ Cw·Vw
  FixPropertyAtom *fix_deltaE   = nullptr;  // optional

  // RNG + radiation helpers
  RanMars *RGen = nullptr;
  double  *raypoint = nullptr;

  // radiation stencils
  int *stencilLength = nullptr;
  int *binStencildx  = nullptr;
  int *binStencilmdx = nullptr;
  int *binStencildy  = nullptr;
  int *binStencilmdy = nullptr;
  int *binStencildz  = nullptr;
  int *binStencilmdz = nullptr;

  // templated heat eval
  template<int HISTFLAG,int CONTACTAREA>
  void post_force_eval(int vflag,int cpl_flag);
};

} // namespace LAMMPS_NS

#endif // LMP_FIX_HEATGRAN_ELECTRIC_H
#endif // FIX_CLASS

