/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine, DCS Computing GmbH, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. See the top-level README and LICENSE files.

    Copyright 2024-     LCT - UGENT, Ghent
    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(heat/gran/electric,FixHeatGranElectric)

#else

#ifndef LMP_FIX_HEATGRAN_ELECTRIC_H
#define LMP_FIX_HEATGRAN_ELECTRIC_H

#include "fix_heat_gran.h"
#include "fix_wall_gran.h"
#include <vector>

namespace LAMMPS_NS {

  class FixHeatGranElectric : public FixHeatGran {
  public:
    FixHeatGranElectric(class LAMMPS *, int, char **);
    ~FixHeatGranElectric();
    virtual void post_create();
    virtual void pre_delete(bool);

    int setmask();
    void init();
    virtual void pre_force(int vflag);
    virtual void post_force(int vflag);
    virtual void end_of_step();

    virtual void cpl_evaluate(class ComputePairGranLocal *);
    void register_compute_pair_local(class ComputePairGranLocal *);
    void unregister_compute_pair_local(class ComputePairGranLocal *);

    virtual void updatePtrs();

  protected:
    struct EdgeCG { int i, j; double G; double fi; };

    int iarg_ = 0;

    template <int,int> void post_force_eval(int,int);
    template <int,int> void electric_end_step_impl();
    template <int,int> void buildEdgeCache();
    template <int,int> void detectedComponentsParallel();
    template <int,int> void iterativeDomainDecompSolve();
    template <int,int> void calculateJouleHeat();
    template <int,int> double conductanceLocalIJ(const int i, const int j, double &fi);

    bool type_is_conductive(int itype) const;
    double sigma_at_typeT(int itype, double T) const;
    double aspot_fraction_for_types(int ti, int tj) const;
    void recoverG0G1(double Gw, double cM, double V0, double V1, double &G0, double &G1);
    int  uf_find(std::vector<int> &p, int x);
    void uf_union(std::vector<int> &p, std::vector<int> &sz, int a, int b);
    void buildLocalIC0();
    void applyPrecond(const std::vector<double> &r, std::vector<double> &z);

    void calculateRadiation();
    void randOnSphere(const double *c, double r, double *ansP, double *ansD);
    void randDir(const double *n, double *d);
    int  trace(int orig_id, int ibin, double *o, double *d, double *buffer3, double *hitp);
    void reflect(int radID, int orig_id, int ibin, double *o, double *d, double flux, double accum_eps, int n, double *buffer3);
    int  nextBin(int ibin, const double *o, const double *d, double *p, int &dx, int &dy, int &dz);
    bool intersectRaySphere(const double *o, const double *d, const double *center, double radius, double &t, double *buffer3);
    void createStencils();
    void addscaled3(const double *v1, const double *v2, double s, double *ans);

    class FixPropertyGlobal *fix_eConductivity_ = NULL;
    double *eConductivity_ = NULL;

    bool store_contact_data_ = false;
    class FixPropertyAtom *fix_conduction_contact_area_ = NULL;
    class FixPropertyAtom *fix_n_conduction_contacts_ = NULL;
    class FixPropertyAtom *fix_wall_heattransfer_coeff_ = NULL;
    class FixPropertyAtom *fix_wall_temperature_ = NULL;
    class FixPropertyAtom *fix_wall_electricPotential_ = NULL;
    double *conduction_contact_area_ = NULL;
    double *n_conduction_contacts_ = NULL;
    double *wall_heattransfer_coeff_ = NULL;
    double *wall_temp_ = NULL;
    double *wall_ep_ = NULL;

    int area_calculation_mode_ = 0;
    double fixed_contact_area_ = 0.;
    int area_correction_flag_ = 0;
    double const* const* deltan_ratio_ = NULL;

    bool electricMode_ = false;
    const double *e_cond = NULL;
    class FixWallGran *wallfix = NULL;

    class FixPropertyAtom *fix_eP = NULL;
    class FixPropertyAtom *fix_eP_old = NULL;
    class FixPropertyAtom *fix_deltaE = NULL;
    class FixPropertyAtom *fix_qj = NULL;
    class FixPropertyAtom *fix_qj_pp = NULL;
    class FixPropertyAtom *fix_qj_pw = NULL;
    class FixPropertyAtom *fix_qr = NULL;
    class FixPropertyAtom *fix_tot_I = NULL;
    class FixPropertyAtom *fix_current = NULL;
    class FixPropertyAtom *fix_compID = NULL;
    class FixPropertyAtom *fix_rhsG = NULL;
    class FixPropertyAtom *fix_diagFlag = NULL;
    class FixPropertyAtom *fix_total_C = NULL;
    class FixPropertyAtom *fix_countCheck = NULL;
    class FixPropertyAtom *fix_cMatrix = NULL;
    class FixPropertyAtom *fix_pwER = NULL;
    class FixPropertyAtom *fix_pwEP = NULL;

    double *eP = NULL;
    double *eP_old = NULL;
    double *qr = NULL;
    double eP0 = 0.;
    bool eP_needs_init_ = false;

    double joule_heat_startTime = 0.;
    bool temp_control_flag_ = false;
    double temp_control_ = 0.;
    int electric_log_stride_ = 1;

    double aspot_fraction_global_ = 1.0;
    class FixPropertyGlobal *fix_aspotFrac_ = NULL;
    const double *aspotFrac_ = NULL;
    double Cmin_branch_ = 0.0;
    const double *sigma0_ = NULL;
    const double *alpha_ = NULL;
    double Tref_ = 293.15;
    double sigma_floor_ = 1e-12;

    double elec_contact_factor_ = 1.0;
    double elec_min_overlap_frac_ = 0.0;
    bool elec_min_overlap_set_ = false;

    std::vector<EdgeCG> edges_;
    std::vector<int> esolve_;
    double V0_global_ = 0.;
    double V1_global_ = 0.;
    bool two_terminal_ = false;
    std::vector<char> active_;

    std::vector<int> uf_p_, uf_sz_;
    std::vector<double> uf_rootA_, uf_rootK_;

    std::vector<double> cg_diag_, cg_x_, cg_r_, cg_z_, cg_p_, cg_Ap_, cg_Vbk_;

    int ic_n_ = 0;
    bool ic_ok_ = false;
    std::vector<int> ic_map_, ic_inv_, ic_mark_;
    std::vector<int> icLp_, icLj_, icCp_, icCi_, icCq_;
    std::vector<double> icLx0_, icLx_, icLd_, icAd_, ic_w_, ic_y_;

    bool radiation_flag_ = false;
    const double *emissivity_ = NULL;
    double TB = 0.;
    double Sigma = 5.670373E-8;
    double cutGhost = 0.10;
    double cutGhostsq = 0.01;
    class RanMars *RGen = NULL;
    double *raypoint = NULL;
    int *stencilLength = NULL;
    int *binStencildx = NULL;
    int *binStencilmdx = NULL;
    int *binStencildy = NULL;
    int *binStencilmdy = NULL;
    int *binStencildz = NULL;
    int *binStencilmdz = NULL;
  };

}

#endif
#endif
