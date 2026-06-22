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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_ELECTRIC,electric,3)
#else

#ifndef COHESION_MODEL_ELECTRIC_H_
#define COHESION_MODEL_ELECTRIC_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include <cmath>
#include <algorithm>
#include "math_extra_liggghts.h"
#include "global_properties.h"
#include "fix_property_atom.h"
#include "neighbor.h"
#include "update.h"

namespace MODEL_PARAMS
{

    //dielectric constant of paricles
    inline static ScalarProperty* createEpsilonParticle(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* epsilonPScalar = MODEL_PARAMS::createScalarProperty(registry, "epsilonP", caller);
      return epsilonPScalar;
    }

    // dielectric constant of gas
    inline static ScalarProperty* createEpsilonGas(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* epsilonGScalar = MODEL_PARAMS::createScalarProperty(registry, "epsilonG", caller);
      return epsilonGScalar;
    }

}

namespace LIGGGHTS {

namespace ContactModels {
    const double epsilon0 = 8.854e-12; // vaccum dielectric constant
    template<>
    class CohesionModel<COHESION_ELECTRIC> : public CohesionModelBase {
    public:
        CohesionModel(LAMMPS *lmp, IContactHistorySetup *hsetup, class ContactModelBase *c) :
            CohesionModelBase(lmp, hsetup, c),

	    epsilonP(0.0),
	    epsilonG(0.0),
	    fix_eField_(0)

            
        {}

        void registerSettings(Settings& settings)
        {
            settings.registerOnOff("tangential_reduce", tangentialReduce_, false);
        }

        inline void postSettings(IContactHistorySetup *hsetup, ContactModelBase *cmb) {}

        void connectToProperties(PropertyRegistry & registry)
        {
            registry.registerProperty("epsilonP", &MODEL_PARAMS::createEpsilonParticle);
            registry.registerProperty("epsilonG", &MODEL_PARAMS::createEpsilonGas);

	    registry.connect("epsilonP", epsilonP, "cohesion_model electric");
            registry.connect("epsilonG", epsilonG, "cohesion_model electric");
    	    
	    //accessing electric field per particle 
	    fix_eField_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("pE","property/atom","vector",0,0,"cohesion_model electric"));
	    if(!fix_eField_)
	    {
		const char* fixarg[12];
		fixarg[0]="pE";
		fixarg[1]="all";
		fixarg[2]="property/atom";
		fixarg[3]="pE";
		fixarg[4]="vector";
		fixarg[5]="no";
		fixarg[6]="yes";
		fixarg[7]="no";
		fixarg[8]="0.";
		fixarg[9]="0.";
		fixarg[10]="0.";
		fix_eField_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg), "cohesion_model electric");
	    }
	    if(!fix_eField_)
		error->one(FLERR, "Electric field not found for the electric cohesion model");

            // error checks on coarsegraining
            //if (force->cg_active())  error->cg(FLERR, "cohesion model electric");
            
            if (force->cg_active()) {
		// Optional: warn once instead of aborting
		if (comm->me == 0)	error->warning(FLERR,"cohesion model electric: running with coarse-graining; ensure cutoffs/parameters are CG-consistent");
	     }
	}


        void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
        {            
	    // Check if the distance is less than 3 mm
            if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;
            
            const double rcut_factor = 1.5;                 // example: 1.5*(ri+rj)
            const double rcut = rcut_factor * (sidata.radi + sidata.radj);
            if (sidata.r >= rcut) return;
            
            //if (sidata.r < 0.003) {
                // Calculate orientation angle Theta
       		int step = update->ntimestep;             
                double v1_x = sidata.v_i[0];
                double v1_y = sidata.v_i[1];
                double v1_z = sidata.v_i[2];

                double v2_x = sidata.v_j[0];
                double v2_y = sidata.v_j[1];
                double v2_z = sidata.v_j[2];
                               
                // Calculate the dot product
                double dotProduct = v1_x * v2_x + v1_y * v2_y + v1_z * v2_z;
		
		// distance vector delta coming from the sidata.delta which is the delx, dely and delz
		double deltaX = sidata.delta[0];
		double deltaY = sidata.delta[1];
		double deltaZ = sidata.delta[2];
		double deltaMag = sqrt(sidata.rsq);

		// magnitude of electric field Emag 
		//double Emag = sqrt(electricX*electricX + electricY*electricY + electricZ*electricZ);               
		double * const * const eField = fix_eField_-> array_atom;
		double * const eField_i = eField[sidata.i];
		double * const eField_j = eField[sidata.j];

		// Calculating the magnitude of electric field strength 
		//double Emag = std::sqrt(eField_i[1] * eField_i[1] + eField_i[2] * eField_i[2] + eField_i[3] * eField_i[3]);
		
		const double Exi = eField_i[0];
		const double Eyi = eField_i[1];
		const double Ezi = eField_i[2];
		double Emag = std::sqrt(Exi*Exi + Eyi*Eyi + Ezi*Ezi);


		//Calculating the angle between the electric field and the center-to-center distance vector
		//double Theta_prod = (electricX* deltaX+ electricY*deltaY + electricZ*deltaZ)/(deltaMag*Emag);
		
		/*double Theta_prod;
		if(Emag!=0 && eField_i[1] !=0 && eField_i[2] !=0 && eField_i[3] !=0)  Theta_prod = (eField_i[1]* deltaX+ eField_i[2]*deltaY + eField_i[3]*deltaZ)/(deltaMag*Emag);
 		else if (Emag==0 && eField_i[1]==0 && eField_i[2]==0 && eField_i[3]==0)  Theta_prod =0;

		double Theta_update = acos(Theta_prod);
                 // Calculate the magnitudes
                double mag1 = std::sqrt(v1_x * v1_x + v1_y * v1_y + v1_z * v1_z);
                double mag2 = std::sqrt(v2_x * v2_x + v2_y * v2_y + v2_z * v2_z);

                // Ensure denominators are not zero to avoid division by zero
                double Theta;
                if (mag2 == 0.0 || mag1 == 0.0) {
                    Theta = 0.0; // Return 0 if either magnitude is zero
                }else{
		    Theta =acos(dotProduct / (mag1 * mag2));
		}*/
		
		const double r = sidata.r;
		if (r <= 0.0) return;

		double cosTheta = 0.0;
		if (Emag > 0.0) {
		    cosTheta = (Exi*sidata.delta[0] + Eyi*sidata.delta[1] + Ezi*sidata.delta[2]) / (Emag * r);
		    cosTheta = std::max(-1.0, std::min(1.0, cosTheta));
		}
		const double Theta = std::acos(cosTheta);


                // declare di-electric of free space (epsilon_0)
                double epsilon0 = 8.85418782e-12;                        //THE VALUE NEEDS TO CHANGE TO THE VALUE OF VACCUM'S DIELECTRIC
			
		// Calculate force using FORCER and FORCETHETA functions     // UPDATE 
                double beta = (epsilonP - epsilonG)/(epsilonP + (2.0 * epsilonG));         //whats epsilon_c
		double C_factor = 10; 			//correction factor
		double d_p = 2 * sidata.radi; 		//diameter of a sphere	
	
                double cteForce = C_factor * 3/16 * M_PI * epsilon0 * 1.0 * epsilonG * pow(beta,2) * pow(Emag, 2) * pow( d_p,2) * pow(d_p/sidata.r, 4);
                double f0RCER1 = cteForce * (3* pow( cos(Theta) , 2.0) -1);
                double f0RCETHETA1 = cteForce * (sin(2*Theta));
	

                if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;               
                const double FviscN = f0RCER1;
                const double FviscT_over_vt = f0RCETHETA1;
                //std::cout << "FORCER: " << f0RCER1<< std::endl;
                //std::cout << "FORCETHETA: " << f0RCETHETA1 << std::endl;
		
	
                // tangential force components
                const double Ft1 = FviscT_over_vt*sidata.vtr1;
                const double Ft2 = FviscT_over_vt*sidata.vtr2;
                const double Ft3 = FviscT_over_vt*sidata.vtr3;

                // torques
                const double tor1 = sidata.en[1] * Ft3 - sidata.en[2] * Ft2;
                const double tor2 = sidata.en[2] * Ft1 - sidata.en[0] * Ft3;
                const double tor3 = sidata.en[0] * Ft2 - sidata.en[1] * Ft1;

                // add to fn, Ft
               if(tangentialReduce_) sidata.Fn += FviscN; 
               //sidata.Ft += ...

               // apply normal and tangential force
               const double fx = (FviscN) * sidata.en[0] + Ft1;
               const double fy = (FviscN) * sidata.en[1] + Ft2;
               const double fz = (FviscN) * sidata.en[2] + Ft3;

               // return resulting forces
               if(sidata.is_wall) {
                const double area_ratio = sidata.area_ratio;
                i_forces.delta_F[0] += fx * area_ratio;
                i_forces.delta_F[1] += fy * area_ratio;
                i_forces.delta_F[2] += fz * area_ratio;
                i_forces.delta_torque[0] += -sidata.cri * tor1 * area_ratio;
                i_forces.delta_torque[1] += -sidata.cri * tor2 * area_ratio;
                i_forces.delta_torque[2] += -sidata.cri * tor3 * area_ratio;
               } else {
                i_forces.delta_F[0] += fx;
                i_forces.delta_F[1] += fy;
                i_forces.delta_F[2] += fz;
                i_forces.delta_torque[0] += -sidata.cri * tor1;
                i_forces.delta_torque[1] += -sidata.cri * tor2;
                i_forces.delta_torque[2] += -sidata.cri * tor3;

                j_forces.delta_F[0] -= fx;
                j_forces.delta_F[1] -= fy;
                j_forces.delta_F[2] -= fz;
                j_forces.delta_torque[0] += -sidata.crj * tor1;
                j_forces.delta_torque[1] += -sidata.crj * tor2;
                j_forces.delta_torque[2] += -sidata.crj * tor3;
               }
            //}
        }

        inline void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
        void beginPass(SurfacesIntersectData&, ForceData&, ForceData&) {}
        void endPass(SurfacesIntersectData&, ForceData&, ForceData&) {}

        void surfacesClose(SurfacesCloseData& scdata, ForceData&, ForceData&)
        {
            if (scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
        }

    private:
        double Emag; // Use pointer for Emag to allow dynamic update
	double epsilonP, epsilonG;
        bool tangentialReduce_;
	FixPropertyAtom *fix_eField_;
    };
}
}

#endif // COHESION_MODEL_ELECTRIC_H_
#endif
