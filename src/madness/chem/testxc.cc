/*
  This file is part of MADNESS.
  
  Copyright (C) 2007,2010 Oak Ridge National Laboratory
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
  
  For more information please contact:
  
  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367
  
  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

#include <madness/world/MADworld.h>
#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>
#include <fstream>
#include "xcfunctional.h"

using namespace madness;

static std::string df_repo_functionals[] = {
"lda_x", 
"lda_c_vwn_rpa", 
"lda_c_vwn", 
"lda_c_pz", 
"lda_c_pw", 
"hyb_gga_xc_b3lyp", 
"gga_xc_hcth_93", 
"gga_xc_hcth_407", 
"gga_xc_hcth_147", 
"gga_xc_hcth_120", 
"gga_xc_edf1", 
"gga_xc_b97_2", 
"gga_xc_b97_1", 
"gga_xc_b97", 
"gga_x_pw91", 
"gga_x_pbe", 
"gga_x_ft97_b", 
"gga_x_b88", 
"gga_c_pw91", 
"gga_c_pbe", 
"gga_c_p86", 
"gga_c_lyp"};


struct xcfunc_data_point
{
  double rhoa, rhob;
  double sigmaaa, sigmaab, sigmabb;
  double zk;
  double vrhoa, vrhob;
  double vsigmaaa, vsigmaab, vsigmabb;
  double v2rhoa2, v2rhoab, v2rhob2;
};

std::vector<xcfunc_data_point> read_test_data(const std::string& dfname,
                                              bool spin_polarized)
{
  std::ifstream fstr(dfname.c_str());
  char buffer[120];
  fstr.getline(buffer, 120);
  fstr.getline(buffer, 120);

  std::string tmpstr;

  std::vector<xcfunc_data_point> dps;
  print(dfname.c_str());
  int i=1;
  while(!fstr.eof())
  {
    xcfunc_data_point dp;

    fstr >> tmpstr; fstr >> dp.rhoa;
    fstr >> tmpstr; fstr >> dp.rhob;
    fstr >> tmpstr; fstr >> dp.sigmaaa;
    fstr >> tmpstr; fstr >> dp.sigmaab;
    fstr >> tmpstr; fstr >> dp.sigmabb;



    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.zk;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.vrhoa;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.vrhob;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.vsigmaaa;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.vsigmaab;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.vsigmabb;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.v2rhoa2;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.v2rhoab;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.v2rhob2;


    for (int iskip = 0; iskip < 13; iskip++)
      fstr.getline(buffer,120);
    if (!spin_polarized)
    {
      if (std::abs(dp.rhoa-dp.rhob) <= 1e-10)
        dps.push_back(dp);
    }
    else
    {
      dps.push_back(dp);
    }
    i=i+1;
    if(i == 26) break;
  }

  return dps;

}

void test_xcfunctional(World& world)
{
    bool spin_polarized ;
    //spin_polarized = true;
    spin_polarized = false;

    int ispin = 0; //alpha=0 beta=1
    int istr = 0;

//vama5     0  "lda_x", 
//vama5     1	"lda_c_vwn_rpa", 
//vama5     2	"lda_c_vwn", 
//vama5     3	"lda_c_pz", 
//vama5     4	"lda_c_pw", 
//vama5     5	"hyb_gga_xc_b3lyp", 
//vama5     6	"gga_xc_hcth_93", 
//vama5     7	"gga_xc_hcth_407", 
//vama5     8	"gga_xc_hcth_147", 
//vama5     9	"gga_xc_hcth_120", 
//vama5    10	"gga_xc_edf1", 
//vama5    11	"gga_xc_b97_2", 
//vama5    12	"gga_xc_b97_1", 
//vama5    13	"gga_xc_b97", 
//vama5    14	"gga_x_pw91", 
//vama5    15	"gga_x_pbe", 

//vama5    16	"gga_x_ft97_b", 
//vama5    17	"gga_x_b88", 
//vama5    18	"gga_c_pw91", 
//vama5    19	"gga_c_pbe", 

//vama5    20	"gga_c_p86", 
//vama5    21	"gga_c_lyp"};

    XCfunctional xcfunc;
    std::string xcfuncstr = df_repo_functionals[istr];
    std::cout << "Testing exchange-correlation functional:  "<< xcfuncstr << std::endl;

    xcfuncstr += " 1.0";
    xcfunc.initialize(xcfuncstr,spin_polarized,world);

    std::string fpath("df_repo/");
    fpath += df_repo_functionals[istr];
    fpath += ".data";
    std::vector<xcfunc_data_point> dps = read_test_data(fpath.c_str(),spin_polarized);

    print("hola");
    Tensor<double> rhoa_t((long)dps.size());
    Tensor<double> rhob_t((long)dps.size());
    Tensor<double> sigmaaa_t((long)dps.size());
    Tensor<double> sigmaab_t((long)dps.size());
    Tensor<double> sigmabb_t((long)dps.size());
    Tensor<double> zk_t((long)dps.size());
    Tensor<double> vrhoa_t((long)dps.size());
    Tensor<double> vrhob_t((long)dps.size());
    Tensor<double> vsigmaaa_t((long)dps.size());
    Tensor<double> vsigmaab_t((long)dps.size());
    Tensor<double> vsigmabb_t((long)dps.size());

    Tensor<double> rhoa1_t((long)dps.size());
    Tensor<double> sigmaaa1_t((long)dps.size());

    std::vector<Tensor<double> > xc_args;
    for (unsigned int idp = 0; idp < dps.size(); idp++)
    {
      rhoa_t(idp) = dps[idp].rhoa;
      rhob_t(idp) = dps[idp].rhob;

      sigmaaa_t(idp) = dps[idp].sigmaaa;
      sigmaab_t(idp) = dps[idp].sigmaab;
      sigmabb_t(idp) = dps[idp].sigmabb;

      zk_t(idp) = dps[idp].zk;

      vrhoa_t(idp) = dps[idp].vrhoa;
      vrhob_t(idp) = dps[idp].vrhob;


      vsigmaaa_t(idp) = dps[idp].vsigmaaa;
      vsigmaab_t(idp) = dps[idp].vsigmaab;
      vsigmabb_t(idp) = dps[idp].vsigmabb;
    }
    if (spin_polarized)
    {
      xc_args.push_back(rhoa_t);
      xc_args.push_back(rhob_t);
      if (xcfunc.is_gga()) {
        xc_args.push_back(sigmaaa_t);
        xc_args.push_back(sigmaab_t);
        xc_args.push_back(sigmabb_t);
      }
    }
    else
    {
      xc_args.push_back(rhoa_t);
      if (xcfunc.is_gga()) {
         xc_args.push_back(sigmaaa_t);
      }
    }

    print("xc_args_size", xc_args.size());
    print("ispin ", ispin);
    print("spin polarized ", spin_polarized);
//vama1  std::cout << "Testing spin-polarized case: " << std::endl << std::endl;

    // xc local vr[0] and semilocal vr[1-3] potential
    std::vector<Tensor<double> > vr =  xcfunc.vxc(xc_args,ispin);

#if 0
    print("\n");
    if(what == 0)
    print("rhoa \t vr \t read\t comp   \n");
    if(what == 1)
    print("rhoa \t vsigmaa \t read\t comp   \n");
    if(what == 2)
    print("rhoa \t vsigmab \t read\t comp   \n");
    if(what == 3)
    print("rhoa \t zk \t read\t comp   \n");
    for (unsigned int idp = 0; idp < dps.size(); idp++)
    {
    if (spin_polarized)
    {
      if(what == 0)
      printf("%25.12e %25.12e  %25.12e %25.12e  \n",
          rhoa_t[idp], vrhoa_t[idp], dps[idp].vrhoa, vr[idp] );
      if(what == 1)
      printf("%i %25.12e %25.12e  %25.12e %25.12e  \n",
          idp, rhoa_t[idp], vsigmaaa_t[idp], dps[idp].vsigmaaa, vr[idp]);
      if(what == 2)
      printf("%i %25.12e %25.12e  %25.12e %25.12e  \n",
          idp, rhoa_t[idp], vsigmaab_t[idp], dps[idp].vsigmaab, vr[idp]);
      if(what == 3) 
      printf("%25.12e %25.12e  %25.12e %25.12e  \n",
          rhoa_t[idp], zk_t[idp], dps[idp].zk, vr[idp]);
    }
    else
    {
      if(what == 0)
      printf("%25.12e %25.12e  %25.12e %25.12e  \n",
          rhoa_t[idp], vrhoa_t[idp], dps[idp].vrhoa, vr[idp] );
      if(what == 1)
      printf("vsga %i %25.12e %25.12e  %25.12e %25.12e  \n",
          idp, rhob_t[idp], vsigmaaa_t[idp], dps[idp].vsigmaaa, vr[idp]);
        //  idp, rhob_t[idp], .5*(vsigmaaa_t[idp] + .5*vsigmaab_t[idp]), dps[idp].vsigmaaa, vr[idp]);
      if(what == 2)
       printf (" \n bad bad bad\n ") ;
      if(what ==3) 
      printf("%25.12e %25.12e  %25.12e %25.12e  \n",
          rhoa_t[idp], zk_t[idp], dps[idp].zk, vr[idp]);
    }
    }
    print("\n\n");
#endif

    if (xcfunc.is_spin_polarized())
    {
        printf("%25s %25s %25s %25s %25s %25s %25s %25s\n","#rhoa","rhob","sigmaaa","sigmaab","sigmabb","vrhoa (input)","vr (output)","vrhoa-vr");
        for (unsigned int idp = 0; idp < dps.size(); idp++)
        {
            printf("%25.12e %25.12e %25.12e %25.12e %25.12e %25.12e %25.12e %25.12e\n",
                    rhoa_t[idp], rhob_t[idp], sigmaaa_t[idp], sigmaab_t[idp], sigmabb_t[idp],
                    dps[idp].vrhoa, vr[0][idp],
                    std::abs(dps[idp].vrhoa - vr[0][idp]));
        }
    }
    else
    {
        printf("%25s %25s  %25s %25s   %25s\n","#rhoa","sigmaaa","vrhoa (input)","vr (output)","vrhoa-vr");
        for (unsigned int idp = 0; idp < dps.size(); idp++)
        {
            printf("%25.12e %25.12e  %25.12e %25.12e   %25.12e\n",
                    rhoa_t[idp], sigmaaa_t[idp], dps[idp].vrhoa, vr[0][idp],
                    std::abs(dps[idp].vrhoa - vr[0][idp]));
        }
    }
    print("\n\n");

}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);

    madness::World world(SafeMPI::COMM_WORLD);
    world.gop.fence();

    test_xcfunctional(world);

    madness::finalize();
    return 0;
}
