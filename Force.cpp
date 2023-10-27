/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Torres                                     *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

//STD
#include<iostream>
#include <list>

// MechSys
#include <mechsys/flbm/Domain.h>

struct UserData
{
    double            rho;
    double             Tf;
    Array<double>     Vel;
    Array<size_t>     Sid;
    std::ofstream oss_ss;       ///< file for stress strain data
};

void Setup (FLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));

    double scale = 1.0;

    // Cells with prescribed velocity (left boundary)
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t i=0; i<dom.Ndim(1); ++i) //y-axis
    for (size_t j=0; j<dom.Ndim(2); ++j) //z-axis
    {
        //double c  = 1;
        double * f = dom.F[0][0][i][j];
        double rho = (f[0]+f[3]+f[4]+f[5]+f[6]+ 2.0*(f[2]+f[8]+f[10]+f[12]+f[14]))/(1.0-scale*dat.Vel[j]/dom.Cs);
        f[1] = f[2] + (2.0/3.0)*rho*scale*dat.Vel[j]/dom.Cs;
        f[7] = f[8] + (1.0/12.0)*rho*scale*dat.Vel[j]/dom.Cs;
        f[9] = f[10] + (1.0/12.0)*rho*scale*dat.Vel[j]/dom.Cs;
        f[11]= f[12] + (1.0/12.0)*rho*scale*dat.Vel[j]/dom.Cs;
        f[13]= f[14] + (1.0/12.0)*rho*scale*dat.Vel[j]/dom.Cs;
        dom.Vel[0][0][i][j] = OrthoSys::O;
        dom.Rho[0][0][i][j] = 0.0;
        for (size_t k=0;k<dom.Nneigh;k++)
        {
            dom.Rho[0][0][i][j] +=  dom.F[0][0][i][j][k];
            dom.Vel[0][0][i][j] +=  dom.F[0][0][i][j][k]*dom.C[k];
        }
        dom.Vel[0][0][i][j] *= dom.Cs/dom.Rho[0][0][i][j];
    }

    // Cells with prescribed density (right boundary)
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t i=0; i<dom.Ndim(1); ++i) //y-axis
    for (size_t j=0; j<dom.Ndim(2); ++j) //z-axis
    {
        double * f = dom.F[0][dom.Ndim(0)-1][i][j];
        double vx = (f[0]+f[3]+f[4]+f[5]+f[6] + 2.0*(f[1]+f[7]+f[9]+f[11]+f[13]))/dat.rho - 1.0;
        f[2] = f[1] - (2.0/3.0)*dat.rho*vx; 
        f[8] = f[7] - (1.0/12.0)*dat.rho*vx;
        f[10]= f[9] - (1.0/12.0)*dat.rho*vx;
        f[12]= f[11] - (1.0/12.0)*dat.rho*vx; 
        f[14]= f[13] - (1.0/12.0)*dat.rho*vx;
        dom.Vel[0][dom.Ndim(0)-1][i][j] = OrthoSys::O;
        dom.Rho[0][dom.Ndim(0)-1][i][j] = 0.0;
        for (size_t k=0;k<dom.Nneigh;k++)
        {
            dom.Rho[0][dom.Ndim(0)-1][i][j] +=  dom.F[0][dom.Ndim(0)-1][i][j][k];
            dom.Vel[0][dom.Ndim(0)-1][i][j] +=  dom.F[0][dom.Ndim(0)-1][i][j][k]*dom.C[k];
        }
        dom.Vel[0][dom.Ndim(0)-1][i][j] *= dom.Cs/dom.Rho[0][dom.Ndim(0)-1][i][j];
    }

    // free slip (top boundary)
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t i=0; i<dom.Ndim(0); ++i) //x-axis
    for (size_t j=0; j<dom.Ndim(1); ++j) //y-axis
    {
        double * f = dom.F[0][i][j][dom.Ndim(2)-1];
        f[6]  = f[5];
        f[8]  = f[10];
        f[9]  = f[7];
        f[12] = f[14];
        f[13] = f[11];
        dom.Vel[0][i][j][dom.Ndim(2)-1] = OrthoSys::O;
        dom.Rho[0][i][j][dom.Ndim(2)-1] = 0.0;
        for (size_t k=0;k<dom.Nneigh;k++)
        {
            dom.Rho[0][i][j][dom.Ndim(2)-1] +=  dom.F[0][i][j][dom.Ndim(2)-1][k];
            dom.Vel[0][i][j][dom.Ndim(2)-1] +=  dom.F[0][i][j][dom.Ndim(2)-1][k]*dom.C[k];
        }
        dom.Vel[0][i][j][dom.Ndim(2)-1] *= dom.Cs/dom.Rho[0][i][j][dom.Ndim(2)-1];

    }
}

void Report (FLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("force.res");
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss << Util::_10_6  <<  "Time" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "M" << std::endl;
    }
    Vec3_t Force = OrthoSys::O;
    double M     = 0.0;
    for (size_t isid=0;isid<dat.Sid.Size();isid++)
    {
        iVec3_t iv;
        FLBM::idx2Pt(dat.Sid[isid],iv,dom.Ndim);
        double rho = dom.Rho[0][iv(0)][iv(1)][iv(2)];
        for (size_t k=0;k<dom.Nneigh;k++)
        {
            double Fvpp     = dom.Feq(dom.Op[k],rho,OrthoSys::O);
            double Fvp      = dom.Feq(k        ,rho,OrthoSys::O);
            double Omega    = dom.F[0][iv(0)][iv(1)][iv(2)][dom.Op[k]] - Fvpp - (dom.F[0][iv(0)][iv(1)][iv(2)][k] - Fvp);
            Force += Omega*dom.C[k];
            M     += Omega*dom.C[k](0)*iv(2);
        }
    }
    dat.oss_ss << dom.Time << Util::_8s << Force(0) << Util::_8s << Force(1) << Util::_8s << Force(2) << Util::_8s << M << std::endl;
}


int main(int argc, char **argv) try
{
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    std::ifstream infile(filename.CStr());
    size_t Nproc = 0.75*omp_get_max_threads();
    if (argc>=3) Nproc = atoi(argv[2]);

    String fileLBM;
    bool   Render = true;
    size_t Step = 1;
    double nu     = 1.0/6.0;
    double Tf     = 1.0e2;
    double dtOut  = 1.0e1;

    infile >> fileLBM;    infile.ignore(200,'\n');
    infile >> Render;     infile.ignore(200,'\n');
    infile >> Step;       infile.ignore(200,'\n');
    infile >> Tf;         infile.ignore(200,'\n');
    infile >> dtOut;      infile.ignore(200,'\n');

    if (!Util::FileExists(fileLBM)) throw new Fatal("Flow field file not found \n");

    hid_t file_id;
    file_id = H5Fopen(fileLBM.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);

    int data[1];
    H5LTread_dataset_int(file_id,"/Nx",data);
    size_t nx = data[0];
    H5LTread_dataset_int(file_id,"/Ny",data);
    size_t ny = data[0];
    H5LTread_dataset_int(file_id,"/Nz",data);
    size_t nz = data[0];

    size_t ncells = nx*ny*nz;

    float * Gamma = new float[  ncells];
    H5LTread_dataset_float(file_id,"/Gamma",Gamma);
    float * Rho   = new float[  ncells];
    H5LTread_dataset_float(file_id,"/Density_0",Rho);
    float * Vel   = new float[3*ncells];
    H5LTread_dataset_float(file_id,"/Velocity_0",Vel);

    FLBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;

    iVec3_t ndims(nx,ny,nz);

    Array<Array<size_t> > lnum(Nproc);

    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        size_t idx = FLBM::Pt2idx(iVec3_t(ix,iy,iz),ndims);
        if (Gamma[idx]==1) 
        {
            Dom.IsSolid[0][ix][iy][iz] = true;
            if (iz==0) lnum[omp_get_thread_num()].Push(idx);
        }
        Vec3_t vel(Vel[3*idx+0],Vel[3*idx+1],Vel[3*idx+2]);
        Dom.Initialize(0,iVec3_t(ix,iy,iz),Rho[idx],vel);
    }

    size_t tnum = 0;
    for (size_t np=0;np<Nproc;np++)
    {
        tnum += lnum[np].Size();
    }

    dat.Sid.Resize(tnum);
    size_t inum = 0;
    for (size_t np=0;np<Nproc;np++)
    {
        for (size_t npp=0;npp<lnum[np].Size();npp++)
        {
            dat.Sid[inum] = lnum[np][npp];
            inum++;
        }
    }


    dat.Vel.Resize(nz);
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t iz=0;iz<Dom.Ndim(2);iz++)
    {
        dat.Vel[iz] = Dom.Vel[0][0][ny/2][iz](0);
    }

    Dom.Step = Step;
    Dom.WriteXDMF("force_initial");
    dat.rho     = 1.0;

    //Solving
    dat.Tf = Tf;
    Dom.Solve(Tf,dtOut,Setup,Report,filekey.CStr(),Render,Nproc);
    dat.oss_ss.close();
//

    return 0;
}
MECHSYS_CATCH


