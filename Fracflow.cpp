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
#include <mechsys/lbm/Domain.h>
#include <mechsys/dem/domain.h>

enum GradCase
{
    Gx,  ///< Gradient in the x direction solely
    Gy,  ///< Gradient in the x direction solely
    Gz   ///< Gradient in the x direction solely
};

struct UserData
{
    Array<Cell *> xmin;
    Array<Cell *> xmax;
    Array<Cell *> ymin;
    Array<Cell *> ymax;
    Array<Cell *> zmin;
    Array<Cell *> zmax;
    double        vmax;
    double        rho;
    std::ofstream oss_ss;       ///< file for stress strain data
};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t i=0;i<dat.xmin.Size();i++)
    {
        Cell * c = dat.xmin[i];
        if(c->IsSolid) continue;
        c->F[1] = 1.0/3.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*c->RhoBC);
        c->F[7] = 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2] +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*c->RhoBC);
        c->F[9] = 1.0/24.0*(-2*c->F[0]+20*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
        c->F[11]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]+20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*c->RhoBC);
        c->F[13]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
    }
    for (size_t i=0;i<dat.xmax.Size();i++)
    {
        Cell * c = dat.xmax[i];
        if(c->IsSolid) continue;
        c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-c->RhoBC));
        c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*c->RhoBC) ;
        c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
    }
    for (size_t i=0;i<dat.ymin.Size();i++)
    {
        Cell * c = dat.ymin[i];
        if(c->IsSolid) continue;
        c->F[3]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[13]- 2*c->F[2]- c->F[4]- 2*c->F[5]- 2*c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
        c->F[7]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[13]- 5*c->F[2]- 4*c->F[4]+ c->F[5]-  5*c->F[6]+ 20*c->F[8]+ 2*c->RhoBC);
        c->F[9]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[10]- 4*c->F[11]- 4*c->F[13]- 5*c->F[2]- 4*c->F[4]- 5*c->F[5]+ c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
        c->F[12]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]+ 20*c->F[11]- 4*c->F[13]+ c->F[2]- 4*c->F[4]-5*c->F[5]+ c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
        c->F[14]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]- 4*c->F[11]+ 20*c->F[13]+ c->F[2]- 4*c->F[4]+ c->F[5]-5*c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
    }
    for (size_t i=0;i<dat.ymax.Size();i++)
    {
        Cell * c = dat.ymax[i];
        if(c->IsSolid) continue;
        c->F[4]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[12]- 4*c->F[14]- 2*c->F[2]- c->F[3]- 2*c->F[5]- 2*c->F[6]- 4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
        c->F[8]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[12]- 4*c->F[14]+ c->F[2]- 4*c->F[3]- 5*c->F[5]+ c->F[6]+  20*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
        c->F[10]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[12]- 4*c->F[14]+ c->F[2]- 4*c->F[3]+ c->F[5]- 5*c->F[6]- 4*c->F[7]+ 20*c->F[9]+ 2*c->RhoBC);
        c->F[11]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[12]- 4*c->F[14]- 5*c->F[2]- 4*c->F[3]+ c->F[5]- 5*c->F[6]-4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
        c->F[13]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[12]+ 20*c->F[14]- 5*c->F[2]- 4*c->F[3]- 5*c->F[5]+ c->F[6]-4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
    }
    for (size_t i=0;i<dat.zmin.Size();i++)
    {
        Cell * c = dat.zmin[i];
        if(c->IsSolid) continue;
        c->F[5] = 1/3.0*(-2*c->F[0] - 2*c->F[1] - 4*c->F[12] - 4*c->F[13] - 2*c->F[2] - 2*c->F[3] - 2*c->F[4] - c->F[6] -  4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
        c->F[7] = 1/24.0*(-2*c->F[0] + c->F[1] - 4*c->F[12] - 4*c->F[13] - 5*c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] +  20*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
        c->F[10] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] - 4*c->F[13] + c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] + 20*c->F[9] + 2*c->RhoBC);
        c->F[11] = 1/24.0*(-2*c->F[0] + c->F[1] + 20*c->F[12] - 4*c->F[13] - 5*c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
        c->F[14] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] + 20*c->F[13] + c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
    }
    for (size_t i=0;i<dat.zmax.Size();i++)
    {
        Cell * c = dat.zmax[i];
        if(c->IsSolid) continue;
        c->F[6]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]- 2*c->F[2]- 2*c->F[3]- 2*c->F[4]- c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
        c->F[8]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]+ c->F[2]- 5*c->F[3]+ c->F[4]- 4*c->F[5]+ 20*c->F[7]+ 2*c->RhoBC);
        c->F[9]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[10]- 4*c->F[11]- 4*c->F[14]- 5*c->F[2]+ c->F[3]- 5*c->F[4]- 4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
        c->F[12]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]+ 20*c->F[11]- 4*c->F[14]+ c->F[2]+ c->F[3]- 5*c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
        c->F[13]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[10]- 4*c->F[11]+ 20*c->F[14]- 5*c->F[2]- 5*c->F[3]+ c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
    }
}

void Report (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("flux.res");
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss << Util::_10_6  <<  "Time" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "M" << std::endl;
    }
    Vec3_t Flux = OrthoSys::O;
    double M    = 0.0;
    size_t nc   = 0;
    for (size_t i=0;i<dom.Lat[0].Ncells;i++)
    {
        Cell * c = dom.Lat[0].Cells[i];
        if (c->IsSolid) continue;
        Vec3_t DF;
        double rho = c->VelDen(DF);
        Flux += rho*DF;
        M += rho;
        nc++;
    }
    Flux/=M;
    dat.oss_ss << dom.Time << Util::_8s << Flux(0) << Util::_8s << Flux(1) << Util::_8s << Flux(2) << Util::_8s << M/nc << std::endl;
}


int main(int argc, char **argv) try
{
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    std::ifstream infile(filename.CStr());
    size_t Nproc = 1;
    if (argc==3) Nproc = atoi(argv[2]);

    String fileLBM;
    String matrixLBM;
    bool   Render = true;
    size_t nx;
    size_t ny;
    size_t nz;
    double nu     = 1.0/6.0;
    double Tf     = 1.0e3;
    double dtOut  = 1.0e1;
    bool   oct    = true;
    double Giso   = 1.0;
    double Gdev   = 3.0;
    double th     = 0.0;
    double DPx;
    double DPy;
    double DPz;

    infile >> fileLBM;    infile.ignore(200,'\n');
    infile >> matrixLBM;  infile.ignore(200,'\n');
    infile >> Render;     infile.ignore(200,'\n');
    infile >> nx;         infile.ignore(200,'\n');
    infile >> ny;         infile.ignore(200,'\n');
    infile >> nz;         infile.ignore(200,'\n');
    infile >> nu;         infile.ignore(200,'\n');
    infile >> Tf;         infile.ignore(200,'\n');
    infile >> dtOut;      infile.ignore(200,'\n');
    infile >> oct;        infile.ignore(200,'\n');
    if (oct)
    {
        infile >> Giso;     infile.ignore(200,'\n');
        infile >> Gdev;     infile.ignore(200,'\n');
        infile >> th;       infile.ignore(200,'\n');
        DPx    = Giso/3.0 + 2.0/3.0*Gdev*sin(M_PI*th/180.0-2.0*M_PI/3.0);
        DPy    = Giso/3.0 + 2.0/3.0*Gdev*sin(M_PI*th/180.0);
        DPz    = Giso/3.0 + 2.0/3.0*Gdev*sin(M_PI*th/180.0+2.0*M_PI/3.0);
    }
    else
    {
        infile >> DPx;      infile.ignore(200,'\n');
        infile >> DPy;      infile.ignore(200,'\n');
        infile >> DPz;      infile.ignore(200,'\n');
    }
    if (!Util::FileExists(fileLBM)) throw new Fatal("Binary map not found \n");

    std::ofstream parfile("param.res");
    parfile << Util::_8s << DPx/(3.0*nx)     << Util::_8s << DPy/(3.0*ny)     << Util::_8s << DPz/(3.0*nz)  << std::endl;
    //parfile << Util::_8s << Inet(0)   << Util::_8s << Inet(1)     << Util::_8s << Inet(2)  << std::endl;
    //parfile << Util::_8s << nx        << Util::_8s << ny          << Util::_8s << nz       << std::endl;
    //parfile << Util::_8s << 1.0-Dom.Lat[0].SolidFraction()        << Util::_8s << poresize << std::endl;
    //parfile << Util::_8s << mean_area << Util::_8s << mean_volume << Util::_8s << mean_sa  << std::endl;
    parfile.close();

    hid_t file_id;
    file_id = H5Fopen(fileLBM.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (!H5LTfind_dataset(file_id,matrixLBM.CStr())) throw new Fatal("The matrix name does not match \n");;

    float * Gamma = new float[nx*ny*nz];
    H5LTread_dataset_float(file_id,matrixLBM.CStr(),Gamma);

    LBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;

    // set initial condition
    double rho = 1.0;
    for (size_t i=0;i<Dom.Lat[0].Ncells;i++)
    {
        if (Gamma[i]<0.5) Dom.Lat[0].Cells[i]->IsSolid = true;
        Dom.Lat[0].Cells[i]->Initialize(rho,OrthoSys::O);
    }


	// set inner obstacle
    for (int i=0;i<nx;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0   ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0   ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,nz-1))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,nz-1))->IsSolid = true;
    }
    for (int i=0;i<ny;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(0   ,i,0   ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0   ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(0   ,i,nz-1))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(nx-1,i,nz-1))->IsSolid = true;
    }
    for (int i=0;i<nz;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(0   ,0   ,i))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(nx-1,0   ,i))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(0   ,ny-1,i))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(nx-1,ny-1,i))->IsSolid = true;
    }

    //Set boundary conditions
    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    {
        dat.zmin.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,0   )));
        dat.zmax.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,nz-1)));
        Dom.Lat[0].GetCell(iVec3_t(i,j,0   ))->RhoBC = rho+0.5*DPz;
        Dom.Lat[0].GetCell(iVec3_t(i,j,nz-1))->RhoBC = rho-0.5*DPz;
    }
    for (int i=0;i<nx;i++)
    for (int j=0;j<nz;j++)
    {
        dat.ymin.Push(Dom.Lat[0].GetCell(iVec3_t(i,0   ,j)));
        dat.ymax.Push(Dom.Lat[0].GetCell(iVec3_t(i,ny-1,j)));
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,j))->RhoBC = rho+0.5*DPy;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,j))->RhoBC = rho-0.5*DPy;
    }
    for (int i=0;i<ny;i++)
    for (int j=0;j<nz;j++)
    {
        dat.xmin.Push(Dom.Lat[0].GetCell(iVec3_t(0   ,i,j)));
        dat.xmax.Push(Dom.Lat[0].GetCell(iVec3_t(nx-1,i,j)));
        Dom.Lat[0].GetCell(iVec3_t(0   ,i,j))->RhoBC = rho+0.5*DPx;
        Dom.Lat[0].GetCell(iVec3_t(nx-1,i,j))->RhoBC = rho-0.5*DPx;
    }

    Dom.WriteXDMF("initial");


    //Solving
    Dom.Solve(Tf,dtOut,Setup,Report,filekey.CStr(),Render,Nproc);
    dat.oss_ss.close();

    Dom.WriteXDMF("final");

    return 0;
}
MECHSYS_CATCH


