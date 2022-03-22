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
};

void Setup (FLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
}

void Report (FLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
}


int main(int argc, char **argv) try
{
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    std::ifstream infile(filename.CStr());
    size_t Nproc = 0.75*omp_get_max_threads();
    bool hor = false;
    if (argc>=3) Nproc = atoi(argv[2]);
    if (argc>=4) hor   = atoi(argv[3]);

    String fileLBM;
    bool   Render = true;
    size_t Step = 1;
    bool   Run  = true;
    size_t nx0;
    size_t ny0;
    size_t nz0;
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
    infile >> Render;     infile.ignore(200,'\n');
    infile >> Step;       infile.ignore(200,'\n');
    infile >> Run;        infile.ignore(200,'\n');
    infile >> nx0;        infile.ignore(200,'\n');
    infile >> ny0;        infile.ignore(200,'\n');
    infile >> nz0;        infile.ignore(200,'\n');
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
    if (!Util::FileExists("after.h5")) throw new Fatal("Binary map not found \n");

    hid_t file_id;
    file_id = H5Fopen("after.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    if (!H5LTfind_dataset(file_id,"data")) throw new Fatal("The matrix name does not match \n");

    hsize_t dims[5];
    H5T_class_t *class_id; 
    size_t *type_size; 

    H5LTget_dataset_info (file_id, "data", dims, class_id, type_size);
    std::cout << " \n Dimensions of the imported array " << dims[3] << " " << dims[2] << " " << dims[1] << std::endl; 
    size_t ncells = dims[1]*dims[2]*dims[3];
    iVec3_t ndims = iVec3_t(dims[3],dims[2],dims[1]);
    if (nx0>=ndims[0]) throw new Fatal("Simulation box out of bounds in x direction \n");
    if (ny0>=ndims[1]) throw new Fatal("Simulation box out of bounds in y direction \n");
    if (nz0>=ndims[2]) throw new Fatal("Simulation box out of bounds in z direction \n");
    if (nx0+nx>=ndims[0])
    {
        nx = ndims[0] - nx0 - 1;
        std::cout << "Simulation box out of bounds in the x direction, domain changed to nx= " << nx << std::endl;
    }
    if (ny0+ny>=ndims[1])
    {
        ny = ndims[1] - ny0 - 1;
        std::cout << "Simulation box out of bounds in the y direction, domain changed to ny= " << ny << std::endl;
    }
    if (nz0+nz>=ndims[2])
    {
        nz = ndims[2] - nz0 - 1;
        std::cout << "Simulation box out of bounds in the z direction, domain changed to nz= " << nz << std::endl;
    }

    float * Gamma = new float[ncells];
    float * Gplot = new float[nx*ny*nz];
    iVec3_t ndimp(nx,ny,nz);
    H5LTread_dataset_float(file_id,"data",Gamma);
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        size_t idx = FLBM::Pt2idx(iVec3_t(ndims[0]-1-nx0-ix,ndims[1]-1-ny0-iy,ndims[2]-1-nz0-iz),ndims);
        size_t idp = FLBM::Pt2idx(iVec3_t(ix,iy,iz),ndimp);
        if (fabs(Gamma[idx]-0.0)<1.0e-12) Gplot[idp] = 0.0;
        if (fabs(Gamma[idx]-1.0)<1.0e-12) Gplot[idp] = 1.0;
        if (fabs(Gamma[idx]-2.0)<1.0e-12) Gplot[idp] = 2.0;
    }

    H5Fclose(file_id);

    String fn(filekey);
    fn.append("_ct.h5");
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dimp[1];
    dimp[0] = nx*ny*nz;
    H5LTmake_dataset_float(file_id,"Gplot",1,dimp,Gplot);

    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"LBM_Mesh\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << nz << " " << ny << " " << nx << "\"/>\n";
    oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
    oss << "       </DataItem>\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << 1.0 << " " << 1.0  << " " << 1.0 << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Gplot\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Gplot\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";

    fn = filekey;
    fn.append("_ct.xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();

    return 0;
}
MECHSYS_CATCH


