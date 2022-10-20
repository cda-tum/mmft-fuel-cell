/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "palabos3D.h"
#include "palabos3D.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor
#define ADESCRIPTOR descriptors::AdvectionDiffusionD3Q7Descriptor

#define DYNAMICS BGKdynamics<T, DESCRIPTOR>(omega)
#define ADYNAMICS AdvectionDiffusionBGKdynamics<T,ADESCRIPTOR>(adOmega)

plint extraLayer = 0.;
plint referenceDirection = 0.;
plint envelopeWidth = 1;
plint extendedEnvelopeWidth = 1;
plint blockSize = 2000000;

T simTime = 0;
plint maxIter = 0;
plint writeInterval = 0;
T epsilon = 0;
bool performOutput = false;

TriangleSet<T>* triangleSet = 0;

#define NMAX 150topDomain

const T pi = (T)4.*std::atan((T)1.);

// Initiate 2D field from STL file
template<typename T>
class domainInitializer3D : public BoxProcessingFunctional3D_S<T> {
public:
    domainInitializer3D(ScalarField3D<T>* write_)
        : write(write_)
    { }
    domainInitializer3D<T>* clone() const {
        return new domainInitializer3D<T>(*this);
    }
    virtual void process(Box3D domain, ScalarField3D<T>& from) {
        Dot3D rel = from.getLocation();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX){
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    write->get(iX+rel.x, iY+rel.y, iZ+rel.z) = from.get(iX,iY,iZ);
                }
            }
        }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const{
        modified[0] = modif::dataStructure;
    }
private:
    ScalarField3D<T>* write;
};

// Initiate 2D field from STL file
template<typename T>
class domainInitializer2D : public BoxProcessingFunctional3D_S<T> {
public:
    domainInitializer2D(ScalarField3D<T>* write_)
        : write(write_)
    { }
    domainInitializer2D<T>* clone() const {
        return new domainInitializer2D<T>(*this);
    }
    virtual void process(Box3D domain, ScalarField3D<T>& to) {
        Dot3D rel = to.getLocation();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX){
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ){
                    to.get(iX,iY,iZ) = write->get(iX+rel.x, iY+rel.y, iZ+rel.z);
                }
            }
        }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const{
        modified[0] = modif::dataStructure;
    }
private:
    ScalarField3D<T>* write;
};

// InstantiateDynamicsFunctional3D
template<typename T, template<typename U> class Descriptor>
class GetFlux : public BoxProcessingFunctional3D_LS<T,Descriptor,int> {
public:
    GetFlux(T LatticeDiffusivity_, T dt_, T dx_, T* fluxP_, T* surfaceP_)
        : LatticeDiffusivity(LatticeDiffusivity_), dt(dt_), dx(dx_), fluxP(fluxP_), surfaceP(surfaceP_)
    {
        totalFLux = 0.0;
        surface = 0.0;
     }
    GetFlux<T,Descriptor>* clone() const {
        return new GetFlux<T,Descriptor>(*this);
    }
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                        ScalarField3D<int>& flagMatrix) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX){
            for (plint iY=domain.y0; iY<=domain.y1; ++iY){
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ){
                    if (flagMatrix.get(iX,iY,iZ) == 2){
                        if (flagMatrix.get(iX-1,iY,iZ) == 1){
                            totalFLux += ((T) 2.0)*lattice.get(iX-1,iY,iZ).computeDensity();
                            surface += (T) 1.0;
                        }
                        if (flagMatrix.get(iX+1,iY,iZ) == 1){
                            totalFLux += ((T) 2.0)*lattice.get(iX+1,iY,iZ).computeDensity();
                            surface += (T) 1.0;
                        }
                        if (flagMatrix.get(iX,iY-1,iZ) == 1){
                            totalFLux += ((T) 2.0)*lattice.get(iX,iY-1,iZ).computeDensity();
                            surface += (T) 1.0;
                        }
                        if (flagMatrix.get(iX,iY+1,iZ) == 1){
                            totalFLux += ((T) 2.0)*lattice.get(iX,iY+1,iZ).computeDensity();
                            surface += (T) 1.0;
                        }
                        if (flagMatrix.get(iX,iY,iZ-1) == 1){
                            totalFLux += ((T) 2.0)*lattice.get(iX,iY,iZ-1).computeDensity();
                            surface += (T) 1.0;
                        }
                        if (flagMatrix.get(iX,iY,iZ+1) == 1){
                            totalFLux += ((T) 2.0)*lattice.get(iX,iY,iZ+1).computeDensity();
                            surface += (T) 1.0;
                        }
                    }
                }
            }
        }
        cout << "Surface is: " << surface << std::endl;
        cout << "The local lattice flux is: " << LatticeDiffusivity*totalFLux << std::endl;
        cout << "The local physical flux is: " << (LatticeDiffusivity*totalFLux)*(dx/dt) << std::endl;
        *fluxP = (LatticeDiffusivity*totalFLux)*(dx/dt);
        *surfaceP = surface;
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const{
        modified[0] = modif::dataStructure;
    }
private:
    T LatticeDiffusivity;
    T dt;
    T dx;
    T totalFLux;
    T surface;
    T* fluxP;
    T* surfaceP;
};

template<typename T, template<typename U> class Descriptor>
class initializeField : public BoxProcessingFunctional3D_LS<T,Descriptor,int> {
public:
    initializeField(T uMax_, plint nz_)
        : uMax(uMax_), nz(nz_)
        { }
    initializeField<T,Descriptor>* clone() const{
        return new initializeField<T,Descriptor>(*this);
    }
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                        ScalarField3D<int>& flagMatrix) {
        T rho = 1.;
        Array<T,3> tempVel (0., 0., 0.);
        Dot3D rel = lattice.getLocation();

        cout << "The proc no. " << global::mpi().getRank() << " is doing the processor" << std::endl;

        for (plint iX=domain.x0; iX<=domain.x1; ++iX){
            for (plint iY=domain.y0; iY<=domain.y1; ++iY){
                for (plint iZ=domain.z0; iZ <=domain.z1; ++iZ){
                    rho = (T)1;
                    if (flagMatrix.get(iX,iY,iZ) == 1){
                        tempVel[0] = T();
                        tempVel[1] = uMax*(((T)iZ+rel.z)/(T)nz);
                        tempVel[2] = T();
                    } else {
                        tempVel[0] = T();
                        tempVel[1] = T();
                        tempVel[2] = T();
                    }
                    iniCellAtEquilibrium(lattice.get(iX,iY,iZ), rho, tempVel);
                }
            }
        }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const{
        for (pluint iBlock=0; iBlock<modified.size(); ++iBlock) {
            modified[iBlock] = modif::staticVariables;
        }
    }
private:
    T uMax;
    plint nz;
};

template<typename T, template<typename U> class adDescriptor>
class initializeConcentration : public BoxProcessingFunctional3D_LS<T,adDescriptor,int> {
public:
    initializeConcentration(T cMax_, plint nz_)
        : cMax(cMax_), nz(nz_)
        { }
    initializeConcentration<T,adDescriptor>* clone() const{
        return new initializeConcentration<T,adDescriptor>(*this);
    }
    virtual void process(Box3D domain, BlockLattice3D<T,adDescriptor>& adLattice,
                                        ScalarField3D<int>& flagMatrix) {
        T concentration = 0.;
        Array<T,3> jEq (0., 0., 0.);
        Dot3D rel = adLattice.getLocation();

        for (plint iX=domain.x0; iX<=domain.x1; ++iX){
            for (plint iY=domain.y0; iY<=domain.y1; ++iY){
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ){
                    if (flagMatrix.get(iX,iY,iZ) == 1){
                        concentration = cMax*(((T)iZ+rel.z)/(T)nz);
                    } else if (flagMatrix.get(iX,iY,iZ) == 3) {
                        concentration = cMax;
                    } else if (flagMatrix.get(iX,iY,iZ) == 2) {
                        concentration = 0.001;
                    } else {
                        concentration = 0.0;
                    }
                    iniCellAtEquilibrium(adLattice.get(iX,iY,iZ), concentration, jEq);
                    if (iX == domain.x1/2){
                        cout << flagMatrix.get(iX,iY,iZ);
                    }
                }
                if (iX == domain.x1/2){
                    cout << std::endl;
                }
            }
        }
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const{
        for (pluint iBlock=0; iBlock<modified.size(); ++iBlock) {
            modified[iBlock] = modif::staticVariables;
        }
    }
private:
    T cMax;
    plint nz;
};

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              MultiBlockLattice3D<T,ADESCRIPTOR>& adLattice,
              T dx, T dt, plint iter)
{
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeDensity(adLattice), "concentration", (T)1.);
    vtkOut.writeData<float>(*computeDensity(lattice), "pressure", (dx*dx)/(dt*dt));
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
}

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);

    if (argc != 4) {
        pcout << "Error the parameters are wrong. The structure must be :\n";
        pcout << "1 : N\n2 : filename\n3 : Outdir\n";
        exit(1);
    }

    const plint maxIter = 200000;
    const plint N = atoi(argv[1]);
    const std::string meshFileName = argv[2];
    const std::string outDir = argv[3];
    global::directories().setOutputDir(outDir);

    triangleSet = new TriangleSet<T>(meshFileName, DBL);
    referenceDirection = 2;
    plint referenceResolution = 20;
    plint margin = 1;
    plint borderWidth = 1;

    pcout << "At least we're getting here." << std::endl;

    DEFscaledMesh<T>* defMesh =
        new DEFscaledMesh<T>(*triangleSet, N, referenceDirection, margin, extraLayer);
    TriangleBoundary3D<T> boundary(*defMesh);
    delete defMesh;

    pcout << "Getting here 2" << std::endl;

    boundary.getMesh().inflate();

    pcout << "Probably failing here" << std::endl;


    T uAveLB = 107.8167e-6; // 2.88e-6
    T uMax = 0.3743636;     // 0.01
    T cMax = 1.61514e-3;    // same
    T Re = 0.03743636;      // 4e-3


    T dx = boundary.getDx();
    T dt = 2.358e-9;
    T nuLB = (uAveLB*N) / Re;
    T Diffusivity = 4.58e-9;
    T LatticeDiffusivity = Diffusivity*dt/(dx*dx);
    T tau = (3.*nuLB+0.5);
    T omega = 1./tau;
    T adOmega = 1.0/(ADESCRIPTOR<T>::invCs2*LatticeDiffusivity + 0.5);

    Array<T,3>location(boundary.getPhysicalLocation());

    bool performOutput = true;

    if (performOutput) {
         pcout << "dx = " << dx << std::endl;
         pcout << "dt = " << dt << std::endl;
         pcout << "nuLB = " << nuLB << std::endl;
         pcout << "tau = " << tau << std::endl;
         pcout << "LatticeDiffusivity = " << LatticeDiffusivity << "\t adOmega = " << adOmega << std::endl;
         pcout << "adTau = " << ADESCRIPTOR<T>::invCs2*LatticeDiffusivity + 0.5 << std::endl;
         pcout << "invCs2 Diffusion: " << ADESCRIPTOR<T>::invCs2 << std::endl;
     }

     if(tau <= 0.5){
         pcout << "tau has invalid value." << std::endl;
         exit(1);
     }

     const int flowType = voxelFlag::inside;
     VoxelizedDomain3D<T> voxelizedDomain (
         boundary, flowType, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize);

     MultiScalarField3D<int> flagMatrix3D((MultiBlock3D&)voxelizedDomain.getVoxelMatrix());
     setToConstant(flagMatrix3D, voxelizedDomain.getVoxelMatrix(),
                     voxelFlag::inside, flagMatrix3D.getBoundingBox(), 1);
     setToConstant(flagMatrix3D, voxelizedDomain.getVoxelMatrix(),
                     voxelFlag::innerBorder, flagMatrix3D.getBoundingBox(), 1);
     setToConstant(flagMatrix3D, voxelizedDomain.getVoxelMatrix(),
                     voxelFlag::outerBorder, flagMatrix3D.getBoundingBox(), 2);

     MultiScalarField3D<int> Flags3D(
         flagMatrix3D.getBoundingBox().getNx(),
         flagMatrix3D.getBoundingBox().getNy(),
         flagMatrix3D.getBoundingBox().getNz());

     Box3D Domain3D = flagMatrix3D.getBoundingBox();

    plint nx = flagMatrix3D.getBoundingBox().getNx();
    plint ny = flagMatrix3D.getBoundingBox().getNy();
    plint nz = flagMatrix3D.getBoundingBox().getNz();

    Box3D topDomain(0,nx-1,0,ny-1,nz-1,nz-1);
    Box3D Inlet(0,nx-1,0,0,1,nz-2);
    //Box3D Outlet(0,nx-1,ny-1,ny-1,1,nz-2);
    Box3D Outlet(0,nx-1,ny-3,ny-1,1,nz-2);
    Box3D Left(0,0,0,ny-1,1,nz-2);
    //Box3D Right(nx-1,nx-1,0,ny-1,1,nz-2);
    Box3D Right(nx-3,nx-1,0,ny-1,1,nz-2);

    Box3D entireDomain(0,nx-1,0,ny-1,0,nz-1);

    pcout << "nx: " << nx << "\t ny: " << ny << "\t nz: " << nz << std::endl;

    int* rank = new int;
    MPI_Comm_rank(MPI_COMM_WORLD,rank);

    int count(nx*ny*nz);
    int* tempArray = new int[count];

    plint numCells = 0;

    ScalarField3D<int>* temp = new ScalarField3D<int> (flagMatrix3D.getBoundingBox().getNx(),
                            flagMatrix3D.getBoundingBox().getNy(),
                            flagMatrix3D.getBoundingBox().getNz());

    applyProcessingFunctional(new domainInitializer3D<int>(temp), Domain3D, flagMatrix3D);

    if (*rank == global::mpi().bossId()){
        for (plint x=0; x<nx; x++){
            for (plint y=0; y<ny; y++) {
                for (plint z=0; z<nz; z++){
                    tempArray[x*ny*nz+y*nz+z] = (int)temp->get(x,y,z);
                }
            }
        }
    }
    global::mpi().barrier();
    pcout << "We broadcast here" << std::endl;
    MPI_Bcast(tempArray, count, MPI_INT, global::mpi().bossId(), MPI_COMM_WORLD);

    if (global::mpi().getRank() != -1){
        for (plint x=0; x<nx; x++){
            for (plint y=0; y<ny; y++) {
                for (plint z=0; z<nz; z++){
                    temp->get(x,y,z) = (int)tempArray[x*ny*nz+y*nz+z];
                }
            }
        }
    }

    applyProcessingFunctional(new domainInitializer2D<int>(temp), Domain3D, Flags3D);

    setToConstant(Flags3D, topDomain, 3);
    setToConstant(Flags3D, Inlet, 1);
    setToConstant(Flags3D, Outlet, 1);
    setToConstant(Flags3D, Left, 1);
    setToConstant(Flags3D, Right, 1);

    delete temp;
    delete tempArray;

    std::unique_ptr<MultiBlockLattice3D<T,DESCRIPTOR>> lattice
        = generateMultiBlockLattice(Domain3D, new DYNAMICS, envelopeWidth);
    std::unique_ptr<MultiBlockLattice3D<T,ADESCRIPTOR>> adLattice
        = generateMultiBlockLattice(Domain3D, new ADYNAMICS, envelopeWidth);
    lattice->toggleInternalStatistics(false);
    adLattice->toggleInternalStatistics(false);

    // Use periodic boundary conditions.
    lattice->periodicity().toggle(0,true);
    adLattice->periodicity().toggle(0,true);
    lattice->periodicity().toggle(1,true);
    adLattice->periodicity().toggle(1,true);

    defineDynamics(*lattice, Flags3D, entireDomain, new BounceBack<T,DESCRIPTOR>((T) 1.), 2);
    defineDynamics(*lattice, topDomain, new VelocityBounceBack<T,DESCRIPTOR>((T) 1., Array<T,3>((T) 0., uMax, (T) 0.)));
    defineDynamics(*adLattice, Flags3D, entireDomain, new AntiBounceBack<T,ADESCRIPTOR>((T) 0., Array<T,3>((T) 0., (T) 0., (T) 0.)), 2);
    defineDynamics(*adLattice, topDomain, new AntiBounceBack<T,ADESCRIPTOR>(cMax, Array<T,3>((T) 0., (T) 0., (T) 0.)));

    applyProcessingFunctional(new initializeField<T,DESCRIPTOR>(uMax,nz), entireDomain, *lattice, Flags3D);
    applyProcessingFunctional(new initializeConcentration<T,ADESCRIPTOR>(cMax,nz), entireDomain, *adLattice, Flags3D);

    plint processorLevel = 1;
    integrateProcessingFunctional (
            new LatticeToPassiveAdvDiff3D<T,DESCRIPTOR,ADESCRIPTOR> (),
            lattice->getBoundingBox(),
            *lattice, *adLattice, processorLevel );

    //RectangleSetup(lattice, parameters,initializeAtEquilibrium(*lattice, lattice->getBoundingBox(), initializeField<T>(lattice->getBoundingBox().getNy(), uMax, flagMatrix2D)); *boundaryCondition);

    plint numCores = global::mpi().getSize();
    pcout << "Number of MPI threads: " << numCores << std::endl;

    global::timer("benchmark").start();
    //global::profiler().turnOn();

    // Loop over main time iteration.
    util::ValueTracer<T> converge(uMax,lattice->getBoundingBox().getNy(),1.0e-3);
    writeVTK(*lattice, *adLattice, dx, dt, (plint) 0);
    for (plint iT=0; iT<maxIter; ++iT) {
        converge.takeValue(getStoredAverageEnergy(*lattice),true);

        if (converge.hasConverged())
        {
            pcout << "Simulation converged after " << iT << " iterations: \n"
                  << (T) (numCells*iT) /
                     global::timer("benchmark").getTime() / 1.e6
                  << " Mega site updates per second.\n"
                  << "The simulation took " << global::timer("benchmark").getTime() << "seconds" << std::endl << std::endl;
            pcout << "Saving VTK file ..." << endl;
            writeVTK(*lattice, *adLattice, dx, dt, iT);
            break;
        }

        if (iT%10==0)
        {
            pcout << "The simulation took " << iT << "steps." << std::endl;
            writeVTK(*lattice, *adLattice, dx, dt, iT);
        }

        // Execute a time iteration.
        adLattice->collideAndStream();
        lattice->collideAndStream();

    }

    T* totalFLux = new T;
    T* totalSurface = new T;

    T* globalFlux = new T;
    T* globalSurface = new T;

    applyProcessingFunctional(new GetFlux<T,ADESCRIPTOR>(LatticeDiffusivity, dt, dx, totalFLux, totalSurface), adLattice->getBoundingBox(), *adLattice, Flags3D);

    MPI_Reduce(totalFLux, globalFlux, global::mpi().getSize(), MPI_DOUBLE, MPI_SUM, global::mpi().bossId(), MPI_COMM_WORLD);
    MPI_Reduce(totalSurface, globalSurface, global::mpi().getSize(), MPI_DOUBLE, MPI_SUM, global::mpi().bossId(), MPI_COMM_WORLD);

    pcout << "The global physical flux is: " << *globalFlux << std::endl;
    pcout << "The global surface is: " << *globalSurface << std::endl;
    pcout << "The global average physical flux is " << (*globalFlux)/(*globalSurface) << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    global::profiler().writeReport();
}
