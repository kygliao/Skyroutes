// Modified Computational Fabrication Assignment #1 By David Levin 2014

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "ofMain.h"
#include "ofApp.h"
#include "CompFab.h"
#include "Mesh.h"
#include "vecmath/vecmath.h"
using namespace std;
using std::stringstream;
using std::cout;
using std::endl;
using std::ends;

/////////////
// GLOBALS //
/////////////

// The grid representation 
CompFab::VoxelGrid *g_voxelGrid;

// Triangles used for initializing labels;
typedef std::vector<CompFab::Triangle> TriangleList;
//TriangleList g_voxelTriangles;
std::vector<TriangleList> g_buildingTriangles;
unsigned int voxelRes;

double detMat(double A11, double A12, double A13,
              double A21, double A22, double A23,
              double A31, double A32, double A33)
{
    return A11*(A22*A33 - A23*A32) - A12*(A21*A33 - A23*A31) + A13*(A21*A32 - A22*A31);
}


////////////////////
// INITIALIZATION //
////////////////////
//TODO:Initialize labels for multiple objects, everything else!

/*
  Ray-Triangle Intersection
  @returns distance if triangle and ray intersect, 0 otherwise
*/
float rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle)
{
    CompFab::Vec3 e1(triangle.m_v2 - triangle.m_v1);
    CompFab::Vec3 e2(triangle.m_v3 - triangle.m_v1);
    CompFab::Vec3 q = ray.m_direction%e2;
    double a = e1*q;
    if(a > -0.000001 and a < 0.000001){
      return 0;
    }

    double f =  1.0/a;
    CompFab::Vec3 s = ray.m_origin - triangle.m_v1;
    double u = f*(s*q); 
    
    if(u < 0.0){
      return 0;
    }
    
    CompFab::Vec3 r = s%e1;
    double v = f*(ray.m_direction*r); 
    if((v < 0.0) or (u + v > 1.0)){
      return 0;
    }
    double t = f*(e2*r);
    if( t <= 0.0 ){
      return 0;
    }
    return (float)t;

}

/*
  Num Intersections with Ray
  @returns number of intersections with surface made by a ray originating at voxel and cast in direction.
*/
void testBuildingIntersect(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir, int *result)
{
    
    unsigned int numHits = 0;
    float dist;
    TriangleList g_triangles;

    CompFab::RayStruct vRay = CompFab::RayStruct(voxelPos, dir);
    for(unsigned int n = 0; n < g_buildingTriangles.size(); n++){
        g_triangles = g_buildingTriangles[n];
        for(unsigned int i = 0; i < g_triangles.size(); i++){
            CompFab::Triangle triangle = g_triangles[i];
            dist = rayTriangleIntersection(vRay, triangle);
            if(dist){
                numHits ++;
            }
        }
        // Intersection with building, return building num
        if(numHits % 2 != 0){
            result[0] = n+1;
            // Check if voxel is on boundary - wavefront
            return;
        }
        numHits = 0;
    }
    // No intersections
    result[0] = 0;
    return;
}

/*
  Load input building mesh and construct the voxel grid
  @set g_voxelGrid
*/
bool loadMesh(char *filename, unsigned int dim)
{
    g_buildingTriangles.clear();
    
    Mesh *tempMesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1,v2,v3;

    unsigned int currI, nextI;
    TriangleList g_objTriangles;
    for(unsigned int objTriI = 0; objTriI <tempMesh->objTriangleMap.size()-1; objTriI++){
        currI = tempMesh->objTriangleMap[objTriI];
        nextI = tempMesh->objTriangleMap[objTriI+1];
        for(unsigned int tri = currI; tri<nextI; ++tri)
        {
            v1 = tempMesh->v[tempMesh->t[tri][0]];
            v2 = tempMesh->v[tempMesh->t[tri][1]];
            v3 = tempMesh->v[tempMesh->t[tri][2]];
            g_objTriangles.push_back(CompFab::Triangle(v1,v2,v3));
        }
        g_buildingTriangles.push_back(g_objTriangles);
        g_objTriangles.clear();
    }

    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*tempMesh, bbMin, bbMax);
    
    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;
    
    if(bbX > bbY && bbX > bbZ)
    {
        spacing = bbX/(double)(dim-2);
    } else if(bbY > bbX && bbY > bbZ) {
        spacing = bbY/(double)(dim-2);
    } else {
        spacing = bbZ/(double)(dim-2);
    }
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    g_voxelGrid = new CompFab::VoxelGrid(bbMin-hspacing, dim, dim, dim, spacing);

    delete tempMesh;
    return true;
   
}

/*
  Convert the voxel representation of the grid into a mesh for saving to file or render.
*/
void triangulateVoxelGrid(const char * outfile, int label)
{
    cout << "Trianglulating\n";

    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;

    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);

    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!g_voxelGrid->getLabel(ii,jj,kk) || 
                    (g_voxelGrid->getLabel(ii,jj,kk) != label)){
                  continue;
                }
                CompFab::Vec3 coord(((double)ii)*spacing, ((double)jj)*spacing, ((double)kk)*spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }
    // Compute the normals
    mout.compute_norm();
    mout.save_obj(outfile);
}

/*
  Create the grid, update voxel labels for input building meshes
  @sets g_voxelgrid
  TODO: do this for multiple objects
*/
void voxelizer(char* filename, char* outfilename, unsigned int voxelres, unsigned int label) 
{

    unsigned int dim = voxelres; //dimension of voxel grid (e.g. 32x32x32)

    // Construct the voxel grid, populate g_buildingTriangles - used in ray tracing for labeling in rest of function
    loadMesh(filename, dim);

    if(label < 1 || label > g_buildingTriangles.size()) {
        std::cout<< "debug label invalid, number of buildings is :" << g_buildingTriangles.size() << "\n";
        exit(0);
    }

    // Assign voxel labels for buildings 
    //TODO: Only works if no triangles at 10,0,0, needs debug
    CompFab::Vec3 direction(10.0,0.0,0.0);

    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;
    CompFab::Vec3 left = g_voxelGrid->m_lowerLeft;
    int bIntersect[2];
    cout << "m_lowerleft" << left.m_x << "," << left.m_y << "," << left.m_z;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    // Iterate over all voxels in g_voxelGrid and test whether they are inside our outside 
    // of the  surface defined by the each input building mesh
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                CompFab::Vec3 vPos(left.m_x + ((double)ii)*spacing, left.m_y + ((double)jj)*spacing, left.m_z +((double)kk)*spacing);
                testBuildingIntersect(vPos, direction, bIntersect);
                if(bIntersect[0]){
                    g_voxelGrid->getLabel(ii,jj,kk) = bIntersect[0];
                }
            }
        }
    }

    // For testing: converts grid representation to a mesh file with voxelized buildings (file will show blocks where getLabel >= 1)
    triangulateVoxelGrid(outfilename, label);
    cout << "Done \n";
}
void initializeWF(){
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;

    std::pair <CompFab::Vec3,CompFab::Vec3> voxel;
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(g_voxelGrid->getLabel(ii,jj,kk)){
                   voxel = std::make_pair(CompFab::Vec3(ii, jj, kk),CompFab::Vec3(ii, jj, kk));
                   g_voxelGrid->wavefront.push(voxel);
                }
            }
        }
    }


}
///////////////
// PROPAGATE //
///////////////
void propagate(){
   //while(!wavefront.empty()){




   //}
}

///////////////
// RENDERING //
///////////////


//////////
// MAIN //
//////////

int main(int argc, char **argv)
{

    if(argc < 4)
    {
        std::cout<<"Usage: InputMeshFilename OutputMeshFilename voxelRes building\n";
        exit(0);
    }
    std::cout<<"Load Mesh file: "<<argv[1]<<"\n";
    std::cout<<"Output Mesh file: "<<argv[2]<<"\n";
    std::cout<<"Grid resolution : "<<argv[3]<<"\n";

    int debugLabel = 1;
    if(argc == 5) {
        std::cout<<"Building label : "<<argv[4]<<"\n";
        debugLabel = atoi(argv[4]);
    }

    voxelRes = atoi(argv[3]);

    // Create the grid, set "getLabel" for initial inputfile 
    // TODO: do this for multiple files
    voxelizer(argv[1], argv[2], voxelRes, debugLabel);
    initializeWF();
    propagate();

	//ofSetupOpenGL(1024,768, OF_WINDOW);			// <-------- setup the GL context
	// this kicks off the running of my app
	// can be OF_WINDOW or OF_FULLSCREEN
	// pass in width and height too:
	//ofRunApp( new ofApp());

}

