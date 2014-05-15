// Modified Computational Fabrication Assignment #1 By David Levin 2014

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cassert>
//#include <vector>
//#include "ofMain.h"
#include "ofApp.h"
//#include "CompFab.h"
//#include "Mesh.h"
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

//The Voronoi Diagram
std::vector<CompFab::Vec3> diagram;

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

// Kai's mesh rendering stuff
Mesh *rMesh;


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
    
    //rMesh = new Mesh(filename, true);
    rMesh = new Mesh(filename, false);
    
    CompFab::Vec3 v1,v2,v3;
    
    unsigned int currI, nextI;
    TriangleList g_objTriangles;
    for(unsigned int objTriI = 0; objTriI <rMesh->objTriangleMap.size()-1; objTriI++){
        currI = rMesh->objTriangleMap[objTriI];
        nextI = rMesh->objTriangleMap[objTriI+1];
        for(unsigned int tri = currI; tri<nextI; ++tri)
        {
            v1 = rMesh->v[rMesh->t[tri][0]];
            v2 = rMesh->v[rMesh->t[tri][1]];
            v3 = rMesh->v[rMesh->t[tri][2]];
            g_objTriangles.push_back(CompFab::Triangle(v1,v2,v3));
        }
        g_buildingTriangles.push_back(g_objTriangles);
        g_objTriangles.clear();
    }
    
    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*rMesh, bbMin, bbMax);
    
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

    g_voxelGrid->initializeLabels();

    return true;
    
}


/*
  Convert the voxel representation of the grid into a mesh for saving to file or render.
*/
void triangulateVoxelGrid(const char * outfile, unsigned int label)
{
    cout << "Trianglulating\n";

    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;

    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    int num_label1 = 0;
    int num_label2 = 0;
    int num_label_both = 0;

    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {

                std::set<unsigned int>::iterator it = g_voxelGrid->getLabels(ii,jj,kk).find(label);
                std::set<unsigned int>::iterator it1 = g_voxelGrid->getLabels(ii,jj,kk).find(1);
                std::set<unsigned int>::iterator it2 = g_voxelGrid->getLabels(ii,jj,kk).find(2);
                bool is_in1 = (it1 != g_voxelGrid->getLabels(ii,jj,kk).end());
                bool is_in2 = (it2 != g_voxelGrid->getLabels(ii,jj,kk).end());
                if(is_in1 && is_in2){
                    num_label_both ++;
                }else if(is_in1){
                    num_label1 ++;
                }else if(is_in2){
                    num_label2 ++;
                }
                bool is_in = (it != g_voxelGrid->getLabels(ii,jj,kk).end());
                if(g_voxelGrid->getLabels(ii,jj,kk).empty() || !is_in){
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
    cout << "FINAL COUNTS: label1:" << num_label1 << " label2:" << num_label2 << " both:" << num_label_both << "\n";
    // Compute the normals
    mout.compute_norm();
    mout.save_obj(outfile);
}

/*
  Create the grid, update voxel labels for input building meshes
  @sets g_voxelgrid
  TODO: do this for multiple objects
*/
void voxelizer(char* filename, unsigned int voxelres, unsigned int label)
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
                    //std::cout << "init label..." << bIntersect[0] << "\n";
                    g_voxelGrid->addLabel(ii,jj,kk, bIntersect[0]);
                }
            }
        }
    }

    // For testing: converts grid representation to a mesh file with voxelized buildings (file will show blocks where getLabel >= 1)
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

                std::set<unsigned int>::iterator it1 = g_voxelGrid->getLabels(ii,jj,kk).find(1);
                bool is_in1 = (it1 != g_voxelGrid->getLabels(ii,jj,kk).end());
                if(!g_voxelGrid->getLabels(ii,jj,kk).empty()){
                   voxel = std::make_pair(CompFab::Vec3(ii, jj, kk),CompFab::Vec3(ii, jj, kk));
                   g_voxelGrid->wavefront.push(voxel);
                }
            }
        }
    }
    std::cout << "WAVEFRONT INIT SIZE: " << g_voxelGrid->wavefront.size() << "\n";


}

bool areCoordinatesValid(double x, double y, double z)
{
    return true;
}

//void printWavefront()
//{
//   std::priority_queue< std::pair<Vec3, Vec3>, std::vector<std::pair<Vec3,Vec3> >, CompareVoxel> copy;
//   copy = g_voxelGrid->wavefront;
//   print "PQ------------------------------------------\n";
//   while(!copy.empty()){
//       std::pair <CompFab::Vec3,CompFab::Vec3> voxelPair = g_voxelGrid->wavefront.top();
//       CompFab::Vec3 voxel = voxelPair.first;
//       CompFab::Vec3 parent = voxelPair.second;
//       std::cout << ", (voxel: " << ii << ";" << jj << ";" << kk << " label:" distance: " << CompFab:distance(voxel, parent); << "\n";
//       g_voxelGrid->wavefront.pop();
//   }
//
//}

///////////////
// PROPAGATE //
///////////////
void propagate(){
   //std::list<CompFab::Vec3> diagram;
    
   while(!g_voxelGrid->wavefront.empty()){
       //printWavefront();
       std::pair <CompFab::Vec3,CompFab::Vec3> voxelPair = g_voxelGrid->wavefront.top();
       CompFab::Vec3 voxel = voxelPair.first;
       CompFab::Vec3 parent = voxelPair.second;
       g_voxelGrid->wavefront.pop();
       
       std::set<unsigned int> curr_labels = g_voxelGrid->getLabels(voxel.m_x, voxel.m_y, voxel.m_z);
       //std::cout << "curr_labels size: " << curr_labels.size() << "\n";
      
       // For boundary case
       if(curr_labels.size() > 1){
            continue; 
       }
       assert(curr_labels.size() == 1);
       int curr_label;
       if (curr_labels.size() == 1){
        curr_label = *curr_labels.begin();
       }
       assert(curr_label != 0);
       ///std::cout << "Curr Voxel: " << voxel.m_x << " " << voxel.m_y << " " << voxel.m_z << " has label " <<  curr_label << "\n";
       for(int i = -1; i <=1; i++)
       {
           for(int j = -1; j <=1; j++)
           {
               for(int k = -1; k <=1; k++)
               {
                   int n_x = (unsigned int)voxel.m_x + i;
                   int n_y = (unsigned int)voxel.m_y + j;
                   int n_z = (unsigned int)voxel.m_z + k;

                   //std::cout << "Voxel neighbor: " << n_x << " " << n_y << " " << n_z << "\n";

                   //If same square as before don't propegate
                   if(!(i == 0 and j == 0 and k == 0))
                   {
                       //Check if coordinates are valid
                       if( (n_x >= 0 and n_x < g_voxelGrid->m_dimX) and (n_y >= 0 and n_y < g_voxelGrid->m_dimY) and (n_z >= 0 and n_z < g_voxelGrid->m_dimZ))
                       {
                           //int n_label = g_voxelGrid->getLabels(n_x, n_y, n_z);
                           //std::cout << "Voxel neighbor original label: " << n_label << "\n";
                           //Check if unlabelled
                           if(g_voxelGrid->getLabels(n_x, n_y, n_z).empty())
                           {
                               //std::cout << "Valid: ";
                               //Assign it the new label
                               g_voxelGrid->addLabel(n_x,n_y,n_z,curr_label);
                               //std::cout << "setting label to " << curr_label << "\n";
                               std::pair <CompFab::Vec3,CompFab::Vec3> newVoxel;
                               newVoxel = std::make_pair(CompFab::Vec3(n_x, n_y, n_z), parent);
                               g_voxelGrid->wavefront.push(newVoxel);
                           } else {
                               if(curr_label != *g_voxelGrid->getLabels(n_x,n_y,n_z).begin())
                               {
                                   diagram.push_back(CompFab::Vec3(n_x, n_y, n_z));
                                   g_voxelGrid->addLabel(n_x,n_y,n_z,curr_label);
                                   //cout << "label count: " << g_voxelGrid->getLabels(n_x,n_y,n_z).size() << "\n";
                               }
                           }
                       } else {
                         //std::cout << "Invalid" << "\n";
                       }
                   }
               }
           }
       }
       
       
   }
   std::cout << "\nfinished propegation\n";
}

//find closest voxel to each point
//TODO check if path goes through building, if so search neighbors to find better path
//do a search of neighbors (check all neighbors, check if on boundary ie have 2 labels, then do BFS
//returns a vector of the voxels on the boundary
vector<CompFab::Vec3> findPath(CompFab::Vec3 a, CompFab::Vec3 b)
{
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    
    //Find closest point to a (start) and b (end)
    //TODO update to check if path intersects a building?
    CompFab::Vec3 start;
    CompFab::Vec3 end;
    double minDistStart = -1.0;
    double minDistEnd = -1.0;
    for(int i = 0; i < diagram.size(); i++) {
        if(minDistStart < 0 || CompFab::distance(diagram[i], a) < minDistStart)
        {
            start = diagram[i];
            minDistStart = CompFab::distance(diagram[i], a);
        }
        if(minDistEnd < 0 || CompFab::distance(diagram[i], b) < minDistEnd)
        {
            end = diagram[i];
            minDistEnd = CompFab::distance(diagram[i], b);
        }
    }
    
    cout << "start: " << start.m_x << " " << start.m_y << " " << start.m_z << "\n";
    cout << "start label count: " << g_voxelGrid->getLabels(start.m_x,start.m_y,start.m_z).size() << "\n";
    cout << "end: " << end.m_x << " " << end.m_y << " " << end.m_z << "\n";
    cout << "end label count: " << g_voxelGrid->getLabels(end.m_x,end.m_y,end.m_z).size() << "\n";
    
    //Find path from start to end
    vector<CompFab::Vec3> path;
    queue<CompFab::Vec3> queue;
    std::map<int, CompFab::Vec3> visited;
    
    queue.push(start);
    visited[(start.m_x*(nx*ny) + start.m_y*(nz) + start.m_z)] = start;
    while(!queue.empty())
    {
        CompFab::Vec3 current = queue.front();
        queue.pop();
        //cout << "checking node: " << current.m_x << " " << current.m_y << " " << current.m_z << "\n";
        
        //found end point
        if((unsigned int)end.m_x == current.m_x && (unsigned int)end.m_y == current.m_y && (unsigned int)end.m_z == current.m_z)
        {
            //cout << "found end\n";
            //loop through parents and add to path
            CompFab::Vec3 node = end;
            while (node.m_x != start.m_x || node.m_y != start.m_y || node.m_z != start.m_z)
            {
                //cout << "adding node to path: " << node.m_x << " " << node.m_y << " " << node.m_z << "\n";
                path.insert(path.begin(), node);
                node = visited.at(node.m_x*(nx*ny) + node.m_y*(nz) + node.m_z); //gets the parent of node
            }
            return path;
        }
        
        //check all neighbors
        for(int i = -1; i <=1; i++)
        {
            for(int j = -1; j <=1; j++)
            {
                for(int k = -1; k <=1; k++)
                {
                    int x = (unsigned int)current.m_x + i;
                    int y = (unsigned int)current.m_y + j;
                    int z = (unsigned int)current.m_z + k;
                    
                    //TODO replace visited with set
                    map<int,CompFab::Vec3>::iterator it = visited.find(x*(nx*ny) + y*(nz) + z);
                    //If not already in visited and not same square as before and on the grid
                    //cout << "nz: " << nz << "z: " << z << "\n";
                    if(it == visited.end() && !(i == 0 && j == 0 && k == 0) && ((x >= 0 && x < nx) and (y >= 0 and y < ny) and (z >= 0 and z < nz)))
                    {
                        //If on boundary, ie has >1 label
                        if(g_voxelGrid->getLabels(x,y,z).size() > 1)
                        {
                            queue.push(CompFab::Vec3(x,y,z));
                            //visited.insert(pair<CompFab::Vec3, CompFab::Vec3>(CompFab::Vec3(x,y,z), current));
                            //visited.insert(std::map<CompFab::Vec3, CompFab::Vec3>::value_type(CompFab::Vec3(x,y,z), current));
                            visited[(x*(nx*ny) + y*(nz) + z)] = current;
                        }
                    }
                }
            }
        }
    }
    
    return path;
}


///////////////
// RENDERING //
///////////////


//////////
// MAIN //
//////////

int main(int argc, char **argv)
{
    std::clock_t prerend;
    prerend = std::clock();
    
    if(argc < 3)
    {
        std::cout<<"Usage: InputMeshFilename voxelRes \n";
        exit(0);
    }
    std::cout<<"Load Mesh file: "<<argv[1]<<"\n";
    std::cout<<"Grid resolution : "<<argv[2]<<"\n";

    int debugLabel = 1;
    if(argc == 4) {
        std::cout<<"Building label : "<<argv[3]<<"\n";
        debugLabel = atoi(argv[3]);
    }

    voxelRes = atoi(argv[2]);

    voxelizer(argv[1], voxelRes, debugLabel);
    initializeWF();
    
    std::clock_t preprop;
    preprop = std::clock();
    propagate();
    
    std::clock_t prepath;
    prepath = std::clock();
    
    cout << "starting path finding" << endl;
    std::vector< CompFab::Vec3 > path = findPath(CompFab::Vec3(0,0,0), CompFab::Vec3(voxelRes,voxelRes/2,(int)(voxelRes*.6)));
    cout << "finished path finding" << endl;
    
    
    
    cout << "render time = " << (preprop - prerend) / (double) CLOCKS_PER_SEC << endl;
    cout << "propagation time = " << (prepath - preprop) / (double) CLOCKS_PER_SEC << endl;
    cout << "path calculation time = " << (std::clock() - prepath) / (double) CLOCKS_PER_SEC << endl;
    cout << "total time = " << (std::clock() - prerend) / (double) CLOCKS_PER_SEC << endl;
    
	ofSetupOpenGL(1024,768, OF_WINDOW);			// <-------- setup the GL context
	// this kicks off the running of my app
	// can be OF_WINDOW or OF_FULLSCREEN
	// pass in width and height too:
	ofRunApp( new ofApp(rMesh, true, path));
    
    
}

