//
//  CompFab.h
//  voxelizer
//
//  Created by David Levin on 2/3/14.
//
//

#ifndef voxelizer_CompFab_h
#define voxelizer_CompFab_h

#define EPSILON 1e-9

#include <cmath>
#include <queue>
#include <set>
#include <vector>
#include <utility>
#include <iostream>

namespace CompFab
{
    
    //Data Types
    typedef struct Vec3Struct
    {
        
        Vec3Struct();
        Vec3Struct(double x, double y, double z);

        union
        {
            double m_pos[3];
            struct { double m_x,m_y,m_z; };
        };
        
        inline double & operator[](unsigned int index) { return m_pos[index]; }
        inline const double & operator[](unsigned int index) const { return m_pos[index]; }
        inline void operator+=(const Vec3Struct &a)
        {
            m_x += a.m_x;
            m_y += a.m_y;
            m_z += a.m_z;
        }
        
        void normalize();
        
    }Vec3;

    //Data Types
    typedef struct Vec3iStruct
    {
        
        Vec3iStruct();
        Vec3iStruct(double x, double y, double z);
        union
        {
            int m_pos[3];
            struct {int m_x,m_y,m_z;};
        };
        
        inline int & operator[](unsigned int index) { return m_pos[index]; }
        inline const int & operator[](unsigned int index) const { return m_pos[index]; }
        
    }Vec3i;

    //Data Types
    typedef struct Vec2fStruct
    {
        
        Vec2fStruct();
        Vec2fStruct(double x, double y);
        
        union
        {
            float m_pos[2];
            struct { float m_x,m_y; };
        };
        
        inline float & operator[](unsigned int index) { return m_pos[index]; }
        inline const float & operator[](unsigned int index) const { return m_pos[index]; }
        
    }Vec2f;

    
    //NOTE: Ray direction must be normalized
    typedef struct RayStruct
    {
        
        RayStruct();
        RayStruct(Vec3 &origin, Vec3 &direction);
        
        Vec3 m_origin;
        Vec3 m_direction;
        
    } Ray;
    
    typedef struct TriangleStruct
    {
        
        TriangleStruct(Vec3 &v1, Vec3 &v2,Vec3 &v3);
        
        Vec3 m_v1, m_v2, m_v3;
        
    }Triangle;
    
    //Some useful operations
    //Compute v1 - v2
    Vec3 operator-(const Vec3 &v1, const Vec3 &v2);
    
    Vec3 operator+(const Vec3 &v1, const Vec3 &v2);
    
    //Cross Product
    Vec3 operator%(const Vec3 &v1, const Vec3 &v2);
    
    //Dot Product
    double operator*(const Vec3 &v1, const Vec3 &v2);
     
    
    double distance(Vec3 &v1, Vec3 &v2);
    class CompareVoxel
    {
    public:
        bool operator()(std::pair<Vec3,Vec3> n1,std::pair<Vec3,Vec3> n2)
        {
          double distn1 = distance(n1.first, n1.second);

          double distn2 = distance(n2.first, n2.second); 
    
          return (distn1>=distn2);
        }
    }; 
    //Grid structure for Voxels
    typedef struct VoxelGridStruct
    {
        //Square voxels only
        VoxelGridStruct(Vec3 lowerLeft, unsigned int dimX, unsigned int dimY, unsigned int dimZ, double spacing);
        ~VoxelGridStruct();

        std::vector< std::set<unsigned int> > voxelLabels;

        inline std::set<unsigned int> & getLabels(unsigned int i, unsigned int j, unsigned int k)
        {
            
            return voxelLabels[k*(m_dimX*m_dimY)+j*m_dimY + i];
        }

        inline void addLabel(unsigned int i, unsigned int j, unsigned int k, unsigned int label)
        {
            //std::cout << "voxelLabels.size " << voxelLabels.size() << "\n";
            //std::cout << "trying to set" << k*(m_dimX*m_dimY)+j*m_dimY + i << "\n";
            voxelLabels[k*(m_dimX*m_dimY)+j*m_dimY + i].insert(label);
        }
        // Initialize the voxelLabels
        inline void initializeLabels(){
            std::set<unsigned int> empty; 
            for(unsigned int ii=0; ii<m_size; ++ii)
            {
                voxelLabels.push_back(empty);
            }
        }
        
        //unsigned int *m_labelArray;
        std::priority_queue< std::pair<Vec3, Vec3>, std::vector<std::pair<Vec3,Vec3> >, CompareVoxel> wavefront;
        unsigned int m_dimX, m_dimY, m_dimZ, m_size;
        double m_spacing;
        Vec3 m_lowerLeft;
        
    } VoxelGrid;
}



#endif
