#pragma once

#include "ofMain.h"
#include "ofxGui.h"
#include "Mesh.h"

class ofApp : public ofBaseApp{
public:
    ofApp(Mesh* m, bool useVoxel, std::vector< CompFab::Vec3 > path);
    void setupMesh();
    void setupVoxel();
    void setupBoundary();
    void setupPath();
    void drawMesh();
    void drawVoxel();
    void drawBoundary();
    void drawPath();
    
    void setup();
    void update();
    void draw();
    
    void keyPressed(int key);
    void keyReleased(int key);
    void mouseMoved(int x, int y);
    void mouseDragged(int x, int y, int button);
    void mousePressed(int x, int y, int button);
    void mouseReleased(int x, int y, int button);
    void windowResized(int w, int h);
    void dragEvent(ofDragInfo dragInfo);
    void gotMessage(ofMessage msg);
    
    ofEasyCam cam;
    Mesh *reMesh;
    ofMesh rendMesh;
    bool showPath;
    std::vector< CompFab::Vec3 > pathCells;
    std::vector< std::vector< ofPoint > > voxObjList;// a vector of voxel objects. Each object is a vector of cell indices
    std::vector< ofPoint > boundaryList;
    std::vector< ofPoint > pathList;
    
    ofPoint lowerleft;
    ofPoint ptStart;
    ofPoint ptDest;
    ofPoint ptVehicle;
    double space;
    
    ofxPanel guiPanel;
    
};
