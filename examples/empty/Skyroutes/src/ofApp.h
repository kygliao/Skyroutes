#pragma once

#include "ofMain.h"
#include "Mesh.h"

class ofApp : public ofBaseApp{
public:
    ofApp(Mesh* m, bool useVoxel);
    void setupMesh();
    void setupVoxel();
    void drawMesh();
    void drawVoxel();
    
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
    bool useVoxel;
    std::vector< std::vector< ofPoint > > voxObjList;
    // a vector of voxel objects. Each object is a vector of cell indices
};
