#include "ofApp.h"
#include <cstdlib>
#include <ctime>


//--------------------------------------------------------------
ofApp::ofApp(Mesh* m, bool vox, vector<CompFab::Vec3> path){
    reMesh = m;
    showPath = vox;
    pathCells = path;
    
    
    ofBaseApp();
}


//--------------------------------------------------------------
void ofApp::setupMesh(){
    
    for (int i=0; i<reMesh->v.size(); i++){
        ofVec3f ofv = ofVec3f(reMesh->v[i][0], reMesh->v[i][1], reMesh->v[i][2]);
        rendMesh.addVertex(ofv);
    }
    
    for (int k=0; k<reMesh->t.size(); k++){
        rendMesh.addIndex(reMesh->t[k][0]);
        rendMesh.addIndex(reMesh->t[k][1]);
        rendMesh.addIndex(reMesh->t[k][2]);
        
    }
    
}

void ofApp::setupVoxel(){
    for (int i = 0; i < reMesh->numObj; i++){
        voxObjList.push_back(vector< ofPoint >());
    }
    
    
    for (int ii = 0; ii < g_voxelGrid->m_dimX; ii++){
        for (int jj = 0; jj < g_voxelGrid->m_dimY; jj++){
            for (int kk = 0; kk < g_voxelGrid->m_dimZ; kk++){
                std::set< unsigned int > label = g_voxelGrid->getLabels(ii,jj,kk);
                if ((int) label.size() == 1){
                    ofPoint pt = ofPoint(ii*space, jj*space, kk*space) + lowerleft;
                    std::set< unsigned int >::iterator it = label.begin();
                    unsigned int index = *it;
                    voxObjList[index].push_back(pt);
                }
                
            }
        }
    }
    
}

void ofApp::setupBoundary(){
    
    
    for (int ii = 0; ii < g_voxelGrid->m_dimX; ii++){
        for (int jj = 0; jj < g_voxelGrid->m_dimY; jj++){
            for (int kk = 0; kk < g_voxelGrid->m_dimZ; kk++){
                std::set< unsigned int > label = g_voxelGrid->getLabels(ii,jj,kk);
                if (label.size() > 1){
                    ofPoint pt = ofPoint(ii*space, jj*space, kk*space) + lowerleft;
                    boundaryList.push_back(pt);
                }
                
            }
        }
    }
    
}

void ofApp::setupPath(){
    for (int i = 0; i < pathCells.size(); i++){
        ofPoint pt(pathCells[i][0]*space, pathCells[i][1]*space, pathCells[i][2]*space);
        pt += lowerleft;
        pathList.push_back(pt);
    }
}
//--------------------------------------------------------------
void ofApp::drawMesh(){
    rendMesh.drawFaces();
    rendMesh.disableColors();
    ofSetColor(120,120,120);
    rendMesh.drawWireframe();
    rendMesh.enableColors();
}

void ofApp::drawVoxel(){
    
    int color = 256 * 3 / reMesh->numObj;
    int div = reMesh->numObj / 3;
    
    for (int i = 1; i < reMesh->numObj; i++){ // look at each voxelized object
        int r = (i/div == 0) ? ((i % div) + 1)*color : 0;
        int g = (i/div == 1) ? ((i % div) + 1)*color : 0;
        int b = (i/div > 1) ? ((i % div) + 1)*color : 0;
        ofSetColor(r, g, b);
        
        cout << "object num " << i << endl;
        for (int j = 0; j < voxObjList[i].size(); j++){ // look at each point in voxelized object
            ofPoint loc = voxObjList[i][j];
            ofDrawBox(loc, space, space, space);
        }
    }
}

void ofApp::drawBoundary(){
    ofColor color(0,200,0, 10);
    ofSetColor(color);
    
    for (int i = 0; i < boundaryList.size(); i++){
        ofPoint loc = boundaryList[i];
        ofDrawBox(loc, space, space, space);
    }
    
}

void ofApp::drawPath(){
    ofSetColor(200, 200, 0, 80);
    
    for (int i = 0; i < pathList.size(); i++){
        ofDrawBox(pathList[i], space, space, space);
    }
}

//--------------------------------------------------------------
void ofApp::setup(){
    
    //change camera angle
    cam.setDistance(30);
    cam.setPosition(3,5,30);
    
    
    space = g_voxelGrid->m_spacing;
    lowerleft = ofPoint(g_voxelGrid->m_lowerLeft[0], g_voxelGrid->m_lowerLeft[1], g_voxelGrid->m_lowerLeft[2]);
    
    setupBoundary();
    setupMesh();
    if(showPath){
        setupPath();
    }
    
    
    // dealing with the end points for the path
    ptStart = lowerleft;
    ptDest = lowerleft + ofPoint(g_voxelGrid->m_dimX*space, g_voxelGrid->m_dimY*space*.5, g_voxelGrid->m_dimZ*space*.6);
    ptVehicle = ptStart;
    
    //guiPanel.setup();
    
}

//--------------------------------------------------------------
void ofApp::update(){
    
}

//--------------------------------------------------------------
void ofApp::draw(){
    
    cam.begin();
    
    ofBackground(231,214,224);
    
    drawMesh();
    
    if(showPath){
        drawPath();
    }
    
    drawBoundary();
    
    ofSetColor(200, 0, 150, 50);
    ofDrawSphere(ptStart, 1);
    ofSetColor(200, 0, 150, 50);
    ofDrawSphere(ptDest, 1);
    
    //drawing the vehicle;
    ptVehicle = ptVehicle + ofPoint(); // TODO: figure out how to update this smoothly
    ofSetColor(0,0,200, 50);
    ofDrawSphere(ptVehicle, 0.5);
    
    cam.end();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
    
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){
    
}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y){
    
}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){
    
}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){
    
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
    
}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){
    
}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){
    
}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 
    
}
