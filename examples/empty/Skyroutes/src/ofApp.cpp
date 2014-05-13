#include "ofApp.h"


int valx;
int valy;
int valz;
int j;



//--------------------------------------------------------------
ofApp::ofApp(Mesh* m, bool vox){
    reMesh = m;
    useVoxel = vox;
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

//--------------------------------------------------------------
void ofApp::setupVoxel(){
    
    ofPoint lowerleft = ofPoint(g_voxelGrid->m_lowerLeft[0], g_voxelGrid->m_lowerLeft[1], g_voxelGrid->m_lowerLeft[2]);
    double space = g_voxelGrid->m_spacing;
    
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
    ofPoint lowerleft = ofPoint(g_voxelGrid->m_lowerLeft[0], g_voxelGrid->m_lowerLeft[1], g_voxelGrid->m_lowerLeft[2]);
    double space = g_voxelGrid->m_spacing;
    
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

void ofApp::drawMesh(){
    ofSetColor(0, 155, 155);
    rendMesh.drawFaces();
    
    ofSetColor(255,0,0);
    rendMesh.drawWireframe();
}

void ofApp::drawVoxel(){
    
    int color = 256 * 3 / reMesh->numObj;
    int div = reMesh->numObj / 3;
    double space = g_voxelGrid->m_spacing;
    
    for (int i = 1; i < reMesh->numObj; i++){ // look at each voxelized object
        //for (int i = reMesh->numObj - 1; i > 0; i--){
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
    double space = g_voxelGrid->m_spacing;
    ofColor color(0,200,0, 10);
    ofSetColor(color);
    
    for (int i = 0; i < boundaryList.size(); i++){
        ofPoint loc = boundaryList[i];
        ofDrawBox(loc, space, space, space);
    }
    
}

//--------------------------------------------------------------
void ofApp::setup(){
    
	//ofSetFrameRate(500);
    
	valx = 500;
	valy = 300;
	valz = 10;
	j = 5;
    
    //change camera angle
    cam.setDistance(15);
    
    //float a = cam.getDistance();
    //cout << "cam dist = " << a << endl;
    
    
    setupBoundary();
    
    if(useVoxel){
        setupVoxel();
    }
    else{
        setupMesh();
    }    
    
}

//--------------------------------------------------------------
void ofApp::update(){
    
}

//--------------------------------------------------------------
void ofApp::draw(){
    
    cam.begin();
    
    ofBackground(50,50,50);
    
    /*
     valz += j;
     if (valz > 500){
     j = -5;
     }
     else if (valz < 0){
     j = 5;
     }
     
     ofSetColor(255,0,255);
     //ofCircle(valx, valy, 100);
     ofDrawSphere(valx, valy, valz, 30);
     
     ofSetColor(100, 100, 255);
     ofDrawBox(200, 400, 0, 40, 300, 50);
     
     ofDrawBox(800, 500, 0, 80, 200, 50);
     */
    
    
    
    if(useVoxel){
        drawVoxel();
    }
    else{
        drawMesh();
    }
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawBoundary();
    //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    
    
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
