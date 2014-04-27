#include "ofApp.h"

int valx;
int valy;
int valz;
int j;

//--------------------------------------------------------------
void ofApp::setup(){

	//ofSetFrameRate(500);

	valx = 500;
	valy = 300;
	valz = 10;
	j = 5;
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){

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