#include "deformer.h"
#include <iostream>

Deformer::Deformer() : mMesh(nullptr),
                       mCholeskySolver(nullptr) {
}

Deformer::~Deformer() {
	clear();
}

void Deformer::clear() {
	if (mCholeskySolver) {
		delete mCholeskySolver;
	}
	mCholeskySolver = nullptr;
	mRoiList.clear();
}

void Deformer::setMesh(Mesh* mesh) {
	mMesh = mesh;
	clear();
	// Record the handle vertices
	for (Vertex* vert : mMesh->vertices()) {
		if (vert->flag() > 0 || vert->isBoundary()) {
			mRoiList.push_back(vert);
		}
	}
	// Build system matrix for deformation
	buildSystemMat();
}


void Deformer::buildSystemMat() {
	/*====== Programming Assignment 2 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Build the matrix of the linear system for 
	/* deformation and do factorization, in order
	/* to reuse and speed up in Deformer::deform().
	/* Handle vertices are maked by Vertex::flag() > 0
	/* Movements of the specified handle are already
	/* recorded in Vertex::position()
	/**********************************************/
	int n_vert = mMesh->vertices().size();
	int n_constraint =  mRoiList.size();

	Eigen::SparseMatrix< double > systemMat;
	Eigen::SparseMatrix< double > A(n_vert + n_constraint, n_vert); 
	std::vector< Eigen::Vector3f > deltas;
	
	Eigen::Vector3f curr_pos;
	Eigen::MatrixXd x (n_vert, 3);
	std::vector< Eigen::Triplet<double> > triplets;

	Eigen::Vector3f p1 = Eigen::Vector3f::Zero();
	Eigen::Vector3f p2 = Eigen::Vector3f::Zero();
	Eigen::Vector3f p3 = Eigen::Vector3f::Zero();

	Eigen::Vector3f hedge_1;
	Eigen::Vector3f hedge_2;
	Eigen::Vector3f hedge_3;
	Eigen::Vector3f hedge_4;

	double cot1;
	double cot2;

	for (int i = 0; i < n_vert ; i ++){
		Vertex* curr_v = mMesh->vertices()[i];
		std::vector< HEdge* > adjHEdgeList;
		OneRingHEdge orhe(curr_v);
		HEdge* curr_e = nullptr;
		// collect neighbours
		while (curr_e = orhe.nextHEdge())
			adjHEdgeList.push_back(curr_e);

		float sum_w = 0.0;
		std::vector < float > ws;
		// calculate weight
		for (int j = 0; j < adjHEdgeList.size(); j++) {
			HEdge* outHEdge = adjHEdgeList[j];
			p1 = outHEdge->end()->position();
			p2 = outHEdge->next()->end()->position();
			p3 = outHEdge->next()->next()->end()->position();
			hedge_1 = p1 - p2;
			hedge_2 = p3 - p2;
			cot1 = hedge_1.dot(hedge_2) / hedge_1.cross(hedge_2).norm();

			HEdge* inHEdge = outHEdge->twin();
			p1 = inHEdge->end()->position();
			p2 = inHEdge->next()->end()->position();
			p3 = inHEdge->next()->next()->end()->position();

			hedge_3 = p1 - p2;
			hedge_4 = p3 - p2;
			cot2 = hedge_3.dot(hedge_4) / hedge_3.cross(hedge_4).norm();
			// sum of cotangent weights 
			sum_w += (cot1 + cot2) / 2;
			ws.push_back((cot1 + cot2) / 2);
		}

		curr_pos = curr_v->position();
		x.row(i) = curr_pos.cast<double>();

		// set diagonal element as 1 
		triplets.push_back(Eigen::Triplet<double>(curr_v->index(), curr_v->index(),1));

		// set non-diagonal element as cotangent weight 
		for (int j = 0 ; j < adjHEdgeList.size(); j++)
			triplets.push_back(Eigen::Triplet<double>(curr_v->index(), adjHEdgeList[j]->end()->index(),-ws[j] / sum_w));		
	}

	// set modelling constraints
	for (int i =0 ; i < n_constraint; i ++){
		int constraint_col = mRoiList[i]->index();
		triplets.push_back(Eigen::Triplet<double>(n_vert + i, constraint_col, 1));				
	}
	A.setFromTriplets(triplets.begin(), triplets.end());
	systemMat = A.transpose() * A;
	Eigen::MatrixXd b = A * x;
	A_ = A;
	b_ = b;

	// Do factorization
	if (systemMat.nonZeros() > 0) {
		mCholeskySolver = new Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >();
		mCholeskySolver->compute(systemMat);
		if (mCholeskySolver->info() != Eigen::Success) {
			// Decomposition failed
			std::cout << "Sparse decomposition failed\n";
		} else {
			std::cout << "Sparse decomposition succeeded\n";
		}
	}
}

void Deformer::deform() {
	if (mCholeskySolver == nullptr) {
		return;
	}
	int n_vert = mMesh->vertices().size();
	int n_constraint =  mRoiList.size();

	// Set modeliing constraint value
	for (int v_idx = 0; v_idx < n_constraint; v_idx++) {
		Vertex* curr_v = mRoiList[v_idx];
		Eigen::Vector3f v_pos = curr_v->position();
		b_.row(n_vert + v_idx) = v_pos.cast<double>();

	}
	Eigen::MatrixXd updated_x; 
	// Solve linear system
	updated_x = mCholeskySolver->solve( A_.transpose() * b_);

	// Update vertices location
	for(int r_idx = 0; r_idx < n_vert; r_idx++) 
		mMesh->vertices()[r_idx]->setPosition(updated_x.row(r_idx).cast<float>());

	std::cout << "Done solving linear system \n";

}
