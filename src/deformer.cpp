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

	for (int i = 0; i < n_vert ; i ++){
		Vertex* curr_v = mMesh->vertices()[i];
		std::vector< HEdge* > adjHEdgeList;
		OneRingHEdge orhe(curr_v);
		HEdge* curr_e = nullptr;
		while (curr_e = orhe.nextHEdge())
			adjHEdgeList.push_back(curr_e);

		curr_pos = curr_v->position();
		x.row(i) = curr_pos.cast<double>();

		triplets.push_back(Eigen::Triplet<double>(curr_v->index(), curr_v->index(),1));
		for (int j = 0 ; j < adjHEdgeList.size(); j++)
			triplets.push_back(Eigen::Triplet<double>(curr_v->index(), adjHEdgeList[j]->end()->index(),-1.0/adjHEdgeList.size()));		
	}

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

	// // Constraint part
	for (int v_idx = 0; v_idx < n_constraint; v_idx++) {
		Vertex* curr_v = mRoiList[v_idx];
		Eigen::Vector3f v_pos = curr_v->position();
		b_.row(n_vert + v_idx) = v_pos.cast<double>();

	}
	Eigen::MatrixXd updated_x; 
	updated_x = mCholeskySolver->solve( A_.transpose() * b_);

	for(int r_idx = 0; r_idx < n_vert; r_idx++) 
		mMesh->vertices()[r_idx]->setPosition(updated_x.row(r_idx).cast<float>());

	std::cout << "Done solving linear system \n";

}
