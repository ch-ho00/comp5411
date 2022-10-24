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
	float weight = 1.0 ;
	int n_vert = (int) mMesh->vertices().size();
	int n_constraint =  mRoiList.size();

	Eigen::SparseMatrix< double > systemMat(3 * n_vert, 3 * n_vert);
	Eigen::SparseMatrix< double > A(3 * (n_vert + n_constraint),3 * n_vert); 
	std::vector< Eigen::Vector3f > deltas;

	/*====== Programming Assignment 2 ======*/
	// initialize triplet for sparse matrix
	std::vector< Eigen::Triplet<double> > triplet;
	// non-diagonal elements
	for (int v_idx = 0; v_idx < n_vert; v_idx++) {
		Vertex* curr_v = mMesh->vertices()[v_idx];
		std::vector< HEdge* > adjHEdgeList;
		OneRingHEdge orhe(curr_v);
		HEdge* curr_e = nullptr;
		while (curr_e = orhe.nextHEdge())
			adjHEdgeList.push_back(curr_e);

		int curIndex = curr_v->index();
		Eigen::Vector3f delta = curr_v->position();
		Eigen::Vector3f tmp_delta;

		for(int i=0; i< 3; i++)
			triplet.push_back(Eigen::Triplet<double>(3 * curIndex + i, 3 * curIndex + i, adjHEdgeList.size() * weight));

		for (int near_idx = 0; near_idx < adjHEdgeList.size(); near_idx++) {
			if (near_idx == 0)
				tmp_delta = delta - adjHEdgeList[near_idx]->end()->position();
			else
				tmp_delta = tmp_delta - adjHEdgeList[near_idx]->end()->position();

			int col_idx = adjHEdgeList[near_idx]->end()->index();
			for (int i=0; i<3; i++) 
				triplet.push_back(Eigen::Triplet<double>(curIndex * 3 + i, col_idx * 3 + i, - weight));
		}
		deltas.push_back(tmp_delta);
	}
	A.setFromTriplets(triplet.begin(), triplet.end());

	// extract b 
	// Laplace operator part
	Eigen::VectorXd b(3 * (n_vert + n_constraint));
	for (int r_idx = 0; r_idx < n_vert; r_idx++) {
		Eigen::Vector3f delta = deltas[r_idx];
		for (int i=0; i<3; i++)
			b[ r_idx * 3 + i] = delta[i];
	}
	// Constraint part
	for (int v_idx = 0; v_idx < n_constraint; v_idx++) {
		Vertex* curr_v = mRoiList[v_idx];
		Eigen::Vector3f v_pos = curr_v->position();
		for (int i=0; i<3; i++)
			b[(n_vert + v_idx) * 3 + i] = v_pos[i];
	}

	systemMat = A.transpose() * A;
	// Please refer to the following link for the usage of sparse linear system solvers in Eigen
	// https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html

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

	// solve linear system
	// Eigen::VectorXd x(3 * n_vert);
	// for (int j = 0; j< 3 * n_vert; j++)
	// 	x[j] = 0.0;

	// solver.compute(A);
	// x = solver.solve(b);

	// // update position
	// for (int v_idx = 0; v_idx < n_vert; v_idx++) {
	// 	Vertex* curr_v = mMesh->vertices()[v_idx];
	// 	Eigen::Vector3f updated;
	// 	for (int i=0; i<3; i++) 
	// 		updated[i] = x[curr_v->index() * 3 + i];
	// 	curr_v->setPosition(updated);
	// }

	}

	/*====== Programming Assignment 2 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* This is the place where the editing techniques 
	/* take place.
	/* Solve for the new vertex positions after the 
	/* specified handles move using the factorized
	/* matrix from Deformer::buildSystemMat(), i.e.,
	/* mCholeskySolver defined in deformer.h
	/**********************************************/

	// Please refer to the following link for the usage of sparse linear system solvers in Eigen
	// https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html


	/*====== Programming Assignment 2 ======*/
}
