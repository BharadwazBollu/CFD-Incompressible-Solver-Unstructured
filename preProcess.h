
for (int i = 0; i < num_bc; i++) {
    num_totBC = num_totBC + num_bcFaces[i];				// Number of total boundary faces
}

vector<double> Vp(num_cells + 1 + num_totBC);			// Volume of each cell
vector<double> xc(num_cells + 1 + num_totBC);			// Cell center X-coordinate
vector<double> yc(num_cells + 1 + num_totBC);			// Cell center Y-coordinate

vector<vector<double> > fcX(num_cells + 1, vector<double> (4, 0));			// Face center X-coordinate
vector<vector<double> > fcY(num_cells + 1, vector<double> (4, 0));			// Face center Y-coordinate
vector<vector<double> > sn1(num_cells + 1, vector<double> (4, 0));			// Surface Normal in X-direction
vector<vector<double> > sn2(num_cells + 1, vector<double> (4, 0));			// Surface Normal in Y-direction
vector<vector<double> > normCoff(num_cells + 1, vector<double> (4, 0));		// Normal coefficient for Diffussion term
vector<vector<double> > crossCoff(num_cells + 1, vector<double> (4, 0));	// Cross  coefficient for Diffussion term
vector<vector<double> > faceWt(num_cells + 1, vector<double> (4, 0));		// Face Weight for each face
vector<vector<int> > snSign(num_cells + 1, vector<int> (4, 0));				// Surface Normal sign 1 = outwards, -1 = inwards 
vector<vector<int> > fout(num_cells + 1, vector<int> (4, 0));				// Outward cell number for corresponding cell
vector<vector<int> > lcf(num_cells + 1, vector<int> (4, 0));				// linking cell to face
vector<vector<int> > lnc(num_cells + 1, vector<int> (4, 0));				// linking node to cell
vector<vector<int> > cellOrder(num_cells + 1, vector<int> (4, 0));			// Cell Order for exporing to Tecplot/ParaView


int k = 1;
for (int j = num_intFaces + 1; j <= num_faces; j++) {
    faces[j][3] = num_cells + k;		// As BC cell number are zeros filling it with fictitious cell number greater than num_cells
    xc[num_cells + k] = 0.5 * (vertices[faces[j][0]][0] + vertices[faces[j][1]][0]);	// Cell center X-coordinate of BC i,e Face Center
    yc[num_cells + k] = 0.5 * (vertices[faces[j][0]][1] + vertices[faces[j][1]][1]);	// Cell center Y-coordinate of BC i,e Face Center
    k = k + 1;
}

for (int i = 1; i <= num_cells; i++) {
    k = 0;
    for (int j = 1; j <= num_faces; j++) {
        if (faces[j][2] == i) {
            lcf[i][k] = j;				// linking all face numbers to corresponding cell 
            xc[i] = xc[i] + 0.5 * 0.25 * (vertices[faces[j][0]][0] + vertices[faces[j][1]][0]);	// Cell center X-coordinate
            yc[i] = yc[i] + 0.5 * 0.25 * (vertices[faces[j][0]][1] + vertices[faces[j][1]][1]);	// Cell center Y-coordinate
            k = k + 1;
        } else if (faces[j][3] == i) {
            lcf[i][k] = j;
            xc[i] = xc[i] + 0.5 * 0.25 * (vertices[faces[j][0]][0] + vertices[faces[j][1]][0]);
            yc[i] = yc[i] + 0.5 * 0.25 * (vertices[faces[j][0]][1] + vertices[faces[j][1]][1]);
            k = k + 1;
        }
    }
}

for (int i = 1; i <= num_cells; i++) {
    Vp[i] = 0;
    for (int k = 0; k < 4; k++) {

        double dx, dy, dot, area_sf, dist, ex, ey, normX, normY, tx, ty;

        sn1[i][k] = vertices[faces[lcf[i][k]][1]][1] - vertices[faces[lcf[i][k]][0]][1];	// X Surface Normal
        sn2[i][k] = vertices[faces[lcf[i][k]][0]][0] - vertices[faces[lcf[i][k]][1]][0];	// Y Surface Normal

        fcX[i][k] = 0.5 * (vertices[faces[lcf[i][k]][0]][0] + vertices[faces[lcf[i][k]][1]][0]) ;	// X Face Center
        fcY[i][k] = 0.5 * (vertices[faces[lcf[i][k]][0]][1] + vertices[faces[lcf[i][k]][1]][1]) ;	// Y Face Center

        dx = fcX[i][k] - xc[i];			// X direction vector for cell to corresponding face center
        dy = fcY[i][k] - yc[i];			// Y direction vector for cell to corresponding face center

        dot = sn1[i][k] * dx + sn2[i][k] * dy;		// Calculating Dot Product of Surface Vector and Direction vector

        if (dot < 0) {								// Checking Whether Surface Normal is inwards < 0 or outwards > 0 
            snSign[i][k] = -1;						// Surface Normal sign for inwards
            sn1[i][k] = -sn1[i][k];					// Flipping sign to make corresponding surface normal outwards 
            sn2[i][k] = -sn2[i][k];					
            fout[i][k] = faces[lcf[i][k]][2];		// Storing corresponding outward cell number for cell for that face number
        } else {
            fout[i][k] = faces[lcf[i][k]][3];
            snSign[i][k] = 1;
        }

		// Calculating Volume of cell using Green Gauss Method
        Vp[i] = Vp[i] + sn1[i][k] * 0.5 * (vertices[faces[lcf[i][k]][0]][0] + vertices[faces[lcf[i][k]][1]][0]);

        area_sf = sqrt(pow(sn1[i][k], 2) + pow(sn2[i][k], 2));		// Magnitude of Surface Normal

        dist = sqrt(pow((xc[faces[lcf[i][k]][3]] - xc[faces[lcf[i][k]][2]]), 2) + pow((yc[faces[lcf[i][k]][3]] - yc[faces[lcf[i][k]][2]]), 2));			// Distance between Cell center and Neighbour Cell center
        
        // Caclculating unit Surface Normals for Over Relaxed Approach
        ex = (xc[faces[lcf[i][k]][3]] - xc[faces[lcf[i][k]][2]]) / dist;
        ey = (yc[faces[lcf[i][k]][3]] - yc[faces[lcf[i][k]][2]]) / dist;

        normX = (pow(sn1[i][k], 2) + pow(sn2[i][k], 2)) * ex / (sn1[i][k] * ex + sn2[i][k] * ey) ;
        normY = (pow(sn1[i][k], 2) + pow(sn2[i][k], 2)) * ey / (sn1[i][k] * ex + sn2[i][k] * ey) ;

		// Caclculating Normal Coefficient
        normCoff[i][k] = sqrt(pow(normX, 2) + pow(normY, 2)) / dist;

        tx = sn1[i][k] - normX;
        ty = sn2[i][k] - normY;

		// Caclculating Cross Coefficient
        crossCoff[i][k] = sqrt(tx * tx + ty * ty) / area_sf;
    }
}

// Calculating Face Weights
for (int i = 1; i <= num_cells; i++) {
    int xx, yy;
    for (int k = 0; k < 4; k++) {
        if (snSign[i][k] > 0) {		// Checking for Outward Cell Number
            xx = 2;
            yy = 3;
        } else {
            xx = 3;
            yy = 2;
        }
        faceWt[i][k] = Vp[faces[lcf[i][k]][yy]] / (Vp[faces[lcf[i][k]][2]] + Vp[faces[lcf[i][k]][3]]);
    }
}

// Calculating Max number of cells a Node Shares for Whole Mesh
int max_cellShare = 0;									// Max number of cells sharing	
for (int j = 1; j <= num_nodes; j++) {
    int zz = 1;
    for (int i = 1; i <= num_cells; i++) {
        for (int k = 0; k < 4; k++) {
            if (faces[lcf[i][k]][0] == j) {				// Checking for Cell Number of corresponding Node
                zz = zz + 1;
            } else if (faces[lcf[i][k]][1] == j) {
                zz = zz + 1;
            }
        }
    }
    if (zz > max_cellShare) {
        max_cellShare = zz;
    }
}
// Substracting 1 as initially zz = 1, and dividing it by 2 as each face has two vertex i,e nodes become twice for a cell of all faces
max_cellShare = (max_cellShare - 1) / 2;

// Allocating Memory for Point Weight(Node), and linking Cell to Node
vector<vector<double> > pointWt(num_nodes + 1, vector<double> (max_cellShare, 0));		// Node Weight
vector<vector<int> > lcn(num_nodes + 1, vector<int> (max_cellShare, 0));				// linking Cell to Node

for (int j = 1; j <= num_nodes; j++) {
    vector<int> pointNode(2 * max_cellShare, 0);
    int zz = 0;
    for (int i = 1; i <= num_cells; i++) {
        for (int k = 0; k < 4; k++) {
            if (faces[lcf[i][k]][0] == j) {
                pointNode[zz] = i;							// Storing node sharing cell numbers
                zz = zz + 1;
            } else if (faces[lcf[i][k]][1] == j) {
                pointNode[zz] = i;
                zz = zz + 1;
            }
        }
    }

    int k = 0;
    int temp = 0;
    for (int zz = 0; zz < 2 * max_cellShare - 1; zz++) {
        if (lcn[j][k] != temp) {
            temp = pointNode[zz];
            k = k + 1;
        }
        lcn[j][k] = pointNode[zz];
    }

    double tmp = 0;
    vector<double> pV(max_cellShare, 0);

    for (int k = 0; k < max_cellShare; k++) {
        if (lcn[j][k] > 0) {
            pV[k] = Vp[lcn[j][k]];
        }
        if (pV[k] > 0) {
            tmp = tmp + 1 / pV[k];
        }
    }

    for (int k = 0; k < max_cellShare; k++) {
        if (pV[k] > 0) {
            pointWt[j][k] = (1 / pV[k]) / tmp;				// Calculating Node Weight = 1/Vp[i]/( Sum of all Vp[i] )
        } else {
            pointWt[j][k] = 0;
        }
    }
}

for (int i = 1; i <= num_cells; i++) {
    int zz = 0;
    for (int j = 1; j <= num_nodes; j++) {
        for (int k = 0; k < max_cellShare; k++) {
            if (lcn[j][k] == i) {
                lnc[i][zz] = j;
                zz = zz + 1;
            }
        }
    }
}

// Calculating Cell Order for Exporting using Sorting
for (int i = 1; i <= num_cells; i++) {
    double angle[4];
    double angleTemp[4];
    for (int k = 0; k < 4; k++) {
        angle[k] = atan2((vertices[lnc[i][k]][1] - yc[i]), (vertices[lnc[i][k]][0] - xc[i]));
        angleTemp[k] = angle[k];
    }
    sort(angleTemp, angleTemp + 4);
    for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
            if (angle[k] == angleTemp[j]) {
                cellOrder[i][j] = lnc[i][k];
            }
        }
    }
}

