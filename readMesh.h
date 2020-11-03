// Reading the Mesh File
if (myfile.is_open()) {
    string cur_line;
    int results[6];
    
    for (int k = 1; k <= 4; k++) {		// Ignoring starting lines
        getline(myfile, cur_line);
    }

    int n = cur_line.length();			// Reading the number nodes line
    char char_array[n + 1];
    strcpy(char_array, cur_line.c_str());
    sscanf(char_array, "(%d (%d %d %x %d %d)(')", &results[0], &results[1], &results[2], &results[3], &results[4], &results[5]);
    num_nodes = results[3] - results[2] + 1;
    cout << " Total num Nodes " << num_nodes << '\n';

	// Allocating memory for vertices ( num nodes = num vertices )
    vertices = (double **) malloc((num_nodes + 1) * sizeof (double *));
    for (int i = 0; i < num_nodes + 1; i++)
        vertices[i] = (double *) malloc(2 * sizeof (double));

    for (int k = 1; k <= 2; k++) {
        getline(myfile, cur_line);
    }

    cout << " Reading Vertices start " << '\n';
    for (int i = 1; i <= num_nodes; i++) {
        getline(myfile, cur_line);
        n = cur_line.length();
        char_array[n + 1];
        strcpy(char_array, cur_line.c_str());
        sscanf(char_array, "%lf %lf", &vertices[i][0], &vertices[i][1]);
    }

    for (int k = 0; k <= 1; k++) {
        getline(myfile, cur_line);
    }

    n = cur_line.length();				// Reading the Total number of cells line
    char_array[n + 1];
    strcpy(char_array, cur_line.c_str());
    sscanf(char_array, "(%d (%d %d %x %d %d)(')", &results[0], &results[1], &results[2], &results[3], &results[4], &results[5]);
    num_cells = results[3] - results[2] + 1;
    cout << " Total num Cells " << num_cells << '\n';

    for (int k = 0; k <= 1; k++) {
        getline(myfile, cur_line);
    }

    n = cur_line.length();				// Reading the Total number of faces line
    char_array[n + 1];
    strcpy(char_array, cur_line.c_str());
    sscanf(char_array, "(%d (%d %d %x %d %d)(')", &results[0], &results[1], &results[2], &results[3], &results[4], &results[5]);
    num_faces = results[3] - results[2] + 1;
    cout << " Total Faces = " << num_faces << '\n';

	// Allocating memory for Faces as we know total num faces
    faces = (int **) malloc((num_faces + 1) * sizeof (int *));
    for (int i = 0; i < num_faces + 1; i++)
        faces[i] = (int *) malloc(4 * sizeof (int));

    for (int k = 1; k <= 2; k++) {
        getline(myfile, cur_line);
    }

    n = cur_line.length();				// Reading the number of Internal faces line i,e without boundaries
    char_array[n + 1];
    strcpy(char_array, cur_line.c_str());
    sscanf(char_array, "(%d (%x %d %x %d %d)(')", &results[0], &results[1], &results[2], &results[3], &results[4], &results[5]);
    num_intFaces = results[3] - results[2] + 1;
    cout << " Total interior Faces = " << num_intFaces << '\n';

    cout << " Reading Interior Faces start " << '\n';
    for (int i = 0; i < num_intFaces; i++) {
        getline(myfile, cur_line);
        n = cur_line.length();
        char_array[n + 1];
        strcpy(char_array, cur_line.c_str());
        sscanf(char_array, "%x %x %x %x", &faces[faceSum][0], &faces[faceSum][1], &faces[faceSum][2], &faces[faceSum][3]);
        faceSum = faceSum + 1;
    }

    // Reading Boundary Conditions

    for (int zz = 0; zz < num_bc; zz++) {
        for (int k = 0; k <= 3; k++) {
            getline(myfile, cur_line);
        }

        n = cur_line.length();			// Reading number of boundary faces
        char_array[n + 1];
        strcpy(char_array, cur_line.c_str());
        sscanf(char_array, "(%d (%x %x %x %d %d)(')", &results[0], &results[1], &results[2], &results[3], &results[4], &results[5]);
        num_bcFaces[zz] = results[3] - results[2] + 1;
        printf(" Total bc%d Faces = %d  \n", zz, num_bcFaces[zz]);		// Printing boundary condition number zz

        printf(" Reading bc%d Faces start \n", zz);
        for (int i = 0; i < num_bcFaces[zz]; i++) {			// Reading boundary condition faces of number zz
            getline(myfile, cur_line);
            n = cur_line.length();
            char_array[n + 1];
            strcpy(char_array, cur_line.c_str());
            sscanf(char_array, "%x %x %x %x", &faces[faceSum][0], &faces[faceSum][1], &faces[faceSum][2], &faces[faceSum][3]);
            faceSum = faceSum + 1;
        }
    }

    myfile.close();						// Reading Mesh completed
}
