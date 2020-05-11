/*
** This program solves the linear elastic equations in a
** random linear elastic material, subject to an applied macroscopic strain,
** using the finite element method. Each pixel in the 3-D digital
** image is a cubic tri-linear finite element, having its own
** elastic moduli tensor. Periodic boundary conditions are maintained.
** In the comments below, (USER) means that this is a section of code that
** the user might have to change for his particular problem. Therefore the
** user is encouraged to search for this string.
*/

//#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>
using namespace std;

//YZ:added 200212
std::ofstream stress_fout("stress.txt"); // store stresses of each CSH after equilibrium of each step
std::ofstream phasechange_fout("phasechange.txt"); // store phase change at which pixel
std::ofstream phasechange_csh_fout("phasechange_csh.txt"); // store phase change at which pixel (only for CSH)
std::ofstream fout("outputfile.txt"); // store the bulk stress due to applied strains



//global variables: 100x100x100 voxels, assume 50 input files, 3 directions
double u[1000001][51][4], gb[1000001][4], b[1000001][4]; //node displacement, gradient, b value
double h[1000001][4], Ah[1000001][4]; // for solving the gb+b=0
double cmod[101][7][7][51], dk[101][9][4][9][4][51]; //modulus matrix

long ib[1000001][28]; //for global node, 1000000 nodes, each node has 27 neighbors (including self)
short pix[1000001], pixstor[1000001], pixt[1000001]; //phase parameters and aging time
int phasechange[1000001]; //phasechange value

double strxx, stryy, strzz, strxz, stryz, strxy; //averaged global stresses
double sxx, syy, szz, sxz, syz, sxy; //averaged global strains
double exx, eyy, ezz, exz, eyz, exy; //applied global strains
double C; //constant C

double str11, str22, str33, str12, str23, str13; //pixel stresses
double s11, s22, s33, s12, s23, s13; //pixel strains
double dxx, dyy, dzz, dxy, dxz, dyz; //to get the strain in one pixel



//Subroutine that sets up microstructural image
void ppixel(int nx, int ny, int nz, int ns, int nphase, char filename[15])
{
	ifstream fin(filename);
	
	std::cout << "  Read in microstructure file: " << filename << std::endl;
	
	for (int k = 1; k <= nz; k++)
	{
		for (int j = 1; j <= ny; j++)
		{
			for (int i = 1; i <= nx; i++)
			{
				int m = nx * ny * (k - 1) + nx * (j - 1) + i;
				fin >> pix[m]; //assign phase labels of 100x100x100 image to array pix[m]
			}
		}
	}

	for (int m = 1; m <= ns; m++)
	{
		//check for wrong phase labels -- less than 1 or greater than nphase
		//YZ: phase 0 is void, phase 1 is water
		/*
		if (pix[m]<1)
		{
			std::cout << "Phase label in pix < 1 -- error at " << m << "\n";
		}
		*/

		if (pix[m] > nphase)
		{
			std::cout << "Phase label in pix > nphase -- error at " << m << "\n";
		}
	}
}



//Subroutine that counts volume fractions
void assig(int ns, int nphase, double (&prob)[101]) //double prob[101] stores volume fraction of each phase
{	
	std::cout << "  Calculate volume fractions. " << std::endl;
	for (int i = 1; i <= nphase; i++)
	{
		prob[i] = 0.0;
	}

	for (int m = 1; m <= ns; m++)
	{
		for (int i = 1; i <= nphase; i++)
		{
			if (pix[m] == i)
			{
				prob[i] = prob[i] + 1;
				break;
			}
		}
	}

	for (int i = 1; i <= nphase; i++)
	{
		prob[i] = prob[i] / (ns * 1.0);

		//TODO: output a file contains volume fractions of nphase 
		//std::ofstream prob_fout("prob.txt");
  		//prob_fout << prob[i] << std::endl;
		//prob_fout.close();
	}
}



//Calculate the strain in one pixel at time qt
void strain(int i, int j, int  k, int nx, int ny, int  nz, int qt)
{
	double dndx[9], dndy[9], dndz[9]; //dndx, dndy, dndz are strain components in 3 directions (x,y,z) of 8 nodes
	//each pixel has 6 independent deformation gradient compuonents, 8 nodes, 3 directions
	double es[7][9][4]; //es is the displacement/deformation gradient that relates nodal displacement to strains, dndx = es * uu
	double uu[9][4]; //uu is the displacement vectors of 8 nodes
	
	//set up single pixel strain matrix
	//dndx, dndy, dndz are the components of the strain matrix in a pixel
	dndx[1] = -0.25; //MANUAL page 170, see figure 1 for node label
	dndx[2] = 0.25;
	dndx[3] = 0.25;
	dndx[4] = -0.25;
	dndx[5] = -0.25;
	dndx[6] = 0.25;
	dndx[7] = 0.25;
	dndx[8] = -0.25;

	dndy[1] = -0.25;
	dndy[2] = -0.25;
	dndy[3] = 0.25;
	dndy[4] = 0.25;
	dndy[5] = -0.25;
	dndy[6] = -0.25;
	dndy[7] = 0.25;
	dndy[8] = 0.25;

	dndz[1] = -0.25;
	dndz[2] = -0.25;
	dndz[3] = -0.25;
	dndz[4] = -0.25;
	dndz[5] = 0.25;
	dndz[6] = 0.25;
	dndz[7] = 0.25;
	dndz[8] = 0.25;

	//build average strain matrix, follows code in femat, 
	//but for average strain over the pixel, not the strain at a point.
	for (int n1 = 1; n1 <= 6; n1++)
	{
		for (int n2 = 1; n2 <= 8; n2++)
		{
			for (int n3 = 1; n3 <= 3; n3++)
			{
				es[n1][n2][n3] = 0.0;
			}
		}
	}

	for (int n = 1; n <= 8; n++)
	{
		// es is deformation gradient, dndx = es * uu, when uu ==1, es = dndx
		es[1][n][1] = dndx[n]; //es11, xx
		es[2][n][2] = dndy[n]; //es22, yy
		es[3][n][3] = dndz[n]; //es33, zz
		es[4][n][1] = dndz[n]; //es41, xz-x
		es[4][n][3] = dndx[n]; //es43, xz-z
		es[5][n][2] = dndz[n]; //es52, yz-y
		es[5][n][3] = dndy[n]; //es53, yz-z
		es[6][n][1] = dndy[n]; //es61, xy-x
		es[6][n][2] = dndx[n]; //es62, xy-y
	}

	//assign nodal displacement uu of a pixel according to global nodel displacement vector u
	int m1 = nx * ny * (k - 1) + nx * (j - 1) + i;
	for (int mm = 1; mm <= 3; mm++)
	{
		//uu[i-th point][direction] displacement = u[m1-voxel][time][direction] displacement
		//ib[m][28] is global vector, stores 27 neigbhor labels of each pixel
		uu[1][mm] = u[m1][qt][mm]; //point-1 displacement = self displacement
		uu[2][mm] = u[ib[m1][3]][qt][mm]; //point-2 displacement = neighbor-3 displacement (see figure 1,table 3 in the MANUAL)
		uu[3][mm] = u[ib[m1][2]][qt][mm]; //point-3 displacement = neighbor-2 displacement
		uu[4][mm] = u[ib[m1][1]][qt][mm]; //point-4 displacement = neighbor-1 displacement
		uu[5][mm] = u[ib[m1][26]][qt][mm]; //point-5 displacement = neighbor-26 displacement
		uu[6][mm] = u[ib[m1][19]][qt][mm]; //point-6 displacement = neighbor-19 displacement
		uu[7][mm] = u[ib[m1][18]][qt][mm]; //point-7 displacement = neighbor-18 displacement
		uu[8][mm] = u[ib[m1][17]][qt][mm]; //point-8 displacement = neighbor-17 displacement
	}

	//Correct for nodal displacements at boundary surfaces, see MANUAL page 111
	//boundary 100 node take displacment from 1 node, should correct for the displacement accross the thickness
	//exx, eyy, ezz, exy, eyz, exz are applied global strains on the system
	if (i == nx) 
	{
		//at x=100 surface, node 2, 3, 6, 7 are on the surface		
		uu[2][1] = uu[2][1] + exx * nx; 
		uu[2][2] = uu[2][2] + exy * nx;
		uu[2][3] = uu[2][3] + exz * nx;
		uu[3][1] = uu[3][1] + exx * nx;
		uu[3][2] = uu[3][2] + exy * nx;
		uu[3][3] = uu[3][3] + exz * nx;
		uu[6][1] = uu[6][1] + exx * nx;
		uu[6][2] = uu[6][2] + exy * nx;
		uu[6][3] = uu[6][3] + exz * nx;
		uu[7][1] = uu[7][1] + exx * nx;
		uu[7][2] = uu[7][2] + exy * nx;
		uu[7][3] = uu[7][3] + exz * nx;
	}

	if (j == ny)
	{
		//at y=100 surface, node 3, 4, 7, 8 are on the surface
		uu[3][1] = uu[3][1] + exy * ny;
		uu[3][2] = uu[3][2] + eyy * ny;
		uu[3][3] = uu[3][3] + eyz * ny;
		uu[4][1] = uu[4][1] + exy * ny;
		uu[4][2] = uu[4][2] + eyy * ny;
		uu[4][3] = uu[4][3] + eyz * ny;
		uu[7][1] = uu[7][1] + exy * ny;
		uu[7][2] = uu[7][2] + eyy * ny;
		uu[7][3] = uu[7][3] + eyz * ny;
		uu[8][1] = uu[8][1] + exy * ny;
		uu[8][2] = uu[8][2] + eyy * ny;
		uu[8][3] = uu[8][3] + eyz * ny;
	}

	if (k == nz)
	{
		//at z=100 surface, node 5, 6, 7, 8 are on the surface
		uu[5][1] = uu[5][1] + exz * nz;
		uu[5][2] = uu[5][2] + eyz * nz;
		uu[5][3] = uu[5][3] + ezz * nz;
		uu[6][1] = uu[6][1] + exz * nz;
		uu[6][2] = uu[6][2] + eyz * nz;
		uu[6][3] = uu[6][3] + ezz * nz;
		uu[7][1] = uu[7][1] + exz * nz;
		uu[7][2] = uu[7][2] + eyz * nz;
		uu[7][3] = uu[7][3] + ezz * nz;
		uu[8][1] = uu[8][1] + exz * nz;
		uu[8][2] = uu[8][2] + eyz * nz;
		uu[8][3] = uu[8][3] + ezz * nz;
	}

	//local strains in a pixel
	dxx = 0.0;
	dyy = 0.0;
	dzz = 0.0;
	dxz = 0.0;
	dyz = 0.0;
	dxy = 0.0;

	//calculate pixel strain = deformation gradient * nodal displacements
	for (int n3 = 1; n3 <= 3; n3++)
	{
		for (int n8 = 1; n8 <= 8; n8++)
		{
			dxx = dxx + es[1][n8][n3] * uu[n8][n3];
			dyy = dyy + es[2][n8][n3] * uu[n8][n3];
			dzz = dzz + es[3][n8][n3] * uu[n8][n3];
			dxz = dxz + es[4][n8][n3] * uu[n8][n3];
			dyz = dyz + es[5][n8][n3] * uu[n8][n3];
			dxy = dxy + es[6][n8][n3] * uu[n8][n3];
		}
	}
}



//Subroutine that sets up the elastic moduli variables,
//the stiffness matrices, dk, 
//the linear term in displacements, b,
//and the constant term, C, that appear in the total energy due to the periodic boundary conditions
/*
** (USER) NOTE: complete elastic modulus matrix is used, so an anisotropic
** matrix could be directly input at any point, since program is written
** to use a general elastic moduli tensor, but is only explicitly
** implemented for isotropic materials.
*/
void femat(int nx, int ny, int nz, int ns, double (&phasemod)[101][51][3], int nphase, int q)
{
	std::cout << "  Compute stiffness matrix dk, and b, C." << std::endl;

	double dndx[9], dndy[9], dndz[9]; //derivatives of shape functions
	double ck[7][7], cmu[7][7]; //bulk and shear matrix
	double g[4][4][4]; //weight vector of Simpson's rule
	double es[7][9][4]; //deformation gradient
	double delta[9][4]; 
	int is[9]; //neighbor labels

	//Set up elastic moduli matrices cmod for each kind of element
	//ck and cmu are the bulk and shear modulus matrices, which need to be weighted by the actual bulk and shear moduli

	//bulk modulus matrix for an element, 6x6 components (6 independent)
	ck[1][1] = 1.0; 
	ck[1][2] = 1.0;
	ck[1][3] = 1.0;
	ck[1][4] = 0.0;
	ck[1][5] = 0.0;
	ck[1][6] = 0.0;

	ck[2][1] = 1.0;
	ck[2][2] = 1.0;
	ck[2][3] = 1.0;
	ck[2][4] = 0.0;
	ck[2][5] = 0.0;
	ck[2][6] = 0.0;

	ck[3][1] = 1.0;
	ck[3][2] = 1.0;
	ck[3][3] = 1.0;
	ck[3][4] = 0.0;
	ck[3][5] = 0.0;
	ck[3][6] = 0.0;

	ck[4][1] = 0.0;
	ck[4][2] = 0.0;
	ck[4][3] = 0.0;
	ck[4][4] = 0.0;
	ck[4][5] = 0.0;
	ck[4][6] = 0.0;
	ck[5][1] = 0.0;
	ck[5][2] = 0.0;
	ck[5][3] = 0.0;
	ck[5][4] = 0.0;
	ck[5][5] = 0.0;
	ck[5][6] = 0.0;

	ck[6][1] = 0.0;
	ck[6][2] = 0.0;
	ck[6][3] = 0.0;
	ck[6][4] = 0.0;
	ck[6][5] = 0.0;
	ck[6][6] = 0.0;

	//shear modulus matrix for an element, 6x6 components
	cmu[1][1] = 4.0 / 3.0; 
	cmu[1][2] = -2.0 / 3.0;
	cmu[1][3] = -2.0 / 3.0;
	cmu[1][4] = 0.0;
	cmu[1][5] = 0.0;
	cmu[1][6] = 0.0;

	cmu[2][1] = -2.0 / 3.0;
	cmu[2][2] = 4.0 / 3.0;
	cmu[2][3] = -2.0 / 3.0;
	cmu[2][4] = 0.0;
	cmu[2][5] = 0.0;
	cmu[2][6] = 0.0;

	cmu[3][1] = -2.0 / 3.0;
	cmu[3][2] = -2.0 / 3.0;
	cmu[3][3] = 4.0 / 3.0;
	cmu[3][4] = 0.0;
	cmu[3][5] = 0.0;
	cmu[3][6] = 0.0;

	cmu[4][1] = 0.0;
	cmu[4][2] = 0.0;
	cmu[4][3] = 0.0;
	cmu[4][4] = 1.0;
	cmu[4][5] = 0.0;
	cmu[4][6] = 0.0;

	cmu[5][1] = 0.0;
	cmu[5][2] = 0.0;
	cmu[5][3] = 0.0;
	cmu[5][4] = 0.0;
	cmu[5][5] = 1.0;
	cmu[5][6] = 0.0;

	cmu[6][1] = 0.0;
	cmu[6][2] = 0.0;
	cmu[6][3] = 0.0;
	cmu[6][4] = 0.0;
	cmu[6][5] = 0.0;
	cmu[6][6] = 1.0;

	//cmod is the elastic moduli tensor (6x6 matrix) with 2 independent material properties to relate strain to stress of in a pixel
	//ck and cmu are bulk and shear modulus matrix, 6x6 components in matrix form (rank is 1)
	//phasemod is bulk and shear modulus single value for each phase-k at time-q
	for (int k = 1; k <= nphase; k++)
	{
		for (int j = 1; j <= 6; j++)
		{
			for (int i = 1; i <= 6; i++)
			{
				cmod[k][i][j][q] = phasemod[k][q][1] * ck[i][j] + phasemod[k][q][2] * cmu[i][j]; //make cmod the complete modulus tensor
				//phasemod [phase-k][time-q][bulk/shear]
				//bulk ck[6][6], shear cmu[6][6]
				//cmod [phase-k][6][6][time-q]
			}
		}
	}

	/*
	** dndx means the derivative of the shape matrix N with respect to x,
	** (see manual, Sec. 2.2), dndy, and dndz are similar.
	*/
	//see MANUAL page 15, set derivatives of shape functions that relate nodal displacements to displacement field in an element	
	for (int k = 1; k <= 3; k++)
	{
		for (int j = 1; j <= 3; j++)
		{
			for (int i = 1; i <= 3; i++)
			{	
				//the loop is just to get 0, 0.5, 1. for x, y, z, no further meaning
				double	x = 1.0 * (i - 1) / 2.0;
				double	y = 1.0 * (j - 1) / 2.0;
				double  z = 1.0 * (k - 1) / 2.0;

				dndx[1] = -(1.0 - y) * (1.0 - z);
				dndx[2] = (1.0 - y) * (1.0 - z);
				dndx[3] = y * (1.0 - z);
				dndx[4] = -y * (1.0 - z);
				dndx[5] = -(1.0 - y) * z;
				dndx[6] = (1.0 - y) * z;
				dndx[7] = y * z;
				dndx[8] = -y * z;

				dndy[1] = -(1.0 - x) * (1.0 - z);
				dndy[2] = -x * (1.0 - z);
				dndy[3] = x * (1.0 - z);
				dndy[4] = (1.0 - x) * (1.0 - z);
				dndy[5] = -(1.0 - x) * z;
				dndy[6] = -x * z;
				dndy[7] = x * z;
				dndy[8] = (1.0 - x) * z;

				dndz[1] = -(1.0 - x) * (1.0 - y);
				dndz[2] = -x * (1.0 - y);
				dndz[3] = -x * y;
				dndz[4] = -(1.0 - x) * y;
				dndz[5] = (1.0 - x) * (1.0 - y);
				dndz[6] = x * (1.0 - y);
				dndz[7] = x * y;
				dndz[8] = (1.0 - x) * y;
			}
		}
	}

	//build deformation gradient that relates nodal displacements to local strain in a pixel
	for (int n1 = 1; n1 <= 6; n1++)
	{
		for (int n2 = 1; n2 <= 8; n2++)
		{
			for (int n3 = 1; n3 <= 3; n3++)
			{
				es[n1][n2][n3] = 0.0;
			}
		}
	}

	for (int n = 1; n <= 8; n++)
	{
		//assign deformation gradient
		es[1][n][1] = dndx[n];
		es[2][n][2] = dndy[n];
		es[3][n][3] = dndz[n];
		es[4][n][1] = dndz[n];
		es[4][n][3] = dndx[n];
		es[5][n][2] = dndz[n];
		es[5][n][3] = dndy[n];
		es[6][n][1] = dndy[n];
		es[6][n][2] = dndx[n];
	}

	/*
	** loop over the nphase kinds of pixels and Simpson's rule quadrature
	** points in order to compute the stiffness matrices dk. Stiffness matrices
	** of trilinear finite elements are quadratic in x, y, and z, so that
	** Simpson's rule quadrature gives exact results.
	*/
	// set up Simpson's integration rule weight vector
	for (int k = 1; k <= 3; k++)
	{
		for (int j = 1; j <= 3; j++)
		{
			for (int i = 1; i <= 3; i++)
			{
				int	nm = 0;
				if (i == 2)
				{
					nm = nm + 1;
				}
				if (j == 2)
				{
					nm = nm + 1;
				}
				if (k == 2)
				{
					nm = nm + 1;
				}
				g[i][j][k] = pow(4.0, nm); //g is weight vector of Simpson's rule
			}
		}
	}

	//dk is the stiffness matrix to relate nodal displacements to nodal forces for each finite element
	//initialize stiffness matrices dk for each phase
	for (int m = 1; m <= nphase; m++)
	{
		for (int l = 1; l <= 3; l++)
		{
			for (int k = 1; k <= 3; k++)
			{
				for (int j = 1; j <= 8; j++)
				{
					for (int i = 1; i <= 8; i++)
					{
						dk[m][i][k][j][l][q] = 0.0; //phase-m, ikjl(8x3x8x3), q-th image
					}
				}
			}
		}
	}

	//calculate stiffness matrix dx for each phase
	for (int ijk = 1; ijk <= nphase; ijk++)
	{
		for (int k = 1; k <= 3; k++)
		{
			for (int j = 1; j <= 3; j++)
			{
				for (int i = 1; i <= 3; i++)
				{
					//Matrix multiply to determine value at (x,y,z), multiply by proper weight,
					//and sum into the stiffness matrix dk
					//see MANUAL page 22, equation 20, integrate es*cmod*es to get dk
					for (int mm = 1; mm <= 3; mm++)
					{
						for (int nn = 1; nn <= 3; nn++)
						{
							for (int ii = 1; ii <= 8; ii++)
							{
								for (int jj = 1; jj <= 8; jj++)
								{
									// Define sum over strain matrices and elastic moduli matrix for stiffness matrix
									double sum = 0.0;
									for (int kk = 1; kk <= 6; kk++)
									{
										for (int ll = 1; ll <= 6; ll++)
										{
											sum = sum + es[kk][ii][mm] * cmod[ijk][kk][ll][q] * es[ll][jj][nn];
										}
									}
									dk[ijk][ii][mm][jj][nn][q] = dk[ijk][ii][mm][jj][nn][q] + g[i][j][k] * sum / 216.;
									//g/216 comes from Simpson's rule
								}
							}
						}
					}
				}
			}
		}
	}


	/*
	** Set up vector for linear term, b, and constant term, C, in the elastic energy. see MANUAL equation 10.
	** This is done using the stiffness matrices dk, and the periodic terms delta in the applied strain at the boundary pixels.
	** It is easier to set b up this way than to analytically write out all the terms involved. 
	*/

	//see MANUAL page 19 (electrical explanation) and page 22 (elastic explanation), and table 7 for delta
	//Usually, En=1/2*u*dk*u, but need to be corrected for some cases, En=1/2*(u+delta)*dk*(u+delta)
	//b=delta*dk, C=1/2*delta*dk*delta, they are used to correct for displacements:
	// 1)aging pixels 
	// 2)CSH pixels 
	// 3)boundaries 

	int qtime;

	for (int m3 = 1; m3 <= 3; m3++)
	{
		for (int m = 1; m <= ns; m++) //for all the 1000000 elements, 3 directions
		{
			b[m][m3] = 0.0;
		}
	}
	C = 0.0;

	/*
	** For all cases, the correspondence between 1-8 finite element node labels and 1-27 neighbor labels is (see Table 4 in manual):
	** 1:ib(m,27), 2:ib(m,3),
	** 3:ib(m,2), 4:ib(m,1),
	** 5:ib(m,26), 6:ib(m,19),
	** 7:ib(m,18), 8:ib(m,17).
	*/
	is[1] = 27;
	is[2] = 3;
	is[3] = 2;
	is[4] = 1;
	is[5] = 26;
	is[6] = 19;
	is[7] = 18;
	is[8] = 17;

	//FOR aging pixels, b matrix is set (correct for energy of newly formed pixels, remove pre-history strains)
	for (int k = 1; k <= nz; k++)
	{
		for (int j = 1; j <= ny; j++)
		{
			for (int i = 1; i <= nx; i++)
			{
				int m = nx * ny * (k - 1) + nx * (j - 1) + i;

				//if phase changes
				if (phasechange[m] == 1) //phasechange[m] is global
				{
					//calculate the local strain in that pixel at previous time step pixt[m], at pixt[m]+1, new phase is formed, and current time is q					
					strain(i, j, k, nx, ny, nz, pixt[m]); 

					//get the b value
					//define delta array
					//dxx, dyy, dzz, dxy, dyz, dxz are strain components of the pixel
					//delta is used to correct for nodal displacements, in this case, for new formed pixels
					for (int i8 = 1; i8 <= 8; i8++)
					{
						for (int i3 = 1; i3 <= 3; i3++)
						{
							delta[i8][i3] = 0.0;
						}

						if (i8 == 2)
						{
							delta[i8][1] = delta[i8][1] - dxx; //remove previous displacement (before formation of the phase) from current step
							delta[i8][2] = delta[i8][2] - dxy;
							delta[i8][3] = delta[i8][3] - dxz;
						}
						if (i8 == 4)
						{
							delta[i8][1] = delta[i8][1] - dxy;
							delta[i8][2] = delta[i8][2] - dyy;
							delta[i8][3] = delta[i8][3] - dyz;
						}
						if (i8 == 5)
						{
							delta[i8][1] = delta[i8][1] - dxz;
							delta[i8][2] = delta[i8][2] - dyz;
							delta[i8][3] = delta[i8][3] - dzz;
						}
						if (i8 == 8)
						{
							delta[i8][1] = delta[i8][1] - dxy - dxz;
							delta[i8][2] = delta[i8][2] - dyy - dyz;
							delta[i8][3] = delta[i8][3] - dyz - dzz;
						}
						if (i8 == 6)
						{
							delta[i8][1] = delta[i8][1] - dxx - dxz;
							delta[i8][2] = delta[i8][2] - dxy - dyz;
							delta[i8][3] = delta[i8][3] - dxz - dzz;
						}
						if (i8 == 3)
						{
							delta[i8][1] = delta[i8][1] - dxx - dxy;
							delta[i8][2] = delta[i8][2] - dxy - dyy;
							delta[i8][3] = delta[i8][3] - dxz - dyz;
						}
						if (i8 == 7)
						{
							delta[i8][1] = delta[i8][1] - dxx - dxy - dxz;
							delta[i8][2] = delta[i8][2] - dxy - dyy - dyz;
							delta[i8][3] = delta[i8][3] - dxz - dyz - dzz;
						}
					}

					//define b and C
					for (int nn = 1; nn <= 3; nn++)
					{
						for (int mm = 1; mm <= 8; mm++)
						{
							double sum = 0.0;
							for (int m3 = 1; m3 <= 3; m3++)
							{
								for (int m8 = 1; m8 <= 8; m8++)
								{
									for (qtime = pixt[m] + 1; qtime <= pixt[m] + 1; qtime++)
									{	//Question:???not a loop, only once
										sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]]; //eq(11), global b array
										C = C + 0.5 * delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]] * delta[mm][nn]; //eq(11), global constant C		 
									}
								}
							}
							b[ib[m][is[mm]]][nn] = b[ib[m][is[mm]]][nn] + sum;
						}
					}
				}
			}
		}
	}

	//FOR viscoelastic
	//If all phases are elastic, delete the following for loop
	for (int k = 1; k <= nz; k++)
	{
		for (int j = 1; j <= ny; j++)
		{
			for (int i = 1; i <= nx; i++)
			{
				int m = nx * ny * (k - 1) + nx * (j - 1) + i;

				//if the pixel is CSH
				if (pix[m] == 12) //phase 12 is CSH
				{
					//Adds up all previous time steps 
					for (qtime = pixt[m] + 1; qtime <= q - 1; qtime++) //Question: why from pixt[m]+1 to q-1?
					{
						strain(i, j, k, nx, ny, nz, qtime);

						for (int i8 = 1; i8 <= 8; i8++)
						{
							for (int i3 = 1; i3 <= 3; i3++)
							{
								delta[i8][i3] = 0.0;
							}

							if (i8 == 2)
							{
								delta[i8][1] = delta[i8][1] - dxx;
								delta[i8][2] = delta[i8][2] - dxy;
								delta[i8][3] = delta[i8][3] - dxz;
							}
							if (i8 == 4)
							{
								delta[i8][1] = delta[i8][1] - dxy;
								delta[i8][2] = delta[i8][2] - dyy;
								delta[i8][3] = delta[i8][3] - dyz;
							}
							if (i8 == 5)
							{
								delta[i8][1] = delta[i8][1] - dxz;
								delta[i8][2] = delta[i8][2] - dyz;
								delta[i8][3] = delta[i8][3] - dzz;
							}
							if (i8 == 8)
							{
								delta[i8][1] = delta[i8][1] - dxy - dxz;
								delta[i8][2] = delta[i8][2] - dyy - dyz;
								delta[i8][3] = delta[i8][3] - dyz - dzz;
							}
							if (i8 == 6)
							{
								delta[i8][1] = delta[i8][1] - dxx - dxz;
								delta[i8][2] = delta[i8][2] - dxy - dyz;
								delta[i8][3] = delta[i8][3] - dxz - dzz;
							}
							if (i8 == 3)
							{
								delta[i8][1] = delta[i8][1] - dxx - dxy;
								delta[i8][2] = delta[i8][2] - dxy - dyy;
								delta[i8][3] = delta[i8][3] - dxz - dyz;
							}
							if (i8 == 7)
							{
								delta[i8][1] = delta[i8][1] - dxx - dxy - dxz;
								delta[i8][2] = delta[i8][2] - dxy - dyy - dyz;
								delta[i8][3] = delta[i8][3] - dxz - dyz - dzz;
							}
						}
						//for csh phase, calculate b array and C constant 
						for (int nn = 1; nn <= 3; nn++)
						{
							for (int mm = 1; mm <= 8; mm++)
							{
								double sum = 0.0;
								for (int m3 = 1; m3 <= 3; m3++)
								{
									for (int m8 = 1; m8 <= 8; m8++)
									{
										sum = sum - delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][q - qtime + 1];
										C = C + 0.5 - delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][q - qtime + 1] * delta[mm][nn];											
									}
								}
								b[ib[m][is[mm]]][nn] = b[ib[m][is[mm]]][nn] + sum;
							}
						}
					}

					//If aging has occurred in that pixel, the strain in the pixel when aging occurs should be considered.
					if (phasechange[m] == 1)
					{

						strain(i, j, k, nx, ny, nz, pixt[m]);

						for (int i8 = 1; i8 <= 8; i8++)
						{

							for (int i3 = 1; i3 <= 3; i3++)
							{
								delta[i8][i3] = 0.0;
							}


							if (i8 == 2)
							{
								delta[i8][1] = delta[i8][1] - dxx;
								delta[i8][2] = delta[i8][2] - dxy;
								delta[i8][3] = delta[i8][3] - dxz;
							}
							if (i8 == 4)
							{
								delta[i8][1] = delta[i8][1] - dxy;
								delta[i8][2] = delta[i8][2] - dyy;
								delta[i8][3] = delta[i8][3] - dyz;
							}
							if (i8 == 5)
							{
								delta[i8][1] = delta[i8][1] - dxz;
								delta[i8][2] = delta[i8][2] - dyz;
								delta[i8][3] = delta[i8][3] - dzz;
							}
							if (i8 == 8)
							{
								delta[i8][1] = delta[i8][1] - dxy - dxz;
								delta[i8][2] = delta[i8][2] - dyy - dyz;
								delta[i8][3] = delta[i8][3] - dyz - dzz;
							}
							if (i8 == 6)
							{
								delta[i8][1] = delta[i8][1] - dxx - dxz;
								delta[i8][2] = delta[i8][2] - dxy - dyz;
								delta[i8][3] = delta[i8][3] - dxz - dzz;
							}
							if (i8 == 3)
							{
								delta[i8][1] = delta[i8][1] - dxx - dxy;
								delta[i8][2] = delta[i8][2] - dxy - dyy;
								delta[i8][3] = delta[i8][3] - dxz - dyz;
							}
							if (i8 == 7)
							{
								delta[i8][1] = delta[i8][1] - dxx - dxy - dxz;
								delta[i8][2] = delta[i8][2] - dxy - dyy - dyz;
								delta[i8][3] = delta[i8][3] - dxz - dyz - dzz;
							}
						}

						for (int nn = 1; nn <= 3; nn++)
						{
							for (int mm = 1; mm <= 8; mm++)
							{
								double sum = 0.0;
								for (int m3 = 1; m3 <= 3; m3++)
								{
									for (int m8 = 1; m8 <= 8; m8++)
									{
										for (qtime = 2; qtime <= q - pixt[m]; qtime++)
										{
											sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime];
											C = C + 0.5 * delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime]* delta[mm][nn];								
										}
									}
								}
								b[ib[m][is[mm]]][nn] = b[ib[m][is[mm]]][nn] + sum;
							}
						}
					}
				}
			}
		}
	}

	//For boundary conditions, see MANUAL table 7
	//x=nx face
	for (int i3 = 1; i3 <= 3; i3++)
	{
		for (int i8 = 1; i8 <= 8; i8++)
		{
			delta[i8][i3] = 0.0;
			if (i8 == 2 || i8 == 3 || i8 == 6 || i8 == 7)
			{
				delta[i8][1] = exx * nx;
				delta[i8][2] = exy * nx;
				delta[i8][3] = exz * nx;
			}

		}
	}
	for (int j = 1; j <= ny - 1; j++)
	{
		for (int k = 1; k <= nz - 1; k++)
		{
			int m = nx * ny * (k - 1) + nx * j; //at x=nx face

			for (int nn = 1; nn <= 3; nn++)
			{
				for (int mm = 1; mm <= 8; mm++)
				{
					double sum = 0.0;
					for (int m3 = 1; m3 <= 3; m3++)
					{
						for (int m8 = 1; m8 <= 8; m8++)
						{
							if (phasechange[m] == 1) //if phase change at a pixel at x=nx surface
							{
								for (qtime = pixt[m] + 1; qtime <= q; qtime++) //from the time step this phase appears to current time step
								{
									sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]]; //sum over history since formed
									C = C + 0.5 * delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]] * delta[mm][nn];									
								}
							}

							if (phasechange[m] == 0) //if phase is always the same chemical
							{
								for (qtime = 1; qtime <= q; qtime++)
								{
									sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime]; //sum over history
									C = C + 0.5 * delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime] * delta[mm][nn];										
								}
							}
						}
					}
					b[ib[m][is[mm]]][nn] = b[ib[m][is[mm]]][nn] + sum;
					//b[m][3]: array b with mx3 elements
					//ib[m][28]: 27 neighbors of pixel m, lable of that particular neighbor
					//is[mm]: is[1-8], e.g. is[8]=17, node 8 of pixel = neighor 17's node 1
				}
			}
		}
	}

	//y=ny face
	for (int i3 = 1; i3 <= 3; i3++)
	{
		for (int i8 = 1; i8 <= 8; i8++)
		{
			delta[i8][i3] = 0.0;
			if (i8 == 3 || i8 == 4 || i8 == 7 || i8 == 8)
			{
				delta[i8][1] = exy * ny;
				delta[i8][2] = eyy * ny;
				delta[i8][3] = eyz * ny;
			}
		}
	}
	for (int i = 1; i <= nx - 1; i++)
	{
		for (int k = 1; k <= nz - 1; k++)
		{
			int m = nx * ny * (k - 1) + nx * (ny - 1) + i; //y=ny face
			for (int nn = 1; nn <= 3; nn++)
			{
				for (int mm = 1; mm <= 8; mm++)
				{
					double sum = 0.0;
					for (int m3 = 1; m3 <= 3; m3++)
					{
						for (int m8 = 1; m8 <= 8; m8++)
						{
							if (phasechange[m] == 1)
							{
								for (qtime = pixt[m] + 1; qtime <= q; qtime++)
								{
									sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]]; //sum over history since the phase is formed
									C = C + 0.5 * delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]] * delta[mm][nn];										
								}
							}

							if (phasechange[m] == 0)
							{

								for (qtime = 1; qtime <= q; qtime++)
								{
									sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime]; //sum over history since the first image
									C = C + 0.5 * delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime] * delta[mm][nn];										
								}
							}
						}
					}
					b[ib[m][is[mm]]][nn] = b[ib[m][is[mm]]][nn] + sum;
				}
			}
		}
	}

	//z=nz face
	for (int i3 = 1; i3 <= 3; i3++)
	{
		for (int i8 = 1; i8 <= 8; i8++)
		{
			delta[i8][i3] = 0.0;
			if (i8 == 5 || i8 == 6 || i8 == 7 || i8 == 8)
			{
				delta[i8][1] = exz * nz;
				delta[i8][2] = eyz * nz;
				delta[i8][3] = ezz * nz;
			}
		}
	}
	for (int i = 1; i <= nx - 1; i++)
	{
		for (int j = 1; j <= ny - 1; j++)
		{
			int m = nx * ny * (nz - 1) + nx * (j - 1) + i; //at z=nz face
			for (int nn = 1; nn <= 3; nn++)
			{
				for (int mm = 1; mm <= 8; mm++)
				{
					double sum = 0.0;
					for (int m3 = 1; m3 <= 3; m3++)
					{
						for (int m8 = 1; m8 <= 8; m8++)
						{
							if (phasechange[m] == 1)
							{
								for (qtime = pixt[m] + 1; qtime <= q; qtime++)
								{
									sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]];
									C = C + 0.5 * delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]] * delta[mm][nn];										
								}
							}

							if (phasechange[m] == 0)
							{
								for (qtime = 1; qtime <= q; qtime++)
								{
									sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime];
									C = C + 0.5 * delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime] * delta[mm][nn];										
								}
							}
						}
					}
					b[ib[m][is[mm]]][nn] = b[ib[m][is[mm]]][nn] + sum;
				}
			}
		}
	}
	
	//x=nx, y=ny edge
	for (int i3 = 1; i3 <= 3; i3++)
	{
		for (int i8 = 1; i8 <= 8; i8++)
		{
			delta[i8][i3] = 0.0;
			if (i8 == 2 || i8 == 6)
			{
				delta[i8][1] = exx * nx;
				delta[i8][2] = exy * nx;
				delta[i8][3] = exz * nx;
			}
			if (i8 == 4 || i8 == 8)
			{
				delta[i8][1] = exy * ny;
				delta[i8][2] = eyy * ny;
				delta[i8][3] = eyz * ny;
			}
			if (i8 == 3 || i8 == 7)
			{
				delta[i8][1] = exy * ny + exx * nx;
				delta[i8][2] = eyy * ny + exy * nx;
				delta[i8][3] = eyz * ny + exz * nx;
			}
		}
	}
	for (int k = 1; k <= nz - 1; k++)
	{
		int m = nx * ny * k; //at x=nx, y=ny edge
		for (int nn = 1; nn <= 3; nn++)
		{
			for (int mm = 1; mm <= 8; mm++)
			{
				double sum = 0.0;
				for (int m3 = 1; m3 <= 3; m3++)
				{
					for (int m8 = 1; m8 <= 8; m8++)
					{
						if (phasechange[m] == 1)
						{
							for (qtime = pixt[m] + 1; qtime <= q; qtime++)
							{
								sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]];
								C = C + 0.5 * delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]] * delta[mm][nn];									
							}
						}

						if (phasechange[m] == 0)
						{
							for (qtime = 1; qtime <= q; qtime++)
							{
								sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime];
								C = C + 0.5 * delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime] * delta[mm][nn];									
							}
						}
					}
				}
				b[ib[m][is[mm]]][nn] = b[ib[m][is[mm]]][nn] + sum;
			}
		}
	}
	
	//x=nx, z=nz edge
	for (int i3 = 1; i3 <= 3; i3++)
	{
		for (int i8 = 1; i8 <= 8; i8++)
		{
			delta[i8][i3] = 0.0;
			if (i8 == 2 || i8 == 3)
			{
				delta[i8][1] = exx * nx;
				delta[i8][2] = exy * nx;
				delta[i8][3] = exz * nx;
			}
			if (i8 == 5 || i8 == 8)
			{
				delta[i8][1] = exz * nz;
				delta[i8][2] = eyz * nz;
				delta[i8][3] = ezz * nz;
			}
			if (i8 == 6 || i8 == 7)
			{
				delta[i8][1] = exz * nz + exx * nx;
				delta[i8][2] = eyz * nz + exy * nx;
				delta[i8][3] = ezz * nz + exz * nx;
			}
		}
	}
	for (int j = 1; j <= ny - 1; j++)
	{
		int m = nx * ny * (nz - 1) + nx * (j - 1) + nx; //at x=nx, z=nz edge
		for (int nn = 1; nn <= 3; nn++)
		{
			for (int mm = 1; mm <= 8; mm++)
			{
				double sum = 0.0;
				for (int m3 = 1; m3 <= 3; m3++)
				{
					for (int m8 = 0; m8 <= 8; m8++)
					{
						if (phasechange[m] == 1)
						{
							for (qtime = pixt[m] + 1; qtime <= q; qtime++)
							{
								sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]];
								C = C + 0.5*delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]] * delta[mm][nn];									
							}
						}

						if (phasechange[m] == 0)
						{
							for (qtime = 1; qtime <= q; qtime++)
							{
								sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime];
								C = C + 0.5*delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime] * delta[mm][nn];									
							}
						}
					}
				}
				b[ib[m][is[mm]]][nn] = b[ib[m][is[mm]]][nn] + sum;
			}
		}
	}
	
	//y=ny z=nz edge
	for (int i3 = 1; i3 <= 3; i3++)
	{
		for (int i8 = 1; i8 <= 8; i8++)
		{
			delta[i8][i3] = 0.0;
			if (i8 == 5 || i8 == 6)
			{
				delta[i8][1] = exz * nz;
				delta[i8][2] = eyz * nz;
				delta[i8][3] = ezz * nz;
			}
			if (i8 == 3 || i8 == 4)
			{
				delta[i8][1] = exy * ny;
				delta[i8][2] = eyy * ny;
				delta[i8][3] = eyz * ny;
			}
			if (i8 == 7 || i8 == 8)
			{
				delta[i8][1] = exy * ny + exz * nz;
				delta[i8][2] = eyy * ny + eyz * nz;
				delta[i8][3] = eyz * ny + ezz * nz;
			}
		}
	}
	for (int i = 1; i <= nx - 1; i++)
	{
		int m = nx * ny * (nz - 1) + nx * (ny - 1) + i; //at y=ny, z=nz edge
		for (int nn = 1; nn <= 3; nn++)
		{
			for (int mm = 1; mm <= 8; mm++)
			{
				double sum = 0.0;
				for (int m3 = 1; m3 <= 3; m3++)
				{
					for (int m8 = 1; m8 <= 8; m8++)
					{
						if (phasechange[m] == 1)
						{
							for (qtime = pixt[m] + 1; qtime <= q; qtime++)
							{
								sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]];
								C = C + 0.5 * delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]] * delta[mm][nn];									
							}
						}

						if (phasechange[m] == 0)
						{
							for (qtime = 1; qtime <= q; qtime++)
							{
								sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime];
								C = C + 0.5*delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime] * delta[mm][nn];									
							}
						}
					}
				}
				b[ib[m][is[mm]]][nn] = b[ib[m][is[mm]]][nn] + sum;
			}
		}
	}

	//x=nx y=ny z=nz corner
	for (int i3 = 1; i3 <= 3; i3++)
	{
		for (int i8 = 1; i8 <= 8; i8++)
		{
			delta[i8][i3] = 0.0;
			if (i8 == 2)
			{
				delta[i8][1] = exx * nx;
				delta[i8][2] = exy * nx;
				delta[i8][3] = exz * nx;
			}
			if (i8 == 4)
			{
				delta[i8][1] = exy * ny;
				delta[i8][2] = eyy * ny;
				delta[i8][3] = eyz * ny;
			}
			if (i8 == 5)
			{
				delta[i8][1] = exz * nz;
				delta[i8][2] = eyz * nz;
				delta[i8][3] = ezz * nz;
			}
			if (i8 == 8)
			{
				delta[i8][1] = exy * ny + exz * nz;
				delta[i8][2] = eyy * ny + eyz * nz;
				delta[i8][3] = eyz * ny + ezz * nz;
			}
			if (i8 == 6)
			{
				delta[i8][1] = exx * nx + exz * nz;
				delta[i8][2] = exy * nx + eyz * nz;
				delta[i8][3] = exz * nx + ezz * nz;
			}
			if (i8 == 3)
			{
				delta[i8][1] = exx * nx + exy * ny;
				delta[i8][2] = exy * nx + eyy * ny;
				delta[i8][3] = exz * nx + eyz * ny;
			}
			if (i8 == 7)
			{
				delta[i8][1] = exx * nx + exy * ny + exz * nz;
				delta[i8][2] = exy * nx + eyy * ny + eyz * nz;
				delta[i8][3] = exz * nx + eyz * ny + ezz * nz;
			}
		}
	}
	int m = nx * ny * nz; //at x=nx, y=ny, z=nz corner
	for (int nn = 1; nn <= 3; nn++)
	{
		for (int mm = 1; mm <= 8; mm++)
		{
			double sum = 0.0;
			for (int m3 = 1; m3 <= 3; m3++)
			{
				for (int m8 = 1; m8 <= 8; m8++)
				{
					if (phasechange[m] == 1)
					{
						for (qtime = pixt[m] + 1; qtime <= q; qtime++)
						{
							sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]];
							C = C + 0.5 * delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime - pixt[m]] * delta[mm][nn];								
						}
					}

					if (phasechange[m] == 0)
					{
						for (qtime = 1; qtime <= q; qtime++)
						{
							sum = sum + delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime];
							C = C + 0.5 * delta[m8][m3] * dk[pix[m]][m8][m3][mm][nn][qtime] * delta[mm][nn];								
						}
					}
				}
			}
			b[ib[m][is[mm]]][nn] = b[ib[m][is[mm]]][nn] + sum;
		}
	}
}



//Subroutine that computes the total energy, utot, and the gradient, gb = Au + b
//minimize gradient of energy gb close to zero, see MANUAL equation 12
void energy(int nx, int ny, int nz, int ns, double& utot, int q)
{	
	std::cout << "      Call energy() to compute total energy and gradient..." << std::endl;
	
	for (int m3 = 1; m3 <= 3; m3++)
	{
		for (int m = 1; m <= ns; m++)
		{
			gb[m][m3] = 0.0;
		}
	}

	/*
	** Do global matrix multiply via small stiffness matrices, gb = A * u
	** The long statement below correctly brings in all the terms from
	** the global matrix A using only the small stiffness matrices.
	*/
	//see MANUAL page 113
	for (int j = 1; j <= 3; j++)
	{
		for (int n = 1; n <= 3; n++)
		{
			for (int m = 1; m <= ns; m++)
			{
				for (int qtime = 1; qtime <= 1; qtime++)
				{
					gb[m][j] = gb[m][j] + u[ib[m][1]][q + 1 - qtime][n] *
						(dk[pix[ib[m][27]]][1][j][4][n][qtime] + dk[pix[ib[m][7]]][2][j][3][n][qtime] +
						dk[pix[ib[m][25]]][5][j][8][n][qtime] + dk[pix[ib[m][15]]][6][j][7][n][qtime]) +
						u[ib[m][2]][q + 1 - qtime][n] * (dk[pix[ib[m][27]]][1][j][3][n][qtime] + dk[pix[ib[m][25]
						]][5][j][7][n][qtime]) + u[ib[m][3]
						][q + 1 - qtime][n] * (dk[pix[ib[m][27]]][1][j][2][n][qtime] + dk[pix[ib[m][5]
						]][4][j][3][n][qtime] + dk[pix[ib[m][13]
						]][8][j][7][n][qtime] + dk[pix[ib[m][25]
						]][5][j][6][n][qtime]) + u[ib[m][4]
						][q + 1 - qtime][n] * (dk[pix[ib[m][5]
						]][4][j][2][n][qtime] + dk[pix[ib[m][13]
						]][8][j][6][n][qtime]) + u[ib[m][5]
						][q + 1 - qtime][n] * (dk[pix[ib[m][6]
						]][3][j][2][n][qtime] + dk[pix[ib[m][5]
						]][4][j][1][n][qtime] + dk[pix[ib[m][14]
						]][7][j][6][n][qtime] + dk[pix[ib[m][13]
						]][8][j][5][n][qtime]) + u[ib[m][6]
						][q + 1 - qtime][n] * (dk[pix[ib[m][6]
						]][3][j][1][n][qtime] + dk[pix[ib[m][14]
						]][7][j][5][n][qtime]) + u[ib[m][7]][q + 1 - qtime][n] * (dk[pix[ib[m][6]
						]][3][j][4][n][qtime] + dk[pix[ib[m][7]
						]][2][j][1][n][qtime] + dk[pix[ib[m][14]
						]][7][j][8][n][qtime] + dk[pix[ib[m][15]
						]][6][j][5][n][qtime]) + u[ib[m][8]
						][q + 1 - qtime][n] * (dk[pix[ib[m][7]
						]][2][j][4][n][qtime] + dk[pix[ib[m][15]
						]][6][j][8][n][qtime]) + u[ib[m][9]
						][q + 1 - qtime][n] * (dk[pix[ib[m][25]
						]][5][j][4][n][qtime] + dk[pix[ib[m][15]
						]][6][j][3][n][qtime]) + u[ib[m][10]
						][q + 1 - qtime][n] * (dk[pix[ib[m][25]
						]][5][j][3][n][qtime]) + u[ib[m][11]
						][q + 1 - qtime][n] * (dk[pix[ib[m][13]
						]][8][j][3][n][qtime] + dk[pix[ib[m][25]
						]][5][j][2][n][qtime]) + u[ib[m][12]
						][q + 1 - qtime][n] * (dk[pix[ib[m][13]
						]][8][j][2][n][qtime]) + u[ib[m][13]
						][q + 1 - qtime][n] * (dk[pix[ib[m][13]
						]][8][j][1][n][qtime] + dk[pix[ib[m][14]
						]][7][j][2][n][qtime]) + u[ib[m][14]
						][q + 1 - qtime][n] * (dk[pix[ib[m][14]
						]][7][j][1][n][qtime]) + u[ib[m][15]
						][q + 1 - qtime][n] * (dk[pix[ib[m][14]
						]][7][j][4][n][qtime] + dk[pix[ib[m][15]
						]][6][j][1][n][qtime]) + u[ib[m][16]
						][q + 1 - qtime][n] * (dk[pix[ib[m][15]
						]][6][j][4][n][qtime]) + u[ib[m][17]
						][q + 1 - qtime][n] * (dk[pix[ib[m][27]
						]][1][j][8][n][qtime] + dk[pix[ib[m][7]
						]][2][j][7][n][qtime]) + u[ib[m][18]
						][q + 1 - qtime][n] * (dk[pix[ib[m][27]
						]][1][j][7][n][qtime]) + u[ib[m][19]
						][q + 1 - qtime][n] * (dk[pix[ib[m][27]
						]][1][j][6][n][qtime] + dk[pix[ib[m][5]
						]][4][j][7][n][qtime]) + u[ib[m][20]
						][q + 1 - qtime][n] * (dk[pix[ib[m][5]
						]][4][j][6][n][qtime]) + u[ib[m][21]
						][q + 1 - qtime][n] * (dk[pix[ib[m][5]
						]][4][j][5][n][qtime] + dk[pix[ib[m][6]
						]][3][j][6][n][qtime]) + u[ib[m][22]
						][q + 1 - qtime][n] * (dk[pix[ib[m][6]
						]][3][j][5][n][qtime]) + u[ib[m][23]
						][q + 1 - qtime][n] * (dk[pix[ib[m][6]
						]][3][j][8][n][qtime] + dk[pix[ib[m][7]
						]][2][j][5][n][qtime]) + u[ib[m][24]
						][q + 1 - qtime][n] * (dk[pix[ib[m][7]
						]][2][j][8][n][qtime]) + u[ib[m][25]
						][q + 1 - qtime][n] * (dk[pix[ib[m][14]
						]][7][j][3][n][qtime] + dk[pix[ib[m][13]
						]][8][j][4][n][qtime] + dk[pix[ib[m][15]
						]][6][j][2][n][qtime] + dk[pix[ib[m][25]
						]][5][j][1][n][qtime]) + u[ib[m][26]
						][q + 1 - qtime][n] * (dk[pix[ib[m][6]
						]][3][j][7][n][qtime] + dk[pix[ib[m][5]
						]][4][j][8][n][qtime] + dk[pix[ib[m][27]
						]][1][j][5][n][qtime] + dk[pix[ib[m][7]
						]][2][j][6][n][qtime]) + u[ib[m][27]
						][q + 1 - qtime][n] * (dk[pix[ib[m][27]
						]][1][j][1][n][qtime] + dk[pix[ib[m][7]
						]][2][j][2][n][qtime] + dk[pix[ib[m][6]
						]][3][j][3][n][qtime] + dk[pix[ib[m][5]
						]][4][j][4][n][qtime] + dk[pix[ib[m][25]
						]][5][j][5][n][qtime] + dk[pix[ib[m][15]
						]][6][j][6][n][qtime] + dk[pix[ib[m][14]
						]][7][j][7][n][qtime] + dk[pix[ib[m][13]
						]][8][j][8][n][qtime]);
				}
			}
		}
	}

	utot = C;

	for (int m3 = 1; m3 <= 3; m3++)
	{
		for (int m = 1; m <= ns; m++)
		{
			if (phasechange[m] == 1)
			{
				utot = utot + 0.5 * (u[m][q][m3] - u[m][pixt[m]][m3])
					* gb[m][m3] + b[m][m3] * (u[m][q][m3] - u[m][pixt[m]][m3]); //utot = 1/2uAu + bu, gb = Au
			}

			if (phasechange[m] == 0)
			{
				utot = utot + 0.5 * u[m][q][m3]
					* gb[m][m3] + b[m][m3] * u[m][q][m3];
			}

			gb[m][m3] = gb[m][m3] + b[m][m3];
		}
	}
}


//Subroutine that use conjugate gradient relaxation method to minimize energy
void dembx(int ns, int &Lstep, double &gg, double(&dk)[101][9][4][9][4][51], double &gtest, int &ldemb, int &kkk, int q)
{
	std::cout << "      Call dembx() to use conjugate gradient relaxation..." << std::endl;
	/*
	** Initialize the conjugate direction vector on first call to dembx only
	** For calls to dembx after the first, we want to continue using the
	** value of h determined in the previous call. Of course, if npoints is
	** greater than 1, this initialization step will be run for every new
	** microstructure used, as kkk is reset to 1 every time the counter micro
	** is increased.
	*/
	double  gamma;
	double lambda;

	int qtime;
	if (kkk == 1)
	{
		for (int m3 = 1; m3 <= 3; m3++)
		{
			for (int m = 1; m <= ns; m++)
			{
				h[m][m3] = gb[m][m3];
			}
		}
	}

	// Lstep counts the number of conjugate gradient steps taken in each call to dembx
	Lstep = 0;

	for (int ijk = 1; ijk <= ldemb; ijk++)
	{
		Lstep = Lstep + 1;

		for (int m3 = 1; m3 <= 3; m3++)
		{
			for (int m = 1; m <= ns; m++)
			{
				Ah[m][m3] = 0.0;
			}
		}
		/*
		** Do global matrix multiply via small stiffness matrices, Ah = A * h
		** The long statement below correctly brings in all the terms from
		** the global matrix A using only the small stiffness matrices dk.
		*/
		//see MANUAL page 115
		for (int j = 1; j <= 3; j++)
		{
			for (int n = 1; n <= 3; n++)
			{
				for (int m = 1; m <= ns; m++)
				{
					qtime = 1;

					Ah[m][j] = Ah[m][j] + h[ib[m][1]][n] * (dk[pix[ib[m][27]]][1][j][4][n][qtime]
						+ dk[pix[ib[m][7]]][2][j][3][n][qtime]
						+ dk[pix[ib[m][25]]][5][j][8][n][qtime] + dk[pix[ib[m][15]]][6][j][7][n][qtime]) +
						h[ib[m][2]][n] * (dk[pix[ib[m][27]]][1][j][3][n][qtime]
						+ dk[pix[ib[m][25]]][5][j][7][n][qtime]) +
						h[ib[m][3]][n] * (dk[pix[ib[m][27]]][1][j][2][n][qtime] + dk[pix[ib[m][5]]][4][j][3][n][qtime] +
						dk[pix[ib[m][13]]][8][j][7][n][qtime] + dk[pix[ib[m][25]]][5][j][6][n][qtime]) +
						h[ib[m][4]][n] * (dk[pix[ib[m][5]]][4][j][2][n][qtime]
						+ dk[pix[ib[m][13]]][8][j][6][n][qtime]) +
						h[ib[m][5]][n] * (dk[pix[ib[m][6]]][3][j][2][n][qtime] + dk[pix[ib[m][5]]][4][j][1][n][qtime] +
						dk[pix[ib[m][14]]][7][j][6][n][qtime] + dk[pix[ib[m][13]]][8][j][5][n][qtime]) +
						h[ib[m][6]][n] * (dk[pix[ib[m][6]]][3][j][1][n][qtime]
						+ dk[pix[ib[m][14]]][7][j][5][n][qtime]) +
						h[ib[m][7]][n] * (dk[pix[ib[m][6]]][3][j][4][n][qtime] + dk[pix[ib[m][7]]][2][j][1][n][qtime] +
						dk[pix[ib[m][14]]][7][j][8][n][qtime] + dk[pix[ib[m][15]]][6][j][5][n][qtime]) +
						h[ib[m][8]][n] * (dk[pix[ib[m][7]]][2][j][4][n][qtime]
						+ dk[pix[ib[m][15]]][6][j][8][n][qtime]) +
						h[ib[m][9]][n] * (dk[pix[ib[m][25]]][5][j][4][n][qtime]
						+ dk[pix[ib[m][15]]][6][j][3][n][qtime]) +
						h[ib[m][10]][n] * (dk[pix[ib[m][25]]][5][j][3][n][qtime]) +
						h[ib[m][11]][n] * (dk[pix[ib[m][13]]][8][j][3][n][qtime]
						+ dk[pix[ib[m][25]]][5][j][2][n][qtime]) +
						h[ib[m][12]][n] * (dk[pix[ib[m][13]]][8][j][2][n][qtime]) +
						h[ib[m][13]][n] * (dk[pix[ib[m][13]]][8][j][1][n][qtime]
						+ dk[pix[ib[m][14]]][7][j][2][n][qtime]) +
						h[ib[m][14]][n] * (dk[pix[ib[m][14]]][7][j][1][n][qtime]) +
						h[ib[m][15]][n] * (dk[pix[ib[m][14]]][7][j][4][n][qtime]
						+ dk[pix[ib[m][15]]][6][j][1][n][qtime]) +
						h[ib[m][16]][n] * (dk[pix[ib[m][15]]][6][j][4][n][qtime]) +
						h[ib[m][17]][n] * (dk[pix[ib[m][27]]][1][j][8][n][qtime]
						+ dk[pix[ib[m][7]]][2][j][7][n][qtime]) +
						h[ib[m][18]][n] * (dk[pix[ib[m][27]]][1][j][7][n][qtime]) +
						h[ib[m][19]][n] * (dk[pix[ib[m][27]]][1][j][6][n][qtime]
						+ dk[pix[ib[m][5]]][4][j][7][n][qtime]) +
						h[ib[m][20]][n] * (dk[pix[ib[m][5]]][4][j][6][n][qtime]) +
						h[ib[m][21]][n] * (dk[pix[ib[m][5]]][4][j][5][n][qtime]
						+ dk[pix[ib[m][6]]][3][j][6][n][qtime]) +
						h[ib[m][22]][n] * (dk[pix[ib[m][6]]][3][j][5][n][qtime]) +
						h[ib[m][23]][n] * (dk[pix[ib[m][6]]][3][j][8][n][qtime]
						+ dk[pix[ib[m][7]]][2][j][5][n][qtime]) +
						h[ib[m][24]][n] * (dk[pix[ib[m][7]]][2][j][8][n][qtime]) +
						h[ib[m][25]][n] * (dk[pix[ib[m][14]]][7][j][3][n][qtime]
						+ dk[pix[ib[m][13]]][8][j][4][n][qtime] +
						dk[pix[ib[m][15]]][6][j][2][n][qtime] + dk[pix[ib[m][25]]][5][j][1][n][qtime]) +
						h[ib[m][26]][n] * (dk[pix[ib[m][6]]][3][j][7][n][qtime]
						+ dk[pix[ib[m][5]]][4][j][8][n][qtime] +
						dk[pix[ib[m][27]]][1][j][5][n][qtime] + dk[pix[ib[m][7]]][2][j][6][n][qtime]) +
						h[ib[m][27]][n] * (dk[pix[ib[m][27]]][1][j][1][n][qtime]
						+ dk[pix[ib[m][7]]][2][j][2][n][qtime] +
						dk[pix[ib[m][6]]][3][j][3][n][qtime] + dk[pix[ib[m][5]]][4][j][4][n][qtime]
						+ dk[pix[ib[m][25]]][5][j][5][n][qtime] +
						dk[pix[ib[m][15]]][6][j][6][n][qtime] + dk[pix[ib[m][14]]][7][j][7][n][qtime] +
						dk[pix[ib[m][13]]][8][j][8][n][qtime]);

				}
			}
		}

		double hAh = 0.0;
		for (int m3 = 1; m3 <= 3; m3++)
		{
			for (int m = 1; m <= ns; m++)
			{
				hAh = hAh + h[m][m3]
					* Ah[m][m3];
			}
		}

		lambda = gg / hAh;
		for (int m3 = 1; m3 <= 3; m3++)
		{
			for (int m = 1; m <= ns; m++)
			{
				u[m][q][m3] = u[m][q][m3] - lambda*h[m][m3];					
				gb[m][m3] = gb[m][m3] - lambda*Ah[m][m3];
			}
		}

		double gglast = gg;
		gg = 0.0;
		for (int m3 = 1; m3 <= 3; m3++)
		{
			for (int m = 1; m <= ns; m++)
			{
				gg = gg + gb[m][m3] * gb[m][m3];
			}
		}
		if (gg<gtest)
		{
			break;
		}

		gamma = gg / gglast;
		for (int m3 = 1; m3 <= 3; m3++)
		{
			for (int m = 1; m <= ns; m++)
			{
				h[m][m3] = gb[m][m3] + gamma*h[m][m3];					
			}
		}
	}
}



// Subroutine that computes the six average stresses and six average strains.
void stress(int nx, int ny, int nz, int ns, int q)
{
	std::cout << "      Call stress() to compute average stresses and strains..." << std::endl;
	
	double dndx[9], dndy[9], dndz[9], es[7][9][4], uu[9][51][4];

	// set up single element strain matrix
	// dndx, dndy, and dndz are the components of the average strain matrix in a pixel
	dndx[1] = -0.25;
	dndx[2] = 0.25;
	dndx[3] = 0.25;
	dndx[4] = -0.25;
	dndx[5] = -0.25;
	dndx[6] = 0.25;
	dndx[7] = 0.25;
	dndx[8] = -0.25;
	dndy[1] = -0.25;
	dndy[2] = -0.25;
	dndy[3] = 0.25;
	dndy[4] = 0.25;
	dndy[5] = -0.25;
	dndy[6] = -0.25;
	dndy[7] = 0.25;
	dndy[8] = 0.25;
	dndz[1] = -0.25;
	dndz[2] = -0.25;
	dndz[3] = -0.25;
	dndz[4] = -0.25;
	dndz[5] = 0.25;
	dndz[6] = 0.25;
	dndz[7] = 0.25;
	dndz[8] = 0.25;
	/*
	** Build averaged strain matrix, follows code in femat, but for average
	** strain over the pixel, not the strain at a point.
	*/
	for (int n1 = 1; n1 <= 6; n1++)
	{
		for (int n2 = 1; n2 <= 8; n2++)
		{
			for (int n3 = 1; n3 <= 3; n3++)
			{
				es[n1][n2][n3] = 0.0;
			}
		}
	}
	for (int n = 1; n <= 8; n++)
	{
		es[1][n][1] = dndx[n];
		es[2][n][2] = dndy[n];
		es[3][n][3] = dndz[n];
		es[4][n][1] = dndz[n];
		es[4][n][3] = dndx[n];
		es[5][n][2] = dndz[n];
		es[5][n][3] = dndy[n];
		es[6][n][1] = dndy[n];
		es[6][n][2] = dndx[n];
	}

	// Compute components of the average stress and strain tensors in each pixel
	strxx = 0.0;
	stryy = 0.0;
	strzz = 0.0;
	strxz = 0.0;
	stryz = 0.0;
	strxy = 0.0;
	sxx = 0.0;
	syy = 0.0;
	szz = 0.0;
	sxz = 0.0;
	syz = 0.0;
	sxy = 0.0;

	int set;

	for (int k = 1; k <= nz; k++)
	{
		for (int j = 1; j <= ny; j++)
		{
			for (int i = 1; i <= nx; i++)
			{
				int m = nx * ny * (k - 1) + nx * (j - 1) + i;
				str11 = 0.0;
				str22 = 0.0;
				str33 = 0.0;
				str13 = 0.0;
				str23 = 0.0;
				str12 = 0.0;

				if (phasechange[m] == 1)
				{
					if (pix[m] == 12) set = pixt[m] + 1; //set = the time step a CSH pixel is formed
					else set = q;

					for (int qtime = set; qtime <= q; qtime++)
					{
						//load in elements of 8-vector using pd. bd. conds.
						for (int mm = 1; mm <= 3; mm++)
						{
							uu[1][qtime][mm] = u[m][qtime][mm] - u[m][pixt[m]][mm];
							uu[2][qtime][mm] = u[ib[m][3]][qtime][mm] - u[ib[m][3]][pixt[m]][mm];								
							uu[3][qtime][mm] = u[ib[m][2]][qtime][mm] - u[ib[m][2]][pixt[m]][mm];								
							uu[4][qtime][mm] = u[ib[m][1]][qtime][mm] - u[ib[m][1]][pixt[m]][mm];								
							uu[5][qtime][mm] = u[ib[m][26]][qtime][mm] - u[ib[m][26]][pixt[m]][mm];								
							uu[6][qtime][mm] = u[ib[m][19]][qtime][mm] - u[ib[m][19]][pixt[m]][mm];								
							uu[7][qtime][mm] = u[ib[m][18]][qtime][mm] - u[ib[m][18]][pixt[m]][mm];								
							uu[8][qtime][mm] = u[ib[m][17]][qtime][mm] - u[ib[m][17]][pixt[m]][mm];						
						}

						// local stresses and strains in a pixel
						for (int n3 = 1; n3 <= 3; n3++)
						{
							for (int n8 = 1; n8 <= 8; n8++)
							{
								for (int n = 1; n <= 6; n++)
								{
									str11 = str11 + cmod[pix[m]][1][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										  										  
									str22 = str22 + cmod[pix[m]][2][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];																				  
									str33 = str33 + cmod[pix[m]][3][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];																				  
									str13 = str13 + cmod[pix[m]][4][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										  										  
									str23 = str23 + cmod[pix[m]][5][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];									  
									str12 = str12 + cmod[pix[m]][6][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										  										  
								}
							}
						}
					}
				}

				if (phasechange[m] == 0)
				{
					if (pix[m] == 12) set = 1;
					else set = q;

					for (int qtime = set; qtime <= q; qtime++)
					{
						//load in elements of 8-vector using pd. bd. conds.
						for (int mm = 1; mm <= 3; mm++)
						{
							uu[1][qtime][mm] = u[m][qtime][mm];
							uu[2][qtime][mm] = u[ib[m][3]][qtime][mm];								
							uu[3][qtime][mm] = u[ib[m][2]][qtime][mm];								
							uu[4][qtime][mm] = u[ib[m][1]][qtime][mm];								
							uu[5][qtime][mm] = u[ib[m][26]][qtime][mm];								
							uu[6][qtime][mm] = u[ib[m][19]][qtime][mm];								
							uu[7][qtime][mm] = u[ib[m][18]][qtime][mm];								
							uu[8][qtime][mm] = u[ib[m][17]][qtime][mm];								
						}

						//correction for period boundary conditions
						if (i == nx)
						{
							uu[2][qtime][1] = uu[2][qtime][1] + exx * nx;								
							uu[2][qtime][2] = uu[2][qtime][2] + exy * nx;								
							uu[2][qtime][3] = uu[2][qtime][3] + exz * nx;								
							uu[3][qtime][1] = uu[3][qtime][1] + exx * nx;								
							uu[3][qtime][2] = uu[3][qtime][2] + exy * nx;								
							uu[3][qtime][3] = uu[3][qtime][3] + exz * nx;								
							uu[6][qtime][1] = uu[6][qtime][1] + exx * nx;								
							uu[6][qtime][2] = uu[6][qtime][2] + exy * nx;								
							uu[6][qtime][3] = uu[6][qtime][3] + exz * nx;								
							uu[7][qtime][1] = uu[7][qtime][1] + exx * nx;								
							uu[7][qtime][2] = uu[7][qtime][2] + exy * nx;								
							uu[7][qtime][3] = uu[7][qtime][3] + exz * nx;								
						}
						if (j == ny)
						{
							uu[3][qtime][1] = uu[3][qtime][1] + exy * ny;								
							uu[3][qtime][2] = uu[3][qtime][2] + eyy * ny;								
							uu[3][qtime][3] = uu[3][qtime][3] + eyz * ny;								
							uu[4][qtime][1] = uu[4][qtime][1] + exy * ny;								
							uu[4][qtime][2] = uu[4][qtime][2] + eyy * ny;								
							uu[4][qtime][3] = uu[4][qtime][3] + eyz * ny;								
							uu[7][qtime][1] = uu[7][qtime][1] + exy * ny;								
							uu[7][qtime][2] = uu[7][qtime][2] + eyy * ny;								
							uu[7][qtime][3] = uu[7][qtime][3] + eyz * ny;								
							uu[8][qtime][1] = uu[8][qtime][1] + exy * ny;								
							uu[8][qtime][2] = uu[8][qtime][2] + eyy * ny;								
							uu[8][qtime][3] = uu[8][qtime][3] + eyz * ny;								
						}
						if (k == nz)
						{
							uu[5][qtime][1] = uu[5][qtime][1] + exz * nz;								
							uu[5][qtime][2] = uu[5][qtime][2] + eyz * nz;								
							uu[5][qtime][3] = uu[5][qtime][3] + ezz * nz;								
							uu[6][qtime][1] = uu[6][qtime][1] + exz * nz;								
							uu[6][qtime][2] = uu[6][qtime][2] + eyz * nz;								
							uu[6][qtime][3] = uu[6][qtime][3] + ezz * nz;								
							uu[7][qtime][1] = uu[7][qtime][1] + exz * nz;								
							uu[7][qtime][2] = uu[7][qtime][2] + eyz * nz;								
							uu[7][qtime][3] = uu[7][qtime][3] + ezz * nz;								
							uu[8][qtime][1] = uu[8][qtime][1] + exz * nz;								
							uu[8][qtime][2] = uu[8][qtime][2] + eyz * nz;								
							uu[8][qtime][3] = uu[8][qtime][3] + ezz * nz;								
						}

						//local stresses and strains in a pixel
						for (int n3 = 1; n3 <= 3; n3++)
						{
							for (int n8 = 1; n8 <= 8; n8++)
							{
								for (int n = 1; n <= 6; n++)
								{
									str11 = str11 + cmod[pix[m]][1][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										 										 
									str22 = str22 + cmod[pix[m]][2][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										 										  
									str33 = str33 + cmod[pix[m]][3][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										 										 
									str13 = str13 + cmod[pix[m]][4][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										 										 
									str23 = str23 + cmod[pix[m]][5][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										 										 
									str12 = str12 + cmod[pix[m]][6][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										 										 
								}
							}
						}
					}
				}

				//sum local strains and stresses into global values
				strxx = strxx + str11;
				stryy = stryy + str22;
				strzz = strzz + str33;
				strxz = strxz + str13;
				stryz = stryz + str23;
				strxy = strxy + str12;
				
			}
		}
	}

	// Volume average of global stresses and strains
	strxx = strxx / (1.0* ns);
	stryy = stryy / (1.0* ns);
	strzz = strzz / (1.0* ns);
	strxz = strxz / (1.0* ns);
	stryz = stryz / (1.0* ns);
	strxy = strxy / (1.0* ns);
	sxx = exx;
	syy = eyy;
	szz = ezz;
	sxz = exz;
	syz = eyz;
	sxy = exy;
}


// Added by Yige, with output of CSH pixel stresses
// Subroutine that computes the six average stresses and six average strains.
void stress2(int nx, int ny, int nz, int ns, int q)
{
	std::cout << "      Call stress2() to compute average stresses and strains..." << std::endl;
	
	double dndx[9], dndy[9], dndz[9], es[7][9][4], uu[9][51][4];
	
	//YZ added to store stresses of CSH
	int counter = 0;

	// set up single element strain matrix
	// dndx, dndy, and dndz are the components of the average strain matrix in a pixel
	dndx[1] = -0.25;
	dndx[2] = 0.25;
	dndx[3] = 0.25;
	dndx[4] = -0.25;
	dndx[5] = -0.25;
	dndx[6] = 0.25;
	dndx[7] = 0.25;
	dndx[8] = -0.25;
	dndy[1] = -0.25;
	dndy[2] = -0.25;
	dndy[3] = 0.25;
	dndy[4] = 0.25;
	dndy[5] = -0.25;
	dndy[6] = -0.25;
	dndy[7] = 0.25;
	dndy[8] = 0.25;
	dndz[1] = -0.25;
	dndz[2] = -0.25;
	dndz[3] = -0.25;
	dndz[4] = -0.25;
	dndz[5] = 0.25;
	dndz[6] = 0.25;
	dndz[7] = 0.25;
	dndz[8] = 0.25;
	/*
	** Build averaged strain matrix, follows code in femat, but for average
	** strain over the pixel, not the strain at a point.
	*/
	for (int n1 = 1; n1 <= 6; n1++)
	{
		for (int n2 = 1; n2 <= 8; n2++)
		{
			for (int n3 = 1; n3 <= 3; n3++)
			{
				es[n1][n2][n3] = 0.0;
			}
		}
	}
	for (int n = 1; n <= 8; n++)
	{
		es[1][n][1] = dndx[n];
		es[2][n][2] = dndy[n];
		es[3][n][3] = dndz[n];
		es[4][n][1] = dndz[n];
		es[4][n][3] = dndx[n];
		es[5][n][2] = dndz[n];
		es[5][n][3] = dndy[n];
		es[6][n][1] = dndy[n];
		es[6][n][2] = dndx[n];
	}

	// Compute components of the average stress and strain tensors in each pixel
	strxx = 0.0;
	stryy = 0.0;
	strzz = 0.0;
	strxz = 0.0;
	stryz = 0.0;
	strxy = 0.0;
	sxx = 0.0;
	syy = 0.0;
	szz = 0.0;
	sxz = 0.0;
	syz = 0.0;
	sxy = 0.0;

	int set;

	for (int k = 1; k <= nz; k++)
	{
		for (int j = 1; j <= ny; j++)
		{
			for (int i = 1; i <= nx; i++)
			{
				int m = nx * ny * (k - 1) + nx * (j - 1) + i;
				str11 = 0.0;
				str22 = 0.0;
				str33 = 0.0;
				str13 = 0.0;
				str23 = 0.0;
				str12 = 0.0;

				if (phasechange[m] == 1)
				{
					if (pix[m] == 12) set = pixt[m] + 1; //set = the time step a CSH pixel is formed
					else set = q;

					for (int qtime = set; qtime <= q; qtime++)
					{
						//load in elements of 8-vector using pd. bd. conds.
						for (int mm = 1; mm <= 3; mm++)
						{
							uu[1][qtime][mm] = u[m][qtime][mm] - u[m][pixt[m]][mm];
							uu[2][qtime][mm] = u[ib[m][3]][qtime][mm] - u[ib[m][3]][pixt[m]][mm];								
							uu[3][qtime][mm] = u[ib[m][2]][qtime][mm] - u[ib[m][2]][pixt[m]][mm];								
							uu[4][qtime][mm] = u[ib[m][1]][qtime][mm] - u[ib[m][1]][pixt[m]][mm];								
							uu[5][qtime][mm] = u[ib[m][26]][qtime][mm] - u[ib[m][26]][pixt[m]][mm];								
							uu[6][qtime][mm] = u[ib[m][19]][qtime][mm] - u[ib[m][19]][pixt[m]][mm];								
							uu[7][qtime][mm] = u[ib[m][18]][qtime][mm] - u[ib[m][18]][pixt[m]][mm];								
							uu[8][qtime][mm] = u[ib[m][17]][qtime][mm] - u[ib[m][17]][pixt[m]][mm];						
						}

						// local stresses and strains in a pixel
						for (int n3 = 1; n3 <= 3; n3++)
						{
							for (int n8 = 1; n8 <= 8; n8++)
							{
								for (int n = 1; n <= 6; n++)
								{
									str11 = str11 + cmod[pix[m]][1][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										  										  
									str22 = str22 + cmod[pix[m]][2][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];																				  
									str33 = str33 + cmod[pix[m]][3][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];																				  
									str13 = str13 + cmod[pix[m]][4][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										  										  
									str23 = str23 + cmod[pix[m]][5][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];									  
									str12 = str12 + cmod[pix[m]][6][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										  										  
								}
							}
						}
					}
				}

				if (phasechange[m] == 0)
				{
					if (pix[m] == 12) set = 1;
					else set = q;

					for (int qtime = set; qtime <= q; qtime++)
					{
						//load in elements of 8-vector using pd. bd. conds.
						for (int mm = 1; mm <= 3; mm++)
						{
							uu[1][qtime][mm] = u[m][qtime][mm];
							uu[2][qtime][mm] = u[ib[m][3]][qtime][mm];								
							uu[3][qtime][mm] = u[ib[m][2]][qtime][mm];								
							uu[4][qtime][mm] = u[ib[m][1]][qtime][mm];								
							uu[5][qtime][mm] = u[ib[m][26]][qtime][mm];								
							uu[6][qtime][mm] = u[ib[m][19]][qtime][mm];								
							uu[7][qtime][mm] = u[ib[m][18]][qtime][mm];								
							uu[8][qtime][mm] = u[ib[m][17]][qtime][mm];								
						}

						//correction for period boundary conditions
						if (i == nx)
						{
							uu[2][qtime][1] = uu[2][qtime][1] + exx * nx;								
							uu[2][qtime][2] = uu[2][qtime][2] + exy * nx;								
							uu[2][qtime][3] = uu[2][qtime][3] + exz * nx;								
							uu[3][qtime][1] = uu[3][qtime][1] + exx * nx;								
							uu[3][qtime][2] = uu[3][qtime][2] + exy * nx;								
							uu[3][qtime][3] = uu[3][qtime][3] + exz * nx;								
							uu[6][qtime][1] = uu[6][qtime][1] + exx * nx;								
							uu[6][qtime][2] = uu[6][qtime][2] + exy * nx;								
							uu[6][qtime][3] = uu[6][qtime][3] + exz * nx;								
							uu[7][qtime][1] = uu[7][qtime][1] + exx * nx;								
							uu[7][qtime][2] = uu[7][qtime][2] + exy * nx;								
							uu[7][qtime][3] = uu[7][qtime][3] + exz * nx;								
						}
						if (j == ny)
						{
							uu[3][qtime][1] = uu[3][qtime][1] + exy * ny;								
							uu[3][qtime][2] = uu[3][qtime][2] + eyy * ny;								
							uu[3][qtime][3] = uu[3][qtime][3] + eyz * ny;								
							uu[4][qtime][1] = uu[4][qtime][1] + exy * ny;								
							uu[4][qtime][2] = uu[4][qtime][2] + eyy * ny;								
							uu[4][qtime][3] = uu[4][qtime][3] + eyz * ny;								
							uu[7][qtime][1] = uu[7][qtime][1] + exy * ny;								
							uu[7][qtime][2] = uu[7][qtime][2] + eyy * ny;								
							uu[7][qtime][3] = uu[7][qtime][3] + eyz * ny;								
							uu[8][qtime][1] = uu[8][qtime][1] + exy * ny;								
							uu[8][qtime][2] = uu[8][qtime][2] + eyy * ny;								
							uu[8][qtime][3] = uu[8][qtime][3] + eyz * ny;								
						}
						if (k == nz)
						{
							uu[5][qtime][1] = uu[5][qtime][1] + exz * nz;								
							uu[5][qtime][2] = uu[5][qtime][2] + eyz * nz;								
							uu[5][qtime][3] = uu[5][qtime][3] + ezz * nz;								
							uu[6][qtime][1] = uu[6][qtime][1] + exz * nz;								
							uu[6][qtime][2] = uu[6][qtime][2] + eyz * nz;								
							uu[6][qtime][3] = uu[6][qtime][3] + ezz * nz;								
							uu[7][qtime][1] = uu[7][qtime][1] + exz * nz;								
							uu[7][qtime][2] = uu[7][qtime][2] + eyz * nz;								
							uu[7][qtime][3] = uu[7][qtime][3] + ezz * nz;								
							uu[8][qtime][1] = uu[8][qtime][1] + exz * nz;								
							uu[8][qtime][2] = uu[8][qtime][2] + eyz * nz;								
							uu[8][qtime][3] = uu[8][qtime][3] + ezz * nz;								
						}

						//local stresses and strains in a pixel
						for (int n3 = 1; n3 <= 3; n3++)
						{
							for (int n8 = 1; n8 <= 8; n8++)
							{
								for (int n = 1; n <= 6; n++)
								{
									str11 = str11 + cmod[pix[m]][1][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										 										 
									str22 = str22 + cmod[pix[m]][2][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										 										  
									str33 = str33 + cmod[pix[m]][3][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										 										 
									str13 = str13 + cmod[pix[m]][4][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										 										 
									str23 = str23 + cmod[pix[m]][5][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										 										 
									str12 = str12 + cmod[pix[m]][6][n][q + 1 - qtime] * es[n][n8][n3] * uu[n8][qtime][n3];										 										 
								}
							}
						}
					}
				}

				if (pix[m] == 12){
					counter = counter + 1;				
					//stress_fout << "Time step " << q << ", CSH voxel " << counter << ", m = " << m << " (" << i << ", " << j << ", "<< k << "), stresses: " 
					//<< str11 << ", " << str22 << ", " << str33 
					//<< ", " << str13 << ", " << str23 << ", " << str12 << std::endl;

					//2020.04.08 for analysis of CSH stress distribution
					stress_fout << str11 << " " << str22 << " " << str33 << " " << str13 << " " << str23 << " " << str12 << std::endl;
				}

				//sum local strains and stresses into global values
				strxx = strxx + str11;
				stryy = stryy + str22;
				strzz = strzz + str33;
				strxz = strxz + str13;
				stryz = stryz + str23;
				strxy = strxy + str12;				
			}
		}
	}

	// Volume average of global stresses and strains
	strxx = strxx / (1.0* ns);
	stryy = stryy / (1.0* ns);
	strzz = strzz / (1.0* ns);
	strxz = strxz / (1.0* ns);
	stryz = stryz / (1.0* ns);
	strxy = strxy / (1.0* ns);
	sxx = exx;
	syy = eyy;
	szz = ezz;
	sxz = exz;
	syz = eyz;
	sxy = exy;
}

// This is a stress relaxation test, apply constant strain to the system
// define Relaxation Modulus function E(t) for viscoelastic materials (CSH) in GPa
// t is the loading time of each specific CSH pixel, i.e. when constant strain is applied
double Fmodulus(double t)
{
	return 11.2 + 11.2 * exp(-0.2 * t);
}



int main()
{
	
	///////////////////////////////////////////////////////////////////////////////
	/*
	** BACKGROUND:
	**
	** This program solves the linear elastic equations in a
	** random linear aging viscoelastic material, subject to an applied macroscopic strain,
	** using the finite element method. Each pixel in the 3-D digital
	** image is a cubic tri-linear finite element, having its own
	** moduli tensor. Periodic boundary conditions are maintained.
	** In the comments below, (USER) means that this is a section of code that
	** the user might have to change for his particular problem. Therefore the
	** user is encouraged to search for this string.
	*/

	/*
	** PROBLEM AND VARIABLE DEFINITION:
	**
	** The problem being solved is the minimization of the energy
	** 1/2 uAu + b u + C, where A is the Hessian matrix composed of the
	** stiffness matrices (dk) for each pixel/element, b is a constant vector
	** and C is a constant that are determined by the applied strain and
	** the periodic boundary conditions, and u is a vector of all the displacements.
	**
	** The solution method used is the conjugate gradient relaxation algorithm.
	**
	** Other variables are: 
	** gb is the gradient = Au+b, 
	** h and Ah are auxiliary variables used in the conjugate gradient algorithm (in dembx),
	** dk(n,i,j) is the stiffness matrix of the n'th phase, 
	** cmod(n,i,j) is the elastic moduli tensor of the n'th phase,
	** pix is a vector that gives the phase label of each pixel,
	** ib is a matrix that gives the labels of the 27 (counting itself) neighbors of a given node,
	** prob is the volume fractions of the various phases,
	** strxx, stryy, strzz, strxz, stryz, strxy are the six Voigt volume averaged total stresses,
	** sxx, syy, szz, sxz, syz, and sxy are the six Voigt volume averaged total strains.
	*/

	/*
	** DIMENSIONS:
	**
	** The vectors u,gb,b,h, and Ah are dimensioned to be the system size,
	** ns=nx*ny*nz, with three components, where the digital image of the
	** microstructure considered is a rectangular paralleliped, nx * ny * nz in size.
	** The arrays pix and ib are also dimensioned to the system size.
	** The array ib has 27 components, for the 27 neighbors of a node.
	** Note that the program is set up at present to have at most 100
	** different phases. This can easily be changed, simply by changing
	** the dimensions of dk, prob, and cmod. The parameter nphase gives the
	** number of phases being considered in the problem.
	** All arrays are passed between subroutines using simple common statements.
	*/

	/*
	** STRONGLY SUGGESTED:  READ THE MANUAL BEFORE USING PROGRAM!!
	*/

	/*
	** (USER) Change these dimensions and in other subroutines at same time.
	** For example, search and replace all occurrences throughout the
	** program of "(8000" by "(64000", to go from a
	** 20 x 20 x 20 system to a 40 x 40 x 40 system.
	*/
	///////////////////////////////////////////////////////////////////////////////

	double phasemod[101][51][3]; //phasemod[phase][time step][bulk/shear single value]
	double prob[101]; //prob[phase], volume fraction of each phase
	int in[28]; //MANUAL table 3: use delta i, delta j, delta k to construct 27 neighbors of 1 element
	int jn[28];
	int kn[28];

	int nphase = 22;//nphase value should not be zero. otherwise change it in ppixel function
	int step = 3; //step number should be equal to number of Thames microstructure images
		
	//set log time: t = 0.05-56, with number of step values
	/*
	double tmin = 0.05;
	double tmax = 56;
	for (int q = 1; q <= step; q++)
	{
		t[q] = tmin*pow(tmax / tmin, (q*1.0 - 1.) / (step*1.0 - 1.)); //0.5*(q*1.0-1.);
	}
	*/
	
	double t[18]; //time step should be set accordingly to Thames microstructure images
	//t[1]=0.0;
	//t[2]=0.1834;
	//t[3]=0.4004;
	//t[4]=0.6571;
	t[1]=1.1; //t[5]=0.9610;
	t[2]=1.3206; //t[6]
	t[3]=1.7461;
	/*
	t[4]=2.2497;
	t[5]=2.8456;
	t[6]=3.5507;
	t[7]=4.3852;
	t[8]=5.3726;
	t[9]=6.5412;
	t[10]=7.9240;
	t[11]=9.5604;
	t[12]=11.4968;
	t[13]=13.7883;
	t[14]=16.5;
	t[15]=19.7090;
	t[16]=23.5063;
	t[17]=28.0;
	*/

	//read in file in string
	string prefilename[22]; //file number should be equal to number of Thames microstructure images
	char filename[15];

	prefilename[1] = "10.img";
	prefilename[2] = "11.img";
	prefilename[3] = "12.img";
	/*
	prefilename[4] = "8.img";
	prefilename[5] = "9.img";
	prefilename[6] = "10.img";
	prefilename[7] = "11.img";
	prefilename[8] = "12.img";
	prefilename[9] = "13.img";
	prefilename[10] = "14.img";
	prefilename[11] = "15.img";
	prefilename[12] = "16.img";
	prefilename[13] = "17.img";
	prefilename[14] = "18.img";
	prefilename[15] = "19.img";
	prefilename[16] = "20.img";
	prefilename[17] = "21.img";

	prefilename[18] = "18.img";
	prefilename[19] = "19.img";
	prefilename[20] = "20.img";
	prefilename[21] = "21.img";
	
	prefilename[22] = "22.img";
	prefilename[23] = "23.img";
	prefilename[24] = "24.img";
	prefilename[25] = "25.img";
	prefilename[26] = "26.img";
	prefilename[27] = "27.img";
	prefilename[28] = "28.img";
	prefilename[29] = "29.img";
	prefilename[30] = "30.img";
	prefilename[31] = "31.img";
	*/

	//composite dimension
	int nx = 100;
	int ny = 100;
	int nz = 100;
	int ns = nx * ny * nz; //total number of sites

	//assign phase initial aging state and aging time	
	for (int k = 1; k <= nz; k++)
	{
		for (int j = 1; j <= ny; j++)
		{
			for (int i = 1; i <= nx; i++)
			{
				int m = nx * ny * (k - 1) + nx * (j - 1) + i;
				phasechange[m] = 0; //initialize aging state, initialize before time loop
				pixt[m] = 0; //initialize aging time
			}
		}
	}

	//assign phase modulus
	/*
	** (USER)
	** The parameter phasemod(i,j) is the bulk (i,1) and shear (i,2) value of the i'th phase.
	** These can be input in terms of Young's moduli E(i,1) and Poisson's ratio nu(i,2).
	** There is a loop below that changes them to bulk and shear moduli.
	** For anisotropic elastic material, one can directly input
	** the elastic moduli tensor cmod in subroutine femat, and skip this part.
	** If you wish to input in terms of bulk (1) and shear (2), then make sure
	** to comment out the loop for conversion.
	*/
	for (int j = 1; j <= step; j++)
	{
		for (int i = 1; i <= nphase; i++)
		{
			phasemod[i][j][1] = 0.;
			phasemod[i][j][2] = 0.;
		}
	}

	//!!!!!!!!!!!!!!!when poisson's ratio is not constant,use laplace transform!!!!!!!
	//These data are from 2005 paper: Modeling the linear elastic properties of Portland cement paste
	//The values are 1-Young's modulus (in GPa) and 2-Poisson's ratio.
	phasemod[0][1][1] = 0.0; //phase 0 is void
	phasemod[0][1][2] = 0.00;
	phasemod[1][1][1] = 0.0001; //phase 1 is water, Origianl 0,0 is wrong, bulk modulus of water should not be zero
	phasemod[1][1][2] = 0.4999921; //use E~0, v~0.5, and K~2.1GPa
	phasemod[2][1][1] = 117.6; //phase 2 is C3S
	phasemod[2][1][2] = 0.314;
	phasemod[3][1][1] = 117.6; //phase 3 is C2S
	phasemod[3][1][2] = 0.314;
	phasemod[4][1][1] = 117.6; //phase 4 is C3A
	phasemod[4][1][2] = 0.314;
	phasemod[5][1][1] = 117.6; //phase 5 is C4AF
	phasemod[5][1][2] = 0.314;

	phasemod[6][1][1] = 44.2; //phase 6 is K2SO4
	phasemod[6][1][2] = 0.269;
	phasemod[7][1][1] = 57.1; //phase 7 is Na2SO4
	phasemod[7][1][2] = 0.281;
	phasemod[8][1][1] = 45.7; //phase 8 is Gypsum (dihydrate)
	phasemod[8][1][2] = 0.33;
	phasemod[9][1][1] = 62.85; //phase 9 is Hemihydrate (may not match THAMES: HEMIANH)
	phasemod[9][1][2] = 0.3025;
	phasemod[10][1][1] = 79.6; //phase 10 is CaCO3
	phasemod[10][1][2] = 0.31;
	phasemod[11][1][1] = 42.3; //phase 11 is CH
	phasemod[11][1][2] = 0.324;

	phasemod[12][1][1] = 22.4; //phase 12 is CSH, average properties
	phasemod[12][1][2] = 0.25;

	phasemod[13][1][1] = 42.3; //phase 13 is AFMC (same as CH)
	phasemod[13][1][2] = 0.324;
	phasemod[14][1][1] = 42.3; //phase 14 is Monosulfate (same as CH)
	phasemod[14][1][2] = 0.324;
	phasemod[15][1][1] = 22.4; //phase 15 is Ettringite (same as CSH)
	phasemod[15][1][2] = 0.25;
	phasemod[16][1][1] = 42.3; //phase 16 is Brucite (Mg(OH)2) where are these values from???
	phasemod[16][1][2] = 0.324;
	phasemod[17][1][1] = 22.4; //phase 17 is Hydrotalcite where are these values from???
	phasemod[17][1][2] = 0.25;

	phasemod[18][1][1] = 42.3; //phase 18 is AFM (same as CH)
	phasemod[18][1][2] = 0.324;
	phasemod[19][1][1] = 117.6; //phase 19 is Lime where are these values from ???
	phasemod[19][1][2] = 0.314;
	phasemod[20][1][1] = 117.6; //phase 20 is Periclase
	phasemod[20][1][2] = 0.314;
	phasemod[21][1][1] = 15.0; //phase 21 is Damage
	phasemod[21][1][2] = 0.4;

	//assign Poisson's ratio of each phase for all time steps
	//assume Poisson's ratio do not change with time for all phases
	for (int i = 2; i <= step; i++)
	{	
		for (int j = 1; j <= nphase; j++)
		{
			phasemod[j][i][2] = phasemod[j][1][2]; 
		}
	}

	//assign Young's modulus of CSH according to the time dependent function
	//if CSH is elastic, delete the following loop
	
	for (int i = 2; i <= step; i++)
	{	
		phasemod[12][i][1] = Fmodulus(t[i])-Fmodulus(t[i-1]); //Fmodulus(t[i]-t[1]); // Fmodulus(t[i] - t[1])-Fmodulus(t[i-1] - t[1]); //TODO, how to update for each CSH pixel that appears at each time step????
	}
	

	//convert (E, v) to (K, G)
	//if input is Young's modulus and Poisson's ration, use this loop. otherwise input Bulk and Shear skip this one.
	for (int i = 1; i <= nphase; i++)
	{
		for (int j = 1; j <= step; j++)
		{
			double save = phasemod[i][j][1];
			phasemod[i][j][1] = phasemod[i][j][1] / 3. / (1. - 2. * phasemod[i][j][2]); //K=E/3/(1-2v)
			phasemod[i][j][2] = save / 2. / (1. + phasemod[i][j][2]); //G=E/2/(1+v)
		}
	}

	//set controlled strain, for shear, 2 should be divided
	/*
	** (USER) Set applied strains
	** Actual shear strain applied in loop is exy, exz, and eyz as given in the statements below.
	** The engineering shear strain, by which the shear modulus is usually defined, is twice of these values.
	*/
	exx = 0.1; // Is it tensile strain???
	eyy = 0.1;
	ezz = 0.1;
	exz = 0.0;
	eyz = 0.0;
	exy = 0.0;

	
	/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
	/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TIME LOOP Starts !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
	/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

	for (int q = 1; q <= step; q++)
	{
		std::cout << '\n';
		std::cout << "Time step: " << q << std::endl;

		//read in microstructure file starting from 1.txt
		strcpy(filename, prefilename[q].c_str());
		ifstream fin(filename);

		/*
		** (USER) gtest is the stopping criterion, 
		** gb is the gradient of the energy,
		** gg = gb * gb is the norm squared of the gradient,
		** gtest = (small number) * ns, 
		** so that when gg < gtest, 
		** the RMS (root mean square) value per pixel of gb is less than sqrt(small number).
		** i.e. gb/sqrt(ns) < (small number)
		*/
		double gtest = pow(10.0, -8) * ns; //modified: 10^(-12)

		//assign each phase and output volume fraction
		//Read in a microstructure in subroutine ppixel, and set up pix[m] with phase label
		ppixel(nx, ny, nz, ns, nphase, filename); 

		
		//Count and output the volume fractions of the different phases
		assig(ns, nphase, prob);

		int counter_csh = 0; //YZ: added to count no. of newly formed CSH
		//check whether aging occurs
		for (int k = 1; k <= nz; k++)
		{
			for (int j = 1; j <= ny; j++)
			{
				for (int i = 1; i <= nx; i++)
				{
					int m = nx * ny * (k - 1) + nx * (j - 1) + i;
					if (q > 1)
					{
						if (pix[m] != pixstor[m])
						{
							
							phasechange[m] = 1;//aging will occur, pixel m is different from last image
							pixt[m] = q - 1; //store the time this pixel changes phase, if pixt[m]=4, pixel m change phase at the 5th image. Each pixel has its formation time 
							
							phasechange_fout << "Time step: " << q << ", Phase change occurs at pixel: " << m << " (" << i << ", " << j << ", " << k << ")"
							<< ", in phase: " << pixstor[m] << ", now is: " << pix[m] << std::endl;

							if (pix[m] == 12) { //if a new CSH voxel is formed, output the location
								counter_csh = counter_csh + 1;
								phasechange_csh_fout << "Time step: " << q << ", Phase change occurs at pixel: " << m << " (" << i << ", " << j << ", " << k << ")"
								<< ", CSH: " << counter_csh
								<< ", in phase: " << pixstor[m] << ", now is: " << pix[m] << std::endl;
							}
						}
					}
					pixstor[m] = pix[m]; //initialize at first time step q=1, stores the previous microstructure
				}
			}
		}


		//construct the neighbor table, ib(m,n)
		//first construct labels: delta i, delta j, delta k (MANUAL table 3)
		in[1] = 0; jn[1] = 1;
		in[2] = 1; jn[2] = 1;
		in[3] = 1; jn[3] = 0;
		in[4] = 1; jn[4] = -1;
		in[5] = 0; jn[5] = -1;
		in[6] = -1; jn[6] = -1;
		in[7] = -1; jn[7] = 0;
		in[8] = -1; jn[8] = 1;	

		for (int n = 1; n <= 8; n++)
		{
			kn[n] = 0; //k1-k8 = 0;
			kn[n + 8] = -1; //k9-k16 = -1;
			kn[n + 16] = 1; //k17-k24 = 1
			in[n + 8] = in[n]; 
			in[n + 16] = in[n];
			jn[n + 8] = jn[n];
			jn[n + 16] = jn[n];
		}

		//27th is the middle self, 25 is the center below, 26 is the center above
		in[25] = 0; jn[25] = 0; kn[25] = -1;
		in[26] = 0; jn[26] = 0; kn[26] = 1;
		in[27] = 0; jn[27] = 0; kn[27] = 0;

		//Now construct neighbor table ib(m,n) according to 1-d labels
		//Matrix ib(m,n) gives the 1-d label of the n'th neighbor (n=1-27) of the node labeled m.
		for (int k = 1; k <= nz; k++)
		{
			for (int j = 1; j <= ny; j++)
			{
				for (int i = 1; i <= nx; i++)
				{
					int m = nx *  ny* (k - 1) + nx * (j - 1) + i;
					for (int n = 1; n <= 27; n++)
					{
						int i1 = i + in[n];
						int j1 = j + jn[n];
						int k1 = k + kn[n];
						//periodic boundary
						if (i1 < 1) i1 = i1 + nx; 
						if (i1 > nx) i1 = i1 - nx;
						if (j1 < 1) j1 = j1 + ny;
						if (j1 > ny) j1 = j1 - ny;
						if (k1 < 1) k1 = k1 + nz;
						if (k1 > nz) k1 = k1 - nz;
						int m1 = nx * ny * (k1 - 1) + nx * (j1 - 1) + i1;
						ib[m][n] = m1; //voxel m's n-th neighbor is voxel m1
					}
				}
			}
		}

		
		//Compute the average stress and strain in each microstructure.
		/*
		** Set up the elastic modulus variables cmod, finite element stiffness matrices dk,
		** the constant, C, and vector, b, required for computing the energy.
		** (USER) If anisotropic elastic moduli tensors are used, these need to be
		** input in subroutine femat.
		*/
		
		femat(nx, ny, nz, ns, phasemod, nphase, q);
		
		
		//Apply a homogeneous macroscopic strain as the initial condition.

		for (int k = 1; k <= nz; k++)
		{
			for (int j = 1; j <= ny; j++)
			{
				for (int i = 1; i <= nx; i++)
				{
					int	m = nx * ny * (k - 1) + nx * (j - 1) + i;
					double x = 1.0 * (i - 1); //YG: 1.0 is the size of 1 pixel, 1um. x,y,z are coordinates of each node
					double y = 1.0 * (j - 1);
					double z = 1.0 * (k - 1);
					//at boundaries, displacements are 0, this will be corrected in strain()
					u[m][q][1] = x * exx + y * exy + z * exz; //u is the nodal displacements vector
					u[m][q][2] = x * exy + y * eyy + z * eyz;
					u[m][q][3] = x * exz + y * eyz + z * ezz;
				}
			}
		}


		//RELAXATION LOOP
		//(USER) 
		//kmax is the maximum number of times dembx will be called,
		//ldemb is the conjugate gradient steps performed during each call. 
		//The total number of conjugate gradient steps allowed for a given elastic computation is kmax * ldemb.
		
		int kmax = 5; //40
		int	ldemb = 10; //50
		int	ltot = 0;
		double utot;
		int Lstep;
		
		//Call energy to get initial energy and initial gradient
		energy(nx, ny, nz, ns, utot, q);

		double gg = 0.0; //gg is the norm squared of the gradient (gg=gb*gb)

		for (int m3 = 1; m3 <= 3; m3++)
		{
			for (int m = 1; m <= ns; m++)
			{
				gg = gg + gb[m][m3] * gb[m][m3];
			}
		}

		for (int kkk = 1; kkk <= kmax; kkk++)
		{
			std::cout << "      kkk loop: " << kkk << std::endl;
			/*
			** call dembx to go into the conjugate gradient solver
			*/
			dembx(ns, Lstep, gg, dk, gtest, ldemb, kkk, q);
			ltot = ltot + Lstep;
			/*
			** Call energy to compute energy after dembx call. If gg < gtest, this
			** will be the final energy. If gg is still larger than gtest, then this
			** will give an intermediate energy with which to check how the
			** relaxation process is coming along.
			*/
			energy(nx, ny, nz, ns, utot, q);

			//If relaxation process is finished, jump out of loop
			if (gg <= gtest) break;
			
			/*
			** If relaxation process will continue, compute and output stresses
			** and strains as an additional aid to judge how the
			** relaxation procedure is progressing.
			*/
			
			stress(nx, ny, nz, ns, q);
		}

		std::cout << "  Total energy is minimized." << std::endl;

		stress2(nx, ny, nz, ns, q);

		//calculate modulus and possion's ratio using calculated stresses and strains
		fout << "Stresses: " << strxx << ", " << stryy << ", " << strzz << ", " << strxy << ", " << stryz << ", " << strxz 
		<< ", Volume fractions of CSH: " << prob[12] << std::endl; 

	}

	return 0;
}
