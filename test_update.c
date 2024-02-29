#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define NROWS 		20
#define NCOLS		20
#define dom_length 	1
#define g 		1.1
#define Iteration 	100000

void update (float *uu, float *vv, float *pp, float *x, float *y, float *xc, float *yc, int nx, int ny, int top, int bottom, int left, int right);
void initialize(int nx, int ny, float *uu, float *vv, float *pp);
void print_array(float *uu, int nx, int ny);
void print_vector(float *uu, int nx);

int main ()
{ 
	float 	u[NROWS-1][NCOLS], 
			v[NROWS][NCOLS-1], 
			p[NROWS][NCOLS];
	float   xc[NROWS], yc[NCOLS], 
		x[NROWS-1], y[NCOLS-1];
		
	float hx = (float) dom_length/(NROWS-2), hy = (float) dom_length/(NCOLS-2);	
	printf ("hx = %f , hy = %f\n", hx , hy);
 	initialize(NROWS, NCOLS, &u[0][0], &v[0][0], &p[0][0]);
 	
 	//Create stretched mesh
 	for (int i = 0; i < NROWS-1; i++)
 	{
 		x[i] = 1 - tanh(g*(1-2*(float)i/(NROWS-2)))/tanh(g);
 		//x[i] = i * hx;
 	}
 	
 	for (int j = 0; j < NCOLS-1; j++)
 	{
	 	y[j] = 1 - tanh(g*(1-2*(float)j/(NCOLS-2)))/tanh(g);
	 	//y[j] = j * hy;
 	}
	
	for (int i = 1; i < NROWS-1; i++)
	{
		xc[i] = (x[i] + x[i-1])/2;
	}
	xc[0] = - xc[1]; xc[NROWS-1] = x[NROWS-2] + (x[NROWS-2] - x[NROWS-3])/2; 
	
	for (int j = 1; j < NCOLS-1; j++)
	{
		yc[j] = (y[j] + y[j-1])/2;
	}
	yc[0] = - yc[1]; yc[NCOLS-1] = y[NCOLS-2] + (y[NCOLS-2] - y[NCOLS-3])/2;
	
	//check the stretched mesh
	print_vector(&x[0], NROWS-1);	print_vector(&xc[0], NROWS);
	
	
	for (int it = 0; it <= Iteration; it ++)
	{
		update (&u[0][0], &v[0][0], &p[0][0], &x[0], &y[0], &xc[0], &yc[0], NROWS, NCOLS, -1, -1, -1, -1);
		float error = 0;
		for (int i = 1; i < NROWS-2; i++)
			{
				for (int j = 1; j < NCOLS-2; j++)
				{
					error = error + fabs((u[i][j] - u[i-1][j])/(xc[i]-xc[i-1]) + (v[i][j] - v[i][j-1])/(yc[j]-yc[j-1]));
				}
			}
		if(it % 200 == 0) printf("\nIteration: %d, 	Error: %f \n", it, error);
		if(it == Iteration) {
		//print_array(&u[0][0], NROWS-1, NCOLS); 
		//print_array(&v[0][0], NROWS, NCOLS-1);
		//print_array(&p[0][0], NROWS, NCOLS);
		}
	}
	float u_final[NROWS-1][NCOLS-1], v_final[NROWS-1][NCOLS-1], p_final[NROWS-1][NCOLS-1];
	for (int i = 0; i < NROWS-1; i++)
			{
				for (int j = 0; j < NCOLS-1; j++)
				{
					u_final[i][j] = 0.5*(u[i][j]+u[i][j+1]);
					v_final[i][j] = 0.5*(v[i][j]+v[i+1][j]);
					p_final[i][j] = 0.25*(p[i][j] + p[i+1][j] + p[i][j+1] + p[i+1][j+1]);
				}
			}
			
			print_array(&u_final[0][0], NROWS-1, NCOLS-1); 
		
		// OUTPUT DATA
			FILE *fout2, *fout3;
			fout2 = fopen("UVP.plt","w+t");
			fout3 = fopen("Central_U.plt","w+t");

			if ( fout2 == NULL )
			{
				printf("\nERROR when opening file\n");
				fclose( fout2 );
			}

			else
			{
				fprintf( fout2, "VARIABLES=\"X\",\"Y\",\"Vel\", \"U\",\"V\",\"P\"\n");
				fprintf( fout2, "ZONE  F=POINT\n");
				fprintf( fout2, "I=%d, J=%d\n", NROWS-1, NCOLS-1 );
				float vel[NROWS-1][NCOLS-1];
				for ( int j = 0 ; j < NCOLS-1 ; j++ )
				{
					for ( int i = 0 ; i < NROWS-1 ; i++ )
					{
						float xpos, ypos;
						xpos = x[i]/2;
						ypos = 1-y[j]/2;
						vel[i][j] = sqrt(pow(u_final[i][j],2) + pow(v_final[i][j],2));
						fprintf( fout2, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, vel[i][j], u_final[i][j], v_final[i][j], p_final[i][j] );
					}
				}
			}
	
			fclose( fout2 );
			float data_u[17] = 
			{1,0.84123,0.78871,0.73722,0.68717,0.23151,0.00332,-0.13641,-0.20581,-0.2109,-0.15662,-0.1015,-0.06434,-0.04775,-0.04192,-0.03717,0};
			float data_y[17] =
			{0,0.0234,0.0312,0.0391,0.0469,0.1484,0.2656,0.3828,0.5,0.5469,0.7187,0.8281,0.8984,0.9297,0.9375,0.9453,1};
			
			// CENTRAL --U
			fprintf(fout3, "VARIABLES=\"U\",\"Y\"\n");
			fprintf(fout3, "ZONE F=POINT\n");
			fprintf(fout3, "I=%d\n", NROWS-1 );
			for ( int j = 0 ; j < NCOLS-1 ; ++j )
			{
				float ypos = y[j]/2;
				if ((NROWS-1) % 2 != 0)fprintf( fout3, "%0.6f \t %0.6f\n", (u_final[(NROWS-1)/2][j]), ypos);
				if ((NROWS-1) % 2 == 0)fprintf( fout3, "%0.6f \t %0.6f\n", (u_final[NROWS/2][j] + u_final[(NROWS/2)-1][j])/2, ypos);
			}
			fprintf(fout3, "VARIABLES=\"U\",\"Y\"\n");
			fprintf(fout3, "ZONE F=POINT\n");
			fprintf(fout3, "I=%d\n", 17 );
			for ( int j = 0 ; j < 17 ; ++j )
			{
				fprintf( fout3, "%0.6f \t %0.6f\n", data_u[j], data_y[j] );
			}
			fclose( fout3 );
	
	return 0;
}



/************************************************************************************************************************************/




void update (float *uu, float *vv, float *pp, float *x, float *y, float *xc, float *yc, int nx, int ny, int top, int bottom, int left, int right)
{
	float 	hx = (float) dom_length/(NROWS-2), hy = (float) dom_length/(NCOLS-2),
			delta = 4.5, dt = 0.001, Re = 100;
	int nux, nvy;
	nux = nx, nvy = ny;
	if (bottom < 0) {nux = nx -1;}
	if (right < 0) {nvy = ny -1;}
	//printf("\n nux: %d, ny: %d\n", nux,ny);
	
	float u[nux][ny], 		v[nx][nvy], 		p[nx][ny];
	float u_new[nux][ny], 	v_new[nx][nvy],		p_new[nx][ny];
	
	for (int i = 0; i < nux; i++) 
	{
        for (int j = 0; j < ny; j++) 
		{
            u[i][j] = 0;
			u[i][j] = *(uu + i*ny + j);
			u_new[i][j] = u[i][j];
        }
    }
	
	for (int i = 0; i < nx; i++) 
	{
        for (int j = 0; j < nvy; j++) 
		{
            v[i][j] = 0;
			v[i][j] = *(vv + i*nvy + j);
			v_new[i][j] = v[i][j];
        }
    }
	
	for (int i = 0; i < nx; i++) 
	{
        for (int j = 0; j < ny; j++) 
		{
            p[i][j] = 0;
			p[i][j] = *(pp + i*ny + j);
			p_new[i][j] = p[i][j];
        }
    }

	float pressure = 0, advection_x = 0, advection_y = 0, diffusion = 0;
	//x-momentum eq. - Interior
	for (int i = 1; i < nux-1; i++)
	{
		for (int j = 1; j < ny-1; j++)
		{
			pressure = -(p[i+1][j] - p[i][j])/(xc[i+1]-xc[i]);
			diffusion = (1/Re)*(
				((u[i+1][j]-u[i][j])/(x[i+1]-x[i]) - (u[i][j]-u[i-1][j])/(x[i]-x[i-1]))/(xc[i+1]-xc[i]) +
				((u[i][j+1]-u[i][j])/(yc[j+1]-yc[j]) - (u[i][j]-u[i][j-1])/(yc[j]-yc[j-1]))/(yc[j]-yc[j-1])
				);
			advection_x = (
				pow((0.5*(u[i+1][j]+u[i][j])),2)-
				pow((0.5*(u[i-1][j]+u[i][j])),2))/(xc[i+1]-xc[i]);  
			advection_y = (
				(0.25*(u[i][j-1]+u[i][j])*(v[i][j-1]+v[i+1][j-1]))-
				(0.25*(u[i][j+1]+u[i][j])*(v[i][j]+v[i+1][j])))/(y[j]-y[j-1]);
			u_new[i][j]= u[i][j]+dt*(pressure + diffusion - (advection_x + advection_y));
		}
	}

	//x-momentum eq. - Boundary
	for (int j=0; j<ny; ++j)
	{
		if(top < 0) 	{u_new[0][j] = 0.0;}
		if(bottom < 0) 	{u_new[nux-1][j] = 0.0;}
	}
	
	for (int i = 0 ; i<nux; ++i)
	{
		if(left < 0) 	{u_new[i][0] = 2 - u_new[i][1];}
		if(right < 0) 	{u_new[i][ny-1] = - u_new[i][ny-2];}
	}
		
	
	//y-momentum eq. - Interior
	for (int i = 1; i < nx-1; i++)
	{
		for (int j = 1; j < nvy-1; j++)
		{
			pressure = -(p[i][j] - p[i][j+1])/fabs(yc[j]-yc[j+1]);
			diffusion = (1/Re)*(
		            	((v[i+1][j]-v[i][j])/(xc[i+1]-xc[i]) - (v[i][j]-v[i-1][j])/(xc[i]-xc[i-1]))/(x[i]-x[i-1]) +
				((v[i][j+1]-v[i][j])/(yc[j+1]-yc[j]) - (v[i][j]-v[i][j-1])/(y[j]-y[j-1]))/(yc[j+1]-yc[j])
		            	);
			advection_y = (
		            pow((0.5*(v[i][j-1]+v[i][j])),2)-
		            pow((0.5*(v[i][j+1]+v[i][j])),2))/fabs(yc[j-1]-yc[j]);
			advection_x = (
		            (0.25*(v[i+1][j]+v[i][j])*(u[i][j+1]+u[i][j]))-
		            (0.25*(v[i-1][j]+v[i][j])*(u[i-1][j]+u[i-1][j+1])))/(x[i]-x[i-1]);
			v_new[i][j]= v[i][j]+dt*(pressure + diffusion - (advection_x + advection_y));
		}
		//printf("\npressure: %f, diffusion: %f, advection_x: %f, advection_y: %f\n", pressure, diffusion, advection_x, advection_y);
	}
	
	//y-momentum eq. - Boundary
	for (int j=0; j<nvy; ++j)
	{
		if(top < 0) 		v_new[0][j] = - v_new[1][j];
		if(bottom < 0)	v_new[nx-1][j] = - v_new[nx-2][j];
	}

	for (int i=0; i<nx; ++i)
	{
		if(left < 0)		v_new[i][0] = 0;
		if(right < 0)	v_new[i][nvy-1] = 0;
	}
	
	//Continuity equation - Interior 
	for (int i = 1; i < nx-1; i++)
	{
		for (int j = 1; j < ny-1; j++)
		{
			p_new[i][j] = p[i][j]-
            dt*delta*(u[i][j]-u[i-1][j])/(xc[i]-xc[i-1]) -
            dt*delta*(v[i][j-1]-v[i][j])/(yc[j]-yc[j-1]);
		}
	}
	
	//Continuity eq. - Boundary
	for (int i=1; i<nx-1; ++i)
	{
		if(left < 0)		p_new[i][0] = p_new[i][1];
		if(right < 0)	p_new[i][ny-1] = p_new[i][ny-2];
	}

	for (int j=0; j<ny; ++j)
	{
		if(top < 0)		p_new[0][j] = p_new[1][j];
		if(bottom < 0)	p_new[nx-1][j] = p_new[nx-2][j];
	}
	
	/*update matrix*/
	for (int i = 0; i < nux; i++) 
	{
        for (int j = 0; j < ny; j++) 
		{
            u[i][j] = u_new[i][j];
			*(uu + i*ny + j) = u[i][j];
        }
    }
	
	for (int i = 0; i < nx; i++) 
	{
        for (int j = 0; j < nvy; j++) 
		{
			v[i][j] = v_new[i][j];
			*(vv + i*nvy + j) = v[i][j];
        }
    }
	
	for (int i = 0; i < nx; i++) 
	{
        for (int j = 0; j < ny; j++) 
		{
			p[i][j] = p_new[i][j];
			*(pp + i*ny + j) = p[i][j];
        }
    }

	/*End update matrix*/
}

void initialize(int nx, int ny, float *uu, float *vv, float *pp)
{
    int i, j;
    for (i = 0; i < nx-1; i++) {
        for (j = 0; j < ny; j++) 
		{
            *(uu + i * ny + j) = 0;
			*(uu + i * ny) = 2;
        }
    }
	
	for (i = 0; i < nx; i++) 
	{
        for (j = 0; j < ny-1; j++) 
		{
            *(vv + i * (ny-1) + j) = 0;
        }
    }
	
	for (i = 0; i < nx; i++) 
	{
        for (j = 0; j < ny; j++) 
		{
            *(pp + i * ny + j) = 0;
        }
    }
	*(pp + (nx-1)*ny + ny -1) = 1;
	
}

void print_array(float *uu, int nx, int ny)
{
	printf("\n");
	for (int i =0; i< nx; i++)
	{
		for (int j = 0; j< ny; j++)
		{
			printf ("%.05f ", uu[i * ny + j]);
		}
	printf("\n");
	}
	printf("\n");
}
void print_vector(float *uu, int nx)
{
	printf("\n");
	for (int i =0; i< nx; i++)
	{
		printf ("%.04f ", uu[i]);
	}
	printf("\n");
}
