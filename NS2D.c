/* Two-dimensional Navier-Stokes equation */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define gnx 32
#define gny 32
#define iloop for(i=1;i<=gnx;i++)
#define i0loop for(i=0;i<=gnx;i++)
#define jloop for(j=1;j<=gny;j++)
#define j0loop for(j=0;j<=gny;j++)
#define ijloop iloop jloop
#define i0jloop i0loop jloop
#define ij0loop iloop j0loop
#define iloopt for(i=1;i<=nxt;i++)
#define i0loopt for(i=0;i<=nxt;i++)
#define jloopt for(j=1;j<=nyt;j++)
#define j0loopt for(j=0;j<=nyt;j++)
#define ijloopt iloopt jloopt
#define i0jloopt i0loopt jloopt
#define ij0loopt iloopt j0loopt

int nx,ny,n_level,p_relax,nt;
double **sor,h,h2,**tu,**tv,**workp,**worku,**workv,**adv_u,
        **adv_v,dt,xleft,xright,yleft,yright,Re;
char bufferu[20],bufferv[20],bufferp[20];

double **dmatrix(long nrl,long nrh,long ncl,long nch){
    double **m;
    long i,nrow=nrh-nrl+1+1,ncol=nch-ncl+1+1;
    m=(double **) malloc((nrow)*sizeof(double*));
    m+=1; m-=nrl;
    m[nrl]=(double *) malloc((nrow*ncol)*sizeof(double));
    m[nrl]+=1; m[nrl]-=ncl;
    for (i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;
    return m;
}

void free_dmatrix(double **m,long nrl,long nrh,long ncl,long nch){
    free(m[nrl]+ncl-1); free(m+nrl-1); 
}

void zero_matrix(double **a,int xl,int xr,int yl,int yr){
    int i,j;
    for (i=xl;i<=xr;i++) {for (j=yl;j<=yr;j++) {a[i][j]=0.0;}} 
}

void mat_copy(double **a,double **b,int xl,int xr,int yl,int yr){
    int i,j;
    for (i=xl;i<=xr;i++) {for (j=yl;j<=yr;j++) {a[i][j]=b[i][j];}} 
}

void mat_add(double **a,double **b,double **c,int xl,int xr,int yl,int yr){
    int i,j;
    for (i=xl;i<=xr;i++) {for (j=yl;j<=yr;j++) {a[i][j]=b[i][j]+c[i][j];}} 
}

void mat_sub(double **a,double **b,double **c,int xl,int xr,int yl,int yr){
    int i,j;
    for (i=xl;i<=xr;i++) {for (j=yl;j<=yr;j++) {a[i][j]=b[i][j]-c[i][j];}} 
}

double mat_max(double **a,int nrl,int nrh,int ncl,int nch){
    int i,j;
    double x = 0.0;
    for (i=nrl;i<=nrh;i++) {for (j=ncl;j<=nch;j++) 
    {if (fabs(a[i][j]) > x) {x = fabs(a[i][j]);}}}
    return x;
}

void div_uv(double **tu,double **tv,double **divuv,int nxt,int nyt){
    int i,j;
    double ht;
    ht = xright / (double) nxt;
ijloopt {divuv[i][j]=(tu[i][j]-tu[i-1][j]+tv[i][j]-tv[i][j-1])/ht;}
}

void grad_p(double **p,double **dpdx,double **dpdy,int nxt,int nyt){
    int i,j;
    double ht;
    ht = xright / (double) nxt;
    i0jloopt {if (i==0) {dpdx[0][j] = 0.0;}
        else if (i==nxt) {dpdx[nxt][j] = 0.0;}
        else {dpdx[i][j] = (p[i+1][j] - p[i][j])/ht;}}
    ij0loopt {if (j==0) {dpdy[i][0] = 0.0;}
        else if (j==nyt) {dpdy[i][nyt] = 0.0;}
        else {dpdy[i][j] = (p[i][j+1] - p[i][j])/ht;}}
}

void laplace_p(double **p,double **lap_p,int nxt,int nyt){
    double **dpdx,**dpdy;
    dpdx = dmatrix(0,nxt,1,nyt);dpdy = dmatrix(1,nxt,0,nyt);
    grad_p(p,dpdx,dpdy,nxt,nyt);div_uv(dpdx,dpdy,lap_p,nxt,nyt);
    free_dmatrix(dpdx,0,nxt,1,nyt);free_dmatrix(dpdy,1,nxt,0,nyt);
}

void prolong(double **u_coarse,double **u_fine,int nxt,int nyt){
    int i,j;
    ijloopt {u_fine[2*i-1][2*j-1]=u_fine[2*i-1][2*j]=
    u_fine[2*i][2*j-1]=u_fine[2*i][2*j]=u_coarse[i][j];} 
}

void restrict(double **u_fine,double **u_coarse,int nxt,int nyt){
    int i,j;
    ijloopt {u_coarse[i][j]=0.25*(u_fine[2*i-1][2*j-1]
    +u_fine[2*i-1][2*j]+u_fine[2*i][2*j-1]+u_fine[2*i][2*j]);} 
}

void residual_p(double **r,double **u,double **f,int nxt,int nyt){
    laplace_p(u,r,nxt,nyt); mat_sub(r,f,r,1,nxt,1,nyt); 
}

void relax_p(double **p,double **f,int nxt,int nyt){
    int i,j,iter;
    double ht,ht2,coef,src;
    ht2 = pow((xright-xleft) / (double) nxt,2);
    for (iter=1; iter<=p_relax; iter++) {ijloopt {src = f[i][j];
            if (i==1) {src -= p[2][j]/ht2; coef = -1.0/ht2;}
            else if (i==nxt) {src -= p[nxt-1][j]/ht2; coef = -1.0/ht2;}
            else {src -= (p[i+1][j] + p[i-1][j])/ht2; coef = -2.0/ht2;}
            if (j==1) {src -= p[i][2]/ht2; coef += -1.0/ht2;}
            else if (j==nyt) {src -= p[i][nyt-1]/ht2; coef += -1.0/ht2;}
            else {src -= (p[i][j+1] + p[i][j-1])/ht2; coef += -2.0/ht2;}
            p[i][j] = src / coef;}}
}

void pressure_update(double **a){
    int i,j;
    double ave = 0.0;
    ijloop {ave = ave + a[i][j];} ave /= (nx+0.0)*(ny+0.0);
    ijloop {a[i][j] -= ave;}
}

void vcycle_uv(double **uf,double **ff,int nxf,int nyf,int ilevel){
    relax_p(uf,ff,nxf,nyf);
    if (ilevel < n_level) {
        int nxc,nyc;
        double **rf,**uc,**fc;
        nxc=nxf / 2; nyc=nyf / 2; rf=dmatrix(1,nxf,1,nyf);
        uc=dmatrix(1,nxc,1,nyc); fc=dmatrix(1,nxc,1,nyc);
        residual_p(rf,uf,ff,nxf,nyf); restrict(rf,fc,nxc,nyc);
        zero_matrix(uc,1,nxc,1,nyc);
        vcycle_uv(uc,fc,nxc,nyc,ilevel + 1);
        prolong(uc,rf,nxc,nyc); mat_add(uf,uf,rf,1,nxf,1,nyf);
        relax_p(uf,ff,nxf,nyf); free_dmatrix(rf,1,nxf,1,nyf);
        free_dmatrix(uc,1,nxc,1,nyc); free_dmatrix(fc,1,nxc,1,nyc);}
}

void MG_Poisson(double **p,double **f)
{
	int i,j,max_it = 2000,it_mg = 1;
	double tol = 1.0e-5,resid = 1.0;
	mat_copy(workv,p,1,nx,1,ny);
	while (it_mg <= max_it && resid >= tol) 
	{
		it_mg++;
		vcycle_uv(p,f,nx,ny,1); pressure_update(p);
		ijloop {sor[i][j] = workv[i][j] - p[i][j];}
		resid=mat_max(sor,1,nx,1,ny);mat_copy(workv,p,1,nx,1,ny);
	}
	printf("Mac iteration = %d  residual = %16.15f \n",it_mg,resid);
}

void source_uv(double **tu,double **tv,double **divuv,int nxt,int nyt){
    int i,j;
    div_uv(tu,tv,divuv,nxt,nyt); ijloopt {divuv[i][j] /= dt;}
}

void augmenuv(double **u,double **v,int nx,int ny){   
    int i,j;
    double bdvel=1.0;
    iloop {u[i][0] = -u[i][1]; u[i][ny+1] = 2.0*bdvel-u[i][ny];}
    jloop {v[0][j] = -v[1][j]; v[nx+1][j] = -v[nx][j];}
}

void Poisson(double **tu,double **tv,double **p){
    source_uv(tu,tv,workp,nx,ny); MG_Poisson(p,workp); 
}

void temp_uv(double **tu,double **tv,double **u,double **v,
        double **adv_u,double **adv_v){
            int i,j;
    for (i=1; i<nx; i++) {jloop {tu[i][j] = u[i][j]+dt*( (u[i+1][j]
    +u[i-1][j]-4.0*u[i][j]+u[i][j+1]+u[i][j-1])/(Re*h2)-adv_u[i][j]);}}
    iloop {for (j=1; j<ny; j++) {tv[i][j] = v[i][j]+dt*( (v[i+1][j]
    +v[i-1][j]-4.0*v[i][j]+v[i][j+1]+v[i][j-1])/(Re*h2)-adv_v[i][j]);}}
}

void advection_uv(double **u,double **v,double **adv_u,double **adv_v){
    int i,j;
    augmenuv(u,v,nx,ny);
    for (i=1; i<nx; i++) { jloop {
    if (u[i][j]>0.0) {adv_u[i][j] = u[i][j]*(u[i][j]-u[i-1][j])/h;}
    else {adv_u[i][j] = u[i][j]*(u[i+1][j]-u[i][j])/h;}            
    if (v[i][j-1]+v[i+1][j-1]+v[i][j]+v[i+1][j]>0.0)
    adv_u[i][j] += 0.25*(v[i][j-1]+v[i+1][j-1]+v[i][j]+v[i+1][j])
                *(u[i][j]-u[i][j-1])/h;
    else adv_u[i][j] += 0.25*(v[i][j-1]+v[i+1][j-1]+v[i][j]+v[i+1][j])
                *(u[i][j+1]-u[i][j])/h;}}
    iloop { for (j=1; j<ny; j++) {
    if (u[i-1][j]+u[i][j]+u[i-1][j+1]+u[i][j+1]>0.0)
    adv_v[i][j] = 0.25*(u[i-1][j]+u[i][j]+u[i-1][j+1]+u[i][j+1])
                *(v[i][j]-v[i-1][j])/h;
    else adv_v[i][j] = 0.25*(u[i-1][j]+u[i][j]+u[i-1][j+1]+u[i][j+1])
                *(v[i+1][j]-v[i][j])/h;
    if (v[i][j]>0.0) {adv_v[i][j] += v[i][j]*(v[i][j]-v[i][j-1])/h;}
    else {adv_v[i][j] += v[i][j]*(v[i][j+1] - v[i][j])/h;}}}
}

void full_step(double **u,double **v,double **nu,double **nv,double **p){
    int i,j;
    advection_uv(u,v,adv_u,adv_v);temp_uv(tu,tv,u,v,adv_u,adv_v);
    Poisson(tu,tv,p); grad_p(p,worku,workv,nx,ny);    
    for (i=1;i<nx;i++){jloop {nu[i][j]=tu[i][j]-dt*worku[i][j];}}
    iloop {for (j=1;j<ny;j++){nv[i][j]=tv[i][j]-dt*workv[i][j];}}
}

void print_data1(double **u,double **v,double **p){
    int i,j;
    FILE *fu,*fv,*fp;
		fu=fopen(bufferu,"a");fv=fopen(bufferv,"a");fp=fopen(bufferp,"a");
		iloop {jloop {fprintf(fu,"  %16.14f",0.5*(u[i][j]+u[i-1][j]));
		fprintf(fv,"  %16.14f",0.5*(v[i][j]+v[i][j-1]));
		fprintf(fp,"  %16.14f",p[i][j]);} fprintf(fu,"\n");fprintf(fv,"\n"); 
		fprintf(fp,"\n");}fclose(fu); fclose(fv); fclose(fp);
}

void initialization(double **p,double **u,double **v){   
    int i,j;
    ijloop p[i][j] = 0.0; ij0loop v[i][j] = 0.0;
    i0jloop u[i][j] = 0.0;
}

int main(){
    int it,max_it,ns,count = 1;
    double **u,**v,**nu,**nv,**p;
    FILE *fu,*fv,*fp;
    p_relax=5;nx=gnx;ny=gny;n_level=(int)(log(ny)/log(2)-0.9);
    xleft=0.0;xright=1.0;yleft=0.0;yright=1.0*gny/gnx*xright;
    h = (xright-xleft)/ (double)nx; h2 = pow(h,2);
    max_it=1;ns=(int)(max_it/10+0.001);Re=100.0;dt=0.1*h*h*Re;
    p = dmatrix(1,nx,1,ny); sor = dmatrix(1,nx,1,ny);
    workp=dmatrix(0,nx+1,0,ny+1);worku=dmatrix(0,nx+1,0,ny+1);
    workv=dmatrix(0,nx+1,0,ny+1);u=dmatrix(-1,nx+1,0,ny+1);
    v = dmatrix(0,nx+1,-1,ny+1);nu = dmatrix(-1,nx+1,0,ny+1); 
    nv = dmatrix(0,nx+1,-1,ny+1);tu = dmatrix(0,nx,1,ny); 
    tv = dmatrix(1,nx,0,ny);adv_u = dmatrix(0,nx,1,ny); 
    adv_v = dmatrix(1,nx,0,ny);zero_matrix(tu,0,nx,1,ny); 
    zero_matrix(tv,1,nx,0,ny);sprintf(bufferu,"u.m"); 
    sprintf(bufferv,"v.m"); sprintf(bufferp,"p.m");
    fu = fopen(bufferu,"w"); fv = fopen(bufferv,"w"); 
    fp = fopen(bufferp,"w");fclose(fu);fclose(fv);fclose(fp);
    initialization(p,u,v);  print_data1(u,v,p);
    mat_copy(nu,u,0,nx,1,ny); mat_copy(nv,v,1,nx,0,ny);
    for (it=1; it<=max_it; it++) {printf("iteration = %d\n",it);
        full_step(u,v,nu,nv,p); mat_copy(u,nu,0,nx,1,ny);
        mat_copy(v,nv,1,nx,0,ny);
        if (it % ns==0) {print_data1(nu,nv,p);
      printf("print out counts %d \n",count);count++;}}
    return 0;
}
