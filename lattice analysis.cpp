#include<stdio.h>
#include<stdlib.h>
#include<math.h>



void resolution(int n,float a[200][200],float b[200],float depr[200]);
void produit(int m,float a[200][200],float b[200],float p[200]);
main(){

int n_elet,n_nodes,i,j,k;
const int dimax=200;
int nodes[dimax][2];
float x_cor[dimax],y_cor[dimax],area[dimax],eyoung[dimax],kg[dimax][dimax]={0};
char fichin[20],fichout[20];

FILE *infile,*outfile;

printf("donner fichier data : ");
scanf("%s",&fichin);
infile=fopen(fichin,"r");

printf("donner fichier resultat : ");
scanf("%s",&fichout);
outfile=fopen(fichout,"w");

fscanf(infile,"%d %d",&n_elet,&n_nodes);
fprintf(outfile,"nombre elements : %d \nnombre nodes : %d\n",n_elet,n_nodes);

for(i=0;i<n_nodes;i++) fscanf(infile,"%f%f",&x_cor[i],&y_cor[i]);

for(i=0;i<n_nodes;i++) fprintf(outfile,"%d\t%f\t%f\n",i+1,x_cor[i],y_cor[i]);

for(i=0;i<n_elet;i++) fscanf(infile,"%d%d%e%e",&nodes[i][0],&nodes[i][1],&eyoung[i],&area[i]);

fprintf(outfile,"\n");


for(i=0;i<n_elet;i++) fprintf(outfile,"element %d : %d \t %d \t %e \t %e\n",i+1,nodes[i][0],nodes[i][1],eyoung[i],area[i]);

fprintf(outfile,"\n");


for(i=0;i<n_elet;i++){
int ig,jg;
int pg[4];


    int node1=nodes[i][0]-1;
    int node2=nodes[i][1]-1;

    float dx =x_cor[node2]-x_cor[node1];
    float dy =y_cor[node2]-y_cor[node1];

    float l= sqrt(dx*dx+dy*dy);
    float c=dx/l;
    float s=dy/l;
    float EAL=area[i]*eyoung[i]/l;

    float matrix[4][4];


    matrix[0][0]=c*c;
    matrix[0][1]=c*s;
    matrix[0][2]=-c*c;
    matrix[0][3]=-c*s;
    matrix[1][0]=matrix[0][1];
    matrix[1][1]=s*s;
    matrix[1][2]=-c*s;
    matrix[1][3]=-s*s;
    matrix[2][0]=matrix[0][2];
    matrix[2][1]=matrix[1][2];
    matrix[2][2]=c*c;
    matrix[2][3]=c*s;
    matrix[3][0]=matrix[0][3];
    matrix[3][1]=matrix[1][3];
    matrix[3][2]=matrix[2][3];
    matrix[3][3]=s*s;

    for(j=0;j<4;j++){
            for(k=0;k<4;k++){fprintf(outfile,"%e\t",EAL*matrix[j][k]);}
            fprintf(outfile,"\n");
    }



fprintf(outfile,"\n\n\n");



pg[0]=2*node1;
pg[1]=2*node1+1;
pg[2]=2*node2;
pg[3]=2*node2+1;

 for(j=0;j<4;j++)
    {
        ig=pg[j];
        for(k=0;k<4;k++)
        {
            jg=pg[k];
            kg[ig][jg]+=EAL*matrix[j][k];
        }

    }

}



for(i=0;i<16;i++){
            for(j=0;j<16;j++){fprintf(outfile,"%e\t",kg[i][j]);}
            fprintf(outfile,"\n");
    }
fprintf(outfile,"\n");

int nbrdep,nbrload;
int nodedep[dimax],nodeload[dimax],id[dimax][2];
float fxy[dimax][2];

for(i=0;i<n_nodes;i++){
    for(j=0;j<2;j++){
        id[i][j]=0;
    }
}
for(i=0;i<n_elet;i++){
    for(j=0;j<2;j++){
        fxy[i][j]=0;
    }
}

fscanf(infile,"%d",&nbrdep);
for(i=0;i<nbrdep;i++){
    fscanf(infile,"%d",&nodedep[i]);
}

for(i=0;i<nbrdep;i++){
    for(j=0;j<2;j++){
        fscanf(infile,"%d",&id[nodedep[i]-1][j]);
    }
}

fprintf(outfile,"le nombre des noeuds avec condition limite de deplacement : %d\n",nbrdep);
fprintf(outfile,"les noeuds avec condition limite de deplacement : \n");
for(i=0;i<nbrdep;i++){
    fprintf(outfile,"%d\t",nodedep[i]);
}
fprintf(outfile,"\n la matrice id : \n");
for(i=0;i<n_nodes;i++){
    for(j=0;j<2;j++){
        fprintf(outfile,"%d\t",id[i][j]);
    }
    fprintf(outfile," \n");
}



fscanf(infile,"%d",&nbrload);
for(i=0;i<nbrload;i++){
    fscanf(infile,"%d",&nodeload[i]);
}

for(i=0;i<nbrload;i++){
    for(j=0;j<2;j++){
        fscanf(infile,"%f",&fxy[nodeload[i]-1][j]);
    }
}


fprintf(outfile,"\n");
fprintf(outfile,"le nombre des noeuds avec une charge : %d\n",nbrload);
fprintf(outfile,"les noeuds avec une charge : \n");
for(i=0;i<nbrload;i++){
    fprintf(outfile,"%d\t",nodeload[i]);
}

fprintf(outfile,"\n la matrice Fxy : \n");
for(i=0;i<n_elet;i++){
    for(j=0;j<2;j++){
        fprintf(outfile,"%e \t",fxy[i][j]);
    }
    fprintf(outfile," \n");
}

int dimkr=0;

for(i=0;i<n_nodes;i++){
    for(j=0;j<2;j++){
        if(id[i][j]==1)
            id[i][j]=0;
        else{
            dimkr++;
            id[i][j]=dimkr;
        }
    }
}

fprintf(outfile,"\n la matrice id modifie : \n");
for(i=0;i<n_nodes;i++){
    for(j=0;j<2;j++){
        fprintf(outfile,"%d\t",id[i][j]);
    }
    fprintf(outfile," \n");
}

int ir,jr;
float kr[dimax][dimax],load[dimax];


    for(int i=0;i<n_nodes;i++)
        for(int ii=0;ii<2;ii++)
        if (id[i][ii]!=0)
            {
            int ir;
            ir=id[i][ii]-1;
            load[ir]=fxy[i][ii];
            for(int j=0;j<n_nodes;j++)
                for(int jj=0;jj<2;jj++)
            if(id[j][jj]!=0)
                {
                int jr;
                jr=id[j][jj]-1;
                kr[ir][jr]=kg[2*i+ii][2*j+jj];
                }
            }

fprintf(outfile,"\n la matrice reduite : \n");
for(i=0;i<dimkr;i++){
    for(j=0;j<dimkr;j++){
        fprintf(outfile,"%e\t",kr[i][j]);
    }
    fprintf(outfile," \n");
}
fprintf(outfile," \n");

fprintf(outfile,"\n le vecteur des forces : \n");
for(i=0;i<dimkr;i++){
    fprintf(outfile,"%e\n",load[i]);
}

fprintf(outfile," \n");




float depr[dimax],depg[dimax];

resolution(dimkr,kr,load,depr);

 fprintf(outfile,"les deplacements sont : \n");
 for(i=0;i<dimkr;i++){
        fprintf(outfile,"%e\t",depr[i]);
        fprintf(outfile,"\n");
    }


for(i=0;i<n_nodes;i++)
    for(j=0;j<2;j++)
        if(id[i][j]==0)
            depg[2*i+j]=0;
        else depg[2*i+j]=depr[id[i][j]-1];

fprintf(outfile,"\n les deplacements globales sont : \n");
 for(i=0;i<n_nodes*2;i++){
        fprintf(outfile,"%e\t",depg[i]);
        fprintf(outfile,"\n");
    }

float R[dimax];
int dimkg=n_nodes*2;
produit(dimkg,kg,depg,R);

fprintf(outfile,"\n les reactions sont : \n");
 for(i=0;i<dimkg;i++){
        fprintf(outfile,"%e\t",R[i]);
        fprintf(outfile,"\n");
    }


fprintf(outfile,"\n \n");
fprintf(outfile,"La barre \t noeuds \t L'effort axial \t  la contrainte \t la deformation (10^3) \n");
fprintf(outfile,"\n");

for(int i=0; i<n_elet ; i++)
{
    int node1=nodes[i][0]-1;
    int node2=nodes[i][1]-1;
    float dx =x_cor[node2]-x_cor[node1];
    float dy =y_cor[node2]-y_cor[node1];
    float l= sqrt(dx*dx+dy*dy);
    float c=dx/l;
    float s=dy/l;
    float EAL=area[i]*eyoung[i]/l;
    float dU=depg[2*(node2)]-depg[2*(node1)];
    float dV=depg[2*(node2)+1]-depg[2*(node1)+1];
    float N=EAL*(dU*c+dV*s);
    float sigma=N/area[i];
    float epsilon=(sigma/eyoung[i])*1000;
    fprintf(outfile,"   %d \t\t  (%d,%d)  \t   %.2f \t\t   %.2f \t          %f\n",i+1, node1+1, node2+1, N, sigma, epsilon);
}
}




void resolution(int n, float a[200][200], float b[200], float x[200])
{
    float sum, r;
    for(int j=0 ; j<n-1 ; j++)
    {
        for(int i=j+1 ; i<n ; i++)
        {
            r=a[i][j]/a[j][j];
            for(int k=j ; k<n ; k++)
            {
                a[i][k]=a[i][k]-r*a[j][k];
            }
            b[i]=b[i]-r*b[j];

        }
    }

    for(int i=n-1 ; i>=0 ; i--)
    {
        sum=0;
        for(int j=i; j<n ; j++)
        {
            sum=sum+a[i][j]*x[j];
        }
        x[i]=(b[i]-sum)/a[i][i];
    }
}


void produit (int m , float M1[200][200] , float M2[200] , float P[200])
{
    int i,k;


    for(i=0;i<2*m;i++ )
    {
        P[i]=0;
            for(k=0;k<2*m;k++)
                P[i] = P[i] + M1[i][k]*M2[k];
    }
}
