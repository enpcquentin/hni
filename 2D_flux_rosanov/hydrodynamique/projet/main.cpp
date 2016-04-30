#include <cmath>
#include<iostream>
#include<ctime>
#include <Imagine/Graphics.h>
using namespace Imagine;


// CAS 2D AVEC BATIMETRIE NON CONSTANTE


const int L = 100;
const int N= 50;                      // compiler avec N pair (nb de segments)
const float dx = float(L)/float(N);          //(il y a N+1 points)
const float H=0.05;
const float h0=H/2.;
const int T=1000; //(nombre de pas de temps)

const float g = 9.81;

const int x0=N/2;

void affiche(float U[N+1][N+1][3],float z[N+1][N+1]){
    Triangle dT[0];

    DoublePoint3* fp=new DoublePoint3 [(N+1)*(N+1)];
    for (int i=0;i<N+1;++i){
        for (int j=0;j<N+1;++j){
            fp[i+(N+1)*j]=DoublePoint3(i,j,30*U[i][j][0]);
        }
    }
    Quad* fq=new Quad[N*N];
    for (int i=0;i<N;++i){
        for (int j=0;j<N;++j){
            fq[i+N*j]=Quad(i+N*j,i+1+N*j,i+N*(j+1),i+1+N*(j+1));
        }
    }
    Mesh fM(fp,(N+1)*(N+1),dT,0,fq,N*N);
    fM.setColor(BLUE);
    noRefreshBegin();
    showMesh(fM);


    // batimétrie
    DoublePoint3* bp=new DoublePoint3 [(N+1)*(N+1)];
    for (int i=0;i<N+1;++i){
        for (int j=0;j<N+1;++j){
            bp[i+(N+1)*j]=DoublePoint3(i,j,z[i][j]);
        }
    }
    Quad* bq=new Quad[N*N];
    for (int i=0;i<N;++i){
        for (int j=0;j<N;++j){
            bq[i+N*j]=Quad(i+N*j,i+1+N*j,i+N*(j+1),i+1+N*(j+1));
        }
    }
    Mesh bM(bp,(N+1)*(N+1),dT,0,bq,N*N);
    bM.setColor(YELLOW);
    showMesh(bM);
    /*
    DoublePoint3 pos(-100,-100,180);
    DoubleVector3 dir (-10,-10,10);
    DoubleVector3 up(0,0,1);
    */
    DoublePoint3 pos(-10,-10,18);
    DoubleVector3 dir (-20,-20,20);
    DoubleVector3 up(0,0,1);
    setCamera(pos,dir,up);
    noRefreshEnd();

}





float val_propres(float ug,float ud,float hg,float hd){
    float lambda1=ug+sqrt(g*hg);
    float lambda2=ud-sqrt(g*hd);
    if (abs(lambda1)<abs(lambda2)){
        return abs(lambda2);
    }
    else
        return abs(lambda1);
}

float maxi(float a,float b){
    if (a>b)
        return a;
    else
        return b;
}




void calcul_flux_x(float G[3],const float Ug[3],const float Ud[3]){
    float vp1,vp2;
    // la troisième valeur propre ne nous intéresse pas car elle est comprise entre les deux valeurs propres ici considérées (et seul le max nous intéresse)
    vp1=maxi(abs(Ug[1]/Ug[0]+sqrt(g*Ug[0])),abs(Ud[1]/Ud[0]+sqrt(g*Ud[0])));
    vp2=maxi(abs(Ug[1]/Ug[0]-sqrt(g*Ug[0])),abs(Ud[1]/Ud[0]-sqrt(g*Ud[0])));
    float vp=maxi(vp1,vp2);
    G[0]=(Ug[1]+Ud[1])/2 - vp* (Ud[0]-Ug[0])/2;
    float GUg1=Ug[1]*Ug[1]/Ug[0]+g*Ug[0]*Ug[0]/2;
    float GUd1=Ud[1]*Ud[1]/Ud[0]+g*Ud[0]*Ud[0]/2;
    G[1]=(GUg1+GUd1)/2 - vp* (Ud[1]-Ug[1])/2;
    if (Ug[1]/Ug[0] + Ud[1]/Ud[0] >0){
        G[2]=(Ug[2]/Ug[0])*G[0];
    }
    else
        G[2]=(Ud[2]/Ud[0])*G[0];

}


void calcul_flux_y(float H[3],const float Ug[3],const float Ud[3]){
    float vp1,vp2;
    // la troisième valeur propre ne nous intéresse pas car elle est comprise entre les deux valeurs propres ici considérées (et seul le max nous intéresse)
    vp1=maxi(abs(Ug[1]/Ug[0]+sqrt(g*Ug[0])),abs(Ud[1]/Ud[0]+sqrt(g*Ud[0])));
    vp2=maxi(abs(Ug[1]/Ug[0]-sqrt(g*Ug[0])),abs(Ud[1]/Ud[0]-sqrt(g*Ud[0])));
    float vp=maxi(vp1,vp2);
    H[0]=(Ug[2]+Ud[2])/2 - vp* (Ud[0]-Ug[0])/2;
    //std::cout << H[0] << "  " << Ug[2] << "  " << Ud[2] << "  " << Ud[0] << "  " << Ug[0] << std::endl;
    float HUg2=Ug[2]*Ug[2]/Ug[0]+g*Ug[0]*Ug[0]/2;
    float HUd2=Ud[2]*Ud[2]/Ud[0]+g*Ud[0]*Ud[0]/2;
    H[2]=(HUg2+HUd2)/2 - vp* (Ud[1]-Ug[1])/2;
    if (Ug[2]/Ug[0] + Ud[2]/Ud[0] >0){
        H[1]=(Ug[1]/Ug[0])*H[0];
    }
    else
        H[1]=(Ud[1]/Ud[0])*H[0];
}


void transfert(float U[N+1][N+1][3],const float Ucopie[N+1][N+1][3]){
    for (int i=0;i<N+1;++i){
        for (int j=0;j<N+1;++j){
            U[i][j][0]=Ucopie[i][j][0];
            U[i][j][1]=Ucopie[i][j][1];
            U[i][j][2]=Ucopie[i][j][2];
        }
    }
}

int main(){
    Imagine::openWindow3D(500,500, "3D");

    float U[N+1][N+1][3];
    float Ucopie[N+1][N+1][3];
    float Upp[3];
    float Upm[3];
    float Ump[3];
    float Umm[3];
    float S;
    float Gc;
    float Hc;
    float Gp[3]; // représente G(i+1/2 , j) à j fixé.
    float Gm[3];
    float Hp[3];
    float Hm[3];
    float z[N+1][N+1];

    // BEGIN  CONDITIONS INITALES

    for (int i=0; i< N+1;i++){
        for (int j=0; j<N+1;++j){
            U[i][j][0]=H;
            U[i][j][1]=0.;
            U[i][j][2]=0.;
            z[i][j]=0.;
        }
    }
    U[N/2][N/2][0]=H+0.1;


    // END  CONDITIONS INITALES

    float* lamb=new float[(N+1)*(N+1)];

    for (int t=0; t<T-1; t++){

        //BEGIN choisir un pas de temps permettant de vérifier la condition de stabilité de Courant-Friedrichs-Lewy

        for (int i=0;i<N+1;i++){
            for (int j=0;j<N+1;++j){
                float u=U[i][j][1]/U[i][j][0];
                float vp1=u+sqrt(g*U[i][j][0]);
                float vp2=u-sqrt(g*U[i][j][0]);

                if (abs(vp1)>abs(vp2))
                    lamb[i+(N+1)*j]=abs(vp1);
                else
                    lamb[i+(N+1)*j]=abs(vp2);
            }
        }
        float max=lamb[0];
        for (int i=0;i<N+1;i++){
            for (int j=0;j<N+1;++j){
                if (lamb[i+j*(N+1)]>max){
                    max=lamb[i+j*(N+1)];
                }
            }
        }

        float dt=dx*0.2/max ;

        //END choisir un pas de temps permettant de vérifier la condition de stabilité de Courant-Friedrichs-Lewy

        for (int j=0;j<N+1;++j){
            //i=0
            Ucopie[0][j][0]=U[N-1][j][0];
            Ucopie[0][j][1]=U[N-1][j][1];
            Ucopie[0][j][2]=U[N-1][j][2];
            //i=N
            Ucopie[N][j][0]=U[1][j][0];
            Ucopie[N][j][1]=U[1][j][1];
            Ucopie[N][j][2]=U[1][j][2];
        }
        for (int i=0;i<N+1;++i){
            //j=0
            Ucopie[i][0][0]=U[i][N-1][0];
            Ucopie[i][0][1]=U[i][N-1][1];
            Ucopie[i][0][2]=U[i][N-1][2];
            //j=N
            Ucopie[i][N][0]=U[i][1][0];
            Ucopie[i][N][1]=U[i][1][1];
            Ucopie[i][N][2]=U[i][1][2];
        }
        for (int k=1;k<2*N-2;++k){
            int i,j;
            if (k<N){
                i=1;
                j=k;
            }
            else{
                i=k-(N-2);
                j=1;
            }
            int compteur=0;
            while ((i!=N) && (j!=N)){
                if (compteur % 2==0){ // un pas selon x

                    //
                    // calcul des flux en i+1/2
                    //
                    float hbarpp=U[i+1][j][0] +z[i+1][j]-maxi(z[i][j],z[i+1][j]);
                    float hpp=maxi(0,hbarpp);

                    Upp[0]=hpp; // Upp= U_(i+1/2 + , j)
                    Upp[1]=hpp*(U[i+1][j][1]/U[i+1][j][0]);
                    Upp[2]=hpp*(U[i+1][j][2]/U[i+1][j][0]);

                    float hbarpm=U[i][j][0] +z[i][j]-maxi(z[i][j],z[i+1][j]);
                    float hpm=maxi(0,hbarpm);

                    Upm[0]=hpm; // Upm= U_(i+1/2 -, j)
                    Upm[1]=hpm*(U[i][j][1]/U[i][j][0]);
                    Upm[2]=hpm*(U[i][j][2]/U[i][j][0]);

                    calcul_flux_x(Gp,Upm,Upp);

                    //
                    // calcul des flux en i-1/2
                    //

                    float hbarmp=U[i][j][0] +z[i][j]-maxi(z[i-1][j],z[i][j]);
                    float hmp=maxi(0,hbarmp);
                    Ump[0]=hmp; // Ump= U_(i-1/2 + , j)
                    Ump[1]=hmp*(U[i][j][1]/U[i][j][0]);
                    Ump[2]=hmp*(U[i][j][2]/U[i][j][0]);

                    float hbarmm=U[i-1][j][0] +z[i-1][j]-maxi(z[i-1][j],z[i][j]);
                    float hmm=maxi(0,hbarmm);
                    Umm[0]=hmm; // Umm= U_(i-1/2 - , j)
                    Umm[1]=hmm*(U[i-1][j][1]/U[i-1][j][0]);
                    Umm[2]=hmm*(U[i-1][j][2]/U[i-1][j][0]);

                    calcul_flux_x(Gm,Umm,Ump);

                    Gc=(-g/2)*(hmp+hpm)*(maxi(z[i+1][j],z[i][j])-maxi(z[i-1][j],z[i][j]));

                    //
                    // calcul du terme source
                    //

                    S=(g/2)*(hpm*hpm - hmp*hmp);

                    Ucopie[i][j][0]=U[i][j][0]-(dt/dx)*(Gp[0]-Gm[0]);
                    //std::cout << t <<" U copie" << Ucopie[i][j][0] << "   " << Ucopie[i][j][1] << " " << Ucopie[i][j][2]<< "  Gp - gm " << Gp[0]-Gm[0] << std::endl;
                    Ucopie[i][j][1]=U[i][j][1]-(dt/dx)*(Gp[1]-Gm[1])+(dt/dx)*Gc+(dt/dx)*S;
                    Ucopie[i][j][2]=U[i][j][2]-(dt/dx)*(Gp[2]-Gm[2]);

                    ++i;
                }
                else {  // un pas selon y

                    //
                    // calcul des flux en j+1/2
                    //
                    float hbarpp=U[i][j+1][0] +z[i][j+1]-maxi(z[i][j],z[i][j+1]);
                    float hpp=maxi(0,hbarpp);

                    Upp[0]=hpp; // Upp= U_(i, j+1/2 + )
                    Upp[1]=hpp*(U[i][j+1][1]/U[i][j+1][0]);
                    Upp[2]=hpp*(U[i][j+1][2]/U[i][j+1][0]);


                    float hbarpm=U[i][j][0] +z[i][j]-maxi(z[i][j],z[i][j+1]);
                    float hpm=maxi(0,hbarpm);

                    Upm[0]=hpm; // Upm= U_(i+1/2 -, j)
                    Upm[1]=hpm*(U[i][j][1]/U[i][j][0]);
                    Upm[2]=hpm*(U[i][j][2]/U[i][j][0]);

                    calcul_flux_y(Hp,Upm,Upp);

                    //
                    // calcul des flux en i-1/2
                    //

                    float hbarmp=U[i][j][0] +z[i][j]-maxi(z[i][j-1],z[i][j]);
                    float hmp=maxi(0,hbarmp);
                    Ump[0]=hmp; // Ump= U_(i,j-1/2 + )
                    Ump[1]=hmp*(U[i][j][1]/U[i][j][0]);
                    Ump[2]=hmp*(U[i][j][2]/U[i][j][0]);

                    float hbarmm=U[i][j-1][0] +z[i][j-1]-maxi(z[i][j-1],z[i][j]);
                    float hmm=maxi(0,hbarmm);
                    Umm[0]=hmm; // Umm= U_(i-1/2 - , j)
                    Umm[1]=hmm*(U[i][j-1][1]/U[i][j-1][0]);
                    Umm[2]=hmm*(U[i][j-1][2]/U[i][j-1][0]);

                    calcul_flux_y(Hm,Umm,Ump);

                    Hc=(-g/2)*(hmp+hpm)*(maxi(z[i][j+1],z[i][j])-maxi(z[i][j-1],z[i][j]));

                    //
                    // calcul du terme source
                    //

                    S=(g/2)*(hpm*hpm - hmp*hmp);

                    Ucopie[i][j][0]=U[i][j][0]-(dt/dx)*(Hp[0]-Hm[0]);
                    Ucopie[i][j][1]=U[i][j][1]-(dt/dx)*(Hp[1]-Hm[1]);
                    Ucopie[i][j][2]=U[i][j][2]-(dt/dx)*(Hp[2]-Hm[2])+(dt/dx)*Hc+(dt/dx)*S;
                    std::cout << t << " " << i << "  " << j <<" U copie" << Ucopie[i][j][0]  << "   " << Ucopie[i][j][1] << " " << Ucopie[i][j][2]<< "  hp - hm " << Hp[0]-Hm[0] << std::endl;

                    ++j;
                }
                ++compteur;
            }
        }
        transfert(U,Ucopie);
        affiche(U,z);
        milliSleep(100);
    }

    endGraphics();
    return 0;
}

















