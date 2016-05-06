#include <cmath>
#include<iostream>
#include<ctime>
#include <Imagine/Graphics.h>
using namespace Imagine;


// CAS 1D AVEC BATIMETRIE NON CONSTANTE
// VOIR DOC RUPTURE DE BARRAGE POUR LE CAS TEST


const int L = 100;
const int N= 100;                      // compiler avec N pair (nb de segments)
const float dx = float(L)/float(N);          //(il y a N+1 points)
const float H=8.;
const float h0=H/2.;
const int T=1000; //(nombre de pas de temps)

const float g = 9.81;

const int x0=N/2;

void affiche(float U[N+1][2],float z[N+1]){
    int largeur = 500/N;
    int saut,n;
    if (largeur == 0){
        largeur=1;
        n=500;
        int saut = N/n;
    }
    else{
        saut=1;
        n=N;
    }
    noRefreshBegin();
    clearWindow();

    for (int i=0; i<n;i++){
        fillRect(i*largeur,500-z[i]*30,largeur,-U[i*saut][0]*30,BLUE);
    }
    for(int j=0;j<N;++j){
        fillRect(j*largeur,500,largeur,-z[j]*30,RED);
    }
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




void calcul_flux(int i,float F[N+1][2],const float Ug[2],const float Ud[2]){
    float vp1,vp2;
    vp1=maxi(abs(Ug[1]/Ug[0]+sqrt(g*Ug[0])),abs(Ud[1]/Ud[0]+sqrt(g*Ud[0])));
    vp2=maxi(abs(Ug[1]/Ug[0]-sqrt(g*Ug[0])),abs(Ud[1]/Ud[0]-sqrt(g*Ud[0])));
    float vp=maxi(vp1,vp2);
    F[i][0]=(Ug[1]+Ud[1])/2 - vp* (Ud[0]-Ug[0])/2;
    float FUg1=Ug[1]*Ug[1]/Ug[0]+g*Ug[0]*Ug[0]/2;
    float FUd1=Ud[1]*Ud[1]/Ud[0]+g*Ud[0]*Ud[0]/2;
    F[i][1]=(FUg1+FUd1)/2 - vp* (Ud[1]-Ug[1])/2;
}


void transfert(float U[N+1][2],const float Ucopie[N+1][2]){
    for (int i=0;i<N+1;++i){
        U[i][0]=Ucopie[i][0];
        U[i][1]=Ucopie[i][1];
    }
}

int main(){
    Imagine::openWindow(500,500);

    float U[N+1][2];
    float Ucopie[N+1][2];
    float Upp[2];
    float Upm[2];
    float Ump[2];
    float Umm[2];
    float S[N+1];
    float Fp[N+1][2];
    float Fm[N+1][2];


    float z[N];

    // BEGIN  CONDITIONS INITALES

/*

    // begin batimétrie triangulaire
    for (int i=0; i< x0;i++){

        z[i]=-(H-h0)*i/x0+H;
    }
    for (int i=x0; i< N+1;i++){

        z[i]=(H-h0)*(i-x0)/(N-x0)+h0;
    }
    for (int i=0;i<N+1;++i){
        U[i][0]=H+0.5-z[i];
        U[i][1]=0.;
    }
    U[N/2][0]=2*H;
    // end batimétrie triangulaire
*/


    // begin batimétrie à bosse
    for (int i=0; i< 40;i++){
        z[i]=3.;
        U[i][0]=H+0.5-z[i];
        U[i][1]=0.;

    }
    for (int i=40; i< 61;i++){
        z[i]=3.+(3./100.)*float(i-40.)*float(i-60.);
        U[i][0]=H+0.5 +abs(i-50);

        U[i][1]=0.;

    }
    for (int i=61; i< N+1;i++){
        z[i]=3.;
        U[i][0]=H+0.5-z[i];
        U[i][1]=0.;

    }
    // end batimétrie à bosse


    for (int i=0;i<N+1;++i){
        std::cout << z[i] << std::endl;
    }
/*

    for (int i=0; i< N+1;i++){
        U[i][0]=H;
        U[i][1]=0.;
        z[i]=0.;
    }
*/


    // END  CONDITIONS INITALES

    float* lamb=new float[N+1];

    for (int t=0; t<T-1; t++){

        //BEGIN choisir un pas de temps permettant de vérifier la condition de stabilité de Courant-Friedrichs-Lewy

        for (int i=0;i<N+1;i++){
            float u=U[i][1]/U[i][0];
            float vp1=u+sqrt(g*U[i][0]);
            float vp2=u-sqrt(g*U[i][0]);

            if (abs(vp1)>abs(vp2))
                lamb[i]=abs(vp1);
            else
                lamb[i]=abs(vp2);
        }
        float max=lamb[0];
        for (int i=0;i<N+1;i++){
            if (lamb[i]>max){
                max=lamb[i];
            }
        }
        float dt=dx*0.2/max ;
        //END choisir un pas de temps permettant de vérifier la condition de stabilité de Courant-Friedrichs-Lewy

        for (int i=0;i<N+1;i++){  // boucle en temps
            if (i!=0 && i!= N){
                //
                // calcul des flux en i+1/2
                //
                float hbarpp=U[i+1][0] +z[i+1]-maxi(z[i],z[i+1]);
                float hpp=maxi(0,hbarpp);

                Upp[0]=hpp;
                Upp[1]=hpp*(U[i+1][1]/U[i+1][0]);

                float hbarpm=U[i][0] +z[i]-maxi(z[i],z[i+1]);
                float hpm=maxi(0,hbarpm);


                Upm[0]=hpm;
                Upm[1]=hpm*(U[i][1]/U[i][0]);


                calcul_flux(i,Fp,Upm,Upp);
                calcul_flux(i-1,Fp,Upm,Upp);


                //
                // calcul des flux en i-1/2
                //

                float hbarmp=U[i][0] +z[i]-maxi(z[i-1],z[i]);
                float hmp=maxi(0,hbarmp);
                Ump[0]=hmp;
                Ump[1]=hmp*(U[i][1]/U[i][0]);

                float hbarmm=U[i-1][0] +z[i-1]-maxi(z[i-1],z[i]);
                float hmm=maxi(0,hbarmm);
                Umm[0]=hmm;
                Umm[1]=hmm*(U[i-1][1]/U[i-1][0]);

                calcul_flux(i,Fm,Umm,Ump);
                calcul_flux(i-1,Fm,Umm,Ump);


                //
                // calcul du terme source
                //

                S[i]=(g/2)*(hpm*hpm - hmp*hmp);
                Ucopie[i][0]=U[i][0]-(dt/dx)*(Fp[i][0]-Fm[i][0]);
                Ucopie[i][1]=U[i][1]-(dt/dx)*(Fp[i-1][1]-Fm[i-1][1])+(dt/dx)*S[i];
            }
            if (i==0){
                Ucopie[i][0]=U[N-1][0];
                Ucopie[i][1]=-U[1][1];
            }
            if (i==N){
                Ucopie[i][0]=U[1][0];
                Ucopie[i][1]=-U[N-1][1];
            }
        }
        transfert(U,Ucopie);
        affiche(U,z);
        milliSleep(50);
    }
    endGraphics();
    return 0;
}
