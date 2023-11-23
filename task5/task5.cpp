#include <iostream>
#include <cmath>
#include <fstream>

#include "../include/tire_constants.h" 
#include "../include/car_constants.h" 

using namespace std;
using namespace TireConstants;
using namespace CarConstants;

class Kesi {
public:
    float X=0.0f, Y=0.0f, theta=0.0f, v_x=0.0f, v_y=0.0f, r=0.0f, omega_f=0.0f,omega_r=0.0f, throttle=0.0f, steering_angle=0.0f, brakes=0.0f;
};

class States_dot {
public:
    float X_dot=0.0f, Y_dot=0.0f, phi_dot=0.0f,vx_dot=0.0f,vy_dot=0.0f,r_dot=0.0f,omega_dot_f=0.0f,omega_dot_r=0.0f; 

    States_dot operator*=(float scalar) {
        return States_dot{ 
            this->X_dot * scalar, 
            this->Y_dot * scalar, 
            this->phi_dot * scalar, 
            this->vx_dot * scalar, 
            this->vy_dot * scalar,
            this->r_dot * scalar,
            this->omega_dot_f * scalar,
            this->omega_dot_r * scalar };
    }

    States_dot operator+(const States_dot& other){
        States_dot result;
        result.X_dot = X_dot + other.X_dot;
        result.Y_dot = Y_dot + other.Y_dot;
        result.phi_dot = phi_dot + other.phi_dot;
        result.vx_dot = vx_dot + other.vx_dot;
        result.vy_dot = vy_dot + other.vy_dot;
        result.r_dot = r_dot + other.r_dot;
        result.omega_dot_f = omega_dot_f + other.omega_dot_f;
        result.omega_dot_r = omega_dot_r + other.omega_dot_r;
        return result;
    }

};

States_dot operator*(float scalar,States_dot other){
    return other*=scalar;
    }

States_dot operator*(States_dot other, float scalar) {
    return other *= scalar;
}

Kesi operator+(const Kesi& kesi, const States_dot& states_dot) {
    Kesi result;
    result.X = kesi.X + states_dot.X_dot;
    result.Y = kesi.Y + states_dot.Y_dot;
    result.theta = kesi.theta + states_dot.phi_dot;
    result.v_x = kesi.v_x + states_dot.vx_dot;
    result.v_y = kesi.v_y + states_dot.vy_dot;
    result.r = kesi.r + states_dot.r_dot;
    result.omega_f = kesi.omega_f+ states_dot.omega_dot_f;
    result.omega_r = kesi.omega_r+ states_dot.omega_dot_r;
    result.throttle = kesi.throttle; 
    result.steering_angle = kesi.steering_angle; 
    result.brakes = kesi.brakes; 
    return result;
}

class MagicTireModel
{
private:
    // float Fz0;
    float By,Cy,Dy,Ey,Ky,SHy,SVy,muy;
    float Bx,Cx,Dx,Ex,Kx,SHx,SVx,mux;
    float Byk,Cyk,Dyk,DVyk,Eyk,Gyk,SHyk,SVyk;
    float Bxa,Cxa,Dxa,Exa,Gxa,SHxa;
    float gamma,gammax,gammay,alphay,kappax,alpha_s,kappa_s;
    float Fz0,Fz,dfz;

    float sgn(float x){
        return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
    }
    
public:
    MagicTireModel(){
        Fz0 = FNOMIN;
        Fz = m*9.8/4;
        dfz = (Fz-Fz0*lFz0)/Fz0*lFz0;
        gamma = 0;

        gammay = gamma*lgammay;
        SHy = (pHy1+pHy2*dfz)*lHy + pHy3*gammay;
        SVy = Fz*((pVy1+pVy2*dfz)*lVy+(pVy3+pVy4*dfz)*gammay)*lmuy;
        Ky = pKy1*Fz0*sin(2*atan(Fz/(pKy2*Fz0*lFz0)))*(1-pKy3*abs(gammay))*lFz0*lKy;
        muy = (pDy1+pDy2*dfz)*(1-pDy3*gammay*gammay)*lmuy;
        Cy = pCy1*lCy;
        Dy = muy*Fz;
        By = Ky/(Cy*Dy);
        Ey = (pEy1+pEy2*dfz);

        gammax = gamma*lgammax;
        SHx = (pHx1+pHx2*dfz)*lHx;
        SVx = Fz*(pVx1+pVx2*dfz)*lVx*lmux;
        Kx = Fz*(pKx1+pKx2*dfz)*exp(pKx3*dfz)*lKx;
        mux = (pDx1+pDx2*dfz)*(1-pDx3*gammax*gammax)*lmux;
        Cx = pCx1*lCx;
        Dx = mux*Fz;
        Bx = Kx/(Cx*Dx);


    }
    
    float r_stat(){ 
        return R_tire - Fz/VERTICAL_STIFFNESS;
        }

    float r_eff(){
        return sin(acos(r_stat()/R_tire))/acos(r_stat()/R_tire)*R_tire;
    } 
    
    float solveFy(float alpha){
        alphay = alpha+SHy;
        float Fy0=Dy*sin(Cy*atan(By*alphay-Ey*(By*alphay-atan(By*alphay))))+SVy;
        return Fy0;
    }
    
    float solveFx(float kappa){
        kappax = kappa + SHx;
        Ex = (pEx1+pEx2*dfz+pEx3*dfz*dfz)*(1-pEx4*sgn(kappax))*lEx;
        float Fx0=Dx*sin(Cx*atan(Bx*kappax-Ex*(Bx*kappax-atan(Bx*kappax))))+SVx;
        // cout<<Bx<<" "<<" "<<Cx<<" "<<Dx<<" "<<Ex<<endl;
        return -Fx0;
    }

    float solveFyCombined(float alpha, float kappa){
        float Fy0=solveFy(alpha);
        Byk = rBy1*cos(atan(rBy2*(alpha - rBy3)))*lyk;
        Cyk = rCy1;
        Dyk = Fy0/(cos(Cyk*atan(Byk*SHyk-Eyk*(Byk*SHyk-atan(Byk*SHyk)))));
        DVyk = muy * Fz*(rVy1+rVy2*dfz+rVy3*gamma)*cos(atan(Byk*SHyk));
        Eyk = rEy1 +rEy2*dfz;
        SHyk = rHy1 + rHy2*dfz;
        kappa_s = kappa + SHyk;
        SVyk = DVyk*sin(rVy5*atan(rVy6*kappa))*lVyk;
        return Dyk*cos(Cyk*atan(Byk*kappa_s-Eyk*(Byk*kappa_s-atan(Byk*kappa_s))))+SVyk;
    }
    
    float solveFxCombined(float alpha, float kappa){
        float Fx0 = solveFx(kappa);
        Bxa = rBx1*cos(atan(rBx2*kappa))*lxa;
        Cxa = rCx1;
        Dxa = Fx0/(cos(Cxa*atan(Bxa*SHxa-Exa*(Bxa*SHxa-atan(Bxa*SHxa)))));
        Exa = rEx1 + rEx2*dfz;
        SHxa = rHx1;
        alpha_s = alpha + SHxa;
        return Dxa*cos(Cxa*atan(Bxa*alpha_s-Exa*(Bxa*alpha_s-atan(Bxa*alpha_s))));
    }
 
    float kappa(float vl, float omega){

            return (omega*r_eff()-vl)/(max(1.0f, vl)); ;
    }

};

class DynaBicycleModel
{
private:
    MagicTireModel Tire;
    Kesi kesi_new;
    Kesi kesi_old;
    States_dot states_dot_dyn;

    float vxf,vxr,vyf,vyr,vlf,vlr,omega;
    float Fdrv,Frrr,Frrf,Fdrag,Fbf,Fbr,Fxf,Fxr,Fyf,Fyr,alpha_f,alpha_r,kappa_f,kappa_r;

public:
    DynaBicycleModel(){}

    States_dot calculateStatesDot(float thr,float ste,float bra, Kesi kesi) {
        kesi_new.throttle = thr;
        kesi_new.steering_angle = ste;
        kesi_new.brakes = bra;

        // Fzf = lr*m*9.8/(4*(lf+lr));
        // Fzr = lf*m*9.8/(4*(lf+lr)); 
        Fdrv = kesi_new.throttle*Cm;
        Fbf = kesi_new.brakes*Cbf*tanh(kesi.v_x);
        Fbr = kesi_new.brakes*Cbr*tanh(kesi.v_x);
        Frrr = Crr*tanh(kesi.v_x);
        Frrf = Crr*tanh(kesi.v_x);
        Fdrag = Cd*kesi.v_x*kesi.v_x;

        vxf = kesi.v_x;
        vxr = kesi.v_x;
        vyf = kesi.v_y + lf * kesi.r;
        vyr = kesi.v_y-lr * kesi.r;
        alpha_f = atan2(vyf,vxf)-kesi_new.steering_angle;
        alpha_r = atan2(vyr,vxr);

        // longitudinal rear speed 
        vlf = vyf*sin(kesi_new.steering_angle)+vxf*cos(kesi_new.steering_angle);
        vlr = vxr;
        kappa_f = Tire.kappa(vlf,kesi.omega_f);
        kappa_r = Tire.kappa(vlr,kesi.omega_r);

        // Fyf = Tire.solveFy(alpha_f);
        // Fyr = Tire.solveFy(alpha_r);
        // Fxf = Tire.solveFx(kappa_f);
        // Fxr = Tire.solveFx(kappa_r);

        Fyf = 2 * Tire.solveFyCombined(alpha_f,kappa_f);
        Fyr = 2 * Tire.solveFyCombined(alpha_r,kappa_r);
        Fxf = 2 * Tire.solveFxCombined(alpha_f,kappa_f);
        Fxr = 2 * Tire.solveFxCombined(alpha_r,kappa_r);



        // cout<<Fdrv<<" "<<Frrr<<" "<<Frrf<<" "<< Fdrag<<" "<< Fbf<<" "<< Fbr<<" "<< alpha_f<<" "<< alpha_r<<" "<< Fyf<<" "<<Fyr<<" "<<endl;

        states_dot_dyn.X_dot = kesi.v_x*cos(kesi.theta)-kesi.v_y*sin(kesi.theta);
        states_dot_dyn.Y_dot = kesi.v_x*sin(kesi.theta)+kesi.v_y*cos(kesi.theta);
        states_dot_dyn.phi_dot = kesi.r;
        states_dot_dyn.vx_dot = 1/m*(m*kesi.v_y*kesi.r
                        //+2*Fdrv-2*Fbf*cos(kesi_new.steering_angle)-2*Fbr
                        +Fxf*cos(kesi_new.steering_angle)+Fxr
                        -Frrr-Frrf*cos(kesi_new.steering_angle)
                        -Fdrag-Fyf*sin(kesi_new.steering_angle)
                        );
        states_dot_dyn.vy_dot = 1/m*(-m*kesi.v_x*kesi.r
                        +Fxf*sin(kesi_new.steering_angle)
                        +Fyf*cos(kesi_new.steering_angle)+Fyr
                        // -(Frrf+2*Fbf)*sin(kesi_new.steering_angle));
                        -Frrf*sin(kesi_new.steering_angle));
                        
        states_dot_dyn.r_dot = 1/Iz*((Fyf*cos(kesi_new.steering_angle)
                        +Fxf*sin(kesi_new.steering_angle)
                        -Frrf*sin(kesi_new.steering_angle))*lf
                        -Fyr*lr);

        states_dot_dyn.omega_dot_f = -(Fxf-Fbf-Frrf)*Tire.r_eff()/Iwz;
        states_dot_dyn.omega_dot_r = (Fdrv-Fbr-Fxr-Frrr)*Tire.r_eff()/Iwz;
        
        return states_dot_dyn;
    }
    

    Kesi rungeKutta(float thr,float ste,float bra, float h) {
        // ofstream outputFile1("output1.txt",std::ios_base::app);

        // outputFile1<<"1 "<<Tire.kappa(vlf,kesi_old.omega_f)<<endl;

        States_dot k1 = calculateStatesDot(thr,ste,bra,kesi_old);
        Kesi kesi_k1 = kesi_old + k1 * (h / 2.0);
        // outputFile1<<"2 "<<Tire.kappa(vlf,kesi_k1.omega_f)<<endl;

        States_dot k2 = calculateStatesDot(thr,ste,bra,kesi_k1);
        Kesi kesi_k2 = kesi_old + k2 * (h / 2.0);
        // outputFile1<<"3 "<<Tire.kappa(vlf,kesi_k2.omega_f)<<endl;

        States_dot k3 = calculateStatesDot(thr,ste,bra,kesi_k2);
        Kesi kesi_k3 = kesi_old + k3 * h;
        // outputFile1<<"4 "<<Tire.kappa(vlf,kesi_k3.omega_f)<<endl;

        States_dot k4 = calculateStatesDot(thr,ste,bra,kesi_k3);
        States_dot k_result = (k1 + 2.0 * k2 + 2.0 * k3  + k4);
        // outputFile1<<"result_dot "<<k_result.omega_dot_f<<" "<<k_result.omega_dot_r<<endl;
        
        Kesi kesi_result = kesi_old + (k1 + 2.0 * k2 + 2.0 * k3  + k4) * (h / 6.0);
        // outputFile1<<"5 "<<Tire.kappa(vlf,kesi_result.omega_f)<<endl;

    return kesi_old + (k1 + 2.0 * k2 + 2.0 * k3  + k4) * (h / 6.0);
    }

    Kesi eulerIntegral(float thr,float ste,float bra, float h) {
        States_dot k = calculateStatesDot(thr,ste,bra,kesi_old);
    return kesi_old + h * k;
    }

    void updatestate(float dt){
        int steps;
        ifstream input("input.txt");
        ofstream outputFile("output.txt");
        input>>steps;
        float thr,ste,bra;

        for (int i = 0; i <= steps; ++i) {
            input >> thr >> ste >> bra;
            for(int j = 0;j < 50; ++j){
 
            //Eular integral
            // kesi_new = eulerIntegral(thr, ste, bra, dt);

            //RungeKutta integral
            // outputFile<<"before vlf "<<vlf<<" "<<kesi_old.omega_f*Tire.r_eff()<<" "<<Tire.kappa(vlf,kesi_old.omega_f)<<"\n";
            
            kesi_new = rungeKutta(thr, ste, bra, dt);
            
            // outputFile<<" Longitudinal FDRV="<<Fdrv<<" FXR="<<Fxr<<" FXF="<<Fxf<<" FRRF="<<Frrf<<" FBF="<<Fbf<<" FBR="<<Fbr<<"\n";
            // outputFile<<" vlf "<<vlf<<" "<<kesi_old.omega_f*Tire.r_eff()<<" "<<Tire.kappa(vlf,kesi_old.omega_f)<<"\n";
            // outputFile<<" Kappa "<<kappa_f<<" "<<kappa_r<<"\n";
            // outputFile<<" OMEGA_DOT "<<states_dot_dyn.omega_dot_f<<" "<<states_dot_dyn.omega_dot_r<<"\n";
            kesi_old = kesi_new;
            }
            outputFile<<i*0.05<<" "<<kesi_new.X<<" "<<kesi_new.Y<<" "<<kesi_new.theta<<" "<<kesi_new.v_x<<" "<<kesi_new.v_y<<" "<<kesi_new.r<<" "<<kesi_new.omega_f<<" "<<kesi_new.omega_r<<"\n"<<endl;
            cout<<kesi_new.X<<" "<<kesi_new.Y<<" "<<kesi_new.theta<<" "<<kesi_new.v_x<<" "<<kesi_new.v_y<<" "<<kesi_new.r<<endl;
        }
    }
};

int main(){
    DynaBicycleModel model1;
    model1.updatestate(0.001);
    // MagicTireModel Tire;
    // std::ofstream outputFile("output_tire.txt");
    // for(float i=-100;i<=100;i+=0.01){
    // // outputFile << i << " " << Tire.solveFy(i) << std::endl;
    // outputFile << i << " " << Tire.solveFx(i) << std::endl;
    // }
    // outputFile.close(); 
}
                              
