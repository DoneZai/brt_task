#include <iostream>
#include <cmath>
#include <fstream>

#include "../include/tire_constants.h" 
#include "../include/car_constants.h" 

using namespace std;
using namespace TireConstants;
using namespace CarConstants;

struct Kesi {
    float X, Y, theta, v_x, v_y, r, omega_f,omega_r, throttle, steering_angle, brakes;
};

class States_dot {
public:
    float X_dot, Y_dot, phi_dot,vx_dot,vy_dot,r_dot,omega_dot_f,omega_dot_r; 

    States_dot operator*=(float scalar) const {
        States_dot result;
        result.X_dot *= scalar;
        result.Y_dot *= scalar;
        result.phi_dot *= scalar;
        result.vx_dot *= scalar;
        result.vy_dot *= scalar;
        result.r_dot *= scalar;
        result.omega_dot_f *= scalar;
        result.omega_dot_r *= scalar;
        return result;
    }

    States_dot operator+(const States_dot& other) const {
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

    States_dot operator*(float scalar)const{
        return *this*=scalar;
        }

};
    States_dot operator*(float scalar,const States_dot& other){
        return other*=scalar;
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
    float gamma,gammax,gammay,alphay,kappax;
    float Fz0,Fz,dfz;
    // float r_eff,r_stat ;
    

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
        // cout<<gammax<<" "<<SHx<<" "<<SVx<<" "<<Kx<<" "<<mux<<" "<<Cx<<" "<<Dx<<" "<<Bx<<" "<<Ex<<" "<<endl;
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
        return Fx0;
    }

    float kappa(float vl, float omega){
        if(vl>r_eff()*omega&&vl!=0){
            cout<<omega<<" braking "<<r_eff()*omega/vl-1<<endl;
            return r_eff()*omega/vl-1;
        }
        if(vl<r_eff()*omega&&omega!=0){
            cout<<vl<<" speeding "<<1-vl/(r_eff()*omega)<<endl;
            return 1-vl/(r_eff()*omega);
        }
        else{
            cout<<"nothing"<<endl;
            return 0;
        }
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
    DynaBicycleModel() :kesi_new({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}),kesi_old({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}),
                        states_dot_dyn({0, 0, 0, 0, 0, 0, 0, 0}){}

    States_dot calculateStatesDot(float thr,float ste,float bra, const Kesi& kesi) {
        kesi_new.throttle = thr;
        kesi_new.steering_angle = ste;
        kesi_new.brakes = bra;

        // Fzf = lr*m*9.8/(4*(lf+lr));
        // Fzr = lf*m*9.8/(4*(lf+lr)); 
        Fdrv = kesi_new.throttle*Cm*kesi_new.throttle;
        Fbf = kesi_new.brakes*Cbf*tanh(kesi.v_x);
        Fbr = kesi_new.brakes*Cbr*tanh(kesi.v_x);
        Frrr = Crr*tanh(kesi.v_x);
        Frrf = Crr*tanh(kesi.v_x);
        Fdrag = Cd*kesi.v_x*kesi.v_x;
        vxf = kesi.v_x;
        vxr = kesi.v_x;
        vyf = kesi.v_y + lf * kesi.r;
        vyr = kesi.v_y-lr*kesi.r;
        alpha_f = atan2(vyf,vxf)-kesi_new.steering_angle;
        alpha_r = atan2(vyr,vxr);
        // longitudinal rear speed 
        vlf = vyf*sin(kesi_new.steering_angle)+vxf*cos(kesi_new.steering_angle);
        vlr = vxr;
        kappa_f = Tire.kappa(vlf,kesi.omega_f);
        kappa_r = Tire.kappa(vlr,kesi.omega_r);

        Fyf = Tire.solveFy(alpha_f);
        Fyr = Tire.solveFy(alpha_r);
        Fxf = Tire.solveFx(kappa_f);
        Fxr = Tire.solveFx(kappa_r);

        // cout<<Fdrv<<" "<<Frrr<<" "<<Frrf<<" "<< Fdrag<<" "<< Fbf<<" "<< Fbr<<" "<< alpha_f<<" "<< alpha_r<<" "<< Fyf<<" "<<Fyr<<" "<<endl;

        states_dot_dyn.X_dot = kesi.v_x*cos(kesi.theta)-kesi.v_y*sin(kesi.theta);
        states_dot_dyn.Y_dot = kesi.v_x*sin(kesi.theta)+kesi.v_y*cos(kesi.theta);
        states_dot_dyn.phi_dot = kesi.r;
        states_dot_dyn.vx_dot = 1/m*(m*kesi.v_y*kesi.r
                        //+2*Fdrv-2*Fbf*cos(kesi_new.steering_angle)-2*Fbr
                        +2*Fxf*cos(kesi_new.steering_angle)+2*Fxr
                        -Frrr-Frrf*cos(kesi_new.steering_angle)
                        -Fdrag-2*Fyf*sin(kesi_new.steering_angle)
                        );
        states_dot_dyn.vy_dot = 1/m*(-m*kesi.v_x*kesi.r
                        +2*Fxf*sin(kesi_new.steering_angle)
                        +2*Fyf*cos(kesi_new.steering_angle)+2*Fyr
                        // -(Frrf+2*Fbf)*sin(kesi_new.steering_angle));
                        -Frrf*sin(kesi_new.steering_angle));
                        
        states_dot_dyn.r_dot = 1/Iz*(2*(Fyf*cos(kesi_new.steering_angle)
                        +Fxf*sin(kesi_new.steering_angle)
                        -Frrf*sin(kesi_new.steering_angle))*lf
                        -2*Fyr*lr);

        // states_dot_dyn.vx_dot = 1/m*(m*kesi.v_y*kesi.r+2*Fdrv
        //                 -Frrr-Frrf*cos(kesi_new.steering_angle)
        //                 -Fdrag-2*Fbf*cos(kesi_new.steering_angle)-2*Fbr
        //                 -2*Fyf*sin(kesi_new.steering_angle)
        //                 );
        // states_dot_dyn.vy_dot = 1/m*(-m*kesi.v_x*kesi.r
        //                 +2*Fyr+2*Fyf*cos(kesi_new.steering_angle)
        //                 -(Frrf+2*Fbf)*sin(kesi_new.steering_angle));
        // states_dot_dyn.r_dot = 1/Iz*((2*Fyf*cos(kesi_new.steering_angle)
        //                 -Frrf*sin(kesi_new.steering_angle)
        //                 -2*Fbf*sin(kesi_new.steering_angle))*lf-2*Fyr*lr);

        states_dot_dyn.omega_dot_f = (Fdrv-Fbf-Fxf)*Tire.r_eff()/Iwz;
        states_dot_dyn.omega_dot_r = (Fdrv-Fbr-Fxr)*Tire.r_eff()/Iwz;
        

        return states_dot_dyn;
    }
    

    Kesi rungeKutta(float thr,float ste,float bra, const Kesi& kesi, float h) {

        States_dot k1 = calculateStatesDot(thr,ste,bra,kesi);
        Kesi kesi_k1 = kesi + k1 * (h / 2.0);

        States_dot k2 = calculateStatesDot(thr,ste,bra,kesi_k1);
        Kesi kesi_k2 = kesi + k2 * (h / 2.0);

        States_dot k3 = calculateStatesDot(thr,ste,bra,kesi_k2);
        Kesi kesi_k3 = kesi + k3 * h;

        States_dot k4 = calculateStatesDot(thr,ste,bra,kesi_k3);

    return kesi + (k1 + k2 *2.0 + k3 *2.0 + k4) * (h / 6.0);
    }

    void updatestate(float dt){
        int steps;
        ifstream input("input.txt");
        ofstream outputFile("output.txt");
        input>>steps;
        float thr,ste,bra;

        for (int i = 0; i < steps; ++i) {
            input >> thr >> ste >> bra;
            for(int j = 0;j < 50; ++j){
            outputFile<<states_dot_dyn.omega_dot_f<<" "<<states_dot_dyn.omega_dot_r<<endl;
 
            //Eular integral
            // states_dot_dyn = calculateStatesDot(thr,ste,bra,kesi_old);
            // kesi_new.X = kesi_old.X + states_dot_dyn.X_dot*dt;
            // kesi_new.Y = kesi_old.Y + states_dot_dyn.Y_dot*dt;
            // kesi_new.theta = kesi_old.theta +states_dot_dyn.phi_dot*dt;
            // kesi_new.v_x = kesi_old.v_x + states_dot_dyn.vx_dot*dt;                
            // kesi_new.v_y = kesi_old.v_y + states_dot_dyn.vy_dot*dt;                 
            // kesi_new.r = kesi_old.r + states_dot_dyn.r_dot*dt;
            // kesi_new.omega_f = kesi_old.omega_f + states_dot_dyn.omega_dot_f*dt;
            // kesi_new.omega_r = kesi_old.omega_r + states_dot_dyn.omega_dot_r*dt;

            kesi_new = rungeKutta(thr, ste, bra,kesi_old, dt);
            outputFile<<i*0.05+j*dt<<" "<<kesi_new.X<<" "<<kesi_new.Y<<" "<<kesi_new.theta<<" "<<kesi_new.v_x<<" "<<kesi_new.v_y<<" "<<kesi_new.r<<" "<<kesi_old.omega_f<<" "<<kesi_old.omega_r<<"\n"<<endl;
            // outputFile<<kesi_new.r<<" "<<kesi_old.omega_f<<" "<<kesi_old.omega_r<<"\n"<<endl;
            // outputFile<<Fxf<<" "<<Fyf<<" "<<Fxr<<" "<<Fyr<<"\n"<<endl;
            // outputFile<<Fdrv<<" "<<Frrr<<" "<<Frrf<<" "<< Fdrag<<" "<< Fbf<<" "<< Fbr<<" "<< alpha_f<<" "<< alpha_r<<" "<< Fyf<<" "<<Fyr<<" "<<endl;
            kesi_old = kesi_new;

            }
        }
        // outputFile.close();
    }
};

int main(){
    DynaBicycleModel model1;
    model1.updatestate(0.001);
    // MagicTireModel Tire;
    // std::ofstream outputFile("output_tire.txt");
    // for(float i=-0.5;i<=0.5;i+=0.01){
    // outputFile << i << " " << Tire.solveFy(i,0) << std::endl;
    // outputFile << i << " " << Tire.solveFx(i) << std::endl;
    // }
    // outputFile.close(); 
}
                              
