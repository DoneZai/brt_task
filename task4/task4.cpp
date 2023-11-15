#include <iostream>
#include <cmath>
#include <fstream>

#include "../include/tire_constants.h" 
#include "car_constants.h" 

using namespace std;
using namespace TireConstants;
using namespace CarConstants;

struct Kesi {
    float X, Y, theta, v_x, v_y, r, throttle, steering_angle, brakes;
};

struct States_dot {
    float X_dot, Y_dot, phi_dot,vx_dot,vy_dot,r_dot; 

    States_dot operator*=(float scalar) const {
        States_dot result;
        result.X_dot *= scalar;
        result.Y_dot *= scalar;
        result.phi_dot *= scalar;
        result.vx_dot *= scalar;
        result.vy_dot *= scalar;
        result.r_dot *= scalar;
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
    float Bx,Cx,Dx,Ex,Kx,SHx,SVx,mux,kappax,kappa;
    float gammay,alphay,gammax,alphax;
    float Fz0,Fz,dfz;
    

    float sgn(float x){
        return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
    }
    
public:
    MagicTireModel(){}
    
    float solveFy(float alpha, float gamma){
        Fz0 = FNOMIN;
        Fz = m*9.8/4;
        dfz = (Fz-Fz0*lFz0)/Fz0*lFz0;
        gammay = gamma*lgammay;
        SHy = (pHy1+pHy2*dfz)*lHy + pHy3*gammay;
        SVy = Fz*((pVy1+pVy2*dfz)*lVy+(pVy3+pVy4*dfz)*gammay)*lmuy;
        Ky = pKy1*Fz0*sin(2*atan(Fz/(pKy2*Fz0*lFz0)))*(1-pKy3*abs(gammay))*lFz0*lKy;
        muy = (pDy1+pDy2*dfz)*(1-pDy3*gammay*gammay)*lmuy;
        Cy = pCy1*lCy;
        Dy = muy*Fz;
        By = Ky/(Cy*Dy);
        Ey = (pEy1+pEy2*dfz);
        alphay = alpha+SHy;
        float Fy0=Dy*sin(Cy*atan(By*alphay-Ey*(By*alphay-atan(By*alphay))))+SVy;
        return Fy0;
    }
    
    float solveFx(float kappa, float gamma){
        Fz0 = FNOMIN;
        Fz = m*9.8/4;
        dfz = (Fz-Fz0*lFz0)/Fz0*lFz0;
        gammax = gamma*lgammax;
        SHx = (pHx1+pHx2*dfz)*lHx;
        SVx = Fz*(pVx1+pVx2*dfz)*lVx*lmux;
        Kx = Fz*(pKx1+pKx2*dfz)*exp(pKx3*dfz)*lKx;
        mux = (pDx1+pDx2*dfz)*(1-pDx3*gammax*gammax)*lmux;
        Cx = pCx1*lCx;
        Dx = mux*Fz;
        Bx = Kx/(Cx*Dx);
        Ex = (pEx1+pEx2*dfz+pEx3*dfz*dfz)*(1-pEx4*sgn(kappax))*lEx;
        kappax = kappa + SHx;
        float Fx0=Dx*sin(Cx*atan(Bx*kappax-Ex*(Bx*kappax-atan(Bx*kappax))))+SVx;
        return Fx0;
    }
};

class DynaBicycleModel
{
private:
    MagicTireModel Tire;
    Kesi kesi_new;
    Kesi kesi_old;
    States_dot states_dot_dyn;
    States_dot states_dot_kin;
    States_dot states_dot_all;

    float vxf,vxr,vyf,vyr,vlf,vlr,omega;
    float Fdrv,Frrr,Frrf,Fdrag,Fbf,Fbr,Ffx,Frx,Ffy,Fry,alpha_f,alpha_r,kappa_f,kappa_r;

public:
    DynaBicycleModel() :kesi_new({0, 0, 0, 0, 0, 0, 0, 0, 0}),kesi_old({0, 0, 0, 0, 0, 0, 0, 0, 0}) {}

    float kappa(float vl, float omega){
        if(vl>R_tire*omega){
            return R_tire*omega/vl-1;
        }
        else{
            return 1-vl/(R_tire*omega);
        }
    }

    States_dot calculateStatesDot(float thr,float ste,float bra, const Kesi& kesi_old) {
        States_dot states_dot_dyn;
            kesi_new.throttle = thr;
            kesi_new.steering_angle = ste;
            kesi_new.brakes = bra;

            // outputFile<<kesi_new.throttle<<" "<<kesi_new.steering_angle<<" "<<kesi_new.brakes<<endl;
            
            // Fzf = lr*m*9.8/(4*(lf+lr));
            // Fzr = lf*m*9.8/(4*(lf+lr)); 
            Fdrv = kesi_new.throttle*Cm*kesi_new.throttle;
            Fbf = kesi_new.brakes*Cbf*tanh(kesi_old.v_x);
            Fbr = kesi_new.brakes*Cbr*tanh(kesi_old.v_x);
            Fbf = kesi_new.brakes*Cbf*tanh(kesi_old.v_x);
            Fbr = kesi_new.brakes*Cbr*tanh(kesi_old.v_x);
            Frrr = Crr*tanh(kesi_old.v_x);
            Frrf = Crr*tanh(kesi_old.v_x);
            Fdrag = Cd*kesi_old.v_x*kesi_old.v_x;
            alpha_f = atan2((kesi_old.v_y+lf*kesi_old.r),kesi_old.v_x)-kesi_new.steering_angle;
            alpha_r = atan2((kesi_old.v_y-lr*kesi_old.r),kesi_old.v_x);
            vxf = kesi_old.v_x;
            vxr = kesi_old.v_x;
            vyf = kesi_old.v_y + lf * kesi_old.r;
            vyr = kesi_old.v_y-lr*kesi_old.r;
            vlf = vyf*sin(kesi_new.steering_angle)+vxr*cos(kesi_new.steering_angle);
            vlr = vyr*sin(kesi_new.steering_angle)+vyr*cos(kesi_new.steering_angle);
            kappa_f = kappa(vlf,omega);
            kappa_r = kappa(vlr,omega);

            Ffy = Tire.solveFy(alpha_f,0);
            Fry = Tire.solveFy(alpha_r,0);
            Ffx = Tire.solveFx(kappa_f,0);
            Frx = Tire.solveFx(kappa_r,0);

            // cout<<Fdrv<<" "<<Frrr<<" "<<Frrf<<" "<< Fdrag<<" "<< Fbf<<" "<< Fbr<<" "<< alpha_f<<" "<< alpha_r<<" "<< Ffy<<" "<<Fry<<" "<<endl;
            // outputFile<<Fdrv<<" "<<Frrr<<" "<<Frrf<<" "<< Fdrag<<" "<< Fbf<<" "<< Fbr<<" "<< alpha_f<<" "<< alpha_r<<" "<< Ffy<<" "<<Fry<<" "<<endl;

            states_dot_dyn.X_dot = kesi_old.v_x*cos(kesi_old.theta)-kesi_old.v_y*sin(kesi_old.theta);
            states_dot_dyn.Y_dot = kesi_old.v_x*sin(kesi_old.theta)+kesi_old.v_y*cos(kesi_old.theta);
            states_dot_dyn.phi_dot = kesi_old.r;
            states_dot_dyn.vx_dot = 1/m*(m*kesi_old.v_y*kesi_old.r
                            //+2*Fdrv-2*Fbf*cos(kesi_new.steering_angle)-2*Fbr
                            -Frrr-Frrf*cos(kesi_new.steering_angle)
                            -Fdrag-2*Ffy*sin(kesi_new.steering_angle)
                            );
            states_dot_dyn.vy_dot = 1/m*(-m*kesi_old.v_x*kesi_old.r
                            +2*Fry+2*Ffy*cos(kesi_new.steering_angle)
                            -(Frrf+2*Fbf)*sin(kesi_new.steering_angle));
            states_dot_dyn.r_dot = 1/Iz*(2*(Ffy*cos(kesi_new.steering_angle)
                            -Frrf*sin(kesi_new.steering_angle)
                            -Fbf*sin(kesi_new.steering_angle))*lf-2*Fry*lr);
            

        return states_dot_dyn;
    }
    

    Kesi rungeKutta(float thr,float ste,float bra,States_dot states_dot, float h) {

        States_dot k1 = calculateStatesDot(thr,ste,bra,kesi_old);
        Kesi kesi_k1 = kesi_old + k1 * (h / 2.0);

        States_dot k2 = calculateStatesDot(thr,ste,bra,kesi_k1);
        Kesi kesi_k2 = kesi_old + k2 * (h / 2.0);

        States_dot k3 = calculateStatesDot(thr,ste,bra,kesi_k2);
        Kesi kesi_k3 = kesi_old + k3 * h;

        States_dot k4 = calculateStatesDot(thr,ste,bra,kesi_k3);

        kesi_new = kesi_old + (k1 + k2 *2.0 + k3 *2.0 + k4) * (h / 6.0);

    return kesi_new;
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
            states_dot_dyn = calculateStatesDot(thr,ste,bra,kesi_old);

            // states_dot_kin = states_dot_dyn;
            // states_dot_kin.vx_dot = 1/m*(2*Fdrv-2*Frrr-2*Frrf*cos(kesi_new.steering_angle)
            //                 -Fdrag*cos(atan(kesi_old.v_y/kesi_old.v_x))-2*Fbf*cos(kesi_new.steering_angle)
            //                 -2*Fbr);
            // states_dot_kin.vy_dot = ((kesi_new.steering_angle-kesi_old.steering_angle)/dt*kesi_old.v_x+kesi_new.steering_angle*states_dot_kin.vx_dot)*(lr/(lr+lf));
            // states_dot_kin.r_dot = ((kesi_new.steering_angle-kesi_old.steering_angle)/dt*kesi_old.v_x+kesi_new.steering_angle*states_dot_kin.vx_dot)*(1/(lr+lf));

            // float lamda = min(max((kesi_old.v_x-3)/(5-3),0.0),1.0);
            // states_dot_all = states_dot_dyn*lamda+states_dot_kin*(1-lamda);
            states_dot_all = states_dot_dyn;

            //Eular integral
            // kesi_new.X = kesi_old.X + states_dot_all.X_dot*dt;
            // kesi_new.Y = kesi_old.Y + states_dot_all.Y_dot*dt;
            // kesi_new.theta = kesi_old.theta +states_dot_all.phi_dot*dt;
            // kesi_new.v_x = kesi_old.v_x + states_dot_all.vx_dot*dt;                
            // kesi_new.v_y = kesi_old.v_y + states_dot_all.vy_dot*dt;                 
            // kesi_new.r = kesi_old.r + states_dot_all.r_dot*dt;

            rungeKutta(thr, ste, bra, states_dot_all, dt);
            if(kesi_new.v_x<0){kesi_new.v_x=0.0001;kesi_new.v_y=0;kesi_new.r=0;} //kesi_new.v_y=0;kesi_new.r=10;
            // cout<<kesi_new.X<<" "<<kesi_new.Y<<" "<<kesi_new.theta<<" "<<kesi_new.v_x<<" "<<kesi_new.v_y<<" "<<kesi_new.r<<"\n"<<endl;
            outputFile<<i*0.05+j*dt<<" "<<kesi_new.X<<" "<<kesi_new.Y<<" "<<kesi_new.theta<<" "<<kesi_new.v_x<<" "<<kesi_new.v_y<<" "<<kesi_new.r<<"\n"<<endl;
            kesi_old = kesi_new;
            
            }
        }
    }
};

int main(){
    DynaBicycleModel model1;
    model1.updatestate(0.001);
    // MagicTireModel Tire;
    // std::ofstream outputFile("output_tire.txt");
    // for(float i=-0.5;i<=0.5;i+=0.01){
    // outputFile << i << " " << Tire.solveFy(i) << std::endl;}
    // outputFile.close(); 
}
                              
