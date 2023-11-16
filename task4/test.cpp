#include <iostream>
#include <cmath>
#include <fstream>


using namespace std;

class Kesi {
public:
    float X=0.0f, Y=0.0f, theta=0.0f, v_x=0.0f, v_y=0.0f, r=0.0f, omega_f=0.0f,omega_r=0.0f, throttle=0.0f, steering_angle=0.0f, brakes=0.0f;
};

class States_dot {
public:
    float X_dot, Y_dot, phi_dot,vx_dot,vy_dot,r_dot,omega_dot_f,omega_dot_r; 

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

Kesi eulerIntegral(Kesi kesi_old, float h) {

    States_dot k;
    k.omega_dot_f = 1;
    k.omega_dot_r = 1;
    k.X_dot = 1;
    k.Y_dot = 1;
    k.phi_dot = 1;
    k.r_dot = 1;
    k.vx_dot = 1;
    k.vy_dot = 1;

return kesi_old + k * h ;
}
States_dot eulerIntegral1(States_dot kk, float h) {
    return h * kk ;
}



int main(){
    Kesi kesi1;
    Kesi kesi2;
    States_dot k;
    k.omega_dot_f = 1;
    k.omega_dot_r = 1;
    k.X_dot = 1;
    k.Y_dot = 1;
    k.phi_dot = 1;
    k.r_dot = 1;
    k.vx_dot = 1;
    k.vy_dot = 1;


    k = eulerIntegral1(k,2);
    cout<<k.X_dot<<" "<<k.Y_dot<<" "<<k.phi_dot<<" "<<k.vx_dot<<" "<<k.vy_dot<<" "<<k.r_dot<<endl;

    cout<<kesi2.X<<" "<<kesi2.Y<<" "<<kesi2.theta<<" "<<kesi2.v_x<<" "<<kesi2.v_y<<" "<<kesi2.r<<endl;
    kesi2 = kesi1 + k;
    cout<<kesi2.X<<" "<<kesi2.Y<<" "<<kesi2.theta<<" "<<kesi2.v_x<<" "<<kesi2.v_y<<" "<<kesi2.r<<endl;
    kesi2 = eulerIntegral(kesi1,0.002);
    cout<<kesi2.X<<" "<<kesi2.Y<<" "<<kesi2.theta<<" "<<kesi2.v_x<<" "<<kesi2.v_y<<" "<<kesi2.r<<endl;
}
                              
