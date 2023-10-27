#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

struct Kesi {
    double X, Y, theta, v_x, v_y, r, throttle, steering_angle, brakes;

};

class DynaBicycleModel
{
private:
    Kesi kesi_new;
    Kesi kesi_old;
    double m=300,Iz=134,lf=0.721,lr=0.823,Cm=3600,Crr=200,Cd=1.53,Cbf=5411,Cbr=2650,Cx=20000;
    double Fdrv,Frrr,Frrf,Fdrag,Fbf,Fbr,Fry,Ffy,alpha_f,alpha_r;

public:
    DynaBicycleModel() :kesi_new({0, 0, 0, 0.00001, 0, 0, 0, 0, 0}),kesi_old({0, 0, 0, 0.00001, 0, 0, 0, 0, 0}) {}

    void updatestate(double dt){
        int steps;
        ifstream input("input.txt");
        input>>steps;
        double controll[steps][3]; 

        for (int i = 0; i < steps; ++i) {
            input >> controll[i][0] >> controll[i][1]>> controll[i][2];
            kesi_new.throttle = controll[i][0];
            kesi_new.steering_angle = controll[i][1];
            kesi_new.brakes = controll[i][2];

            Fdrv = kesi_new.throttle*Cm*kesi_new.throttle;
            Frrr = Crr*tanh(kesi_old.v_x);
            Frrf = Crr*tanh(kesi_old.v_x);
            Fdrag = Cd*kesi_old.v_x*kesi_old.v_x;
            Fbf = kesi_new.brakes*Cbf*tanh(kesi_old.v_x);
            Fbr = kesi_new.brakes*Cbr*tanh(kesi_old.v_x);
            // alpha_f = (kesi_old.v_y+lf*kesi_old.r)/kesi_old.v_x-kesi_new.steering_angle;
            // alpha_r = (kesi_old.v_y-lr*kesi_old.r)/kesi_old.v_x;
            alpha_f = atan(((kesi_old.v_y+lf*kesi_old.r)*cos(kesi_new.steering_angle)-kesi_old.v_x*sin(kesi_new.steering_angle))/((kesi_old.v_y+lf*kesi_old.r)*sin(kesi_new.steering_angle)+kesi_old.v_x*cos(kesi_new.steering_angle)));
            alpha_r = atan((kesi_old.v_y-lr*kesi_old.r)/kesi_old.v_x);
            Ffy = Cx*alpha_f;
            Fry = Cx*alpha_r;
            // cout<<Fdrv<<" "<<Frrr<<" "<<Frrf<<" "<< Fdrag<<" "<< Fbf<<" "<< Fbr<<" "<< alpha_f<<" "<< alpha_r<<" "<< Ffy<<" "<<Fry<<" "<<endl;
            // cout<<kesi_new.X<<" "<<kesi_new.Y<<" "<<kesi_new.theta<<" "<<kesi_new.v_x<<" "<<kesi_new.v_y<<" "<<kesi_new.r<<"\n"<<endl;

            kesi_new.X = kesi_old.X + (kesi_old.v_x*cos(kesi_old.theta)-kesi_old.v_y*sin(kesi_old.theta))*dt;
            kesi_new.Y = kesi_old.Y + (kesi_old.v_x*sin(kesi_old.theta)+kesi_old.v_y*cos(kesi_old.theta))*dt;
            kesi_new.theta = kesi_old.theta + kesi_old.r*dt;
            kesi_new.v_x =kesi_old.v_x + 1/m*(m*kesi_old.v_y*kesi_old.r+2*Fdrv-2*Frrr-2*Frrf*cos(kesi_new.steering_angle)
                            -Fdrag*cos(atan(kesi_old.v_y/kesi_old.v_x))-2*Fbf*cos(kesi_new.steering_angle)
                            -2*Fbr-2*Ffy*sin(kesi_new.steering_angle))*dt;
            if(kesi_old.v_x<0){kesi_old.v_x=0.00001;}
            kesi_new.v_y =kesi_old.v_y + 1/m*(+m*kesi_old.v_x*kesi_old.r+2*Fry+2*Ffy*cos(kesi_new.steering_angle)
                            -2*(Frrf+Fbf)*sin(kesi_new.steering_angle)
                            -Fdrag*sin(atan(kesi_old.v_y/kesi_old.v_x)))*dt;
            kesi_new.r = kesi_old.r + 1/Iz*((Ffy*cos(kesi_new.steering_angle)-Frrf*sin(kesi_new.steering_angle)
                            -Fbf*sin(kesi_new.steering_angle))*lf-Fry*lr)*dt;
            cout<<kesi_new.X<<" "<<kesi_new.Y<<" "<<kesi_new.theta<<" "<<kesi_new.v_x<<" "<<kesi_new.v_y<<" "<<kesi_new.r<<"\n"<<endl;
            kesi_old = kesi_new;
        }
    }
} ;

int main(){
    DynaBicycleModel model1;
    model1.updatestate(0.5);
}

// Ffy() = 0
//         Fbf() = 0                                         
//         Frrf() = 0
//         Fry() = 0
//         Frrf() = 0
//         Flateral() = 36
//         Ftransversal() = 0
// x(0.5) = {0 0 0 0.0597003 0.00599001 0 }
//         Ffy() = 0
//         Fbf() = 0
//         Frrf() = 11.9259
//         Fry() = 2000
//         Frrf() = 11.9259
//         Flateral() = 120.202
//         Ftransversal() = 1998.81
// x(1) = {0.0298501 0.002995 0 3.38039 0.339171 -6.14499 }
//         Ffy() = -19605.4
//         Fbf() = 0
//         Frrf() = 199.537
//         Fry() = 20223.7
//         Frrf() = 199.537
//         Flateral() = 1685.8
//         Ftransversal() = 694.676
// x(1.5) = {1.72005 0.17258 -3.0725 6.40407 0.642551 -73.5692 }
//         Ffy() = -28983.7
//         Fbf() = 5410.97
//         Frrf() = 199.999
//         Fry() = 29330.3
//         Frrf() = 199.999
//         Flateral() = -8523.39
//         Ftransversal() = 340.354
// x(2) = {-1.45217 -0.368999 -39.8571 20.621 0.642551 -163.64 }
//         Ffy() = -27936.8
//         Fbf() = 5411
//         Frrf() = 200
//         Fry() = 28391.4
//         Frrf() = 200
//         Flateral() = -9111.28
//         Ftransversal() = 434.391
// x(2.5) = {-6.89718 -9.13039 -121.677 35.8238 0.642551 -250.827 }


// cout<<Fdrv<<" "<<Frrr<<" "<<Frrf<<" "<< Fdrag<<" "<< Fbf<<" "<< Fbr<<" "<< alpha_f<<" "<< alpha_r<<" "<< Ffy<<" "<<Fry<<" "<<endl;
// 36 0.002 0.002 1.53e-10 0 0 -0.1 0 -2000 0 
// 5e-06 0 0 0.785553 -6.63336 -5.35372

// 144 131.176 131.176 0.944153 0 0 -13.458 -2.83527 -269159 -56705.4 
// 0.392781 -3.31668 -2.67686 107.72 -1090.51 -551.754

// 144 200 200 17753.5 0 0 -13.9166 -5.90808 -278333 -118162 
// -292.13 459.971 -278.554 301044 -32095.6 -934.003

// 0 200 200 1.3866e+11 5411 2650 -0.108851 -0.104061 -2177.03 -2081.22 
// -89386.4 -121919 -745.556 -2.14508e+08 -1.1612e+08 -933.469

// 0 -200 -200 7.04009e+16 -5411 -2650 0.541337 0.54133 10826.7 10826.6 
// 1.0683e+08 -5.87965e+07 -1212.29 -1.03132e+14 1.17335e+14 -937.589