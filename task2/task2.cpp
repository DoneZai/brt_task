#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

struct Kesi {
    double X, Y, theta, v_x, v_y, r, throttle, steering_angle, brakes;
};

struct states_dot {
    double X_dot, Y_dot, phi_dot,vx_dot,vy_dot,r_dot; 

    states_dot operator*(double scalar) const {
        states_dot result;
        result.X_dot = X_dot * scalar;
        result.Y_dot = Y_dot * scalar;
        result.phi_dot = phi_dot * scalar;
        result.vx_dot = vx_dot * scalar;
        result.vy_dot = vy_dot * scalar;
        result.r_dot = r_dot * scalar;
        return result;
    }

    states_dot operator+(const states_dot& other) const {
        states_dot result;
        result.X_dot = X_dot + other.X_dot;
        result.Y_dot = Y_dot + other.Y_dot;
        result.phi_dot = phi_dot + other.phi_dot;
        result.vx_dot = vx_dot + other.vx_dot;
        result.vy_dot = vy_dot + other.vy_dot;
        result.r_dot = r_dot + other.r_dot;
        return result;
    }
};

class DynaBicycleModel
{
private:
    Kesi kesi_new;
    Kesi kesi_old;
    states_dot states_dot_dyn;
    states_dot states_dot_kin;
    states_dot states_dot_all;

    double m=300,Iz=134,lf=0.721,lr=0.823,Cm=3600,Crr=200,Cd=1.53,Cbf=5411,Cbr=2650,Cx=20000;
    double Fdrv,Frrr,Frrf,Fdrag,Fbf,Fbr,Fry,Ffy,alpha_f,alpha_r;

public:
    DynaBicycleModel() :kesi_new({0, 0, 0, 0, 0, 0, 0, 0, 0}),kesi_old({0.001, 0.001, 0.001, 20, 0.001, 0.001, 0, 0, 0}) {}

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

            alpha_f = atan2((kesi_old.v_y+lf*kesi_old.r),kesi_old.v_x)-kesi_new.steering_angle;
            alpha_r = atan2((kesi_old.v_y-lr*kesi_old.r),kesi_old.v_x);
            // if(kesi_new.v_x<0.1){alpha_f = 0}
            // if (alpha_f>0.087){alpha_f = 0.087;}
            // if (alpha_r>0.087){alpha_r = 0.087;}
            // if (alpha_f<-0.087){alpha_f = -0.087;}
            // if (alpha_r<-0.087){alpha_r = -0.087;} //Ограничить угол скольжения шины линейной областью magic formula

            Ffy = -Cx*alpha_f;
            Fry = -Cx*alpha_r;
            cout<<Fdrv<<" "<<Frrr<<" "<<Frrf<<" "<< Fdrag<<" "<< Fbf<<" "<< Fbr<<" "<< alpha_f<<" "<< alpha_r<<" "<< Ffy<<" "<<Fry<<" "<<endl;

            states_dot_dyn.X_dot = kesi_old.v_x*cos(kesi_old.theta)-kesi_old.v_y*sin(kesi_old.theta);
            states_dot_dyn.Y_dot = kesi_old.v_x*sin(kesi_old.theta)+kesi_old.v_y*cos(kesi_old.theta);
            states_dot_dyn.phi_dot = kesi_old.r;
            states_dot_dyn.vx_dot = 1/m*(m*kesi_old.v_y*kesi_old.r+2*Fdrv-2*Frrr-2*Frrf*cos(kesi_new.steering_angle)
                            -Fdrag*cos(atan(kesi_old.v_y/kesi_old.v_x))-2*Fbf*cos(kesi_new.steering_angle)
                            -2*Fbr-2*Ffy*sin(kesi_new.steering_angle));
            states_dot_dyn.vy_dot = 1/m*(-m*kesi_old.v_x*kesi_old.r+2*Fry+2*Ffy*cos(kesi_new.steering_angle)
                            -2*(Frrf+Fbf)*sin(kesi_new.steering_angle)
                            -Fdrag*sin(atan(kesi_old.v_y/kesi_old.v_x)));
            states_dot_dyn.r_dot = 1/Iz*((Ffy*cos(kesi_new.steering_angle)-Frrf*sin(kesi_new.steering_angle)
                            -Fbf*sin(kesi_new.steering_angle))*lf-Fry*lr);

            states_dot_dyn.vx_dot = 1/m*(m*kesi_old.v_y*kesi_old.r+2*Fdrv-2*Frrr-2*Frrf-Fdrag-2*Fbf-2*Fbr-2*Ffy*sin(kesi_new.steering_angle));
            states_dot_dyn.vy_dot = 1/m*(-m*kesi_old.v_x*kesi_old.r+2*Fry+2*Ffy*cos(kesi_new.steering_angle));
            states_dot_dyn.r_dot = 1/Iz*((2*Ffy*cos(kesi_new.steering_angle))*lf-2*Fry*lr);


            // states_dot_kin = states_dot_dyn;
            // states_dot_kin.vx_dot = 1/m*(2*Fdrv-2*Frrr-2*Frrf*cos(kesi_new.steering_angle)
            //                 -Fdrag*cos(atan(kesi_old.v_y/kesi_old.v_x))-2*Fbf*cos(kesi_new.steering_angle)
            //                 -2*Fbr);
            // states_dot_kin.vy_dot = ((kesi_new.steering_angle-kesi_old.steering_angle)/dt*kesi_old.v_x+kesi_new.steering_angle*states_dot_kin.vx_dot)*(lr/(lr+lf));
            // states_dot_kin.r_dot = ((kesi_new.steering_angle-kesi_old.steering_angle)/dt*kesi_old.v_x+kesi_new.steering_angle*states_dot_kin.vx_dot)*(1/(lr+lf));

            // double lamda = min(max((kesi_old.v_x-3)/(5-3),0.0),1.0);
            // states_dot_all = states_dot_dyn*lamda+states_dot_kin*(1-lamda);
            states_dot_all = states_dot_dyn;

            // kesi_new.X = kesi_old.X + (kesi_old.v_x*cos(kesi_old.theta)-kesi_old.v_y*sin(kesi_old.theta))*dt;
            // kesi_new.Y = kesi_old.X + (kesi_old.v_x*sin(kesi_old.theta)+kesi_old.v_y*cos(kesi_old.theta))*dt;
            // kesi_new.theta = kesi_old.theta + kesi_old.r*dt;
            // kesi_new.v_x =kesi_old.v_x + 1/m*(m*kesi_old.v_y*kesi_old.r+2*Fdrv-2*Frrr-2*Frrf*cos(kesi_new.steering_angle)
            //                 -Fdrag*cos(atan(kesi_old.v_y/kesi_old.v_x))-2*Fbf*cos(kesi_new.steering_angle)
            //                 -2*Fbr-2*Ffy*sin(kesi_new.steering_angle))*dt;
            // kesi_new.v_y =kesi_old.v_y + 1/m*(-m*kesi_old.v_x*kesi_old.r+2*Fry+2*Ffy*cos(kesi_new.steering_angle)
            //                 -2*(Frrf+Fbf)*sin(kesi_new.steering_angle)
            //                 -Fdrag*sin(atan(kesi_old.v_y/kesi_old.v_x)))*dt;
            // kesi_new.r = kesi_old.r + 1/Iz*((Ffy*cos(kesi_new.steering_angle)-Frrf*sin(kesi_new.steering_angle)
            //                 -Fbf*sin(kesi_new.steering_angle))*lf-Fry*lr)*dt;
            kesi_new.X = kesi_old.X + states_dot_all.X_dot*dt;
            kesi_new.Y = kesi_old.Y + states_dot_all.Y_dot*dt;
            kesi_new.theta = kesi_old.theta +states_dot_all.phi_dot*dt;
            kesi_new.v_x = kesi_old.v_x + states_dot_all.vx_dot*dt;                
            kesi_new.v_y = kesi_old.v_y + states_dot_all.vy_dot*dt;                 
            kesi_new.r = kesi_old.r + states_dot_all.r_dot*dt;
              
            if(kesi_new.v_x<0){kesi_new.v_x=0.00001;kesi_new.v_y=0;kesi_new.r=0;} //kesi_new.v_y=0;kesi_new.r=10;
            cout<<kesi_new.X<<" "<<kesi_new.Y<<" "<<kesi_new.theta<<" "<<kesi_new.v_x<<" "<<kesi_new.v_y<<" "<<kesi_new.r<<"\n"<<endl;
            kesi_old = kesi_new;
        }
    }
} ;

int main(){
    DynaBicycleModel model1;
    model1.updatestate(0.5);
}
                              
