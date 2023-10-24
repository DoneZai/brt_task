#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

struct Kesi {
    double X, Y, theta, velocity, L, v_dot, steering_angle;
};

class KineBicycleModel
{
private:
    Kesi kesi;
public:
    KineBicycleModel() :kesi({0, 0, 0, 0, 1.5, 0, 0}) {}

    void updatestate(double dt){
        int steps;

        ifstream input("input.txt");
        input>>steps;
        double controll[steps][2]; 

        for (int i = 0; i < steps; ++i) {
            input >> controll[i][0] >> controll[i][1];
            kesi.v_dot = controll[i][0];
            kesi.steering_angle = controll[i][1];
            kesi.X += kesi.velocity*cos(kesi.theta)*dt;
            kesi.Y += kesi.velocity*sin(kesi.theta)*dt;
            kesi.theta += kesi.velocity*tan(kesi.steering_angle)/kesi.L*dt;
            kesi.velocity += kesi.v_dot*dt;
            cout<<kesi.X<<" "<<kesi.Y<<" "<<kesi.theta<<" "<<kesi.velocity<<"\n"<<endl;
        }
    }
} ;

int main(){
    KineBicycleModel model1;
    model1.updatestate(0.5);
}
