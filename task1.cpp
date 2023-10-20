#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

class KineBicycleModel
{
private:
    double X, Y, theta, velocity, v_dot, steering_angle, L;
public:
    KineBicycleModel(){
        X=0;
        Y=0; 
        theta=0;
        velocity=0;
        L=1.5;
    }

    void updateState(double dt){
        int steps;
        double v_dot, steering_angle;

        ifstream input("input.txt");
        input>>steps;
        double controll[steps][2]; 

        for (int i = 0; i < steps; ++i) {
            input >> controll[i][0] >> controll[i][1];
        }

        for (int i = 0; i < steps; ++i) {
            v_dot = controll[i][0];
            steering_angle = controll[i][1];
            X += velocity*cos(theta)*dt;
            Y += velocity*cos(theta)*dt;
            theta += velocity*tan(steering_angle)/L;
            velocity += v_dot*dt;
            cout<<X<<" "<<Y<<" "<<theta<<" "<<velocity<<"\n"<<endl;
        }
    }
} ;

int main(){
    KineBicycleModel model1;
    model1.updateState(0.5);
}
