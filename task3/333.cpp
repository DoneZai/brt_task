#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;


float sgn(float v) {
    return (v > 0) - (v < 0);
}


class SimpleDynamicBycicleModel {
private:
    struct x {
        float X;
        float Y;
        float yaw;
        float vx;
        float vy;
        float r;
    };

    struct u {
        float throttle;
        float steeringAngle;
        float brakes;
    };

    x old = { 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0 };
    u instant;
    float t = 0.0f;

    const float dt = 0.001f;
    const float m = 300.0f;
    const float lf = 0.721f;
    const float lr = 0.823f;
    const float Cm = 3600.0f;
    const float Crr = 200.0f;
    const float Cd = 1.53f;
    const float Cbf = 5411.0f;
    const float Cbr = 2650.0f;
    const float Cx = 20000.0f;
    const float Iz = 134.0f;
    float ar = 0.0f;
    float af = 0.0f;

    float Fdrv() {
        return instant.throttle * Cm * instant.throttle;
    }

    float Frrr() {
        return Crr * tanh(old.vx);
    }

    float Frrf() {
        return Frrr();
    }

    float Fdrag() {
        return Cd * old.vx * old.vx;
    }

    float Fbf() {
        return instant.brakes * Cbf * tanh(old.vx);
    }

    float Fbr() {
        return instant.brakes * Cbr * tanh(old.vx);
    }

    float Fry() {
        //if (ar != ar) {
        //    return 0.0f;
        //}
        //return Cx * ar;
        return magicFy(ar, 0, m / 4);
    }

    float Ffy() {
        //if (af != af) {
        //    return 0.0f;
        //}
        //return Cx * af;
        return magicFy(af, 0, m / 4);
    }

float magicFy(float alpha, float gamma, float Fz) {
    float pDy1 = +2.716E+000;
    float pDy2 = -5.444E-001;
    float lgammay = 0.75;
    float gammay = gamma * lgammay;
    float pDy3 = +5.190E+000;
    float lmuy = 1;
    float Fz0 = +8.000e+002;
    float lFz0 = 1;
    float dfz = (Fz - Fz0 * lFz0) / (Fz0 * lFz0);
    float muy = (pDy1 + pDy2 * dfz) * (1 - pDy3 * gammay * gammay) * lmuy;

    float pCy1 = +1.434E+000;
    float lCy = 1;
    float Cy = pCy1 * lCy;

    // float Fz;
    float Dy = muy * Fz;

    float pKy1 = -5.322E+001;
    float pKy2 = +2.060E+000;
    float pKy3 = +8.336E-001;
    float lKy = 1;
    float Ky = pKy1 * Fz0 * sin(2 * atan2(Fz, (pKy2 * Fz0 * lFz0))) * (1 - pKy3 * abs(gammay)) * lFz0 * lKy;
    float By = Ky / (Cy * Dy);

    float PHy1 = +0.000e+000;
    float PHy2 = +0.000e+000;
    float PHy3 = -2.030E-002;
    float lHy = 1;
    float SHy = (PHy1 + PHy2 * dfz) * lHy + PHy3 * gammay;

    float pEy1 = -4.869E-001;
    float pEy2 = -1.487E+000;
    float pEy3 = +6.282E-002;
    float pEy4 = +1.154E+000;
    float alphay = alpha + SHy;
    float lEy = 0.55;
    float Ey = (pEy1 + pEy2 * dfz) * (1 - (pEy3 + pEy4 * gammay) * sgn(alphay)) * lEy;

    float pVy1 = +0.000e+000;
    float pVy2 = +0.000e+000;
    float pVy3 = -2.713E+000;
    float pVy4 = -1.517E+000;
    float lVy = 1;
    float SVy = Fz * ((pVy1 + pVy2 * dfz) * lVy + (pVy3 + pVy4 * dfz) * gammay) * lmuy;

    float Fy0 = Dy * sin(Cy * atan(By * alphay - Ey * (By * alphay - atan(By * alphay)))) + SVy;
    return Fy0;
}

float magicFx(float alpha, float gamma, float Fz) {
    float pDy1 = +2.716E+000;
    float pDy2 = -5.444E-001;
    float lgammay = 0.75;
    float gammay = gamma * lgammay;
    float pDy3 = +5.190E+000;
    float lmuy = 1;
    float Fz0 = +8.000e+002;
    float lFz0 = 1;
    float dfz = (Fz - Fz0 * lFz0) / (Fz0 * lFz0);
    float muy = (pDy1 + pDy2 * dfz) * (1 - pDy3 * gammay * gammay) * lmuy;

    float pCy1 = +1.434E+000;
    float lCy = 1;
    float Cy = pCy1 * lCy;

    // float Fz;
    float Dy = muy * Fz;

    float pKy1 = -5.322E+001;
    float pKy2 = +2.060E+000;
    float pKy3 = +8.336E-001;
    float lKy = 1;
    float Ky = pKy1 * Fz0 * sin(2 * atan2(Fz, (pKy2 * Fz0 * lFz0))) * (1 - pKy3 * abs(gammay)) * lFz0 * lKy;
    float By = Ky / (Cy * Dy);

    float PHy1 = +0.000e+000;
    float PHy2 = +0.000e+000;
    float PHy3 = -2.030E-002;
    float lHy = 1;
    float SHy = (PHy1 + PHy2 * dfz) * lHy + PHy3 * gammay;

    float pEy1 = -4.869E-001;
    float pEy2 = -1.487E+000;
    float pEy3 = +6.282E-002;
    float pEy4 = +1.154E+000;
    float alphay = alpha + SHy;
    float lEy = 0.55;
    float Ey = (pEy1 + pEy2 * dfz) * (1 - (pEy3 + pEy4 * gammay) * sgn(alphay)) * lEy;

    float pVy1 = +0.000e+000;
    float pVy2 = +0.000e+000;
    float pVy3 = -2.713E+000;
    float pVy4 = -1.517E+000;
    float lVy = 1;
    float SVy = Fz * ((pVy1 + pVy2 * dfz) * lVy + (pVy3 + pVy4 * dfz) * gammay) * lmuy;

    float Fy0 = Dy * sin(Cy * atan(By * alphay - Ey * (By * alphay - atan(By * alphay)))) + SVy;
    return Fy0;
}

    float Flateral() {
        float k = cos(atan(old.vy / old.vx));
        float Faero = 0.0f;
        if (k == k) {
            Faero = Fdrag() * k;
        }
        return -Frrr() - Fbr() + Fdrv()
            - Faero
            - Fbf() * cos(instant.steeringAngle) - Frrf() * cos(instant.steeringAngle) - Ffy() * sin(instant.steeringAngle);
    }

    float Ftransversal() {
        float k = sin(atan(old.vy / old.vx));
        float Faero = 0.0f;
        if (k == k) {
            Faero = Fdrag() * k;
        }
        return Fry() - Faero - (Fbf() + Frrf()) * sin(instant.steeringAngle) + Ffy() * cos(instant.steeringAngle);
    }

    float Ftotal() {
        cout << "\tFlateral() = " << Flateral() << endl;
        cout << "\tFtransversal() = " << Ftransversal() << endl;
        return sqrt(pow(Flateral(), 2) + pow(Ftransversal(), 2));
    }


    float L() {  // Angular momentum
        //cout << "\tFfy() = " << Ffy() << endl;
        //cout << "\tFbf() = " << Fbf() << endl;
        //cout << "\tFrrf() = " << Frrf() << endl;
        //cout << "\tFry() = " << Fry() << endl;
        //cout << "\tFrrf() = " << Frrf() << endl;
        return -lr * Fry() + lf * (Ffy() * cos(instant.steeringAngle) - Fbf() * sin(instant.steeringAngle) - Frrf() * sin(instant.steeringAngle));
    }

public:
    float getX() const { return old.X; }
    float getY() const { return old.Y; }
    float getyaw() const { return old.yaw; }
    float getvx() const { return old.vx; }
    float getvy() const { return old.vy; }
    float getr() const { return old.r; }
    float gett() const { return t; }

    void update(float throttle, float steeringAngle, float brakes) {
        instant.throttle = throttle;
        instant.steeringAngle = steeringAngle;
        instant.brakes = brakes;

        float X_p1 = old.X + (old.vx * cos(old.yaw) + old.vy * sin(old.yaw)) * dt;
        float Y_p1 = old.Y + (old.vx * sin(old.yaw) + old.vy * cos(old.yaw)) * dt;
        float yaw_p1 = old.yaw + old.r * dt;
        float vx_p1 = old.vx + (Flateral() / m + old.vy * old.r) * dt;
        float vy_p1 = old.vy + (Ftransversal() / m - old.vx * old.r) * dt;
        float r_p1 = L() / Iz;

        t += dt;
        old = { X_p1, Y_p1, yaw_p1, vx_p1, vy_p1, r_p1 };
    }
};


ostream& operator<<(ostream& s, const SimpleDynamicBycicleModel& sdbm) {
    return (s << sdbm.getX() <<
        " " << sdbm.getY() << " " << sdbm.getyaw() <<
        " " << sdbm.getvx() << " " << sdbm.getvy() << " " << sdbm.getr() );
}


int main()
{
    int iterations_by_one_step = 50;
    SimpleDynamicBycicleModel A;
    std::ifstream in("input.txt");
    ofstream outputFile("output.txt");
    if (in.is_open())
    {
        int n;
        in >> n;
        for (int i = 0; i < n; i++) {
            float a, sa, br;
            in >> a >> sa >> br;
            for (int j = 0; j < iterations_by_one_step; j++) {
                A.update(a, sa, br);
                // cout << A << endl;
                outputFile<<A<<"\n"<<endl;
                
            }
        }
    }
    in.close();
    outputFile.close();
}