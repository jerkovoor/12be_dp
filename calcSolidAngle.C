#include <fstream>
#include <sstream>

void calcSolidAngle() {
    FILE* file = fopen("csRing01.out", "w");

    Float_t cmEnergy = 5.0;

    TH1F *h = new TH1F("lab_ik", "lab_ik", 360, 0, 180);
    TH1F *h2 = new TH1F("cm_fk", "cm_fk", 360, 0, 180);

    Float_t z = 85.;


    Float_t angle;
    Float_t sum = 0;

    Float_t 
    Float_t elementSizeAngle = 0.1;
    
    for(Float_t dangle = 0.; dangle < 0.733; dangle += elementSize) {
        
    }

    for(Float_t dx = 99; dx <= 149; dx += elementSize) {
        for(Float_t dy = 10.8; dy <= 60.8; dy += elementSize) {
            Float_t dr = sqrt(dx*dx + dy*dy + z*z);
            angle = acos(z/dr)/3.14159*180;
            sum += elementSize*elementSize*z/(dr*dr*dr);
            h->Fill(angle);
            h2->Fill(180. - 2.*angle);
        }
    }
    sum *= 8.;

    if(sum == 0.) continue;

    Float_t change_bin_content = sum;

    printf("%f %f %f %f\n", cmEnergy, change_bin_content, h->GetMean(), h2->GetMean());
    fprintf(file, "%f %f %f %f\n", cmEnergy, change_bin_content, h->GetMean(), h2->GetMean());

    delete h;
    delete h2;

    fflush(file);
    fclose(file);
}
