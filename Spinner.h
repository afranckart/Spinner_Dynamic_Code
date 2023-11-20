//
// Définition de la class spinner
//

#ifndef UNTITLED8_SPINNER_H
#define UNTITLED8_SPINNER_H

#endif //UNTITLED8_SPINNER_H

#include <vector>
#include <iostream>
#include <string>

class Spinner
{
private:
    int m_index;
    //rotation
    double m_theta;

    std::vector<int> m_site;

public:
    Spinner(int,std::vector<int>,int);
    ~Spinner();
    void print(std::string);
    int orientation();
    void update_orientation(int);
    std::vector<int> charge();

};


/* struct pour passer les différent couplages */
struct couplage {
    float J; 
};
