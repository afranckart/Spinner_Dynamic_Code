#include "Spinner.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


Spinner::Spinner(int i_index,std::vector<int> i_site,int i_theta)
{
    m_index=i_index;
    m_theta=i_theta;
    m_site=i_site;
}
Spinner::~Spinner()
{

}

void Spinner::print (std::string add)
{
    std::ostringstream fileName;
    fileName<<add;
    std::ofstream position;
    position.open(fileName.str(),std::ios::app);
    position<<m_index<<"\t";
    for(int i=0; i<m_site.size();i++)
    {
        position<<m_site[i]<<"\t";
    }
    position<<m_theta<<std::endl;
    position.close();
}

int Spinner::orientation()
{
    return m_theta;
}
void Spinner::update_orientation(int phi)
{
    m_theta=phi;
}

std::vector <int> Spinner::charge()
{
    return m_site;
}



