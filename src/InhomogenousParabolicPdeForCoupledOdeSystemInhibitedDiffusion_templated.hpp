#ifndef InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusion_templated_HPP_
#define InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusion_templated_HPP_


#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"
#include <string>
#include "StateVariableRegister.hpp"
//#include "ChemicalDomainFieldForCellCoupling.hpp"
#include "ChemicalDomainField_templated.hpp"
#include "AbstractDomainField_templated.hpp"

// class defining the pde system neded when running chemical reaction-diffusion systems 

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusionTemplated : public InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{

public:

    InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusionTemplated()
        :   InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>()
            
    {
    }

    InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusionTemplated(std::vector<double> diffusionVector, std::vector<std::string> nameVector = std::vector<std::string>())
        :   InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(diffusionVector,nameVector)
  
    {        
    }

    InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusionTemplated(std::vector<double> diffusionVector, StateVariableRegister* p_stateVariableRegister)
        :   InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(diffusionVector,p_stateVariableRegister)

    {        
    }

    InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusionTemplated(AbstractDomainFieldTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* p_domainField)
        :   InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(p_domainField)
    {
    }

    virtual ~InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusionTemplated()
    {   
    }


    virtual std::string GetPdeOdeSystemType()
    {
        return "InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusionTemplated";
    }

    virtual double DiffusionFunction(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex, Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)
    {
        // virtual function, use constant value as base case
        if(this->mIsDomainDiffusionVector == true)
        {
            // use the pdeIndex and stateVaribale register to search for a species-domain pair to retrive value?
            return InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetDiffusionRateConstantByIndex(pdeIndex);
        }
        else if(this->mIsDomainDiffusionField == true)
        {
            return this->mpDomainField -> GetDiffusionValueBasedOnPoint(rX,pdeIndex);
        }
        else
        {
            return 0.0;
        }
    }

    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex, Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)
    {
        // get a pdeIndex and a point in space,  assume isotropic diffusion
        assert(pdeIndex < PROBLEM_DIM);

        c_matrix<double, SPACE_DIM, SPACE_DIM> diffusion_term;
        if(std::pow((rX[0]-50),2)+std::pow((rX[1]-50),2) <= 400)
        {

            diffusion_term = (0.01*DiffusionFunction(rX, pdeIndex,pElement))*identity_matrix<double>(SPACE_DIM);
        }
        else
        {
            diffusion_term = DiffusionFunction(rX, pdeIndex,pElement)*identity_matrix<double>(SPACE_DIM);
        }


        return diffusion_term;
    }
};


#endif 