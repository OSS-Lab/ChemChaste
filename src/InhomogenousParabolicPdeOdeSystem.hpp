#ifndef INHOMOGENOUSPARABOLICPDEODESYSTEM_HPP_
#define INHOMOGENOUSPARABOLICPDEODESYSTEM_HPP_

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"
#include <string>
#include "StateVariableRegister.hpp"
#include "AbstractDomainField.hpp"
#include "ChemicalDomainField.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class InhomogenousParabolicPdeOdeSystem : public AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{

protected:

    // register for controlling the index values for variables to different ode systems
    // if a state variable isn't in this register then the pde doens't know about the variable, ghost state
    StateVariableRegister* mpStateVariableRegister;

    // domain for spatially dependent diffusion and ODEs
    AbstractDomainField* mpDomainField;
    // switch for domain field diffusion data structure
    bool mIsDomainDiffusionField = false;

    // constant diffusion rate for each of the variabels
    std::vector<double> mDiffusionRateConstantVector;

    // switch for non-constant diffusion rate
    bool mIsDomainDiffusionVector = false;


public:

    InhomogenousParabolicPdeOdeSystem()
        :   AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>()
            
    {
    }

    InhomogenousParabolicPdeOdeSystem(std::vector<double> diffusionVector)
        :   AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(),
            mDiffusionRateConstantVector(diffusionVector)
    {        
        mIsDomainDiffusionVector = true;   
    }

    InhomogenousParabolicPdeOdeSystem(AbstractDomainField* p_domainField)
        :   AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(),
            mpDomainField(p_domainField)
    {
        mIsDomainDiffusionField = true;
        if(mpDomainField->GetFieldType()=="ChemicalDomainField")
        {
            SetStateVariableRegister(dynamic_cast<ChemicalDomainField*>(mpDomainField) ->GetStateVariableVector());
        }

    }

    virtual ~InhomogenousParabolicPdeOdeSystem()
    {   
    }

    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rX, unsigned index)
    {
        return 1.0;
    }

    double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX, c_vector<double,PROBLEM_DIM>& rU, std::vector<double>& rOdeSolution, unsigned pdeIndex)
    {
        // the interpolated ode solution at a point rX, will be 0 if the surrounding nodes do not process the variable
        return rOdeSolution[pdeIndex];
    }

    void SetDiffusionRateConstantVector(std::vector<double> diffusionRateConstantVector)
    {
        mIsDomainDiffusionVector = true;
        
        mDiffusionRateConstantVector = diffusionRateConstantVector;
    }

    std::vector<double> GetDiffusionRateConstantVector()
    {
        return mDiffusionRateConstantVector;
    }

    void SetDiffusionRateConstantByIndex(double diffusionValue, unsigned index)
    {
        mDiffusionRateConstantVector[index] = diffusionValue;
    }

    double GetDiffusionRateConstantByIndex(unsigned index)
    {
        assert(mIsDomainDiffusionVector == true);
        assert(index < PROBLEM_DIM);
        return mDiffusionRateConstantVector[index];
    }


    virtual double DiffusionFunction(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex, Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)
    {
        // virtual function, use constant value as base case
        if(mIsDomainDiffusionVector == true)
        {
            // use the pdeIndex and stateVaribale register to search for a species-domain pair to retrive value?
            return GetDiffusionRateConstantByIndex(pdeIndex);
        }
        else if(mIsDomainDiffusionField == true)
        {
            return mpDomainField -> GetDiffusionValueBasedOnPoint(rX,pdeIndex);
        }
        else
        {
            return 0.0;
        }
    }




    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex, Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)
    {
        // get a pdeIndex and a point in space,  assume isotropic diffusion
        assert(pdeIndex < PROBLEM_DIM);


        c_matrix<double, SPACE_DIM, SPACE_DIM> diffusion_term;

        diffusion_term = DiffusionFunction(rX, pdeIndex,pElement)*identity_matrix<double>(SPACE_DIM);

        return diffusion_term;
    }


    void SetStateVariableRegister(StateVariableRegister* p_stateVariableRegister)
    {
        mpStateVariableRegister = p_stateVariableRegister;
    };

    StateVariableRegister* GetStateVariableRegister()
    {
        return mpStateVariableRegister;
    }

    

    void SetIsDomainDiffusionVector(bool isDomainDiffusionVector)
    {
        mIsDomainDiffusionVector = isDomainDiffusionVector;
    }

    bool GetIsDomainDiffusionVector()
    {
        return mIsDomainDiffusionVector;
    }

    void SetIsDomainDiffusionField(bool isDomainDiffusionField)
    {
        mIsDomainDiffusionField = isDomainDiffusionField;
    }

    bool GetIsDomainDiffusionField()
    {
        return mIsDomainDiffusionField;
    }

    AbstractDomainField* GetDomainField()
    {
        return mpDomainField;
    }

    void SetDomainField(AbstractDomainField* p_domainField)
    {
        mIsDomainDiffusionField = true;
        mpDomainField = p_domainField;
    }


    std::string GetPdeOdeSystemType()
    {
        return "InhomogenousParabolicPdeOdeSystem";
    }

};


#endif 