#ifndef InhomogenousParabolicPdeForCoupledOdeSystem_templated_HPP_
#define InhomogenousParabolicPdeForCoupledOdeSystem_templated_HPP_


#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"
#include <string>
#include "StateVariableRegister.hpp"
//#include "ChemicalDomainFieldForCellCoupling.hpp"
#include "ChemicalDomainField_templated.hpp"
#include "AbstractDomainField_templated.hpp"

// class defining the pde system neded when running chemical reaction-diffusion systems 

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class InhomogenousParabolicPdeForCoupledOdeSystemTemplated : public AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{

protected:

    // register for controlling the index values for variables to different ode systems
    // if a state variable isn't in this register then the pde doens't know about the variable, ghost state
    StateVariableRegister* mpStateVariableRegister;

    // domain for spatially dependent diffusion and ODEs
    AbstractDomainFieldTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* mpDomainField;
    // switch for domain field diffusion data structure
    bool mIsDomainDiffusionField = false;

    // constant diffusion rate for each of the variabels
    std::vector<double> mDiffusionRateConstantVector;

    // switch for non-constant diffusion rate
    bool mIsDomainDiffusionVector = false;


public:

    InhomogenousParabolicPdeForCoupledOdeSystemTemplated()
        :   AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>()
            
    {
    }

    InhomogenousParabolicPdeForCoupledOdeSystemTemplated(std::vector<double> diffusionVector, std::vector<std::string> nameVector = std::vector<std::string>())
        :   AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(),
            mDiffusionRateConstantVector(diffusionVector)
    {        
        mIsDomainDiffusionVector = true;

        if(nameVector.empty())
        {
            std::cout<<"empty name vector"<<std::endl;
            for(unsigned i=0; i<PROBLEM_DIM; i++)
            {
                nameVector.push_back(std::to_string(i));
            } 
        }

        StateVariableRegister* p_stateVariableRegister = new StateVariableRegister(nameVector);
          
        SetStateVariableRegister(p_stateVariableRegister);
    }

    InhomogenousParabolicPdeForCoupledOdeSystemTemplated(std::vector<double> diffusionVector, StateVariableRegister* p_stateVariableRegister)
        :   AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(),
            mDiffusionRateConstantVector(diffusionVector),
            mpStateVariableRegister(p_stateVariableRegister)
    {        
        mIsDomainDiffusionVector = true;
    }

    InhomogenousParabolicPdeForCoupledOdeSystemTemplated(AbstractDomainFieldTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* p_domainField)
        :   AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(),
            mpDomainField(p_domainField)
    {

        mIsDomainDiffusionField = true;
        if(mpDomainField->GetFieldType()=="ChemicalDomainFieldTemplated")
        //if(dynamic_cast<ChemicalDomainFieldTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*>(mpDomainField)!= nullptr)
        {

            SetStateVariableRegister(dynamic_cast<ChemicalDomainFieldTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*>(mpDomainField) ->GetStateVariableVector());
        }

        //if(dynamic_cast<ChemicalDomainFieldForCellCoupling<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*>(mpDomainField)!= nullptr)
        if(mpDomainField->GetFieldType()=="ChemicalDomainFieldForCellCoupling")
        {

            SetStateVariableRegister(mpDomainField ->GetStateVariableVector());
            //SetStateVariableRegister(dynamic_cast<ChemicalDomainFieldForCellCoupling<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*>(mpDomainField) ->GetStateVariableVector());
        }
        
    }

    virtual ~InhomogenousParabolicPdeForCoupledOdeSystemTemplated()
    {   
    }

    virtual double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rX, unsigned index)
    {
        return 1.0;
    }
    /*
    virtual double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX, c_vector<double,PROBLEM_DIM>& rU, std::vector<double>& rOdeSolution, unsigned pdeIndex)
    {
        // the interpolated ode solution at a point rX, will be 0 if the surrounding nodes do not process the variable
        return rOdeSolution[pdeIndex];
    }
    */
    virtual double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX, c_vector<double,PROBLEM_DIM>& rU, std::vector<double>& rOdeSolution, unsigned pdeIndex)
    {

        std::string odeLabel = mpDomainField ->ReturnNodeOdeLabelAtPosition(rX.rGetLocation());

        std::vector<std::string> node_labels = mpDomainField ->GetNodeLabels();

        std::vector<double> rUvec(PROBLEM_DIM,0.0);
        for(unsigned i=0;i<PROBLEM_DIM;i++)
        {
            rUvec[i] = rU(i);
        }

        
        AbstractInhomogenousOdeSystemForCoupledPdeSystem* p_system = dynamic_cast<ChemicalDomainFieldTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*>(mpDomainField)->GetOdeSystem()[0];
        AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem* p_system_cast = dynamic_cast<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*>(p_system);
        std::vector<AbstractReaction*>  p_reactionVector = p_system_cast ->GetReactionSystem() -> GetReactionVector();
        
        if(std::pow((rX[0]-50),2)+std::pow((rX[1]-50),2) > 400)
        {

            
        }
        else
        {
           
        
        for(unsigned i=0; i<node_labels.size();i++)
        {
            if(node_labels[i] == odeLabel)
            {
                // found the correct chemical ode
                dynamic_cast<ChemicalDomainFieldTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*>(mpDomainField)->GetOdeSystem()[i] -> EvaluateYDerivatives(1.0, rUvec, rOdeSolution);

                return rUvec[pdeIndex];
            }
        }
        }
        return rUvec[pdeIndex];
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

    virtual std::string GetPdeOdeSystemType()
    {
        return "InhomogenousParabolicPdeForCoupledOdeSystemTemplated";
    }


    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex, Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)
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

    AbstractDomainFieldTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* GetDomainField()
    {
        return mpDomainField;
    }

    void SetDomainField(AbstractDomainFieldTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* p_domainField)
    {
        mIsDomainDiffusionField = true;
        mpDomainField = p_domainField;
    }
};


#endif 