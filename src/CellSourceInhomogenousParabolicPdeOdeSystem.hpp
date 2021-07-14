#ifndef AVERAGEDSOURCEINHOMOGENOUSPARABOLICPDEODESYSTEM_HPP_
#define AVERAGEDSOURCEINHOMOGENOUSPARABOLICPDEODESYSTEM_HPP_

#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"
#include <math.h> 
#include <string>
#include "StateVariableRegister.hpp"
#include "AbstractDomainField.hpp"
#include "ChemicalDomainField.hpp"
#include "ChemicalDomainFieldForCellCoupling.hpp"

#include "ApoptoticCellProperty.hpp"
#include "ChemicalCellProperty.hpp"
#include "MembraneCellProperty.hpp"
#include "TransportCellProperty.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class CellSourceInhomogenousParabolicPdeOdeSystem : public InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:
    using InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetPdeOdeSystemType;
    using InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeSourceTerm;
    using InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::DiffusionFunction;

protected:

    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& mrCellPopulation;

    /** Vector of averaged cell densities on elements of the coarse mesh. */
    std::vector<double> mCellDensityOnCoarseElements;

    bool mIsDiffusionCellDensityScaled = true;

public:

    CellSourceInhomogenousParabolicPdeOdeSystem(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
        :   InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(),
        mrCellPopulation(rCellPopulation)
            
    {
    }

    CellSourceInhomogenousParabolicPdeOdeSystem(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                    std::vector<double> diffusionVector)
        :   InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(diffusionVector),
        mrCellPopulation(rCellPopulation)
    {        
        this->mIsDomainDiffusionVector = true;   
    }

    CellSourceInhomogenousParabolicPdeOdeSystem(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                    AbstractDomainField* p_domainField)
        :   InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(p_domainField),
        mrCellPopulation(rCellPopulation)
    {
    }

    virtual ~CellSourceInhomogenousParabolicPdeOdeSystem()
    {   
    }


    virtual std::string GetPdeOdeSystemType()
    {
        return "AveragedSourceInhomogenousParabolicPdeOdeSystem";
    }

    virtual double DiffusionFunction(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex, Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)
    {
        // virtual function, use constant value as base case
        double diffusion_rate=0.0;
        if(mIsDomainDiffusionVector == true)
        {
            // use the pdeIndex and stateVaribale register to search for a species-domain pair to retrive value?
            diffusion_rate = GetDiffusionRateConstantByIndex(pdeIndex);
        }
        else if(mIsDomainDiffusionField == true)
        {
            diffusion_rate = mpDomainField -> GetDiffusionValueBasedOnPoint(rX,pdeIndex);
        }

        if(mIsDiffusionCellDensityScaled && pElement != NULL)
        {
            diffusion_rate = diffusion_rate * ElementDensityScaleFunction(mCellDensityOnCoarseElements[pElement->GetIndex()])
        }
       
    }

    virtual double ElementDensityScaleFunction(double cell_density)
    {
        // function to scale the base diffusion in the sub domain through consideration of the presence of cells
        // modeling the role of extracellular polymeric substance etc. Default to expenential decrease in diffusion
        // rate with increaing cell density in the FeMesh element
        return exp(-cell_density);
    }

    const AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rGetCellPopulation() const
    {
        return mrCellPopulation;
    }


    void SetupSourceTerms(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rFeMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap) // must be called before solve
    {
        // calculate cell density and store as mCellDensityOnCoarseElements for each of the elements of the FeMesh
        CalculateCellDensity(rFeMesh,pCellPdeElementMap);
    }

    void CalculateCellDensity(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rFeMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap)
    {
        // calculates the cell density of non-apoptotic cells on the pde FeMesh

        // Allocate memory
        mCellDensityOnCoarseElements.resize(rFeMesh.GetNumElements());
        for (unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
        {
            mCellDensityOnCoarseElements[elem_index] = 0.0;
        }

        // Loop over cells, find which coarse element it is in, and add 1 to mSourceTermOnCoarseElements[elem_index]
        for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = mrCellPopulation.Begin();
            cell_iter != mrCellPopulation.End();
            ++cell_iter)
        {
            unsigned elem_index = 0;
            const ChastePoint<SPACE_DIM>& r_position_of_cell = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

            if (pCellPdeElementMap != nullptr)
            {
                elem_index = (*pCellPdeElementMap)[*cell_iter];
            }
            else
            {
                elem_index = rFeMesh.GetContainingElementIndex(r_position_of_cell);
            }

            // Update element map if cell has moved
            //bool cell_is_apoptotic = cell_iter->template HasCellProperty<ApoptoticCellProperty>();


            bool cell_is_apoptotic = cell_iter->template rGetCellPropertyCollection().HasProperty<ApoptoticCellProperty>();

            if (!cell_is_apoptotic)
            {
                mCellDensityOnCoarseElements[elem_index] += 1.0;
            }
        }

        // Then divide each entry of mSourceTermOnCoarseElements by the element's area
        c_matrix<double, ELEMENT_DIM,SPACE_DIM> jacobian;
        double det;
        for (unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
        {
            rCoarseMesh.GetElement(elem_index)->CalculateJacobian(jacobian, det);
            mCellDensityOnCoarseElements[elem_index] /= rCoarseMesh.GetElement(elem_index)->GetVolume(det);
        }
    }

    double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX, c_vector<double,PROBLEM_DIM>& rU, Element<ELEMENT_DIM,SPACE_DIM>* pElement, std::vector<double>& rOdeSolution, unsigned pdeIndex)
    {
        // term to utilise the transport property of the cells to determine what the constant source/sink term with respect ot the 
        // pde FeMesh, term will act upon the cell data with the opposite polarity (remove if acts as positive source to pde mesh etc)

        // add to the ode solution interpolated from the node
        double ode_contribution_at_rX = InhomogenousParabolicPdeForCoupledOdeSystemTemplated::ComputeSourceTerm(rX, rU,rOdeSolution,pdeIndex); 
        // input args different? pElement need to add as argument Element<ELEMENT_DIM,SPACE_DIM>* pElement=nullptr?
        
        double cell_contribution_at_rX = 0;


        return ode_contribution_at_rX + cell_contribution_at_rX;
    }



};


#endif 