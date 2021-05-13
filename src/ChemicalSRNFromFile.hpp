#ifndef CHEMICALSRNFROMFILE_HPP_
#define CHEMICALSRNFROMFILE_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <iostream>

// reaction includes
#include "AbstractReactionSystem.hpp"
#include "AbstractReactionSystemFromFile.hpp"
#include "AbstractReaction.hpp"
#include "AbstractChemistry.hpp"
#include "ReactionTypeDatabase.hpp"
#include "ChemicalSrnModel.hpp"
 
class ChemicalSRNFromFile
{
protected:

    std::string mSrnReactionFilename;

    ChemicalSrnModel* mpChemicalSrnModel;

public:

    ChemicalSRNFromFile(std::string);

    virtual ~ChemicalSRNFromFile()
    {
    };

    void SetChemicalSrnModel(ChemicalSrnModel*);

    ChemicalSrnModel* GetChemicalSrnModel();

    std::string GetReactionFileName();

};

ChemicalSRNFromFile::ChemicalSRNFromFile(std::string reactionFilename) 
    :   mSrnReactionFilename(reactionFilename)
{
    AbstractReactionSystemFromFile* p_file_reaction_system = new  AbstractReactionSystemFromFile(reactionFilename);

    AbstractChemistry* this_cell_chemistry = p_file_reaction_system->GetSystemChemistry();
    unsigned numberOfChemicals = this_cell_chemistry->GetNumberChemicals();

    ChemicalSrnModel* pChemicalSrnModel = new ChemicalSrnModel(p_file_reaction_system);
    SetChemicalSrnModel(pChemicalSrnModel);

    mpChemicalSrnModel -> SetReactionSystem(p_file_reaction_system);
    mpChemicalSrnModel -> Initialise();
}

void ChemicalSRNFromFile::SetChemicalSrnModel(ChemicalSrnModel* pChemicalSrnModel)
{
    mpChemicalSrnModel = pChemicalSrnModel;
}

ChemicalSrnModel* ChemicalSRNFromFile::GetChemicalSrnModel()
{
    return mpChemicalSrnModel;
}

std::string ChemicalSRNFromFile::GetReactionFileName()
{
    return mSrnReactionFilename;
}



#endif