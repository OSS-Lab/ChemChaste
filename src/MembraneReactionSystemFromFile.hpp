#ifndef MEMBRANEREACTIONSYSTEMFROMFILE_HPP_
#define MEMBRANEREACTIONSYSTEMFROMFILE_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <iostream>
#include <fstream>

// reaction includes
#include "AbstractMembraneReactionSystem.hpp"
#include "AbstractMembraneReaction.hpp"
#include "AbstractChemistry.hpp"
#include "ReactionTypeDatabase.hpp"


class MembraneReactionSystemFomFile : public AbstractMembraneReactionSystem
{
protected:

    //inherit:
        // AbstractChemistry* mpBulkChemistry;

        // AbstractChemistry* mpCellChemistry;

        // std::vector<AbstractTransportReaction*> mpReactionVector;

        // unsigned mNumberOfReactions;

    std::string mInputFileName;
    std::string mDomain;
    std::string mStringDelimiter = "";
    std::string mDataDelimiter = ";";
    std::string mForwardDelimiter = "->";
    std::string mReverseDelimiter = "<-";
    std::string mReverDelimiter = "<->";
    std::string mSpeciesSeparator = " + ";
    std::string mTypeDelimiter = " : ";
    std::string mReactionDelimiter = " | ";

    bool mIsReversedReaction = false;
    
    std::vector<bool> mIsCoupledReaction = std::vector<bool>();

    std::vector<std::vector<unsigned>> mNumberSubstrates = std::vector<std::vector<unsigned>>(); 

    std::vector<std::vector<unsigned>> mNumberProducts = std::vector<std::vector<unsigned>>(); 

public:

    MembraneReactionSystemFomFile(std::string);

    virtual ~MembraneReactionSystemFomFile()
    {
    };

    virtual void SetFileDeliminator();

    void FormReactionSystemObjectFromTuple(std::vector<std::tuple<std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string>>);

    std::vector<AbstractChemical*> FromChemicalNameVectorToAbstractChemicalVector(std::vector<std::string>);

    void ParseSystemChemistry(std::vector<std::string>, std::vector<std::string>);

    std::vector<std::tuple<std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string>> ReactionSystemFromFile();

    bool TestReversibility(std::string);

    std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>> ParseReactionString(std::string);

    std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>> AddCoupledReactions(std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>>, std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>>);

    std::string ParseReactionInformation(std::string);

    std::string GetFileName();


    AbstractMembraneReactionSystem* GetReactionSystem();

    std::vector<std::vector<unsigned>> GetNumberOfSubstrates();

    std::vector<std::vector<unsigned>> GetNumberOfProducts();

    void SetNumberOfSubstrates(std::vector<std::vector<unsigned>>);

    void SetNumberOfProducts(std::vector<std::vector<unsigned>>);

    std::vector<unsigned> ReturnSubstratesForReactionIndex(unsigned);

    std::vector<unsigned> ReturnProductsForReactionIndex(unsigned);

};

MembraneReactionSystemFomFile::MembraneReactionSystemFomFile(std::string InputFileName)
    : mInputFileName(InputFileName)
{
    
    std::vector<AbstractMembraneReaction*> mpReactionVector = std::vector<AbstractMembraneReaction*>();
    

    SetFileDeliminator();
    std::vector<std::tuple<std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string>> system = ReactionSystemFromFile();
    
    std::vector<std::string> bulkChemicalNames;

    std::vector<std::string> cellChemicalNames;

    unsigned numberOfBulkSubstrates=0;
    unsigned numberOfCellSubstrates=0;
    unsigned numberOfBulkProducts=0;
    unsigned numberOfCellPRoducts=0;
    
    // for each reaction in the system, parse the bulk and cell chemistry
    for(unsigned reaction=0; reaction<system.size(); reaction++)
    {
        numberOfBulkSubstrates = ReturnSubstratesForReactionIndex(reaction)[0];
        numberOfBulkProducts = ReturnProductsForReactionIndex(reaction)[0];
        for(unsigned species=0; species<numberOfBulkSubstrates; species++)
        {
            // for the substrates
            bulkChemicalNames.push_back(std::get<2>(system[reaction])[species]);
        }
        for(unsigned species=0; species<numberOfBulkProducts; species++)
        {
            // for the products
            bulkChemicalNames.push_back(std::get<3>(system[reaction])[species]);
        }
        if(mIsCoupledReaction[reaction]==true)
        {
            // also has separate cell chemistry
            numberOfCellSubstrates = ReturnSubstratesForReactionIndex(reaction)[1];
            numberOfCellProducts = ReturnProductsForReactionIndex(reaction)[1];
            for(unsigned species=0; species<numberOfCellSubstrates; species++)
            {
                // for the substrates
                cellChemicalNames.push_back(std::get<2>(system[reaction])[species]);
            }
            for(unsigned species=0; species<numberOfCellProducts; species++)
            {
                // for the products
                cellChemicalNames.push_back(std::get<3>(system[reaction])[species]);
            }
        }
    }
    
    // convert chemical names into the chemistries
    ParseSystemChemistry(bulkChemicalNames,cellChemicalNames);
    
    // form the vector of membrane reaction systems and initilaise throuigh referencing the membrane tablet
    FormReactionSystemObjectFromTuple(system);

}

void MembraneReactionSystemFomFile::SetFileDeliminator()
{
    mStringDelimiter = "";
    mDataDelimiter = ";";
    mForwardDelimiter = "->";
    mReverseDelimiter = "<-";
    mReverDelimiter = "<->";
    mSpeciesSeparator = " + ";
    mTypeDelimiter = " : ";
    mReactionDelimiter = " | ";
}

void MembraneReactionSystemFomFile::FormReactionSystemObjectFromTuple(std::vector<std::tuple<std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string>> system)
{
    // for each reaction whose data is in the tuple, form the corresponding reaction class
    // denote in the reaction class the function necessary to parse reaction information

    SetNumberOfReactions(system.size());

    unsigned numberOfBulkSubstrates=0;
    unsigned numberOfCellSubstrates=0;
    unsigned numberOfBulkProducts=0;
    unsigned numberOfCellPRoducts=0;

    for( unsigned reaction =0; reaction<mNumberOfReactions; reaction++)
    {
        AbstractMembraneReaction* p_reaction = new AbstractMembraneReaction();

        std::vector<std::string> bulk_substrate_name = std::vector<std::string>();
        std::vector<unsigned> bulk_substrate_stoich = std::vector<unsigned>();
        std::vector<std::string> bulk_product_name = std::vector<std::string>();
        std::vector<unsigned> bulk_product_stoich = std::vector<unsigned>();
        std::vector<std::string> cell_substrate_name = std::vector<std::string>();
        std::vector<unsigned> cell_substrate_stoich = std::vector<unsigned>();
        std::vector<std::string> cell_product_name = std::vector<std::string>();
        std::vector<unsigned> cell_product_stoich = std::vector<unsigned>();

        numberOfBulkSubstrates = ReturnSubstratesForReactionIndex(reaction)[0];
        numberOfBulkProducts = ReturnProductsForReactionIndex(reaction)[0];

        for(unsigned i=0; i<numberOfBulkSubstrates; i++)
        {
            bulk_substrate_name.push_back(std::get<2>(system[reaction])[i]);
            bulk_substrate_stoich.push_back(std::get<4>(system[reaction])[i]);
        }
        for(unsigned i=0; i<numberOfBulkProducts;i++)
        {
            bulk_product_name.push_back(std::get<3>(system[reaction])[i]);
            bulk_product_stoich.push_back(std::get<5>(system[reaction])[i]);
        }

        std::vector<AbstractChemical*> bulk_substrates_chemical_vector = FromChemicalNameVectorToAbstractChemicalVector(bulk_substrate_name);
        std::vector<AbstractChemical*> bulk_products_chemical_vector = FromChemicalNameVectorToAbstractChemicalVector(bulk_product_name);
        std::vector<AbstractChemical*> cell_substrates_chemical_vector = std::vector<AbstractChemical*>();
        std::vector<AbstractChemical*> cell_products_chemical_vector = std::vector<AbstractChemical*>();

        if(mIsCoupledReaction[reaction]==true)
        {

            numberOfCellSubstrates = ReturnSubstratesForReactionIndex(reaction)[1];
            numberOfCellProducts = ReturnProductsForReactionIndex(reaction)[1];

            for(unsigned i=0; i<numberOfCellSubstrates;i++)
            {
                cell_substrate_name.push_back(std::get<2>(system[reaction])[i+numberOfBulkSubstrates]);
                cell_substrate_stoich.push_back(std::get<4>(system[reaction])[i+numberOfBulkSubstrates]);
            }
            for(unsigned i=0; i<numberOfCellProducts;i++)
            {
                cell_product_name.push_back(std::get<3>(system[reaction])[i+numberOfBulkProducts]);
                cell_product_stoich.push_back(std::get<5>(system[reaction])[i+numberOfBulkProducts]);
            }
            
            cell_substrates_chemical_vector = FromChemicalNameVectorToAbstractChemicalVector(cell_substrate_name);
            cell_products_chemical_vector = FromChemicalNameVectorToAbstractChemicalVector(cell_product_name);
        }
        

        MembraneTablet(p_reaction,std::get<0>(system[reaction]),bulk_substrates_chemical_vector,bulk_products_chemical_vector,cell_substrates_chemical_vector,cell_products_chemical_vector,bulk_substrate_stoich,bulk_product_stoich,cell_substrate_stoich,cell_product_stoich,std::get<6>(system[reaction]),std::get<1>(system[reaction]),mpBulkChemistry,mpCellChemistry);
        
        mpReactionVector.push_back(p_reaction);
    }
}


std::vector<AbstractChemical*> MembraneReactionSystemFomFile::FromChemicalNameVectorToAbstractChemicalVector(std::vector<std::string> nameVector)
{
    std::vector<AbstractChemical*> chemicalVector = std::vector<AbstractChemical*>();
    for(unsigned i=0; i<nameVector.size(); i++)
    {
        chemicalVector.push_back(new AbstractChemical(nameVector[i]));
    }

    return chemicalVector;
}

void MembraneReactionSystemFomFile::ParseSystemChemistry(std::vector<std::string> bulk_names, std::vector<std::string> cell_names)
{
    // take inthe names of the species in both the bulk and cell and form the necessary chemistries
    AbstractChemistry* bulkChemistry = new AbstractChemistry();
    AbstractChemistry* cellChemistry = new AbstractChemistry();

    for(unsigned i =0; i<bulk_names.size(); i++)
    {
        AbstractChemical* candidate_chemical = new AbstractChemical(bulk_names[i]);
        bulkChemistry -> AddChemical(candidate_chemical);
    }

    for(unsigned i =0; i<cell_names.size(); i++)
    {
        AbstractChemical* candidate_chemical = new AbstractChemical(cell_names[i]);
        cellChemistry -> AddChemical(candidate_chemical);
    }

    SetBulkChemistry(bulkChemistry);
    SetCellChemistry(cellChemistry);
}

std::vector<std::tuple<std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string>> MembraneReactionSystemFomFile::ReactionSystemFromFile()
{
    // read a reaction from the file
    std::vector<std::tuple<std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string>>  system;  

    std::ifstream inputFile(mInputFileName);

    std::string line;

    std::vector<std::string> reactionString;

    unsigned numberOfReactions=0;

    if(inputFile.is_open())
    {
        // open the reaction file
        while (getline(inputFile,line))
        {
            // for each non-empty reation file line, parse the reactions into
            // data structures on line by line basis
            if(!line.empty())
            {
                if(line.at(0)=='#')
                {
                    // line starts with escape character

                    //std::cout<<"Escape line: "<<line<<std::endl;
                }
                else
                {
                    // check if two sub-reactions are defined (either side of the membrane)
                    size_t coupled_reactions_delimiter_loc = line.find(mReactionDelimiter);
                    if(coupled_reactions_delimiter_loc != std::string::npos)
                    {
                        mIsCoupledReaction.push_back(true);
                        // coupled reactions are defined
                        // assume the next line is a reacton line
                        std::string mStringDelimiter = "";

                        // point at which the reaction type text and the reaction string deliniate
                        size_t separate_point_type_reaction = line.find(mTypeDelimiter);
                        // point at which the reaction string and the reaction information deliniate
                        size_t separate_point_reaction_info = line.find(mDataDelimiter);
                        
                        // the reaction string to be parsed
                        std::string reactionInfo = line.substr(separate_point_reaction_info+2,std::string::npos);

                        // determine the reaction type to add to reaction system
                        std::string reactionType = line.substr(0,separate_point_type_reaction);

                        line.erase(separate_point_reaction_info, std::string::npos);

                        // first reaction substring
                        std::string first_reaction = line.substr(separate_point_type_reaction+3,coupled_reactions_delimiter_loc);
                        

                        // first reaction substring
                        std::string second_reaction = line.substr(coupled_reactions_delimiter_loc+1,std::string::npos);

                        // don't want to parse info in this function/class
                        
                        // test for reversibility in reaction string
                        bool IsReversible = TestReversibility(first_reaction);
                        // if the first reaction is reverisble then the second must also be reverisble as coupled reactions
        
                        // species, stoich
                        std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>> FirstReactionStringTuple = ParseReactionString(first_reaction);
                        std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>> SecondReactionStringTuple = ParseReactionString(second_reaction);

                        mNumberSubstrates.push_back({(unsigned int)std::get<0>(FirstReactionStringTuple)[0].size(),(unsigned int)std::get<0>(SecondReactionStringTuple)[0].size()}); 
                        mNumberProducts.push_back({(unsigned int)std::get<0>(FirstReactionStringTuple)[1].size(),(unsigned int)std::get<0>(SecondReactionStringTuple)[1].size()}); 

                        // concatenate the reaciton tuples; add the dimensions of each separate reaction tuples ot the reactionInfo for parsing within the reaction type class
                        std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>> ReactionStringTuple  = AddCoupledReactions(FirstReactionStringTuple, SecondReactionStringTuple);

                        // reaction type, IsReversible, Substrates, Products, stoichSubstrates, stoichProducts, reaction information
                        // std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string 
                        system.push_back(std::make_tuple(reactionType, IsReversible, std::get<0>(ReactionStringTuple)[0], std::get<0>(ReactionStringTuple)[1], std::get<1>(ReactionStringTuple)[0], std::get<1>(ReactionStringTuple)[1], reactionInfo));
                        
                        
                        
                        numberOfReactions +=1;
                    }else{

                        // assume the next line is a reaction line, single reaction implies reacts woith membrane
                        mIsCoupledReaction.push_back(false);
                        // point at which the reaction type text and the reaction string deliniate
                        size_t separate_point_type_reaction = line.find(mTypeDelimiter);
                        // point at which the reaction string and the reaction information deliniate
                        size_t separate_point_reaction_info = line.find(mDataDelimiter);
                        
                        // the reaction string to be parsed
                        std::string reactionInfo = line.substr(separate_point_reaction_info+2,std::string::npos);

                        // determine the reaction type to add to reaction system
                        std::string reactionType = line.substr(0,separate_point_type_reaction);

                        line.erase(separate_point_reaction_info, std::string::npos);
                        std::string react = line.substr(separate_point_type_reaction+3,std::string::npos);
                        
                        // don't want to parse info in this function
                        
                        
                        // test for reversibility in reaction string
                        bool IsReversible = TestReversibility(react);
        
                        // species, stoich
                        std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>> ReactionStringTuple = ParseReactionString(react);
                        mNumberSubstrates.push_back({(unsigned int)std::get<0>(ReactionStringTuple)[0].size()}); 
                        mNumberProducts.push_back({(unsigned int)std::get<0>(ReactionStringTuple)[1].size()}); 

                        // reaction type, IsReversible, Substrates, Products, stoichSubstrates, stoichProducts, reaction information
                        // std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string 
                        system.push_back(std::make_tuple(reactionType, IsReversible, std::get<0>(ReactionStringTuple)[0], std::get<0>(ReactionStringTuple)[1], std::get<1>(ReactionStringTuple)[0], std::get<1>(ReactionStringTuple)[1], reactionInfo));
                        numberOfReactions +=1;

                    }
                }
            }
        }
        inputFile.close();
    }
    else
    {
        std::cout<<"Error filename not found: "<<mInputFileName<<std::endl;
    }
    

    return system;
}

bool MembraneReactionSystemFomFile::TestReversibility(std::string line)
{
    // test for reversibility in reaction string
    bool IsReversible =  false;
    if(line.find(mReverDelimiter) != std::string::npos){
        // set reversible switch
        IsReversible =  true;
    }

    return IsReversible;
}

std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>> MembraneReactionSystemFomFile::ParseReactionString(std::string line)
{
    
    std::string tempString;
    std::vector<std::vector<std::string>> species;
    std::vector<std::vector<unsigned>> stoich;
    
    // parse into complexes, two complexes per reaction (head and tail of reaction arrows)
    std::vector<std::string> complexes;

    bool IsReversible = TestReversibility(line);
    std::string delim;
    if(IsReversible)
    {
        delim = mReverDelimiter;
    }
    else
    {
        // try to find each of the irreversible delimiters ("->","<-")
        // if present then position != npos
        if (line.find(mForwardDelimiter) != std::string::npos) {
            delim = mForwardDelimiter;
        }
        else
        {
            // how does the reverse work?
            delim = mReverseDelimiter;
        }
    }
    

    // split the reaction string on the deliminator, whether reversible or irreversible
    size_t posC = line.find(delim);


    complexes.push_back(line.substr(0,posC-1));
    complexes.push_back(line.substr(posC+3,std::string::npos));

    // parse complexes based on additon of new molecules
    for(unsigned complexNumber=0; complexNumber<2; complexNumber++)
    {
        std::string str=complexes[complexNumber];
        // posiitons of characters in the reaction string to parse
        size_t posSnew = 0;
        size_t posSold = 0;
        
        std::vector<std::string> reactants;
        std::vector<unsigned> stoichVector;

        // remove potential whitespace from zeroth character, seen in case <->
        if(isspace(str.c_str()[0])){
            str.erase(str.begin());
        }

        while(posSnew != std::string::npos){

            // update position of character pointer based on species separator " + " default
            posSnew = str.find(mSpeciesSeparator,posSold);
            
            std::string strT = str.substr(posSold,posSnew-posSold);
        
            posSold=posSnew+3;
            // determine the stoich ratio value of the species

            unsigned stoichValue=1;

            if(isdigit(strT.c_str()[0])){
                stoichValue=std::stoul(strT.c_str());
                
                unsigned i=0;
                while(isdigit(strT.c_str()[i])){i++;}
                    tempString=strT.substr(i,std::string::npos);
            }else{tempString=strT;}

            // remove whitespace and memory container size to ensure like for like comparisons are consistent
            tempString.erase(remove_if(tempString.begin(), tempString.end(), isspace), tempString.end());

            if(stoichValue==0){
                // skip the push_back of the zero entry
            }
            else
            {
                reactants.push_back(tempString);
                stoichVector.push_back(stoichValue);
            }
            
        }
        species.push_back(reactants);
        stoich.push_back(stoichVector);

    }

    return std::make_tuple(species, stoich);

    // return reaction type, stoich, species, IsReversible, reaction info (info is dependent on reaction type)
}

std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>> MembraneReactionSystemFomFile::AddCoupledReactions(std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>> FirstReactionStringTuple, std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>> SecondReactionStringTuple)
{
    // concatenate the two reaction informations together, bulk then cell

    std::vector<std::vector<std::string>> species = std::get<0>(FirstReactionStringTuple); //[0][] substrates, [1][] products
    std::vector<std::vector<unsigned>> stoich = std::get<1>(FirstReactionStringTuple); //[0][] substrates, [1][] products
    std::vector<std::vector<std::string>> speciesCell = std::get<0>(SecondReactionStringTuple); //[0][] substrates, [1][] products
    std::vector<std::vector<unsigned>> stoichCell = std::get<1>(SecondReactionStringTuple); //[0][] substrates, [1][] products

    species[0].insert(species[0].end(), speciesCell[0].begin(), speciesCell[0].end()); // substrates
    species[1].insert(species[1].end(), speciesCell[1].begin(), speciesCell[1].end()); // products

    stoich[0].insert(stoich[0].end(), stoichCell[0].begin(), stoichCell[0].end()); // substrates
    stoich[1].insert(stoich[1].end(), stoichCell[1].begin(), stoichCell[1].end()); // products

    return std::make_tuple(species, stoich);
}


std::string MembraneReactionSystemFomFile::ParseReactionInformation(std::string line)
{
    // at the moment leave empty, populated in the individual reaction classes
    return line;
}

std::string MembraneReactionSystemFomFile::GetFileName()
{
    return mInputFileName;
}

std::vector<std::vector<unsigned>> MembraneReactionSystemFomFile::GetNumberOfSubstrates()
{
    return mNumberSubstrates;
}

std::vector<std::vector<unsigned>> MembraneReactionSystemFomFile::GetNumberOfProducts()
{
    return mNumberProducts;
}

void MembraneReactionSystemFomFile::SetNumberOfSubstrates(std::vector<std::vector<unsigned>> numberSubstrates)
{
    mNumberSubstrates = numberSubstrates;
}

void MembraneReactionSystemFomFile::SetNumberOfProducts(std::vector<std::vector<unsigned>> numberProducts)
{
    mNumberProducts = numberProducts;
}

std::vector<unsigned> MembraneReactionSystemFomFile::ReturnSubstratesForReactionIndex(unsigned index)
{
    if(index<mNumberSubstrates.size())
    {
        return mNumberSubstrates[index];
    }
    return std::vector<unsigned>();
}

std::vector<unsigned> MembraneReactionSystemFomFile::ReturnProductsForReactionIndex(unsigned index)
{
    if(index<mNumberProducts.size())
    {
        return mNumberProducts[index];
    }
    return std::vector<unsigned>();
}


#endif
