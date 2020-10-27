#ifndef ABSTRACTREACTIONSYSTEMFROMFILE_HPP_
#define ABSTRACTREACTIONSYSTEMFROMFILE_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <iostream>

// reaction includes
#include "AbstractReactionSystem.hpp"
#include "AbstractReaction.hpp"
#include "AbstractChemistry.hpp"
#include "ReactionTypeDatabase.hpp"

// class containg the emthods to read a system of reactions from a text file. A set of file 
// deliminators for p[asring the differnet information string are provided defaults. For each reaction
// the chemical species participating and their stoichiometry as substrates and/or products is 
// determined for the fiel row string. The reaction type is read form the file and the 
// ReactionTypeDatabase with the ReactionTablet function is used to typecast the desired reaction 
// type which inherit from AbstractReaction. The full system chemistry is determined.

class AbstractReactionSystemFromFile : public AbstractReactionSystem
{
protected:

    //inherit:
        // mpSystemChemistry;
        // mNumberOfReactions;
        // std::vector<AbstractReaction*> mpReactionVector;
    std::string mInputFileName;

    std::string mDomain;

    std::string mStringDelimiter = "";
    std::string mDataDelimiter = ";";
    std::string mIrreverDelimiter = "->";
    std::string mReverDelimiter = "<->";
    std::string mSpeciesSeparator = " + ";
    std::string mTypeDelimiter = " : ";


public:

    AbstractReactionSystemFromFile(std::string);

    virtual ~AbstractReactionSystemFromFile()
    {
    };

    virtual void SetFileDeliminator();

    void FormReactionSystemObjectFromTuple(std::vector<std::tuple<std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string>>);

    std::vector<AbstractChemical*> FromChemicalNameVectorToAbstractChemicalVector(std::vector<std::string>);

    void ParseSystemChemistry(std::vector<std::string>);

    std::vector<std::tuple<std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string>> ReactionSystemFromFile();

    bool TestReversibility(std::string);

    std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>> ParseReactionString(std::string);

    std::string ParseReactionInformation(std::string);

    std::string GetFileName();

    std::string GetDomain();

    void SetDomain(std::string);

    AbstractReactionSystem* GetReactionSystem();

};

AbstractReactionSystemFromFile::AbstractReactionSystemFromFile(std::string InputFileName)
    : mInputFileName(InputFileName)
{

    std::vector<AbstractReaction*> mpReactionVector = std::vector<AbstractReaction*>();

    SetFileDeliminator();

    std::vector<std::tuple<std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string>> system = ReactionSystemFromFile();
    std::vector<std::string> chemicalsNames;
    

    for(unsigned reaction=0; reaction<system.size(); reaction++)
    {

        for(unsigned species=0; species<std::get<2>(system[reaction]).size(); species++)
        {
            chemicalsNames.push_back(std::get<2>(system[reaction])[species]);
        }
        for(unsigned species=0; species<std::get<3>(system[reaction]).size(); species++)
        {
            chemicalsNames.push_back(std::get<3>(system[reaction])[species]);
        }
    }

    ParseSystemChemistry(chemicalsNames);

    FormReactionSystemObjectFromTuple(system);

}

void AbstractReactionSystemFromFile::SetFileDeliminator()
{
    mStringDelimiter = "";
    mDataDelimiter = ";";
    mIrreverDelimiter = "->";
    mReverDelimiter = "<->";
    mSpeciesSeparator = " + ";
    mTypeDelimiter = " : ";
}

void AbstractReactionSystemFromFile::FormReactionSystemObjectFromTuple(std::vector<std::tuple<std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string>> system_tuple)
{
    // for each reaction whose data is in the tuple, form the corresponding reaction class
    // denote in the reaction class the function necessary to parse reaction information

    SetNumberOfReactions(system_tuple.size());

    for( unsigned reaction =0; reaction<mNumberOfReactions; reaction++)
    {
        AbstractReaction* p_reaction = new AbstractReaction();

        std::vector<AbstractChemical*> substrates_chemical_vector = FromChemicalNameVectorToAbstractChemicalVector(std::get<2>(system_tuple[reaction]));
        std::vector<AbstractChemical*> products_chemical_vector = FromChemicalNameVectorToAbstractChemicalVector(std::get<3>(system_tuple[reaction]));

        ReactionTablet(p_reaction,std::get<0>(system_tuple[reaction]),substrates_chemical_vector,products_chemical_vector,std::get<4>(system_tuple[reaction]),std::get<5>(system_tuple[reaction]),std::get<6>(system_tuple[reaction]),std::get<1>(system_tuple[reaction]),mpSystemChemistry);
        
        mpReactionVector.push_back(p_reaction);
    }

}

std::vector<AbstractChemical*> AbstractReactionSystemFromFile::FromChemicalNameVectorToAbstractChemicalVector(std::vector<std::string> nameVector)
{
    std::vector<AbstractChemical*> chemicalVector = std::vector<AbstractChemical*>();
    for(unsigned i=0; i<nameVector.size(); i++)
    {
        chemicalVector.push_back(new AbstractChemical(nameVector[i]));
    }

    return chemicalVector;
}

void AbstractReactionSystemFromFile::ParseSystemChemistry(std::vector<std::string> species_names)
{
    AbstractChemistry* mpSystemChemistry = new AbstractChemistry();
    for(unsigned i =0; i<species_names.size(); i++)
    {
        AbstractChemical* candidate_chemical = new AbstractChemical(species_names[i]);
        mpSystemChemistry -> AddChemical(candidate_chemical);
    }

    SetSystemChemistry(mpSystemChemistry);
}

std::vector<std::tuple<std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string>> AbstractReactionSystemFromFile::ReactionSystemFromFile()
{
    // read a reaction from the file
    std::vector<std::tuple<std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string>>  system;  

    std::ifstream inputFile(mInputFileName);

    std::string line;

    std::vector<std::string> reactionString;
    
    bool DomainFound = false;

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
                if(!DomainFound)
                {
                    if(line.find("Domain : ") != std::string::npos)
                    {
                        SetDomain(line.substr(line.find("Domain : ")+1,std::string::npos));
                        DomainFound = true;
                        
                    }
                }else{
                
                    // assume the enxt line is a reacton line
                    // point at which the reaction type text and the reaction string deliniate
                    size_t separate_point_type_reaction = line.find(mTypeDelimiter);
                    // point at which the reaction string and the reaction information deliniate
                    size_t separate_point_reaction_info = line.find(mDataDelimiter);
                    // the reaction string to be parsed
                    std::string reactionInfo = line.substr(separate_point_reaction_info+1,std::string::npos);
        
                    // determine the reaction type to add to reaction system
                    std::string reactionType = line.substr(0,separate_point_type_reaction);
                    line.erase(separate_point_reaction_info, std::string::npos);
                    std::string react = line.substr(separate_point_type_reaction+3,std::string::npos);
                    
                    // don't want to parse info in this function
                    // test for reversibility in reaction string
                    bool IsReversible = TestReversibility(react);
                    // species, stoich
                    std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>> ReactionStringTuple = ParseReactionString(react);
                    // reaction type, IsReversible, Substrates, Products, stoichSubstrates, stoichProducts, reaction information
                    // std::string, bool, std::vector<std::string>, std::vector<std::string>, std::vector<unsigned>, std::vector<unsigned>, std::string 
                    system.push_back(std::make_tuple(reactionType, IsReversible, std::get<0>(ReactionStringTuple)[0], std::get<0>(ReactionStringTuple)[1], std::get<1>(ReactionStringTuple)[0], std::get<1>(ReactionStringTuple)[1], reactionInfo));
                    numberOfReactions +=1;
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

bool AbstractReactionSystemFromFile::TestReversibility(std::string line)
{
    // test for reversibility in reaction string
    bool IsReversible =  false;
    if(line.find(mReverDelimiter) != std::string::npos){
        // set reversible switch
        IsReversible =  true;
    }

    return IsReversible;
}

std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<unsigned>>> AbstractReactionSystemFromFile::ParseReactionString(std::string line)
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
        delim = mIrreverDelimiter;
    }
    

    // split the reaction string on the deliminator, whether reversible or irreversible
    size_t posC = line.find(delim);


    complexes.push_back(line.substr(0,posC-1));
    complexes.push_back(line.substr(posC+3,std::string::npos));

    // parse complexes based on additon of new molecuels
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

std::string AbstractReactionSystemFromFile::ParseReactionInformation(std::string line)
{
    // at the moment leave empty
    return line;
}

std::string AbstractReactionSystemFromFile::GetFileName()
{
    return mInputFileName;
}

std::string AbstractReactionSystemFromFile::GetDomain()
{
    return mDomain;
}

void AbstractReactionSystemFromFile::SetDomain(std::string domain)
{
    mDomain = domain;
}

#endif