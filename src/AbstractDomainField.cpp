#include "AbstractDomainField.hpp"
#include "RandomNumberGenerator.hpp"

AbstractDomainField::AbstractDomainField(   std::string domainLabelFilename, 
                                            std::string domainKeyFilename, 
                                            std::string odeLabelFilename, 
                                            std::string odeKeyFilename, 
                                            std::string diffusionFilename)
    :   mDomainLabelFilename(domainLabelFilename),
        mDomainKeyFilename(domainKeyFilename),
        mOdeLabelFilename(odeLabelFilename),
        mOdeKeyFilename(odeKeyFilename),
        mDiffusionFilename(diffusionFilename)
        
        
{
    // read the domain labels and label keys form files
    SetupAndInitialiseLabelDomain();

    // scale the input file dimensions to the chaste rectilinear cartesian grid
    MapToChasteCartesianDomain();

    // translate the ode labels to ode nodes, form the ode vector

    // hard code the honey comb mesh scale factor and dimensions
    std::vector<double> meshScaleXY;
    meshScaleXY.push_back(1.0);
    meshScaleXY.push_back(0.866025);

    std::vector<unsigned> meshDimensions;
    double meshElementX = 0.0;
    double meshElementY = 0.0;
    
    // determine the mesh dimensions up to the maximum allowed through the scaling
    // cover the whole domain with the mesh
    unsigned count=0;
    while( meshElementX < mCartesianCellDimensions[0]*mCartesianCellScaleXY[0])
    {
        meshElementX  = count*meshScaleXY[0];
        count=count+1;
    }
    meshDimensions.push_back(count-1);

    count=0;
    while( meshElementY < mCartesianCellDimensions[1]*mCartesianCellScaleXY[1])
    {
        meshElementY  = count*meshScaleXY[1];
        count=count+1;
    }
    meshDimensions.push_back(count-1);

    // use the aforementioned mesh dimension to produce a new honeycomb mesh, as mutable mesh type

    //HoneycombMeshGenerator generator(meshDimensions[0], meshDimensions[1], 0);
    HoneycombMeshGenerator* p_generator = new HoneycombMeshGenerator(meshDimensions[0], meshDimensions[1], 0);
    //MutableMesh<2,2>* p_mesh = generator.GetMesh();
    MutableMesh<2,2>* p_mesh = p_generator -> GetMesh();
    // store the mesh properties as the domain field
    SetMeshDimensions(meshDimensions);
    SetMeshScale(meshScaleXY);
    SetDomainMesh(p_mesh);
    SetMeshGenerator(p_generator);

    // read in the ode information, set the number of ode systems and 
    SetupAndInitialiseOdes();
    MapOdeLabelsToCartesianOdeDomain();

    // read in the domain information, set the number of domains
    SetupAndInitialiseDiffusionDatabase();

    // store the node labels as std::string type
    SetNodeLabels(LabelNodesWithOdes());
    SetNodeDomains(LabelNodesWithDomains());
}

void AbstractDomainField::SetupAndInitialiseLabelDomain()
{
    // read in and store the domain grid labels for mapping to mesh nodes and chaste rectilinear grid
    ReadDomainLabels();
    ReadDomainKey();
    
    std::vector<unsigned> cartesianCellDimensions;
    std::vector<double> cartesianCellScaleXY;

    // assume the domain scaling values are unity, may change with different input grids to increase fidelity
    cartesianCellDimensions.push_back(mDomainLabels[0].size());
    cartesianCellDimensions.push_back(mDomainLabels.size());
    cartesianCellScaleXY.push_back(1.0);
    cartesianCellScaleXY.push_back(1.0);

    // store these values
    SetCartesianCellDimensions(cartesianCellDimensions);
    SetCartesianCellScale(cartesianCellScaleXY);

}

void AbstractDomainField::SetupAndInitialiseOdes()
{
    // read in and store the ode labels from the input grid, take in the ode label keys 
    ReadOdeLabels();
    ReadOdeKey();

    std::vector<unsigned> cartesianOdeDimensions;
    std::vector<double> cartesianOdeScaleXY;

    // assume the domain scaling values are unity, may change with different input grids to increase fidelity
    cartesianOdeDimensions.push_back(mOdeLabels[0].size());
    cartesianOdeDimensions.push_back(mOdeLabels.size());
    cartesianOdeScaleXY.push_back(1.0);
    cartesianOdeScaleXY.push_back(1.0);

    // store these values
    SetCartesianOdeDimensions(cartesianOdeDimensions);
    SetCartesianOdeScale(cartesianOdeScaleXY);

}

void AbstractDomainField::SetupAndInitialiseDiffusionDatabase()
{
    // to be overriden for further meshing and diffusion schemes
    ReadDiffusionDatabase();
}

std::vector<std::string> AbstractDomainField::LabelNodesWithDomains()
{
    // retireve the node positions and compare to domain labels in order to label the nodes
    // nodes are defined bottom-left of domain to top-right in a serialised manner
    MutableMesh<2,2>* p_mesh = GetDomainMesh();

    unsigned numberOfNodes = p_mesh -> GetNumNodes();

    // form the serialised node label array
    std::vector<std::string> serialisedNodeDomains(numberOfNodes,"");
    
    for (AbstractTetrahedralMesh<2,2>::NodeIterator iter = p_mesh ->GetNodeIteratorBegin();
             iter != p_mesh ->GetNodeIteratorEnd();
             ++iter)
    {
        // for each node, retireve index and position, label with correct ode label
        unsigned node_index = iter ->GetIndex();
        // retireve the position of the node of index node_index
        const c_vector<double,2>& position = iter->rGetLocation();
        // search through the position ode labels for the correct label to give the node
        serialisedNodeDomains.at(node_index) = ReturnDomainLabelAtPosition(position);
    }

    return serialisedNodeDomains;
}

std::vector<std::string> AbstractDomainField::LabelNodesWithOdes()
{
    // retireve the node positions and compare to domain labels in order to label the nodes
    // nodes are defined bottom-left of domain to top-right in a serialised manner
    MutableMesh<2,2>* p_mesh = GetDomainMesh();

    unsigned numberOfNodes = p_mesh -> GetNumNodes();

    // form the serialised node label array
    std::vector<std::string> serialisedNodeLabels(numberOfNodes,"");
    
    for (AbstractTetrahedralMesh<2,2>::NodeIterator iter = p_mesh ->GetNodeIteratorBegin();
             iter != p_mesh ->GetNodeIteratorEnd();
             ++iter)
    {
        // for each node, retireve index and position, label with correct ode label
        unsigned node_index = iter ->GetIndex();
        // retireve the position of the node of index node_index
        const c_vector<double,2>& position = iter->rGetLocation();
        // search through the position ode labels for the correct label to give the node
        serialisedNodeLabels.at(node_index) = ReturnNodeOdeLabelAtPosition(position);
    }

    return serialisedNodeLabels;
}

void AbstractDomainField::ParseInitialConditionsFromFile(std::string initialConditionsFilename, bool IsPerturbInitialConditions)
{
    // read the input intial conditions file as a matrix of strings, first column state name, 
    // 2nd column subdomain, 3rd column value.

    // if the state is not within the state register then the state doesn't evolve over time, remaining constant
    std::vector<std::vector<std::string>> fullInitialConditionsAsStrings = ReadMatrix(initialConditionsFilename);
    std::vector<std::vector<std::string>> initialConditionsAsStrings = ReturnStateDataWithinRegister(fullInitialConditionsAsStrings);

    // use state name against stateVariableRegister to remove any state value not in the system
    // replace any state in register but not in intitial conditions with value of 0.0
    
    // for the vector in order of StateVariableRegister
    unsigned numberOfStateVariables = mStateVariableVector->GetNumberOfStateVariables();
    unsigned numberOfNodes = GetDomainMesh() -> GetNumNodes();
    
 
    // determine the nodal initial condition, test for perturbation
    std::vector<double> init_conds(numberOfStateVariables*numberOfNodes,0.0);
    std::vector<std::string> node_domain_labels = GetNodeDomains();
    for (unsigned node_index=0; node_index<numberOfNodes; node_index++)
    {   
        // each node set as being a vector of the read in state condition values

        // test for the domain of the node, retireve the condition for the specific domain
        for(unsigned pdeDim=0; pdeDim<numberOfStateVariables; pdeDim++)
        {   
            // determine which state the pde dimension refers to and then retrieve the condition if present
            std::string stateName = mStateVariableVector -> RetrieveStateVariableName(pdeDim);
       
            bool IsFoundState =false;
            // run through input data, line by line
            for(unsigned inputState=0; inputState<initialConditionsAsStrings.size();inputState++)
            {
                // test whether the state variable is specified within the input data
                if(stateName==initialConditionsAsStrings[inputState][0])
                {
                    // the state is in the input date but now test for the correct domain
                    if(initialConditionsAsStrings[inputState].size()>2)
                    {
                        
                        // then sub domain is also specified, test for sub domain
                        if(node_domain_labels[node_index]==initialConditionsAsStrings[inputState][1]||ReturnDomainKeyFromDomainLabel(node_domain_labels[node_index])==initialConditionsAsStrings[inputState][1])
                        {
                            // sub domain has been found, state found, fill data record
                            IsFoundState = true;
                           
                            // store the condition, serialised for nodes

                            if(PerturbInitialConditionTest(initialConditionsAsStrings[inputState])||IsPerturbInitialConditions)
                            {
                            
                                init_conds[numberOfStateVariables*node_index + pdeDim] =fabs(std::stod(initialConditionsAsStrings[inputState][2]) + RandomNumberGenerator::Instance()->ranf());
                            }
                            else
                            {
                                init_conds[numberOfStateVariables*node_index + pdeDim] =std::stod(initialConditionsAsStrings[inputState][2]);
                            }
                            break;
                        }
                        // subdomain specified, but not found.  State not present in particular subdomain
                    }
                    else
                    {
                        //state data found, but no sub domain is specied, assume present for all sub domains if sub domains are indeed present
                        IsFoundState = true;
                        // store the condition, serialised for nodes
                        if(IsPerturbInitialConditions)
                        {
                            init_conds[numberOfStateVariables*node_index + pdeDim] =fabs(std::stod(initialConditionsAsStrings[inputState][2]) + RandomNumberGenerator::Instance()->ranf());
                        }
                        else
                        {
                            init_conds[numberOfStateVariables*node_index + pdeDim] =std::stod(initialConditionsAsStrings[inputState][2]);
                        }
                        break;
                    }   
                }
                // data record isn't correct for the pde state name, if at end of data record state not found then state not present in input data
            }  
            if(!IsFoundState)
            {
                // state in system but not in the input data, default to 0.0
                // initialse the missing condition, serialised for nodes    
                if(IsPerturbInitialConditions)
                {
                    // initialise to a small perturbation, otherwise leave as intitialised 0.0
                    init_conds[numberOfStateVariables*node_index + pdeDim] =fabs(RandomNumberGenerator::Instance()->ranf());
                }
            }
        }
        // each node pde is set
    }

    SetInitialNodeConditions(init_conds);
}

bool AbstractDomainField::PerturbInitialConditionTest(std::vector<std::string> inputVector)
{
    // determine whether or not the initial conditions is to be perturbed
    // perturbation on more specific basis
        
    if(inputVector.size()>3)
    {
        // then have the perturbation bool selections
        if(inputVector[3]=="true"||inputVector[3]=="True"||inputVector[3]=="TRUE")
        {
            return true;
            //IsPerturbSpecificConditions[mStateVariableVector -> RetrieveStateVariableIndex(initialConditionsAsStrings[i][0])] = true;
        }
        else
        {
            return false;
        }
        
    }
    else
    {
        return false;
    }
}



std::vector<std::vector<std::string>> AbstractDomainField::ReturnStateDataWithinRegister(std::vector<std::vector<std::string>> inputMatrix)
{
    // take an input matrix read from a file wherein the first column denotes a state
    // check against the StateVariableRegister for states that are in the system and remove input
    // that are not in the system
    std::string inputName;
    std::vector<std::vector<std::string>> newMatrix;

    for(unsigned i=0; i<inputMatrix.size();i++)
    {
        inputName = inputMatrix[i][0];

        if(mStateVariableVector -> IsStateVariablePresent(inputName))
        {
            newMatrix.push_back(inputMatrix[i]);
        }
    }
    // newMatrix conatins only the data records of states in the system
    return newMatrix;
}

void AbstractDomainField::ParseBoundaryConditionsFromFile(std::string boundaryConditionsFilename)
{
    // read the input boundary conditions file as a matrix of strings, first column state name, 
    // 2nd column boundary condition type, 3rd column value.
    std::vector<std::vector<std::string>> fullBoundaryConditionsAsStrings = ReadMatrix(boundaryConditionsFilename);
    std::vector<std::vector<std::string>> boundaryConditionsAsStrings = ReturnStateDataWithinRegister(fullBoundaryConditionsAsStrings);
    // use state name against stateVariableRegister to remove any state value not in the system
    // replace any state in register but not in boundary conditions with value of 0.0, type Derichlet

    unsigned numberOfStateVariables = GetStateVariableVector() ->GetNumberOfStateVariables();

    std::vector<std::string> boundaryConditionTypes(numberOfStateVariables,"Dirichlet");
    std::vector<double> boundaryConditionValues(numberOfStateVariables,0.0);

    for(unsigned pdeDim=0; pdeDim<numberOfStateVariables; pdeDim++)
    {   
        // determine which state the pde dimension refers to and then retrieve the condition if present
        std::string stateName = mStateVariableVector -> RetrieveStateVariableName(pdeDim);

        for(unsigned inputState=0; inputState<boundaryConditionsAsStrings.size();inputState++)
        {
            // test whether the state variable is specified within the input data
            if(stateName==boundaryConditionsAsStrings[inputState][0])
            {
                // state found, read data
                if(boundaryConditionsAsStrings[inputState].size()>2)
                {
                    // type value is specified also
                    boundaryConditionTypes[pdeDim] = boundaryConditionsAsStrings[inputState][1];
                    boundaryConditionValues[pdeDim] = std::stod(boundaryConditionsAsStrings[inputState][2]);
                }
                else
                {
                    // type value not specified, assume Dirichlet
                    boundaryConditionValues[pdeDim] = std::stod(boundaryConditionsAsStrings[inputState][1]);
                }
                break;
            }
        }
        // if state isn't found, leave as default, 0.0 Dirichlet
    }

    SetBoundaryConditionTypes(boundaryConditionTypes);

    SetBoundaryConditionValues(boundaryConditionValues);
}


void AbstractDomainField::PrintDiffusionDomain()
{
    std::cout<<"Diffusion domain:"<<std::endl;
    printMatrix(mDomainLabels);
}

void AbstractDomainField::PrintMappedDiffusionDomain()
{
    std::cout<<"Mapped diffusion domain:"<<std::endl;
    printMatrix(mCartesianChasteDomain);
}

void AbstractDomainField::PrintDomainLabelKeys()
{
    std::cout<<"Domain label keys:"<<std::endl;
    printMatrix(mDomainKeys);
}

void AbstractDomainField::PrintODEDomain()
{
    std::cout<<"ODE domain:"<<std::endl;
    printMatrix(mOdeLabels);
}

void AbstractDomainField::PrintMappedODEDomain()
{
    std::cout<<"Mapped ODE domain:"<<std::endl;
    printMatrix(mCartesianOdeDomain);
}

void AbstractDomainField::PrintODELabelKeys()
{
    std::cout<<"ODE label keys:"<<std::endl;
    printMatrix(mOdeKeys);
}

void AbstractDomainField::PrintDiffusionDatabase()
{
    std::cout<<"Diffusion database:"<<std::endl;
    printMatrix(mDiffusionDatabase);
}




std::string AbstractDomainField::ReturnNodeOdeLabelAtPosition(const c_vector<double,2>& position)
{
    // search through the node labels to provide the specific label at the given position

    // scale the dimension to provide the X,Y indices of the grid element containing the node
    unsigned labelXIndex =0;
    unsigned labelYIndex =0;

    // for the y position

    for(unsigned i=1; i<(mCartesianOdeDimensions[1]+1); i++) // top to bottom
    {
        if(position[1] <= i*mCartesianOdeScaleXY[1])
        {
            labelYIndex = i-1;
            break;
        }
    }

    // for the x position

    for(unsigned i=1; i<(mCartesianOdeDimensions[0]+1); i++)
    {
        if(position[0] <= i*mCartesianOdeScaleXY[0])
        {
            labelXIndex = i-1;
            break;
        }
    }

    // return the ode label for the grid element at the X,Y index containing the node
    return mOdeLabels[labelYIndex][labelXIndex];
}

std::string AbstractDomainField::ReturnDomainLabelAtPosition(const c_vector<double,2>& position)
{
    // search through the domain labels to provide the specific label at the given position

    // scale the dimension to provide the X,Y indices of the grid element containing the position under inspection
    unsigned labelXIndex =0;
    unsigned labelYIndex =0;
    // for the y position

    for(unsigned i=1; i<(mCartesianChasteDimensions[1]+1); i++) // top to bottomw
    {
        if(position[1] <= i*mCartesianChasteScaleXY[1])
        {
            labelYIndex = i-1;
            break;
        }
    }

    // for the x position

    for(unsigned i=1; i<(mCartesianChasteDimensions[0]+1); i++)
    {
        if(position[0] <= i*mCartesianChasteScaleXY[0])
        {
            labelXIndex = i-1;
            break;
        }
    }

    // return the domain label for the grid element at the X,Y index containing the position
    return mCartesianChasteDomain[labelYIndex][labelXIndex];
}

void AbstractDomainField::MapOdeLabelsToCartesianOdeDomain()
{
    // map the labels of the input grid file to the scaled chaste cartesian grid.
    
    std::vector<double> cartesianChasteScaleXY = GetCartesianChasteScale();

    // determine the number of steps in a mapped chaste grid based on the the number of input domain dimensions and the
    // ration of scaling parameters, round to the nearest integer value.
    unsigned numberOfXSteps = std::round((mCartesianOdeDimensions[0]*mCartesianOdeScaleXY[0])/cartesianChasteScaleXY[0]); 
    unsigned numberOfYSteps = std::round((mCartesianOdeDimensions[1]*mCartesianOdeScaleXY[1])/cartesianChasteScaleXY[1]); 

    std::vector<std::vector<std::string>> cartesianOdeDomain;

    // define the boundary with the edge steps included for easier comparisons of (double) spatial values
    for(unsigned yindex=0; yindex<numberOfYSteps; yindex++)
    {
        std::vector<std::string> rowVector;
        for(unsigned xindex=0; xindex<numberOfXSteps; xindex++)
        {
            // push back the ode label of the file domain onto the mapped domain, here everything is in terms of std::string label
            rowVector.push_back(ReturnMappedOdeLabel(xindex,yindex,mCartesianChasteScaleXY));

        }
        cartesianOdeDomain.push_back(rowVector);
    }

    //printMatrix(mOdeLabels);
    //printMatrix(cartesianOdeDomain);

    SetCartesianOdeDomain(cartesianOdeDomain);
}

void AbstractDomainField::MapToChasteCartesianDomain()
{
    
    std::vector<double> cartesianChasteScaleXY;
    cartesianChasteScaleXY.push_back(0.25);
    cartesianChasteScaleXY.push_back(0.144338);

    SetCartesianChasteScale(cartesianChasteScaleXY);

    unsigned numberOfXSteps = std::round((mCartesianCellDimensions[0]*mCartesianCellScaleXY[0])/cartesianChasteScaleXY[0]); 
    unsigned numberOfYSteps = std::round((mCartesianCellDimensions[1]*mCartesianCellScaleXY[1])/cartesianChasteScaleXY[1]); 
 
    std::vector<std::vector<std::string>> cartesianChasteDomain;
    std::vector<unsigned> cartesianChasteDimensions;
    cartesianChasteDimensions.push_back(numberOfXSteps);
    cartesianChasteDimensions.push_back(numberOfYSteps);

    // define the boundary with the edge steps included for easier comparisons of (double) spatial values
    for(unsigned yindex=0; yindex<numberOfYSteps; yindex++)
    {
        std::vector<std::string> rowVector;
        for(unsigned xindex=0; xindex<numberOfXSteps; xindex++)
        {
            // push back the domain label of the file domain onto the mapped domain, here everything is in terms of std::string label
            rowVector.push_back(ReturnMappedDomainLabel(xindex,yindex,mCartesianChasteScaleXY));

        }
        cartesianChasteDomain.push_back(rowVector);
    }


    SetCartesianChasteDomain(cartesianChasteDomain);
    SetCartesianChasteDimensions(cartesianChasteDimensions);


}

double AbstractDomainField::GetDiffusionValueBasedOnPoint(const ChastePoint<2>& chastePoint, unsigned stateIndex)
{
    
    // take in the point and state index then determine the state name and domain label, then determine the diffusion value
    if(stateIndex<mStateVariableVector ->GetNumberOfStateVariables())
    {
        if(mDiffusionDatabase[0].size() ==1 )
        {
            // only a diffusion value vector
            return std::stod(mDiffusionDatabase[stateIndex][0]);
        }
        else
        {
            // retrive state name from stateVariableRegister
            std::string stateName = mStateVariableVector -> RetrieveStateVariableName(stateIndex);

            // retrieve label from domainLabels
            std::string domainLabel = ReturnDomainLabelAtPosition(chastePoint.rGetLocation());

            return ReturnDiffusionValueFromStateNameAndDomainLabel(stateName,domainLabel);
        }
    }
    else
    {
        // state hasn't been found therefore return a diffusion value of 0.0
        std::cout<<"Error: AbstractDomainField::GetDiffusionValueBasedOnPoint: State not in state variable"<<std::endl;
        return 0.0;
    }
}

double AbstractDomainField::ReturnDiffusionValueFromStateNameAndDomainLabel(std::string stateName, std::string domainLabel)
{
    // look for diffusion value (from diffusiondatabase) based on state name and domain label

    unsigned numberOfDiffusiveEntries = mDiffusionDatabase.size();
    bool IsStateFound=false;
    if(mDiffusionDatabase[0].size() ==2)
    {
        // take to be only state name, diffusion value structure
        // find stateName
        for(unsigned diffusive_index=0; diffusive_index<numberOfDiffusiveEntries; diffusive_index++)
        {
            if(stateName == mDiffusionDatabase[diffusive_index][0])
            {
                IsStateFound = true;
                return std::stod(mDiffusionDatabase[diffusive_index][1]);
                break;
                
            }
        }
        // if the state isn't found in the diffusion database then state cannot diffuse, return value of 0.0
        if(!IsStateFound)
        {
            return 0.0;
        }
    }
    else
    {
        // find stateName, test for domain name
        for(unsigned diffusive_index=0; diffusive_index<numberOfDiffusiveEntries; diffusive_index++)
        {
            if(stateName == mDiffusionDatabase[diffusive_index][0])
            {
                // found state, then test for domain
                if(domainLabel == mDiffusionDatabase[diffusive_index][1])
                {
                    // state and domain found, return diffusivity
                    IsStateFound = true;
                    return std::stod(mDiffusionDatabase[diffusive_index][2]);
                    break;
                }
            }
        }
        // if the state isn't found in the diffusion database then state cannot diffuse, return value of 0.0
        if(!IsStateFound)
        {
            return 0.0;
        }
    }
    return 0.0;
}

std::string AbstractDomainField::ReturnMappedDomainLabel(unsigned xindex, unsigned yindex, std::vector<double> scaleFactor)
{
    // compare the xindex and yindex of the new chaste domain to the label domain
    // return the label of the chaste domain, work from the top left to bottom right

    unsigned labelXIndex =0;
    unsigned labelYIndex =0;

    double positionX = xindex*scaleFactor[0];
    double positionY = yindex*scaleFactor[1];
    
    // for the y position

    for(unsigned i=1; i<(mCartesianCellDimensions[1]+1); i++) // top to bottom
    {
        if(positionY <= i*mCartesianCellScaleXY[1])
        {
            labelYIndex = i-1;
            break;
        }
    }

    // for the x position

    for(unsigned i=1; i<(mCartesianCellDimensions[0]+1); i++)
    {
        if(positionX <= i*mCartesianCellScaleXY[0])
        {
            labelXIndex = i-1;
            break;
        }
    }

    return mDomainLabels[labelYIndex][labelXIndex];
}

std::string AbstractDomainField::ReturnMappedOdeLabel(unsigned xindex, unsigned yindex,std::vector<double> scaleFactor)
{
    // compare the xindex and yindex of the new chaste domain to the label ode
    // return the label of the chaste domain, work from the top left to bottom right

    unsigned labelXIndex =0;
    unsigned labelYIndex =0;

    std::string label;
    double positionX = xindex*scaleFactor[0];
    double positionY = yindex*scaleFactor[1];
    
    // for the y position

    for(unsigned i=1; i<(mCartesianOdeDimensions[1]+1); i++) // top to bottomw
    {
        if(positionY <= i*mCartesianOdeScaleXY[1])
        {
            labelYIndex = i-1;
            break;
        }
    }

    // for the x position

    for(unsigned i=1; i<(mCartesianOdeDimensions[0]+1); i++)
    {
        if(positionX <= i*mCartesianOdeScaleXY[0])
        {
            labelXIndex = i-1;
            break;
        }
    }

    return mOdeLabels[labelYIndex][labelXIndex];
}

std::vector<std::vector<std::string>> AbstractDomainField::ReadMatrix(std::string filename)
{
    // parse a matrix file (.csv) line by line, ignore escape line,s containing file information
    // that is lines starting with '#' 
    
    std::string line;
    std::ifstream inputFile(filename);

    // read all data types as std::string therefore return the matrix of strings for personalised
    // methods down stream
    std::vector<std::vector<std::string>> outputMatrix = std::vector<std::vector<std::string>>();

    // check file exists and is openable
    if(inputFile.is_open()){
        // open the matrix file
        while (getline(inputFile,line)){
            // while the file still has lines not read.
            // read line left to right, top to bottom.
            if(!line.empty()){
                if(line.at(0)=='#')
                {
                    //std::cout<<"Escape line: "<<line<<std::endl;
                }
                else
                {
                    outputMatrix.push_back(parseMatrixLineString(line));
                }   
            }
        }
        inputFile.close();

        return outputMatrix;
    }else{
        std::cout<<"Error: Unable to open file: "<<filename<<std::endl;
        return outputMatrix;
    }
}

std::vector<std::string> AbstractDomainField::parseMatrixLineString(std::string line)
{
    // for a line string in the matrix read, parse into vector data entries based on delimiters ','
    std::vector<std::string> rowVector = std::vector<std::string>();

    // delimiter, may be modified by further methods
    std::string delim = ",";
    std::string matrixCell;

    // determine the position of the delimiter
    size_t posSnew=line.find(delim);

    bool IsEndOfLine = false;
    
    while(!IsEndOfLine)
    {
        // while not at the end of the file, sample sub strings from the posiiton of the delimiter
        if(posSnew == std::string::npos)
        {
            IsEndOfLine = true;
        }
        
        // sample substring from begining of the string to the delimiter positioon, store as data entry
        matrixCell = line.substr(0,posSnew);

        // remove the sampled entry from the string
        line = line.substr(posSnew+1,std::string::npos);

        rowVector.push_back(matrixCell);

        // update delimiter position
        posSnew=line.find(delim);
    }
    return rowVector;
}

void AbstractDomainField::ReadDomainLabels()
{
    // read and store the information contained within the domain labels file
    SetDomainLabels(ReadMatrix(mDomainLabelFilename));

    SetDomainLabelVector(ReturnUnique(GetDomainLabels()));
}

void AbstractDomainField::ReadDomainKey()
{
    // read the label keys for the domain, providing more information as to the type of the sun domain
    std::vector<std::vector<std::string>> domainKeys;

    domainKeys = ReadMatrix(mDomainKeyFilename);

    SetDomainKeys(domainKeys);
}

void AbstractDomainField::ReadOdeLabels()
{
    // read and store the ode labels, form a vector of the unique labels to store as a reference vector
    std::vector<std::vector<std::string>> odeLabels;
    
    std::vector<std::string> odeLabelVector;

    odeLabels = ReadMatrix(mOdeLabelFilename);

    //printMatrix(odeLabels);

    odeLabelVector = ReturnUnique(odeLabels);

    SetOdeLabelVector(mOdeLabelVector);
    SetOdeLabels(odeLabels);
}

void AbstractDomainField::ReadOdeKey()
{
    // read the ode keys from a .csv file, providing the label used within the domain and
    // the respective ode to be used
    std::vector<std::vector<std::string>> odeKeys;

    odeKeys = ReadMatrix(mOdeKeyFilename);

    //printMatrix(odeKeys);
    SetOdeKeys(odeKeys);
}

void AbstractDomainField::ReadDiffusionDatabase()
{
    // read and store the diffusion data as a matrix of strings 
    std::vector<std::vector<std::string>> diffusionDatabase;

    diffusionDatabase = ReadMatrix(mDiffusionFilename);

    //printMatrix(diffusionDatabase);

    // vector to store the state variables
    std::vector<std::string> stateVector;
    std::vector<std::string> domainVector;
    for(unsigned i=0; i<diffusionDatabase.size(); i++)
    {
        // the first element of the read in database is the identifying names
        // of the state variables that diffuse
        stateVector.push_back(diffusionDatabase[i][0]);
        if(diffusionDatabase[0].size()>2)
        {
            domainVector.push_back(diffusionDatabase[i][1]);
        }
    }

    SetNumberOfDomains(ReturnUnique(domainVector).size());
    // for a new state variable registers with the unique identifying names from the diffusion file
    StateVariableRegister*   p_stateVariableVector = new StateVariableRegister(ReturnUnique(stateVector));
    unsigned numberOfStates = p_stateVariableVector -> GetNumberOfStateVariables();
    

    // store the diffusion information
    SetStateVariableVector(p_stateVariableVector);

    SetNumberOfStates(numberOfStates);

    SetDiffusionDatabase(diffusionDatabase);
}



std::vector<std::string> AbstractDomainField::ReturnUnique(std::vector<std::string> candidateVector)
{
    // serach through a vector of strings and store the first occurance of each unique string
    std::vector<std::string> resultUnique = std::vector<std::string>();
    resultUnique.push_back(candidateVector.at(0));
    unsigned uniqueCount =1;
    bool IsFound = false;
    for(unsigned i=1; i<candidateVector.size(); i++)
    {
        IsFound = false;
        for(unsigned j=0; j<uniqueCount; j++)
        {

            if(candidateVector[i] == resultUnique[j])
            {
                IsFound =true;
                break;
            }
        }
        if(!IsFound)
        {
            resultUnique.push_back(candidateVector[i]);
            uniqueCount++;
        }
    }
    return resultUnique;
}

std::vector<std::string> AbstractDomainField::ReturnUnique(std::vector<std::vector<std::string>> candidateMatrix)
{
    // serach through a matrix of strings and store the first occurance of each unique string
    std::vector<std::string> resultUnique = std::vector<std::string>();
    resultUnique.push_back(candidateMatrix[0][0]);
    unsigned uniqueCount =1;
    bool IsFound = false;

    for(unsigned j=0; j<candidateMatrix.size(); j++)
    {
        for(unsigned i=0; i<candidateMatrix[j].size(); i++)
        {
            IsFound = false;
            for(unsigned k=0; k<uniqueCount; k++)
            {
                if(candidateMatrix[j][i] == resultUnique[k])
                {
                    IsFound =true;
                    break;
                }
            }
            if(!IsFound)
            {
                resultUnique.push_back(candidateMatrix[j][i]);
                uniqueCount++;
            }
        }
    }
    return resultUnique;
}

void AbstractDomainField::printMatrix(std::vector<std::vector<std::string>> matrix)
{
  for (unsigned int i = 0; i < matrix.size(); i++)
  {
    for (unsigned int j = 0; j < matrix[i].size(); j++)
    {
      std::cout << matrix[i][j]<< ' ';
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;

    return;
};

void AbstractDomainField::printVector(std::vector<std::string> vec)
{
  for (unsigned int i = 0; i < vec.size(); i++)
  {
    std::cout << vec[i]<< ' ';
  }
  std::cout<<std::endl;

    return;
};



void AbstractDomainField::SetNumberOfDomains(unsigned numberOfDomains)
{
    mNumberOfDomains = numberOfDomains;
}

unsigned AbstractDomainField::GetNumberOfDomains()
{
    return mNumberOfDomains;
}

void AbstractDomainField::SetDomainLabels(std::vector<std::vector<std::string>> domainLabels)
{
    mDomainLabels = domainLabels;
}

std::vector<std::vector<std::string>> AbstractDomainField::GetDomainLabels()
{
    return mDomainLabels;
}

void AbstractDomainField::SetDomainLabelVector(std::vector<std::string> domainLabelVector)
{
    mDomainLabelVector = domainLabelVector;
}

std::string AbstractDomainField::GetDomainLabelByIndex(unsigned index)
{
    if(index<mNumberOfDomains)
    {
        return mDomainLabelVector[index];
    }
    else
    {
        std::cout<<"Error: AbstractDomainField::GetDomainLabelByIndex: index out of bounds"<<std::endl;
        return "Error";
    }
    
}

std::vector<std::string> AbstractDomainField::GetDomainLabelVector()
{
    return mDomainLabelVector;
}

void AbstractDomainField::SetDomainKeys(std::vector<std::vector<std::string>> domainKeys)
{
    mDomainKeys = domainKeys;
}

std::vector<std::vector<std::string>> AbstractDomainField::GetDomainKeys()
{
    return mDomainKeys;
}

void AbstractDomainField::SetNumberOfOdeSystems(unsigned numberSystems)
{
    mNumberOfOdeSystems = numberSystems;
}

unsigned AbstractDomainField::GetNumberOfOdeSystems()
{
    return mNumberOfOdeSystems;
}

void AbstractDomainField::SetOdeLabels(std::vector<std::vector<std::string>> odeLabels)
{
    mOdeLabels = odeLabels;
}

std::vector<std::vector<std::string>> AbstractDomainField::GetOdeLabels()
{
    return mOdeLabels;
}

void AbstractDomainField::SetOdeLabelVector(std::vector<std::string> labelVector)
{
    mOdeLabelVector = labelVector;
}

std::vector<std::string> AbstractDomainField::GetOdeLabelVector()
{
    return mOdeLabelVector;
}

void AbstractDomainField::SetOdeKeys(std::vector<std::vector<std::string>> odeKeys)
{
    mOdeKeys = odeKeys;
}

std::vector<std::vector<std::string>> AbstractDomainField::GetOdeKeys()
{
    return mOdeKeys;
}

void AbstractDomainField::SetStateVariableVector(StateVariableRegister* stateVector)
{
    mStateVariableVector = stateVector;
}

StateVariableRegister* AbstractDomainField::GetStateVariableVector()
{
    return mStateVariableVector;
}

void AbstractDomainField::SetNumberOfStates(unsigned numberStates)
{
    mNumberOfStates = numberStates;
} 

unsigned AbstractDomainField::GetNumberOfStates()
{
    return mNumberOfStates;
}

void AbstractDomainField::SetCartesianCellDimensions(std::vector<unsigned> cellDims)
{
    mCartesianCellDimensions = cellDims;
}

std::vector<unsigned> AbstractDomainField::GetCartesianCellDimensions()
{
    return mCartesianCellDimensions;
}

void AbstractDomainField::SetCartesianOdeDimensions(std::vector<unsigned> odeDims)
{
    mCartesianOdeDimensions = odeDims;
}

std::vector<unsigned> AbstractDomainField::GetCartesianOdeDimensions()
{
    return mCartesianOdeDimensions;
}

void AbstractDomainField::SetCartesianCellScale(std::vector<double> cellScale)
{
    mCartesianCellScaleXY = cellScale;
}

std::vector<double> AbstractDomainField::GetCartesianCellScale()
{
    return mCartesianCellScaleXY;
}

void AbstractDomainField::SetCartesianOdeScale(std::vector<double> odeScale)
{
    mCartesianOdeScaleXY = odeScale;
}

std::vector<double> AbstractDomainField::GetCartesianOdeScale()
{
    return mCartesianOdeScaleXY;
}

void AbstractDomainField::SetCartesianChasteScale(std::vector<double> chasteScale)
{
    mCartesianChasteScaleXY = chasteScale;
}

std::vector<double> AbstractDomainField::GetCartesianChasteScale()
{
    return mCartesianChasteScaleXY;
}

std::vector<std::vector<std::string>> AbstractDomainField::GetCartesianChasteDomain()
{
    return mCartesianChasteDomain;
}

std::vector<std::vector<std::string>> AbstractDomainField::GetCartesianOdeDomain()
{
    return mCartesianOdeDomain;
}

void AbstractDomainField::SetCartesianChasteDomain(std::vector<std::vector<std::string>> cartesianChasteDomain)
{
    mCartesianChasteDomain = cartesianChasteDomain;
}

void AbstractDomainField::SetCartesianOdeDomain(std::vector<std::vector<std::string>> cartesianOdeDomain)
{
    mCartesianOdeDomain = cartesianOdeDomain;
}

void AbstractDomainField::SetCartesianChasteDimensions(std::vector<unsigned> dimensions)
{
    mCartesianChasteDimensions = dimensions;
}

std::vector<unsigned> AbstractDomainField::GetCartesianChasteDimensions()
{
    return mCartesianChasteDimensions;
}

void AbstractDomainField::SetMeshDimensions(std::vector<unsigned> meshDimensions)
{
    mMeshDimensions = meshDimensions;
}

std::vector<unsigned> AbstractDomainField::GetMeshDimensions()
{
    return mMeshDimensions;
}

void AbstractDomainField::SetDomainMesh(MutableMesh<2,2>* p_mesh)
{
    mpMesh = p_mesh;
}

void AbstractDomainField::SetMeshGenerator(HoneycombMeshGenerator* p_generator)
{
    mpGenerator = p_generator;
}

MutableMesh<2,2>* AbstractDomainField::GetDomainMesh()
{
    return mpMesh;
}

HoneycombMeshGenerator* AbstractDomainField::GetMeshGenerator()
{
    return mpGenerator;
}

void AbstractDomainField::SetMeshScale(std::vector<double> meshScale)
{
    mMeshScale = meshScale;
}

std::vector<double> AbstractDomainField::GetMeshScale()
{
    return mMeshScale;
}

void AbstractDomainField::SetDiffusionDatabase(std::vector<std::vector<std::string>> diffusionDatabase)
{
    mDiffusionDatabase = diffusionDatabase;
}

std::vector<std::vector<std::string>> AbstractDomainField::GetDiffusionDatabase()
{
    return mDiffusionDatabase;
}

void AbstractDomainField::SetNodeLabels(std::vector<std::string> nodeLabels)
{
    mLabelledNodes = nodeLabels;
}

std::vector<std::string> AbstractDomainField::GetNodeLabels()
{
    return mLabelledNodes;
}

void AbstractDomainField::SetNodeDomains(std::vector<std::string> nodeDomains)
{
    mDomainNodes = nodeDomains;
}

std::vector<std::string> AbstractDomainField::GetNodeDomains()
{
    return mDomainNodes;
}

std::vector<double> AbstractDomainField::GetInitialNodeConditions()
{
    return mInitalNodeConditions;
}

void AbstractDomainField::SetInitialNodeConditions(std::vector<double> initialNodeConditions)
{
    mInitalNodeConditions = initialNodeConditions;
}

std::vector<std::string> AbstractDomainField::GetBoundaryConditionTypes()
{
    return mBoundaryConditionTypes;
}

std::vector<double> AbstractDomainField::GetBoundaryConditionValues()
{
    return mBoundaryConditionValues;
}


void AbstractDomainField::SetBoundaryConditionTypes(std::vector<std::string> boundaryConditionsTypes)
{
    mBoundaryConditionTypes = boundaryConditionsTypes;
}

void AbstractDomainField::SetBoundaryConditionValues(std::vector<double> boundaryConditionValues)
{
    mBoundaryConditionValues = boundaryConditionValues;
}


std::string AbstractDomainField::ReturnKeyValueFromNodeLabel(std::string nodeLabel)
{
    bool IsFound=false;

    for(unsigned key_index=0; key_index<mOdeKeys.size(); key_index++)
    {
        if(mOdeKeys[key_index][0] == nodeLabel)
        {
            return mOdeKeys[key_index][1];
            IsFound=true;
            break;
        }
    }
    if(!IsFound)
    {
        std::cout<<"Error: AbstractDomainField::ReturnKeyValueFromNodeLabel, node label not found"<<std::endl;
        return "Null";
    }
    return "Null";
}


std::string AbstractDomainField::ReturnDomainKeyFromDomainLabel(std::string domainLabel)
{
    bool IsFound=false;

    for(unsigned key_index=0; key_index<mDomainKeys.size(); key_index++)
    {
        
        if(mDomainKeys[key_index][0] == domainLabel)
        {
            return mDomainKeys[key_index][1];
            IsFound=true;
            break;
        }
    }
    if(!IsFound)
    {
        std::cout<<"Error: AbstractDomainField::ReturnDomainKeyFromDomainLabel, domain label not found"<<std::endl;
        return "Null";
    }
    return "Null";
}