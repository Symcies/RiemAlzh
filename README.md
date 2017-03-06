# RiemAlzh

XXX is a software for the statistical analysis of longitudinal data.
It essentially computes a reconstruction of a mean scenario of evolution out of series of short-term longitudinal data.

It comes with two different types of data :
- Scalar data
- Data distributed on a network

## Requirements
#### tinyxml2 version XX
*XML-related library*: it is used to read the xml parameter files
It is part of the git submodules
Just type 'git submodule init' then 'git submodule update'

#### Armadillo version XX
*Linear Algebra Library*: it is used to work with the vectors and the matrix
Use brew to install it quickly on mac

#### Google test (In the future)
*Testing Library*
Used for unit and fonctional testing

## Installation
1. Open a terminal and type 'git clone https://github.com/Symcies/RiemAlzh'
2. mkdir build WW && cd WW
3. cmake --options
4. make --options

## Run the software
The software is to be run with three files: 'algorithm_settings.xml', 'data_settings.xml' and 'model_settings.xml'
Open a terminal, go to the build directory and type './XX model_settings.xml algorithm_settings.xml data_settings.xml'

## Parameter files

#### model_settings.xml

#### algorithm_settings.xml

#### data_settings.xml

## Examples
The examples are to be run in the terminal "sh launch_simulation.sh"

#### Scalar Inputs
*Multivariate Model*

#### Network Inputs
*Network Mode*

## Documentation & Articles

#### Software Usage Documentation

#### Software Development Documentation

#### Model & Math Documentation


## Contacts
http://www.aramislab.fr/
